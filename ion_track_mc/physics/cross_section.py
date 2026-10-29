"""Cross sections: how likely each collision process is, at each energy.

A cross section sigma(E) [nm^2] is the effective target area one water
molecule presents to the projectile for a given process. The code reads them
from tables (computed with CDW-EIS, Miller-Green and classical models, see
data/README_data_format.md) and interpolates between tabulated energies.

This module holds:
  - `CrossSectionSet`, the interface the transport engine talks to;
  - the table loaders, which convert to internal units once, on load;
  - `TabulatedCrossSectionSet`, the implementation used everywhere.

See physics/README_physics.md §2.
"""
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np

from ion_track_mc import constants
from ion_track_mc.physics.channels import Channel


class CrossSectionSet(ABC):
    """Contract for the cross-section set of a (projectile, charge, medium).

    The transport engine only ever calls these four methods, so a different
    source of cross sections (another model, another medium) can be plugged
    in without touching the engine.

    Implementations MUST raise (never silently extrapolate) when E is outside
    the tabulated validity range [E_min, E_max]: outside the tables the
    physics is simply unknown to the code.
    """

    @abstractmethod
    def channels(self, q: int) -> list[Channel]:
        """Channels available for charge state q."""
        ...

    @abstractmethod
    def sigma(self, channel: Channel, E: float, q: int) -> float:
        """Cross section [nm^2] for `channel` at energy E [eV], charge q."""
        ...

    @abstractmethod
    def sigma_total(self, E: float, q: int) -> float:
        """sigma_tot = sum of sigma over all open channels [nm^2]."""
        ...

    @abstractmethod
    def energy_range(self, q: int) -> tuple[float, float]:
        """Return (E_min, E_max) in eV of the validity range for charge q."""
        ...


def _parse_two_column_table(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse a two-column whitespace-separated table, raw values (no unit conversion).

    File format (see data/README_data_format.md):
        header line     -> optional, auto-detected: the first non-blank,
                            non-comment line is skipped if it does NOT parse
                            as two floats (some files carry a text header
                            line, others do not).
        blank lines     -> skipped
        lines starting with '#' -> skipped (comments)
        remaining lines -> "<x>  <y>", whitespace-separated.

    Returns (x, y), raw values as tabulated, in file order.
    """
    xs: list[float] = []
    ys: list[float] = []
    header_checked = False
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if not header_checked:
            header_checked = True
            try:
                x, y = (float(tok) for tok in stripped.split())
            except ValueError:
                continue  # this line is a text header; skip it
            xs.append(x)
            ys.append(y)
            continue
        x, y = (float(tok) for tok in stripped.split())
        xs.append(x)
        ys.append(y)
    return np.array(xs), np.array(ys)


def load_cross_section_table(
    path: Path,
    energy_unit: str = constants.DEFAULT_ENERGY_UNIT,
    cross_section_unit: str = constants.DEFAULT_CROSS_SECTION_UNIT,
) -> tuple[np.ndarray, np.ndarray]:
    """Parse a two-column TCS file and convert to internal units.

    `energy_unit` / `cross_section_unit` select the conversion factor from
    `constants.ENERGY_UNIT_TO_EV` / `constants.CROSS_SECTION_UNIT_TO_NM2`;
    they default to the units the shipped tables use (eV and 1e-16 cm^2).

    Returns (E_eV, sigma_nm2), both ascending, converted once on load.
    """
    E_raw, sigma_raw = _parse_two_column_table(path)

    # Interpolation needs a strictly increasing energy grid.
    if np.any(np.diff(E_raw) <= 0):
        raise ValueError(f"{path}: energies must be strictly increasing")

    # Convert once here, so everything downstream works in eV and nm^2.
    E_eV = E_raw * constants.ENERGY_UNIT_TO_EV[energy_unit]
    sigma_nm2 = sigma_raw * constants.CROSS_SECTION_UNIT_TO_NM2[cross_section_unit]
    return E_eV, sigma_nm2


def load_mean_energy_table(
    path: Path,
    energy_unit: str = constants.DEFAULT_ENERGY_UNIT,
    value_unit: str = constants.DEFAULT_ENERGY_UNIT,
) -> tuple[np.ndarray, np.ndarray]:
    """Parse a two-column mean-energy file (Edepmeantot / Eemeantot) and convert.

    Same file format as `load_cross_section_table`, but the value column is
    itself an energy (mean deposited or mean emitted energy per collision),
    so it is converted via `constants.ENERGY_UNIT_TO_EV` rather than the
    cross-section unit table.

    Returns (E_eV, value_eV), both ascending, converted once on load.
    """
    E_raw, value_raw = _parse_two_column_table(path)

    if np.any(np.diff(E_raw) <= 0):
        raise ValueError(f"{path}: energies must be strictly increasing")

    E_eV = E_raw * constants.ENERGY_UNIT_TO_EV[energy_unit]
    value_eV = value_raw * constants.ENERGY_UNIT_TO_EV[value_unit]
    return E_eV, value_eV


def interp_log_log(E: float, E_table: np.ndarray, value_table: np.ndarray) -> float:
    """Linear interpolation in log-log space: a straight line between
    neighbouring table points in (log E, log value).

    Cross sections and mean energies vary roughly as power laws of E over
    decades, so a straight line in log-log follows them far better than
    linear interpolation would. Same rule as TabulatedCrossSectionSet.sigma;
    physics/energy_loss.py uses it for the mean-energy tables.

    Raises ValueError if E is outside [E_table[0], E_table[-1]] (no
    extrapolation) or if any value in value_table is <= 0
    (log-log interpolation is undefined there -- e.g. never call this on
    the loss channel's all-zero Eemeantot file; loss has no secondary
    electron in this model and nothing should query that table).
    """
    E_min, E_max = E_table[0], E_table[-1]
    if E < E_min or E > E_max:
        raise ValueError(f"E={E} eV outside tabulated range [{E_min}, {E_max}] eV")
    if np.any(value_table <= 0):
        raise ValueError("interp_log_log requires strictly positive values for log-log interpolation")
    log_E = np.log(E_table)
    log_V = np.log(value_table)
    # Interpolate log(value) linearly in log(E), then undo the log.
    return float(np.exp(np.interp(np.log(E), log_E, log_V)))


class TabulatedCrossSectionSet(CrossSectionSet):
    """Cross sections interpolated from tabulated (energy, sigma) grids.

    The caller supplies one (E_eV, sigma_nm2) table per (charge, channel),
    usually built with `load_cross_section_table` (see the notebooks). The
    channels open for a charge state are exactly the ones a table was given
    for: H+ has ionization, excitation, capture and elastic; H0 has
    ionization, excitation, loss and elastic.

    Interpolation: linear in log-log space (see `interp_log_log`).
    """

    def __init__(self, tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]]):
        self._tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]] = {}
        self._log_tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]] = {}
        for key, (E_eV, sigma_nm2) in tables.items():
            self._tables[key] = (E_eV, sigma_nm2)
            # Logs are taken once here, not on every sigma() call.
            self._log_tables[key] = (np.log(E_eV), np.log(sigma_nm2))
        # Sorted channel list per charge state, built once here instead of on
        # every channels() call (it is queried several times per collision).
        self._channels: dict[int, tuple[Channel, ...]] = {}
        for q in {charge for (charge, _channel) in self._tables}:
            self._channels[q] = tuple(sorted(
                {channel for (charge, channel) in self._tables if charge == q},
                key=lambda c: c.value,
            ))

    def channels(self, q: int) -> list[Channel]:
        # A fresh list each call, so a caller modifying it cannot corrupt the cache.
        return list(self._channels.get(q, ()))

    def energy_range(self, q: int) -> tuple[float, float]:
        chans = self.channels(q)
        if not chans:
            raise ValueError(f"no cross-section tables loaded for q={q}")
        # Valid range = where EVERY channel is tabulated: the overlap of all
        # the channels' energy grids.
        e_mins = [self._tables[(q, c)][0][0] for c in chans]
        e_maxs = [self._tables[(q, c)][0][-1] for c in chans]
        return max(e_mins), min(e_maxs)

    def sigma(self, channel: Channel, E: float, q: int) -> float:
        key = (q, channel)
        if key not in self._tables:
            raise KeyError(f"channel {channel} not open for q={q}")
        E_eV, _ = self._tables[key]
        E_min, E_max = E_eV[0], E_eV[-1]
        if E < E_min or E > E_max:
            raise ValueError(
                f"E={E} eV outside tabulated range [{E_min}, {E_max}] eV "
                f"for channel={channel}, q={q}"
            )
        log_E, log_sigma = self._log_tables[key]
        return float(np.exp(np.interp(np.log(E), log_E, log_sigma)))

    def sigma_total(self, E: float, q: int) -> float:
        return sum(self.sigma(c, E, q) for c in self.channels(q))
