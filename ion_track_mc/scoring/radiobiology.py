"""Radiation-quality quantities: LET and single-event microdosimetry.

LET (linear energy transfer) is the energy the projectile loses per unit
path length. Along a track it varies from place to place, so two averages
are used (ICRU definitions; scoring/README_scoring.md §3):
  - track-averaged LET_T: total energy lost / total path length;
  - dose-averaged LET_D: each piece of track weighted by the energy it
    deposits, so high-LET regions count more. Biological effect tracks
    LET_D more closely than LET_T.

Microdosimetry (§4) looks instead at the energy epsilon one particle leaves
in a small sphere (a "site", e.g. a cell nucleus or a chromatin domain):
  - lineal energy y = epsilon / l_bar, with l_bar the mean chord length;
  - specific energy z = epsilon / m, the dose in that one site.

Not implemented: the distributions f(y), d(y) over many events, and
biological models built on them (e.g. MKM).

Units note: LET in eV/nm and keV/um are numerically IDENTICAL
(1 eV/nm = 1000 eV / 1000 nm = 1 keV/um), so no conversion is applied --
values returned here in eV/nm can be read directly as keV/um, the
conventional LET reporting unit.
"""
import numpy as np

from ion_track_mc import constants
from ion_track_mc.physics.medium import Medium
from ion_track_mc.scoring.scorer import Scorer


def track_averaged_LET_from_sums(total_s: float, total_dE: float) -> float:
    """LET_T = total_dE / total_s  [eV/nm], from sums already accumulated.

    Lets `analysis.beam.run_beam` add up the sums track by track and get a
    beam-averaged LET at the end, without keeping every collision."""
    if total_s == 0.0:
        raise ValueError("total path length is zero: track-averaged LET is undefined")
    return total_dE / total_s


def dose_averaged_LET_from_sums(total_segment_E2_over_length: float, total_dE: float) -> float:
    """LET_D = total_segment_E2_over_length / total_dE  [eV/nm], from sums.

    The numerator is sum_k(E_k^2 / l_k) over track segments (see
    `LET_segment_sums` and `dose_averaged_LET`). Like
    `track_averaged_LET_from_sums`, it works on running sums over many
    tracks."""
    if total_dE == 0.0:
        raise ValueError("total delta_E is zero: dose-averaged LET is undefined")
    return total_segment_E2_over_length / total_dE


def LET_segment_sums(
    collision_log: list[tuple[float, float]], segment_length_nm: float
) -> tuple[float, float]:
    """Cut one track into segments of length Delta along its path and return
    (sum_k E_k^2 / l_k, sum_k E_k), the two sums dose-averaged LET needs.

    Collision i sits at path position x_i = s_1 + ... + s_i and belongs to
    segment k = floor(x_i / Delta); E_k is the energy lost in segment k.
    All segments have l_k = Delta, except that the leftover piece at the end
    of the track (shorter than Delta) is merged into the last full segment --
    a tiny last segment would bring back the 1/l blow-up this is avoiding.
    A track shorter than Delta is a single segment of its own length.
    Segments with no collisions have E_k = 0 and add nothing.
    """
    if segment_length_nm <= 0.0:
        raise ValueError("segment_length_nm must be positive")
    total_path = sum(s for s, _ in collision_log)
    # Number of full segments; at least one, even for a very short track.
    n_segments = max(1, int(total_path // segment_length_nm))
    lengths = [segment_length_nm] * n_segments
    # The last segment absorbs the leftover piece at the end of the track.
    lengths[-1] = total_path - (n_segments - 1) * segment_length_nm

    # Walk along the track and add each collision's energy loss to the
    # segment its position falls in.
    energies = [0.0] * n_segments
    x = 0.0
    for s, dE in collision_log:
        x += s
        k = min(int(x // segment_length_nm), n_segments - 1)
        energies[k] += dE

    # L_k * E_k = E_k^2 / l_k: each segment's LET weighted by its energy.
    sum_E2_over_length = sum(E * E / l for E, l in zip(energies, lengths) if E > 0.0)
    return sum_E2_over_length, sum(energies)


def track_averaged_LET(collision_log: list[tuple[float, float]]) -> float:
    """LET_T = sum(delta_E_i) / sum(s_i)  [eV/nm]: total energy lost over
    total path length, i.e. the local LET averaged with path length as the
    weight (ICRU track-averaged LET).

    `collision_log` is `TrackEngine.collision_log`: (step_length_nm,
    delta_E_eV) pairs for the primary's own collisions.
    """
    if not collision_log:
        raise ValueError("empty collision_log: nothing to average")
    total_s = sum(s for s, _ in collision_log)
    total_dE = sum(dE for _, dE in collision_log)
    return track_averaged_LET_from_sums(total_s, total_dE)


def dose_averaged_LET(collision_log: list[tuple[float, float]], segment_length_nm: float) -> float:
    """LET_D = sum_k(L_k * E_k) / sum_k(E_k) = sum_k(E_k^2 / l_k) / sum_k(E_k)
    [eV/nm], L_k = E_k / l_k -- the dose-weighted mean of the local LET over
    track segments of length `segment_length_nm` (see `LET_segment_sums`).
    Weights high-LET segments more heavily than track-averaged LET.

    Why segments and not single collisions: with L_i = delta_E_i / s_i per
    collision, the steps s_i are exponentially distributed, E[1/s] diverges,
    and the estimator never converges -- rare tiny steps dominate it. A fixed
    segment length removes the 1/s. The result depends on the segment
    length by construction (it sets the scale LET is measured at): it tends
    to LET_t as the segment covers the whole track, and back towards the
    divergent single-collision value as it shrinks to the mean free path.
    Pick it much longer than the mean free path (~nm) and much shorter than
    the distance over which LET changes. See scoring/README_scoring.md §3.
    """
    if not collision_log:
        raise ValueError("empty collision_log: nothing to average")
    total_E2_over_length, total_dE = LET_segment_sums(collision_log, segment_length_nm)
    return dose_averaged_LET_from_sums(total_E2_over_length, total_dE)


def mean_chord_length_sphere(diameter_nm: float) -> float:
    """Mean chord length of a sphere, l_bar = (2/3) * diameter.

    Cauchy's formula: for any convex body crossed by straight lines from
    random directions, the mean chord is l_bar = 4V/S. For a sphere,
    4 (pi d^3/6) / (pi d^2) = 2d/3. ICRU uses it to define lineal energy."""
    return (2.0 / 3.0) * diameter_nm


def lineal_energy(epsilon_eV: float, diameter_nm: float) -> float:
    """y = epsilon / l_bar  [eV/nm], the single-event lineal energy for
    energy `epsilon` deposited in one spherical site of `diameter_nm`.

    y is the microdosimetric counterpart of LET: the energy actually left
    in the site, per mean chord length. 1 eV/nm = 1 keV/um."""
    return epsilon_eV / mean_chord_length_sphere(diameter_nm)


def specific_energy(epsilon_eV: float, diameter_nm: float, medium: Medium) -> float:
    """z = epsilon / m  [Gy], the single-event specific energy for energy
    `epsilon` deposited in one spherical site of `diameter_nm` in `medium`."""
    radius_nm = diameter_nm / 2.0
    volume_nm3 = (4.0 / 3.0) * np.pi * radius_nm ** 3
    volume_cm3 = volume_nm3 / constants.CM3_TO_NM3
    mass_kg = medium.density_g_cm3 * volume_cm3 * 1.0e-3   # g -> kg
    energy_J = epsilon_eV * constants.EV_TO_JOULE
    return energy_J / mass_kg


class SphericalSiteScorer(Scorer):
    """Adds up all energy deposited inside one spherical site.

    After a track, `epsilon_eV` is the energy imparted to the site by that
    event, and `y()` / `z(medium)` give its lineal and specific energy.
    Call `reset()` before the next independent event.
    """

    def __init__(self, center_nm: np.ndarray, diameter_nm: float):
        self.center_nm = np.asarray(center_nm, dtype=float)
        self.diameter_nm = diameter_nm
        self.radius_nm = diameter_nm / 2.0
        self.epsilon_eV = 0.0

    def score(self, position: np.ndarray, energy: float) -> None:
        if np.linalg.norm(position - self.center_nm) <= self.radius_nm:
            self.epsilon_eV += energy

    def reset(self) -> None:
        self.epsilon_eV = 0.0

    def y(self) -> float:
        return lineal_energy(self.epsilon_eV, self.diameter_nm)

    def z(self, medium: Medium) -> float:
        return specific_energy(self.epsilon_eV, self.diameter_nm, medium)
