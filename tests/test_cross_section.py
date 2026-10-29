"""Unit tests for the tabulated cross-section loader and interpolation.

Uses the real data sample shipped in
ion_track_mc/data/cross_sections/H1+_+_H2O_ionization_TCS_CDW-EIS.dat
(collision energy in eV, cross section in units of 1e-16 cm^2 -- see
constants.DEFAULT_ENERGY_UNIT / DEFAULT_CROSS_SECTION_UNIT).
"""
from pathlib import Path

import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import (
    TabulatedCrossSectionSet,
    interp_log_log,
    load_cross_section_table,
    load_mean_energy_table,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "ion_track_mc" / "data" / "cross_sections"
SAMPLE_FILE = DATA_DIR / "H1+_+_H2O_ionization_TCS_CDW-EIS.dat"


def test_sample_file_exists():
    assert SAMPLE_FILE.exists()


def test_load_converts_units_correctly():
    E_eV, sigma_nm2 = load_cross_section_table(SAMPLE_FILE)

    # First raw row: 8.00E+03 eV, 1.52664862E-01 (in units of 1e-16 cm^2).
    expected_E0 = 8.00e3 * constants.ENERGY_UNIT_TO_EV["eV"]
    expected_sigma0 = 1.52664862e-01 * constants.CROSS_SECTION_UNIT_TO_NM2["1e-16cm2"]
    assert E_eV[0] == pytest.approx(expected_E0, rel=1e-6)
    assert sigma_nm2[0] == pytest.approx(expected_sigma0, rel=1e-6)

    # 40 data rows (headerless file -- see README_data_format.md).
    assert len(E_eV) == 40
    assert len(sigma_nm2) == 40


def test_load_energies_strictly_increasing():
    E_eV, _ = load_cross_section_table(SAMPLE_FILE)
    assert np.all(np.diff(E_eV) > 0)


def test_load_respects_explicit_units():
    # Same file, re-interpreted as keV/u + cm^2 (the old convention), must
    # give a different, internally-consistent result -- i.e. unit selection
    # is a real parameter, not hardcoded.
    E_eV_default, sigma_nm2_default = load_cross_section_table(SAMPLE_FILE)
    E_eV_other, sigma_nm2_other = load_cross_section_table(
        SAMPLE_FILE, energy_unit="keV/u", cross_section_unit="cm2"
    )
    assert not np.allclose(E_eV_default, E_eV_other)
    assert not np.allclose(sigma_nm2_default, sigma_nm2_other)


def _make_xs_set() -> TabulatedCrossSectionSet:
    table = load_cross_section_table(SAMPLE_FILE)
    return TabulatedCrossSectionSet({(1, Channel.IONIZATION): table})


def test_channels_open_for_q():
    xs = _make_xs_set()
    assert xs.channels(1) == [Channel.IONIZATION]
    assert xs.channels(0) == []


def test_sigma_matches_table_at_grid_point():
    xs = _make_xs_set()
    E_eV, sigma_nm2 = load_cross_section_table(SAMPLE_FILE)
    # Exact grid point: interpolation should reproduce the tabulated value.
    got = xs.sigma(Channel.IONIZATION, E_eV[5], q=1)
    assert got == pytest.approx(sigma_nm2[5], rel=1e-10)


def test_sigma_total_equals_sigma_for_single_channel():
    xs = _make_xs_set()
    E_eV, _ = load_cross_section_table(SAMPLE_FILE)
    E_mid = E_eV[10]
    assert xs.sigma_total(E_mid, q=1) == pytest.approx(
        xs.sigma(Channel.IONIZATION, E_mid, q=1)
    )


def test_sigma_raises_below_range():
    xs = _make_xs_set()
    E_eV, _ = load_cross_section_table(SAMPLE_FILE)
    with pytest.raises(ValueError):
        xs.sigma(Channel.IONIZATION, E_eV[0] * 0.5, q=1)


def test_sigma_raises_above_range():
    xs = _make_xs_set()
    E_eV, _ = load_cross_section_table(SAMPLE_FILE)
    with pytest.raises(ValueError):
        xs.sigma(Channel.IONIZATION, E_eV[-1] * 2.0, q=1)


def test_sigma_raises_for_unavailable_channel():
    xs = _make_xs_set()
    with pytest.raises(KeyError):
        xs.sigma(Channel.EXCITATION, 1000.0, q=1)


def test_energy_range():
    xs = _make_xs_set()
    E_eV, _ = load_cross_section_table(SAMPLE_FILE)
    e_min, e_max = xs.energy_range(1)
    assert e_min == pytest.approx(E_eV[0])
    assert e_max == pytest.approx(E_eV[-1])


def test_interpolation_is_log_log_linear_between_points():
    xs = _make_xs_set()
    E_eV, sigma_nm2 = load_cross_section_table(SAMPLE_FILE)
    # Midpoint in log-E between two grid points must satisfy the log-log
    # linear interpolation formula exactly.
    E_lo, E_hi = E_eV[3], E_eV[4]
    s_lo, s_hi = sigma_nm2[3], sigma_nm2[4]
    log_E_mid = 0.5 * (np.log(E_lo) + np.log(E_hi))
    E_mid = np.exp(log_E_mid)
    expected_log_sigma = 0.5 * (np.log(s_lo) + np.log(s_hi))
    got = xs.sigma(Channel.IONIZATION, E_mid, q=1)
    assert np.log(got) == pytest.approx(expected_log_sigma, rel=1e-10)


# --- Mean-energy tables (Edepmeantot / Eemeantot) --------------------------
# Only the loader is exercised here; the engine's use of these tables is
# covered by test_energy_loss.py.

EDEP_FILE = DATA_DIR / "H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat"


def test_mean_energy_file_exists():
    assert EDEP_FILE.exists()


def test_load_mean_energy_table_converts_units():
    E_eV, value_eV = load_mean_energy_table(EDEP_FILE)
    # First raw row: 8.000...E+03 eV, 1.425434242818631070e+01 eV.
    assert E_eV[0] == pytest.approx(8.0e3, rel=1e-6)
    assert value_eV[0] == pytest.approx(14.25434242818631070, rel=1e-6)
    assert len(E_eV) == len(value_eV) == 40


def test_load_mean_energy_table_energies_strictly_increasing():
    E_eV, _ = load_mean_energy_table(EDEP_FILE)
    assert np.all(np.diff(E_eV) > 0)


def test_load_loss_constant_energy_tables():
    # User-decided policy: loss deposits a fixed 13.6 eV (H ground-state
    # binding energy) and emits 0 eV, tabulated as flat files on the same
    # energy grid as the loss TCS table so the loader needs no special case.
    E_edep, dep = load_mean_energy_table(DATA_DIR / "H0_+_H2O_loss_Edepmeantot_constant.dat")
    E_eemean, eemean = load_mean_energy_table(DATA_DIR / "H0_+_H2O_loss_Eemeantot_constant.dat")
    E_tcs, _ = load_cross_section_table(DATA_DIR / "H0_+_H2O_loss_TCS_Miller-Green.dat")
    assert np.allclose(E_edep, E_tcs)
    assert np.allclose(E_eemean, E_tcs)
    assert np.all(dep == pytest.approx(13.6))
    assert np.all(eemean == pytest.approx(0.0))


def test_electron_impact_tables_read_as_plain_ev_not_kev_per_u():
    # "keV/u" (energy per nucleon) presupposes a projectile with a nucleon
    # count -- meaningless for an electron. e- tables must load with the
    # default "eV" unit (a no-op conversion), not "keV/u" (which would
    # inflate these by a further x1000): the BEB ionization threshold is a
    # known ~13.1 eV, not ~13.1 keV, and the Olivero excitation threshold is
    # a known ~7.4 eV, not ~7.4 keV.
    E_ion, _ = load_cross_section_table(DATA_DIR / "e-_+_H2O_ionization_TCS_BEB.dat")
    E_exc, _ = load_cross_section_table(DATA_DIR / "e-_+_H2O_excitation_TCS_Olivero.dat")
    assert E_ion[0] == pytest.approx(13.1, rel=1e-3)
    assert E_exc[0] == pytest.approx(7.4, rel=1e-3)


def test_load_electron_impact_excitation_edepmeantot_real():
    # e- excitation Edepmeantot is real now (Olivero model).
    # NOTE: its energy grid does NOT match the excitation TCS file's grid
    # (7.4 eV threshold) -- it matches the *ionization* TCS/BEB grid instead
    # (13.1 eV threshold, 101 points, same upper bound 8.7 MeV). Documented
    # here rather than assumed silently; a future consumer combining this
    # with the excitation TCS table needs its own interpolation, not a
    # pointwise zip.
    E_exc_dep, dep_exc = load_mean_energy_table(
        DATA_DIR / "e-_+_H2O_excitation_Edepmeantot_Olivero.dat"
    )
    E_ion_tcs, _ = load_cross_section_table(
        DATA_DIR / "e-_+_H2O_ionization_TCS_BEB.dat"
    )
    assert len(E_exc_dep) == 101
    assert np.allclose(E_exc_dep, E_ion_tcs)
    assert np.all(dep_exc > 0.0)
    assert np.all(np.diff(E_exc_dep) > 0.0)


def test_load_electron_impact_ionization_beb_tables():
    # e- ionization TCS/Edepmeantot/Eemeantot are real (BEB model),
    # replacing the earlier placeholders -- all three share one energy grid.
    E_tcs, sigma = load_cross_section_table(
        DATA_DIR / "e-_+_H2O_ionization_TCS_BEB.dat"
    )
    E_dep, dep_ion = load_mean_energy_table(
        DATA_DIR / "e-_+_H2O_ionization_Edepmeantot_BEB.dat"
    )
    E_eemean, eemean_ion = load_mean_energy_table(
        DATA_DIR / "e-_+_H2O_ionization_Eemeantot_BEB.dat"
    )
    assert np.allclose(E_dep, E_tcs)
    assert np.allclose(E_eemean, E_tcs)
    assert len(E_tcs) == 101
    assert np.all(sigma > 0.0)
    assert np.all(dep_ion > 0.0)
    assert np.all(eemean_ion >= 0.0)


# --- Full (q=1, q=0) channel set, using the real CDW-EIS/Miller-Green/Classical
# files - Real TCS tables for all five (q=1: ionization/excitation/capture/
# elastic; q=0: ionization/excitation/loss/elastic) transport-relevant
# channels (q=1 -> ionization+excitation+capture, q=0 ->
# ionization+excitation+loss, both plus elastic per §4.1/§4.3:
# elastic contributes to sigma_tot/path sampling AND to real energy deposit
# + direction change, see physics/energy_loss.py::elastic_scatter).

_FULL_TABLE_FILES = {
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_TCS_CDW-EIS.dat",
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_TCS_Miller-Green.dat",
    (1, Channel.CAPTURE): "H1+_+_H2O_capture_TCS_CDW-EIS.dat",
    (1, Channel.ELASTIC): "H1+_+_H2O_elastic_TCS_Classical.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_TCS_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_TCS_Miller-Green.dat",
    (0, Channel.LOSS): "H0_+_H2O_loss_TCS_Miller-Green.dat",
    (0, Channel.ELASTIC): "H0_+_H2O_elastic_TCS_Classical.dat",
}


def _make_full_xs_set() -> TabulatedCrossSectionSet:
    tables = {
        key: load_cross_section_table(DATA_DIR / fname)
        for key, fname in _FULL_TABLE_FILES.items()
    }
    return TabulatedCrossSectionSet(tables)


def test_full_set_all_files_exist():
    for fname in _FULL_TABLE_FILES.values():
        assert (DATA_DIR / fname).exists()


def test_full_set_channels_per_charge_state():
    xs = _make_full_xs_set()
    assert xs.channels(1) == [
        Channel.CAPTURE, Channel.ELASTIC, Channel.EXCITATION, Channel.IONIZATION,
    ]
    assert xs.channels(0) == [
        Channel.ELASTIC, Channel.EXCITATION, Channel.IONIZATION, Channel.LOSS,
    ]


def test_full_set_sigma_total_sums_all_open_channels():
    xs = _make_full_xs_set()
    E_mid = 5.0e5  # within the tabulated [8 keV, 300 MeV] range
    expected_q1 = sum(xs.sigma(c, E_mid, q=1) for c in xs.channels(1))
    assert xs.sigma_total(E_mid, q=1) == pytest.approx(expected_q1)
    expected_q0 = sum(xs.sigma(c, E_mid, q=0) for c in xs.channels(0))
    assert xs.sigma_total(E_mid, q=0) == pytest.approx(expected_q0)


def test_full_set_all_sigmas_positive_across_grid():
    xs = _make_full_xs_set()
    for q in (0, 1):
        e_min, e_max = xs.energy_range(q)
        probe_Es = np.geomspace(e_min, e_max, 25)
        for c in xs.channels(q):
            for E in probe_Es:
                assert xs.sigma(c, E, q) > 0.0


# --- Log-log interpolation (shared helper for energy-loss physics) ----------


def test_interp_log_log_reproduces_table_value_at_grid_point():
    E_table = np.array([10.0, 100.0, 1000.0, 10000.0])
    V_table = np.array([1.0, 2.0, 5.0, 20.0])
    assert interp_log_log(100.0, E_table, V_table) == pytest.approx(2.0)


def test_interp_log_log_matches_manual_midpoint():
    E_table = np.array([10.0, 100.0, 1000.0, 10000.0])
    V_table = np.array([1.0, 2.0, 5.0, 20.0])
    E_mid = np.sqrt(100.0 * 1000.0)  # geometric mean = midpoint in log-E
    expected = np.exp(0.5 * (np.log(2.0) + np.log(5.0)))
    assert interp_log_log(E_mid, E_table, V_table) == pytest.approx(expected)


def test_interp_log_log_raises_outside_range():
    E_table = np.array([10.0, 100.0])
    V_table = np.array([1.0, 2.0])
    with pytest.raises(ValueError):
        interp_log_log(5.0, E_table, V_table)
    with pytest.raises(ValueError):
        interp_log_log(200.0, E_table, V_table)


def test_interp_log_log_raises_on_nonpositive_values():
    E_table = np.array([10.0, 100.0])
    V_table = np.array([1.0, 0.0])
    with pytest.raises(ValueError):
        interp_log_log(50.0, E_table, V_table)


def test_interp_log_log_matches_real_excitation_table_at_first_grid_point():
    E_tab, V_tab = load_mean_energy_table(
        DATA_DIR / "H1+_+_H2O_excitation_Edepmeantot_Miller-Green.dat"
    )
    assert interp_log_log(E_tab[0], E_tab, V_tab) == pytest.approx(V_tab[0], rel=1e-10)


def test_channels_cache_not_corrupted_by_caller_mutation(full_xs):
    """`channels(q)` is cached at construction; mutating a returned list
    must not change what later calls return."""
    first = full_xs.channels(1)
    first.append(Channel.LOSS)
    assert Channel.LOSS not in full_xs.channels(1)
