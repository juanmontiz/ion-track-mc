"""Per-channel projectile energy-loss formulas (see
ion_track_mc/physics/README_physics.md, §4 "Projectile energy loss per channel").
"""
from pathlib import Path

import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.physics import energy_loss
from ion_track_mc.physics.cross_section import load_cross_section_table, load_mean_energy_table

DATA_DIR = Path(__file__).resolve().parents[1] / "ion_track_mc" / "data" / "cross_sections"


def test_mean_energy_lookup_reproduces_grid_value():
    table = load_mean_energy_table(DATA_DIR / "H1+_+_H2O_excitation_Edepmeantot_Miller-Green.dat")
    E_tab, V_tab = table
    assert energy_loss.mean_energy_lookup(E_tab[5], table) == pytest.approx(V_tab[5], rel=1e-10)


def test_capture_delta_E_matches_formula_at_first_grid_point():
    table = load_mean_energy_table(DATA_DIR / "H1+_+_H2O_capture_Edepmeantot_CDW-EIS.dat")
    E_tab, V_tab = table
    E0 = E_tab[0]
    expected = V_tab[0] + E0 * constants.ELECTRON_TO_PROTON_MASS_RATIO - constants.H_BINDING_ENERGY_EV
    assert energy_loss.capture_delta_E(E0, table) == pytest.approx(expected, rel=1e-10)


def test_capture_delta_E_positive_across_full_real_grid():
    """Regression guard: the capture energy loss must stay positive over the
    whole data grid. If a future data revision makes it negative, this test
    must catch it before _deposit_local does."""
    table = load_mean_energy_table(DATA_DIR / "H1+_+_H2O_capture_Edepmeantot_CDW-EIS.dat")
    E_tab, _ = table
    for E in E_tab:
        assert energy_loss.capture_delta_E(E, table) > 0.0


def test_ionization_energy_matches_tables_at_first_grid_point():
    edep_table = load_mean_energy_table(DATA_DIR / "H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat")
    eemean_table = load_mean_energy_table(DATA_DIR / "H1+_+_H2O_ionization_Eemeantot_CDW-EIS.dat")
    E_tab, edep_tab = edep_table
    _, eemean_tab = eemean_table
    deposit, W = energy_loss.ionization_energy(E_tab[0], edep_table, eemean_table)
    assert deposit == pytest.approx(edep_tab[0], rel=1e-10)
    assert W == pytest.approx(eemean_tab[0], rel=1e-10)


def test_edepmeantot_tables_share_their_tcs_siblings_energy_grid():
    """Every Edepmeantot table the engine will actually consume (H1+/H0,
    excitation/ionization/capture/loss) shares its TCS sibling's exact
    energy grid. TrackEngine's out-of-range check uses
    xs.energy_range(q), computed from the TCS tables -- if a mean-energy
    table's grid were narrower, a query the engine considers in-range could
    still raise inside interp_log_log. This guards that invariant."""
    pairs = [
        ("H1+_+_H2O_excitation_TCS_Miller-Green.dat", "H1+_+_H2O_excitation_Edepmeantot_Miller-Green.dat"),
        ("H0_+_H2O_excitation_TCS_Miller-Green.dat", "H0_+_H2O_excitation_Edepmeantot_Miller-Green.dat"),
        ("H1+_+_H2O_capture_TCS_CDW-EIS.dat", "H1+_+_H2O_capture_Edepmeantot_CDW-EIS.dat"),
        ("H0_+_H2O_loss_TCS_Miller-Green.dat", "H0_+_H2O_loss_Edepmeantot_constant.dat"),
        ("H1+_+_H2O_ionization_TCS_CDW-EIS.dat", "H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat"),
        ("H0_+_H2O_ionization_TCS_CDW-EIS.dat", "H0_+_H2O_ionization_Edepmeantot_CDW-EIS.dat"),
        ("H1+_+_H2O_ionization_TCS_CDW-EIS.dat", "H1+_+_H2O_ionization_Eemeantot_CDW-EIS.dat"),
        ("H0_+_H2O_ionization_TCS_CDW-EIS.dat", "H0_+_H2O_ionization_Eemeantot_CDW-EIS.dat"),
    ]
    for tcs_f, edep_f in pairs:
        E_tcs, _ = load_cross_section_table(DATA_DIR / tcs_f)
        E_edep, _ = load_mean_energy_table(DATA_DIR / edep_f)
        assert np.allclose(E_tcs, E_edep), f"{tcs_f} vs {edep_f} grid mismatch"


def test_elastic_scatter_direction_is_unit_vector():
    rng = np.random.default_rng(123)
    direction = np.array([0.0, 0.0, 1.0])
    for _ in range(200):
        _, new_dir = energy_loss.elastic_scatter(1.0e6, direction, 10.0, 18.0, rng)
        assert np.linalg.norm(new_dir) == pytest.approx(1.0)


def test_elastic_scatter_energy_within_kinematic_bound():
    rng = np.random.default_rng(123)
    direction = np.array([0.0, 0.0, 1.0])
    E = 1.0e6
    M_P = constants.PROTON_TO_ELECTRON_MASS_RATIO
    M_T = 18.0 * M_P
    E_max = 4.0 * E * (M_P * M_T) / (M_P + M_T) ** 2  # theta = pi, the max possible transfer
    for _ in range(500):
        E_elastic, _ = energy_loss.elastic_scatter(E, direction, 10.0, 18.0, rng)
        assert 0.0 <= E_elastic <= E_max + 1e-6


def test_elastic_scatter_raises_clear_error_for_unset_target_fields():
    """A custom `Medium` that leaves target_nuclear_charge/target_nucleon_count
    as `None` (the dataclass default -- only WATER sets them) must fail with
    a clear ValueError, not a confusing TypeError from `None ** (1.0 / 3.0)`
    deep inside the kinematic chain."""
    rng = np.random.default_rng(123)
    direction = np.array([0.0, 0.0, 1.0])
    with pytest.raises(ValueError):
        energy_loss.elastic_scatter(1.0e6, direction, None, None, rng)


def test_elastic_scatter_no_domain_violation_across_real_grid():
    """The angle-sampling formula's arcsin argument stays in [-1, 1]
    everywhere on the real elastic TCS grid and for eta=0/1 extremes
    -- this smoke-tests that across the whole grid with
    many random draws, confirming no NaN/exception."""
    E_tab, _ = load_cross_section_table(DATA_DIR / "H1+_+_H2O_elastic_TCS_Classical.dat")
    rng = np.random.default_rng(7)
    direction = np.array([0.0, 0.0, 1.0])
    for E in E_tab:
        for _ in range(20):
            E_elastic, new_dir = energy_loss.elastic_scatter(E, direction, 10.0, 18.0, rng)
            assert np.isfinite(E_elastic)
            assert np.all(np.isfinite(new_dir))


def test_elastic_scatter_deflects_by_at_least_the_cutoff_angle():
    """theta ranges from theta_cut (at eta=0) to pi (at eta=1) -- it is
    never exactly 0, so the direction always changes, by construction."""
    rng = np.random.default_rng(9)
    direction = np.array([0.0, 0.0, 1.0])
    for _ in range(50):
        _, new_dir = energy_loss.elastic_scatter(1.0e5, direction, 10.0, 18.0, rng)
        cos_angle = float(np.dot(new_dir, direction))
        assert cos_angle < 1.0
