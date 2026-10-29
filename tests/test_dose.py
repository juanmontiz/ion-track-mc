import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.scoring.dose import RadialDoseScorer, VoxelDoseScorer, make_bragg_curve_scorer
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

E_START_EV = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT


def test_score_bins_into_correct_voxel(water):
    scorer = VoxelDoseScorer(water, x_edges=np.array([0.0, 1.0, 2.0]),
                              y_edges=np.array([0.0, 1.0]), z_edges=np.array([0.0, 1.0]))
    scorer.score(np.array([1.5, 0.5, 0.5]), 100.0)
    assert scorer.energy_grid_eV[1, 0, 0] == pytest.approx(100.0)
    assert scorer.energy_grid_eV[0, 0, 0] == 0.0


def test_score_outside_grid_is_dropped(water):
    scorer = VoxelDoseScorer(water, x_edges=np.array([0.0, 1.0]),
                              y_edges=np.array([0.0, 1.0]), z_edges=np.array([0.0, 1.0]))
    scorer.score(np.array([100.0, 100.0, 100.0]), 50.0)
    assert scorer.total_deposited_eV() == 0.0


def test_dose_Gy_matches_manual_calculation(water):
    # Single 1x1x1 nm^3 voxel, water density 1 g/cm^3, deposit exactly 1000 eV.
    scorer = VoxelDoseScorer(water, x_edges=np.array([0.0, 1.0]),
                              y_edges=np.array([0.0, 1.0]), z_edges=np.array([0.0, 1.0]))
    scorer.score(np.array([0.5, 0.5, 0.5]), 1000.0)

    energy_J = 1000.0 * constants.EV_TO_JOULE
    volume_cm3 = 1.0 / constants.CM3_TO_NM3  # 1 nm^3 -> cm^3
    mass_kg = water.density_g_cm3 * volume_cm3 * 1.0e-3
    expected_dose = energy_J / mass_kg

    assert scorer.dose_Gy()[0, 0, 0] == pytest.approx(expected_dose)


def test_voxel_volumes_shape_and_value(water):
    scorer = VoxelDoseScorer(water, x_edges=np.array([0.0, 2.0]),
                              y_edges=np.array([0.0, 3.0]), z_edges=np.array([0.0, 4.0, 8.0]))
    vols = scorer.voxel_volumes_nm3()
    assert vols.shape == (1, 1, 2)
    assert vols[0, 0, 0] == pytest.approx(2.0 * 3.0 * 4.0)


def test_bragg_curve_scorer_full_run_energy_conservation(full_xs, water, full_edep, full_eemean):
    """Grid sized generously (huge z-range and transverse width) so nothing
    escapes -- turns this into a pure energy-conservation check, independent
    of whatever the actual track range happens to be."""
    z_edges = np.linspace(-1.0e9, 1.0e9, 21)
    scorer = make_bragg_curve_scorer(water, z_edges=z_edges, transverse_half_width_nm=1.0e9)

    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, scorers=[scorer], seed=5)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    e0 = p.energy
    spent = engine.transport(p)

    assert scorer.total_deposited_eV() == pytest.approx(e0 - spent.energy)
    assert np.all(scorer.dose_Gy() >= 0.0)


def test_bragg_curve_shape_is_depth_binned(full_xs, water):
    """Sanity check the Bragg-curve convenience constructor actually bins by
    z: two well-separated bins, deposits go to the bin matching their depth."""
    z_edges = np.array([0.0, 500.0, 1000.0])
    scorer = make_bragg_curve_scorer(water, z_edges=z_edges, transverse_half_width_nm=1.0e6)
    scorer.score(np.array([0.0, 0.0, 100.0]), 10.0)
    scorer.score(np.array([0.0, 0.0, 900.0]), 20.0)
    curve = scorer.dose_Gy()[0, 0, :]
    assert curve[0] > 0.0
    assert curve[1] > 0.0
    # Same eV in same-sized bins -> dose ratio should match energy ratio.
    assert curve[1] / curve[0] == pytest.approx(20.0 / 10.0)


def test_radial_scorer_bins_by_perpendicular_distance_to_z_axis(water):
    scorer = RadialDoseScorer(water, r_edges=np.array([0.0, 5.0, 10.0]), length_nm=100.0)
    scorer.score(np.array([2.0, 0.0, 50.0]), 10.0)   # r=2 -> shell 0
    scorer.score(np.array([7.0, 0.0, 50.0]), 20.0)   # r=7 -> shell 1
    scorer.score(np.array([100.0, 0.0, 50.0]), 30.0)  # r=100 -> outside, dropped
    assert scorer.energy_shell_eV[0] == pytest.approx(10.0)
    assert scorer.energy_shell_eV[1] == pytest.approx(20.0)
    assert scorer.total_deposited_eV() == pytest.approx(30.0)


def test_radial_scorer_ignores_position_along_axis(water):
    """z-position (along the default z-axis) must not affect binning -- only
    perpendicular distance matters."""
    scorer = RadialDoseScorer(water, r_edges=np.array([0.0, 5.0]), length_nm=100.0)
    scorer.score(np.array([2.0, 0.0, -1000.0]), 10.0)
    scorer.score(np.array([2.0, 0.0, 1000.0]), 10.0)
    assert scorer.energy_shell_eV[0] == pytest.approx(20.0)


def test_radial_scorer_dose_Gy_matches_manual_calculation(water):
    scorer = RadialDoseScorer(water, r_edges=np.array([0.0, 3.0]), length_nm=10.0)
    scorer.score(np.array([1.0, 0.0, 5.0]), 500.0)

    volume_nm3 = np.pi * (3.0 ** 2 - 0.0 ** 2) * 10.0
    volume_cm3 = volume_nm3 / constants.CM3_TO_NM3
    mass_kg = water.density_g_cm3 * volume_cm3 * 1.0e-3
    expected = (500.0 * constants.EV_TO_JOULE) / mass_kg
    assert scorer.dose_Gy()[0] == pytest.approx(expected)


def test_radial_scorer_full_run_energy_conservation(full_xs, water, full_edep, full_eemean):
    """Generously sized radial grid (huge outer radius, huge length) so
    nothing escapes -- pure energy-conservation check."""
    scorer = RadialDoseScorer(water, r_edges=np.array([0.0, 1.0e9]), length_nm=1.0e9)
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, scorers=[scorer], seed=5)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    e0 = p.energy
    spent = engine.transport(p)
    assert scorer.total_deposited_eV() == pytest.approx(e0 - spent.energy)
