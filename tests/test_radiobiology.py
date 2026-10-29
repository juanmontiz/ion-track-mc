import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.scoring.radiobiology import (
    SphericalSiteScorer,
    dose_averaged_LET,
    dose_averaged_LET_from_sums,
    lineal_energy,
    mean_chord_length_sphere,
    specific_energy,
    track_averaged_LET,
    track_averaged_LET_from_sums,
)
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

E_START_EV = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT


def test_track_averaged_LET_known_value():
    log = [(10.0, 50.0), (20.0, 50.0)]
    assert track_averaged_LET(log) == pytest.approx(100.0 / 30.0)


def test_dose_averaged_LET_known_value_with_merged_trailing_segment():
    """Path positions x = 30, 60, 90, 115 nm; Delta = 50 nm -> 2 segments,
    the trailing 15 nm merged into the last one (lengths 50, 65).
    Collision at 115 nm (would be segment 2) goes to segment 1.
    E = [10, 20+30+40] = [10, 90]; L_D = (10^2/50 + 90^2/65) / 100."""
    log = [(30.0, 10.0), (30.0, 20.0), (30.0, 30.0), (25.0, 40.0)]
    expected = (10.0**2 / 50.0 + 90.0**2 / 65.0) / 100.0
    assert dose_averaged_LET(log, segment_length_nm=50.0) == pytest.approx(expected)


def test_dose_averaged_LET_known_value_with_empty_segments():
    """x = 10, 210 nm; Delta = 50 -> 4 segments, lengths 50, 50, 50, 60
    (trailing 10 nm merged). E = [5, 0, 0, 7]; empty segments add nothing
    to the numerator."""
    log = [(10.0, 5.0), (200.0, 7.0)]
    expected = (5.0**2 / 50.0 + 7.0**2 / 60.0) / 12.0
    assert dose_averaged_LET(log, segment_length_nm=50.0) == pytest.approx(expected)


def test_dose_averaged_LET_equals_track_averaged_when_segment_covers_track():
    """One segment spanning the whole track: L_D = E^2/X / E = E/X = L_T."""
    log = [(10.0, 50.0), (20.0, 50.0), (5.0, 3.0)]
    assert dose_averaged_LET(log, segment_length_nm=35.0) == pytest.approx(track_averaged_LET(log))
    assert dose_averaged_LET(log, segment_length_nm=1.0e6) == pytest.approx(track_averaged_LET(log))


def test_dose_averaged_LET_converges_with_exponential_steps():
    """The old single-collision estimator sum(dE^2/s)/sum(dE) has infinite
    expectation for exponential s (E[1/s] diverges) and does not settle as
    collisions accumulate. Segment-averaged L_D must: 1e5 and 1e6
    collisions of the same synthetic process agree within 2%."""
    rng = np.random.default_rng(0)

    def synthetic_log(n):
        s = rng.exponential(2.0, n)          # mean free path 2 nm
        dE = rng.exponential(50.0, n)        # mean loss 50 eV
        return list(zip(s.tolist(), dE.tolist()))

    small = dose_averaged_LET(synthetic_log(100_000), segment_length_nm=100.0)
    large = dose_averaged_LET(synthetic_log(1_000_000), segment_length_nm=100.0)
    assert small == pytest.approx(large, rel=0.02)


def test_dose_averaged_LET_rejects_nonpositive_segment_length():
    with pytest.raises(ValueError):
        dose_averaged_LET([(10.0, 5.0)], segment_length_nm=0.0)


def test_dose_averaged_geq_track_averaged():
    """ICRU property (Cauchy-Schwarz): dose-averaged LET >= track-averaged LET."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        log = [(float(rng.uniform(1, 50)), float(rng.uniform(0, 200))) for _ in range(10)]
        if sum(dE for _, dE in log) == 0.0:
            continue
        for delta in (5.0, 20.0, 100.0):
            assert dose_averaged_LET(log, segment_length_nm=delta) >= track_averaged_LET(log) - 1e-9


def test_LET_raises_on_empty_log():
    with pytest.raises(ValueError):
        track_averaged_LET([])
    with pytest.raises(ValueError):
        dose_averaged_LET([], segment_length_nm=10.0)


def test_LET_from_sums_matches_list_based_functions():
    log = [(10.0, 50.0), (20.0, 50.0), (5.0, 0.0)]  # includes a zero-dE collision
    total_s = sum(s for s, _ in log)
    total_dE = sum(dE for _, dE in log)
    # Delta = 15: x = 10, 30, 35 -> 2 segments, lengths 15, 20; E = [50, 50]
    total_segment_E2_over_length = 50.0**2 / 15.0 + 50.0**2 / 20.0

    assert track_averaged_LET_from_sums(total_s, total_dE) == pytest.approx(track_averaged_LET(log))
    assert dose_averaged_LET_from_sums(total_segment_E2_over_length, total_dE) == pytest.approx(
        dose_averaged_LET(log, segment_length_nm=15.0))


def test_LET_from_sums_raises_on_zero_denominator():
    with pytest.raises(ValueError):
        track_averaged_LET_from_sums(0.0, 0.0)
    with pytest.raises(ValueError):
        dose_averaged_LET_from_sums(0.0, 0.0)


def test_dose_averaged_LET_raises_when_total_dE_zero():
    with pytest.raises(ValueError):
        dose_averaged_LET([(10.0, 0.0), (20.0, 0.0)], segment_length_nm=10.0)


def test_mean_chord_length_sphere():
    assert mean_chord_length_sphere(3.0) == pytest.approx(2.0)


def test_lineal_energy():
    assert lineal_energy(10.0, diameter_nm=3.0) == pytest.approx(5.0)


def test_specific_energy_matches_manual_calculation(water):
    diameter_nm = 10.0
    epsilon_eV = 1000.0
    radius_nm = diameter_nm / 2.0
    volume_nm3 = (4.0 / 3.0) * np.pi * radius_nm ** 3
    volume_cm3 = volume_nm3 / constants.CM3_TO_NM3
    mass_kg = water.density_g_cm3 * volume_cm3 * 1.0e-3
    expected = (epsilon_eV * constants.EV_TO_JOULE) / mass_kg
    assert specific_energy(epsilon_eV, diameter_nm, water) == pytest.approx(expected)


def test_spherical_site_scorer_only_counts_inside():
    site = SphericalSiteScorer(center_nm=np.array([0.0, 0.0, 0.0]), diameter_nm=10.0)
    site.score(np.array([1.0, 0.0, 0.0]), 100.0)   # inside (r=1 < 5)
    site.score(np.array([100.0, 0.0, 0.0]), 999.0)  # outside
    assert site.epsilon_eV == pytest.approx(100.0)


def test_spherical_site_scorer_y_and_z(water):
    site = SphericalSiteScorer(center_nm=np.array([0.0, 0.0, 0.0]), diameter_nm=10.0)
    site.score(np.array([0.0, 0.0, 0.0]), 500.0)
    assert site.y() == pytest.approx(lineal_energy(500.0, 10.0))
    assert site.z(water) == pytest.approx(specific_energy(500.0, 10.0, water))


def test_spherical_site_scorer_reset():
    site = SphericalSiteScorer(center_nm=np.array([0.0, 0.0, 0.0]), diameter_nm=10.0)
    site.score(np.array([0.0, 0.0, 0.0]), 500.0)
    site.reset()
    assert site.epsilon_eV == 0.0


def test_collision_log_reflects_primary_energy_loss_exactly(full_xs, water, full_edep, full_eemean):
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=9)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    e0 = p.energy
    spent = engine.transport(p)

    total_from_log = sum(dE for _, dE in engine.collision_log)
    assert total_from_log == pytest.approx(e0 - spent.energy)


def test_track_averaged_LET_from_full_run_is_positive_and_finite(full_xs, water, full_edep, full_eemean):
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=9)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    engine.transport(p)

    LET_t = track_averaged_LET(engine.collision_log)
    LET_d = dose_averaged_LET(engine.collision_log, segment_length_nm=100.0)
    assert 0.0 < LET_t < np.inf
    assert LET_d >= LET_t - 1e-9


def test_spherical_site_scorer_integrates_with_engine(full_xs, water, full_edep, full_eemean):
    """A huge site centered on the origin should capture ~all deposited
    energy from a track that starts there (sanity, not a precision claim)."""
    site = SphericalSiteScorer(center_nm=np.array([0.0, 0.0, 0.0]), diameter_nm=2.0e9)
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, scorers=[site], seed=9)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    e0 = p.energy
    spent = engine.transport(p)
    assert site.epsilon_eV == pytest.approx(e0 - spent.energy)
