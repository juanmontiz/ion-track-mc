import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.analysis.beam import run_beam
from ion_track_mc.scoring.dose import make_bragg_curve_scorer

# Kept close to the tabulated floor (~10 keV/u) so each track only takes a
# few dozen collisions -- full-energy (100s of keV/u+) tracks run into the
# tens of thousands of collisions (pure-Python loop), too slow for a test
# suite at the particle counts needed for decent statistics.
E_BEAM = 15.0e3 * constants.PROTON_NUCLEON_COUNT  # 15 keV/u


def test_monoenergetic_beam_all_start_at_E_beam(full_xs, water, full_edep, full_eemean):
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=20, E_beam_eV=E_BEAM, E_sigma_eV=0.0, seed=1)
    assert np.all(result.E_start_eV == E_BEAM)
    assert result.n_rejected_out_of_range == 0
    assert result.E_start_eV.shape == (20,)
    assert result.E_final_eV.shape == (20,)
    assert result.q_final.shape == (20,)


def test_gaussian_beam_energies_within_tabulated_range(full_xs, water, full_edep, full_eemean):
    e_min, e_max = full_xs.energy_range(1)
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=50, E_beam_eV=E_BEAM,
                       E_sigma_eV=0.05 * E_BEAM, seed=2)
    assert np.all(result.E_start_eV >= e_min)
    assert np.all(result.E_start_eV <= e_max)


def test_gaussian_beam_energy_statistics_roughly_match(full_xs, water, full_edep, full_eemean):
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=800, E_beam_eV=E_BEAM,
                       E_sigma_eV=0.02 * E_BEAM, seed=3)
    assert result.E_start_eV.mean() == pytest.approx(E_BEAM, rel=0.02)
    assert result.E_start_eV.std() == pytest.approx(0.02 * E_BEAM, rel=0.25)


def test_final_energies_lower_than_start(full_xs, water, full_edep, full_eemean):
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=10, E_beam_eV=E_BEAM, seed=4)
    assert np.all(result.E_final_eV <= result.E_start_eV)


def test_beam_energy_conservation_aggregate(full_xs, water, full_edep, full_eemean):
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=15, E_beam_eV=E_BEAM, seed=5)
    total_lost = float((result.E_start_eV - result.E_final_eV).sum())
    total_deposited = sum(d.energy for d in result.engine.deposits)
    assert total_deposited == pytest.approx(total_lost)


def test_beam_feeds_shared_scorer(full_xs, water, full_edep, full_eemean):
    scorer = make_bragg_curve_scorer(water, z_edges=np.linspace(-1e9, 1e9, 11),
                                      transverse_half_width_nm=1e9)
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=10, E_beam_eV=E_BEAM,
                       scorers=[scorer], seed=6)
    total_lost = float((result.E_start_eV - result.E_final_eV).sum())
    assert scorer.total_deposited_eV() == pytest.approx(total_lost)


def test_impossible_gaussian_raises(full_xs, water, full_edep, full_eemean):
    e_min, _ = full_xs.energy_range(1)
    with pytest.raises(ValueError):
        run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                  n_particles=5, E_beam_eV=e_min * 0.5,
                  E_sigma_eV=1.0, max_resample_attempts=3, seed=7)


def test_beam_reproducible_with_same_seed(full_xs, water, full_edep, full_eemean):
    r1 = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                   n_particles=10, E_beam_eV=E_BEAM,
                   E_sigma_eV=0.05 * E_BEAM, seed=42)
    r2 = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                   n_particles=10, E_beam_eV=E_BEAM,
                   E_sigma_eV=0.05 * E_BEAM, seed=42)
    assert np.array_equal(r1.E_start_eV, r2.E_start_eV)
    assert np.array_equal(r1.E_final_eV, r2.E_final_eV)


def test_beam_averaged_LET_total_delta_E_matches_energy_loss(full_xs, water, full_edep, full_eemean):
    """total_delta_E_eV is summed from every particle's own collision_log,
    which (per test_collision_log_reflects_primary_energy_loss_exactly in
    test_radiobiology.py) equals that particle's own E_start - E_final
    exactly -- so summed over the beam it must match the aggregate loss."""
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=15, E_beam_eV=E_BEAM, seed=8)
    total_lost = float((result.E_start_eV - result.E_final_eV).sum())
    assert result.total_delta_E_eV == pytest.approx(total_lost)


def test_beam_averaged_LET_methods(full_xs, water, full_edep, full_eemean):
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=15, E_beam_eV=E_BEAM, seed=9, let_segment_length_nm=20.0)
    LET_t = result.track_averaged_LET()
    LET_d = result.dose_averaged_LET()
    assert LET_t == pytest.approx(result.total_delta_E_eV / result.total_path_length_nm)
    assert LET_d >= LET_t - 1e-9  # Cauchy-Schwarz, same property as the single-track case


def test_beam_averaged_LET_differs_from_last_particle_only(full_xs, water, full_edep, full_eemean):
    """Guard against silently regressing to 'just use engine.collision_log'
    (which only ever holds the LAST transported particle) -- the beam
    aggregate must generally differ from a single track's LET."""
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=15, E_beam_eV=E_BEAM, seed=10)
    from ion_track_mc.scoring.radiobiology import track_averaged_LET
    last_particle_LET = track_averaged_LET(result.engine.collision_log)
    assert result.total_path_length_nm > sum(s for s, _ in result.engine.collision_log)
    # Not asserting inequality of the LET values themselves (could coincide by
    # chance); the path-length check above is the real regression guard.


def test_beam_dose_averaged_LET_requires_segment_length(full_xs, water, full_edep, full_eemean):
    """L_D depends on the segment length, so a beam run without one can't
    report it -- it must say so rather than return a number."""
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=3, E_beam_eV=E_BEAM, seed=11)
    result.track_averaged_LET()  # still available
    with pytest.raises(ValueError, match="let_segment_length_nm"):
        result.dose_averaged_LET()


def test_beam_dose_averaged_LET_segments_each_track_separately(full_xs, water, full_edep, full_eemean):
    """Beam L_D = (sum over tracks of each track's own segment sum) / total dE.
    For a 1-particle beam that must equal the single-track function on the
    same collision_log exactly -- guards the running-sum bookkeeping."""
    from ion_track_mc.scoring.radiobiology import dose_averaged_LET
    result = run_beam(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean,
                       n_particles=1, E_beam_eV=E_BEAM, seed=12, let_segment_length_nm=20.0)
    assert result.dose_averaged_LET() == pytest.approx(
        dose_averaged_LET(result.engine.collision_log, segment_length_nm=20.0))
