"""Correctness anchor: emergent S(E) must reproduce
the analytic expectation of the same tables/rules.

Both S_MC and S_analytic now use real, energy-dependent
Edepmeantot/Eemeantot-driven formulas (physics/energy_loss.py) for every
channel except elastic, which stays at 0.0 by policy (matches the reference
TILDA-V model excluding elastic from its electronic-SP sum -- see
ion_track_mc/physics/README_physics.md). This is still a self-consistency check (MC sampling
machinery vs. the exact expectation of the same rules), not yet a real
PSTAR/ICRU comparison -- see analysis/stopping_power.py module docstring.
"""
import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.analysis.stopping_power import (
    analytic_equilibrium_stopping_power,
    analytic_stopping_power,
    emergent_equilibrium_stopping_power,
    emergent_stopping_power,
    equilibrium_charge_fractions,
    mean_delta_E,
    sample_equilibrium_charge_fractions,
    sample_equilibrium_stopping_power,
)
from ion_track_mc.physics.channels import Channel
from ion_track_mc.transport.engine import TrackEngine

N_SAMPLES = 20000
SEED = 2024

# Probe energies in the well-tabulated mid-range (avoid grid edges where
# interpolation kinks add noise): 100, 500, 2000 keV/u.
PROBE_ENERGIES_KEV_U = [100.0, 500.0, 2000.0]


@pytest.mark.parametrize("E_kev_u", PROBE_ENERGIES_KEV_U)
def test_emergent_stopping_power_matches_analytic(full_xs, full_edep, full_eemean, water, E_kev_u):
    E = E_kev_u * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=SEED)

    S_mc = emergent_stopping_power(engine, E, q=1, n_samples=N_SAMPLES)
    S_analytic = analytic_stopping_power(water, full_xs, E, q=1, edep_tables=full_edep, eemean_tables=full_eemean)

    rel_error = abs(S_mc - S_analytic) / S_analytic
    assert rel_error < 0.08, (
        f"S_MC={S_mc:.6g} eV/nm vs S_analytic={S_analytic:.6g} eV/nm "
        f"at E={E_kev_u} keV/u: rel. error {rel_error:.2%}"
    )


def test_analytic_stopping_power_is_positive_across_grid(full_xs, full_edep, full_eemean, water):
    e_min, e_max = full_xs.energy_range(1)
    for E in np.geomspace(e_min, e_max, 15):
        assert analytic_stopping_power(water, full_xs, E, q=1, edep_tables=full_edep, eemean_tables=full_eemean) > 0.0


def test_elastic_contributes_zero_to_analytic_S_by_design(full_edep, full_eemean):
    """Unlike capture/loss (now real formulas, see below), elastic's zero
    contribution is not a placeholder to revisit -- it matches the reference
    model (Quinto et al. 2017 Eq. 16) excluding elastic from the
    electronic-SP sum on purpose."""
    E = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    assert mean_delta_E(Channel.ELASTIC, E, 1, full_edep, full_eemean) == 0.0


def test_capture_and_loss_contribute_real_nonzero_energy_now(full_edep, full_eemean):
    """Capture and loss deposit real energy, from the Edepmeantot tables
    (ion_track_mc/physics/README_physics.md)."""
    E = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    assert mean_delta_E(Channel.CAPTURE, E, 1, full_edep, full_eemean) > 0.0
    assert mean_delta_E(Channel.LOSS, E, 0, full_edep, full_eemean) > 0.0


# --- Charge-state-equilibrium blending (Eqs. III.3/III.4/III.8, Alcocer-Avila
# 2021 PhD thesis).
# `analytic_stopping_power(..., q=1, ...)` alone assumes the beam stays H+ the
# whole track; a real beam spends part of its path as H0 via capture/loss
# cycling, which matters most at low energy where sigma_capture is non-negligible.

@pytest.mark.parametrize("E_kev_u", PROBE_ENERGIES_KEV_U)
def test_equilibrium_charge_fractions_sum_to_one(full_xs, E_kev_u):
    E = E_kev_u * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    f_Hplus, f_H0 = equilibrium_charge_fractions(full_xs, E)
    assert f_Hplus + f_H0 == pytest.approx(1.0)
    assert 0.0 < f_Hplus < 1.0
    assert 0.0 < f_H0 < 1.0


def test_equilibrium_charge_fractions_match_cross_section_ratio(full_xs):
    # Eqs. III.3/III.4: f_H+ = sigma_L / (sigma_L + sigma_C), f_H0 = sigma_C / (sigma_L + sigma_C).
    E = 100.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    sigma_C = full_xs.sigma(Channel.CAPTURE, E, 1)
    sigma_L = full_xs.sigma(Channel.LOSS, E, 0)
    f_Hplus, f_H0 = equilibrium_charge_fractions(full_xs, E)
    assert f_Hplus == pytest.approx(sigma_L / (sigma_L + sigma_C))
    assert f_H0 == pytest.approx(sigma_C / (sigma_L + sigma_C))


def test_equilibrium_charge_fraction_favors_Hplus_at_high_energy(full_xs):
    # Capture cross section falls off fast with energy (harder to capture an
    # electron at high projectile velocity) -- f_H+ should climb toward 1.
    E_low = 100.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    E_high = 2000.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    f_Hplus_low, _ = equilibrium_charge_fractions(full_xs, E_low)
    f_Hplus_high, _ = equilibrium_charge_fractions(full_xs, E_high)
    assert f_Hplus_high > f_Hplus_low


@pytest.mark.parametrize("E_kev_u", PROBE_ENERGIES_KEV_U)
def test_analytic_equilibrium_stopping_power_is_weighted_average(full_xs, full_edep, full_eemean, water, E_kev_u):
    # Eq. III.8: SP = f_H+ * SP_H+ + f_H0 * SP_H0 must land between the two
    # pure-charge-state values (or equal one of them if a fraction is ~0).
    E = E_kev_u * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    S1 = analytic_stopping_power(water, full_xs, E, q=1, edep_tables=full_edep, eemean_tables=full_eemean)
    S0 = analytic_stopping_power(water, full_xs, E, q=0, edep_tables=full_edep, eemean_tables=full_eemean)
    S_eff = analytic_equilibrium_stopping_power(water, full_xs, E, edep_tables=full_edep, eemean_tables=full_eemean)
    assert min(S1, S0) <= S_eff <= max(S1, S0)


def test_analytic_equilibrium_stopping_power_differs_from_q1_only(full_xs, full_edep, full_eemean, water):
    # At 100 keV/u, f_H0 ~ 15% (see equilibrium-fraction test above) -- large
    # enough that blending must move S away from the pure-q=1 value.
    E = 100.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    S1 = analytic_stopping_power(water, full_xs, E, q=1, edep_tables=full_edep, eemean_tables=full_eemean)
    S_eff = analytic_equilibrium_stopping_power(water, full_xs, E, edep_tables=full_edep, eemean_tables=full_eemean)
    assert S_eff != pytest.approx(S1, rel=1e-6)


# --- Emergent (Monte Carlo) equilibrium stopping power: two different ways
# to make Eq. III.8 "emergent" instead of purely closed-form.
#
# (A) emergent_equilibrium_stopping_power blends two independently-sampled,
#     fixed-q `emergent_stopping_power` campaigns using the SAME closed-form
#     f_H+/f_H0 weights `analytic_equilibrium_stopping_power` uses -- cheap,
#     but doesn't test anything the existing fixed-q MC-vs-analytic tests
#     don't already cover.
#
# (B) sample_equilibrium_charge_fractions / sample_equilibrium_stopping_power
#     run ONE stationary-mode collision chain (TILDA-V's own term, Alcocer-
#     Avila 2021 thesis Ch. III) at fixed E, letting q flip via real
#     CAPTURE/LOSS events exactly as TrackEngine._interact applies them --
#     the charge mix genuinely emerges from the sampled dynamics instead of
#     being blended in by a closed-form weight. This validates the actual
#     sampling machinery (sampler.select_channel + _interact's q-flipping)
#     against the closed-form Eqs. III.3/III.4 prediction, not just the
#     blending arithmetic.
#
# Equilibration is only practical to observe within a reasonable sample
# budget at low energy -- P(capture | q=1) is ~6% at 100 keV/u but
# ~0.09%/~0.002% at 500/2000 keV/u (captures become too rare to see enough
# flips), so the chain-based tests (B) only run at 100 keV/u.

CHAIN_E_KEV_U = 100.0
N_CHAIN_COLLISIONS = 30000


@pytest.mark.parametrize("E_kev_u", PROBE_ENERGIES_KEV_U)
def test_emergent_equilibrium_stopping_power_matches_analytic(full_xs, full_edep, full_eemean, water, E_kev_u):
    E = E_kev_u * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=SEED)

    S_eff_mc = emergent_equilibrium_stopping_power(engine, E, n_samples=N_SAMPLES)
    S_eff_analytic = analytic_equilibrium_stopping_power(water, full_xs, E, edep_tables=full_edep, eemean_tables=full_eemean)

    rel_error = abs(S_eff_mc - S_eff_analytic) / S_eff_analytic
    assert rel_error < 0.08, (
        f"S_eff_MC={S_eff_mc:.6g} eV/nm vs S_eff_analytic={S_eff_analytic:.6g} eV/nm "
        f"at E={E_kev_u} keV/u: rel. error {rel_error:.2%}"
    )


def test_sample_equilibrium_charge_fractions_matches_closed_form(full_xs, full_edep, full_eemean, water):
    E = CHAIN_E_KEV_U * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=SEED)

    f_Hplus_mc, f_H0_mc = sample_equilibrium_charge_fractions(engine, E, q0=1, n_collisions=N_CHAIN_COLLISIONS)
    f_Hplus_analytic, f_H0_analytic = equilibrium_charge_fractions(full_xs, E)

    assert f_Hplus_mc == pytest.approx(f_Hplus_analytic, abs=0.03)
    assert f_H0_mc == pytest.approx(f_H0_analytic, abs=0.03)


def test_sample_equilibrium_charge_fractions_independent_of_starting_state(full_xs, full_edep, full_eemean, water):
    # A real equilibrium shouldn't depend on which charge state the chain
    # started in, given enough collisions to forget the initial condition.
    E = CHAIN_E_KEV_U * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    engine_from_1 = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=SEED)
    engine_from_0 = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=SEED + 1)

    f_Hplus_from_1, _ = sample_equilibrium_charge_fractions(engine_from_1, E, q0=1, n_collisions=N_CHAIN_COLLISIONS)
    f_Hplus_from_0, _ = sample_equilibrium_charge_fractions(engine_from_0, E, q0=0, n_collisions=N_CHAIN_COLLISIONS)

    assert f_Hplus_from_1 == pytest.approx(f_Hplus_from_0, abs=0.03)


def test_sample_equilibrium_stopping_power_matches_analytic(full_xs, full_edep, full_eemean, water):
    E = CHAIN_E_KEV_U * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=SEED)

    S_eff_mc = sample_equilibrium_stopping_power(engine, E, q0=1, n_collisions=N_CHAIN_COLLISIONS)
    S_eff_analytic = analytic_equilibrium_stopping_power(water, full_xs, E, edep_tables=full_edep, eemean_tables=full_eemean)

    rel_error = abs(S_eff_mc - S_eff_analytic) / S_eff_analytic
    assert rel_error < 0.15, (
        f"S_eff_MC={S_eff_mc:.6g} eV/nm vs S_eff_analytic={S_eff_analytic:.6g} eV/nm "
        f"at E={CHAIN_E_KEV_U} keV/u: rel. error {rel_error:.2%}"
    )
