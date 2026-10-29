"""Stopping power S(E): mean energy lost per unit path [eV/nm].

Stopping power is the most basic check of a track-structure code: it sets
the range of the ion and the shape of the Bragg curve. This module computes
it in two independent ways from the SAME tables
(analysis/README_analysis.md §2):

  - Monte Carlo ("emergent"): run many single collisions at a fixed energy
    E and divide the total energy lost by the total path flown,
    S_MC = sum(delta_E_i) / sum(s_i);
  - analytic: the exact expectation of those same rules,
    S_analytic = N * sum_c sigma_c(E) * <delta_E_c>.

Agreement between the two shows the sampling machinery is correct. It is a
self-consistency check, not a comparison with measurements; notebook 01
compares against NIST PSTAR and experimental data separately.

Charge-state equilibrium. A real proton beam keeps capturing and losing
electrons, so at each energy it spends a fraction f_H+ of its path as H+
and f_H0 as H0. Balancing capture against loss gives
f_H+ = sigma_L / (sigma_L + sigma_C), and the stopping power to compare with
experiment is the blend f_H+ S(E, q=1) + f_H0 S(E, q=0)
(`analytic_equilibrium_stopping_power`; Alcocer-Avila 2021, Eqs.
III.3/III.4/III.8). `sample_equilibrium_*` reproduce the same fractions from
a sampled chain of collisions in which the charge state really changes.

Elastic scattering: a deliberate, permanent difference. The Monte Carlo
includes the small energy elastic collisions transfer to the nuclei; the
analytic formula leaves it out, because by convention the "electronic"
stopping power counts only collisions with electrons (as in the reference
TILDA-V model, Quinto et al. 2017, Eq. 16). The gap is tiny, about
0.008%-0.08% of S over 8 keV-2 MeV, and is not a bug.
"""
import numpy as np

from ion_track_mc.physics import energy_loss
from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import CrossSectionSet
from ion_track_mc.physics.medium import Medium
from ion_track_mc.transport import sampler
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle


def sample_single_collision(engine: TrackEngine, E: float, q: int) -> tuple[float, Channel, float]:
    """One collision event at fixed (E, q): (step length [nm], channel, deposited delta_E [eV]).

    Reuses `TrackEngine._interact` so the energy-loss dispatch logic lives in
    exactly one place (transport/engine.py); this module never re-implements it.
    """
    p = Particle(kind="proton", position=np.zeros(3), direction=np.array([0.0, 0.0, 1.0]),
                 energy=E, q=q)
    s = sampler.sample_step(E, q, engine.medium.N, engine.xs, engine.rng)
    channel = sampler.select_channel(E, q, engine.xs, engine.rng)
    e_before = p.energy
    engine._interact(p, channel)
    delta_E = e_before - p.energy
    return s, channel, delta_E


def emergent_stopping_power(engine: TrackEngine, E: float, q: int, n_samples: int) -> float:
    """S_MC(E) [eV/nm] = sum(delta_E_i) / sum(s_i) over `n_samples` i.i.d.
    single-collision events at fixed (E, q).

    Each event starts from the same (E, q): energy loss is NOT carried from
    one collision to the next, so this measures S at exactly E."""
    total_s = 0.0
    total_dE = 0.0
    for _ in range(n_samples):
        s, _channel, delta_E = sample_single_collision(engine, E, q)
        total_s += s
        total_dE += delta_E
    return total_dE / total_s


def mean_delta_E(
    channel: Channel,
    E: float,
    q: int,
    edep_tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]],
    eemean_tables: dict[int, tuple[np.ndarray, np.ndarray]],
) -> float:
    """<delta_E_c> at energy E, matching TrackEngine._interact's rules
    exactly for every channel except ELASTIC, which is deliberately
    excluded here by policy -- see this module's docstring."""
    if channel == Channel.EXCITATION:
        return energy_loss.mean_energy_lookup(E, edep_tables[(q, channel)])
    if channel == Channel.IONIZATION:
        edep, W = energy_loss.ionization_energy(E, edep_tables[(q, channel)], eemean_tables[q])
        return edep + W
    if channel == Channel.CAPTURE:
        return energy_loss.capture_delta_E(E, edep_tables[(q, channel)])
    if channel == Channel.LOSS:
        return energy_loss.mean_energy_lookup(E, edep_tables[(q, channel)])
    if channel == Channel.ELASTIC:
        return 0.0  # by convention, not missing data: the electronic
        # stopping power excludes collisions with nuclei, as in the
        # reference model (Quinto et al. 2017, Eq. 16, which sums only the
        # electronic processes).
    raise ValueError(f"unhandled channel: {channel}")


def analytic_stopping_power(
    medium: Medium,
    xs: CrossSectionSet,
    E: float,
    q: int,
    edep_tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]],
    eemean_tables: dict[int, tuple[np.ndarray, np.ndarray]],
) -> float:
    """S_analytic(E) [eV/nm] = N * sum_c sigma_c(E, q) * <delta_E_c>.

    N sigma_c is the number of collisions of type c per nm; each one costs
    <delta_E_c> on average, so the sum is the energy lost per nm."""
    return medium.N * sum(
        xs.sigma(c, E, q) * mean_delta_E(c, E, q, edep_tables, eemean_tables)
        for c in xs.channels(q)
    )


def equilibrium_charge_fractions(xs: CrossSectionSet, E: float) -> tuple[float, float]:
    """(f_H+, f_H0) at collision energy E, Eqs. III.3/III.4 (Alcocer-Avila
    2021 PhD thesis, Université de Bordeaux,
    citing [216]): f_H+ = sigma_L / (sigma_L + sigma_C), f_H0 = sigma_C /
    (sigma_L + sigma_C), where sigma_C/sigma_L are the total electron-
    capture (q=1) / electron-loss (q=0) cross sections at E. This is the
    same two-process balance TILDA-V uses to weight its SP(E) blend.

    Why: in equilibrium the rate of H+ -> H0 (capture, f_H+ sigma_C) equals
    the rate of H0 -> H+ (loss, f_H0 sigma_L); with f_H+ + f_H0 = 1 this
    gives the formulas above (analysis/README_analysis.md §2.2)."""
    sigma_C = xs.sigma(Channel.CAPTURE, E, 1)
    sigma_L = xs.sigma(Channel.LOSS, E, 0)
    total = sigma_C + sigma_L
    return sigma_L / total, sigma_C / total


def analytic_equilibrium_stopping_power(
    medium: Medium,
    xs: CrossSectionSet,
    E: float,
    edep_tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]],
    eemean_tables: dict[int, tuple[np.ndarray, np.ndarray]],
) -> float:
    """S(E) [eV/nm] = f_H+ * S(E, q=1) + f_H0 * S(E, q=0), Eq. III.8 -- the
    charge-state-equilibrium-blended stopping power TILDA-V actually reports
    (Fig. 6 in Quinto et al. 2017), as opposed to `analytic_stopping_power`
    with a fixed `q`, which implicitly assumes the beam never changes charge
    state. This is the physically appropriate quantity to compare against
    experimentally measured "proton" stopping-power data, since a real beam
    spends part of its path as H0 via capture/loss cycling -- most
    significant at low energy, where sigma_C is not yet negligible."""
    f_Hplus, f_H0 = equilibrium_charge_fractions(xs, E)
    S_Hplus = analytic_stopping_power(medium, xs, E, 1, edep_tables, eemean_tables)
    S_H0 = analytic_stopping_power(medium, xs, E, 0, edep_tables, eemean_tables)
    return f_Hplus * S_Hplus + f_H0 * S_H0


def emergent_equilibrium_stopping_power(engine: TrackEngine, E: float, n_samples: int) -> float:
    """S_eff(E) [eV/nm] via Eq. III.8's blend, but using the Monte Carlo
    `emergent_stopping_power` for each charge state instead of the
    closed-form `analytic_stopping_power` -- i.e. the same closed-form
    f_H+/f_H0 weights (Eqs. III.3/III.4) applied to two independently
    sampled, fixed-q MC campaigns. Cheap, but doesn't validate anything the
    per-charge-state MC-vs-analytic tests don't already cover, since the
    charge mix itself is still supplied by the closed form, not sampled.
    See `sample_equilibrium_stopping_power` for a genuinely single-campaign
    alternative where the charge mix emerges from the sampled dynamics."""
    f_Hplus, f_H0 = equilibrium_charge_fractions(engine.xs, E)
    S_Hplus = emergent_stopping_power(engine, E, q=1, n_samples=n_samples)
    S_H0 = emergent_stopping_power(engine, E, q=0, n_samples=n_samples)
    return f_Hplus * S_Hplus + f_H0 * S_H0


def sample_equilibrium_charge_chain(
    engine: TrackEngine, E: float, q0: int, n_collisions: int,
) -> list[tuple[int, float, float]]:
    """`n_collisions` sequential collisions in TILDA-V's own "stationary
    mode" (Alcocer-Avila 2021 PhD thesis, Chapter III "Charged particle
    transport simulations": each projectile is followed collision by
    collision at a FIXED incident energy E -- no energy loss carried
    between draws, only the charge state evolves), starting from charge
    state `q0`. Reuses `TrackEngine._interact` so the same real
    CAPTURE/LOSS q-flipping and energy-deposit rules apply as in an actual
    track (including elastic's real, non-mean-lookup recoil energy -- the
    same deliberate divergence from `mean_delta_E`/`analytic_stopping_power`
    documented at the top of this module for `emergent_stopping_power`).

    Returns a list of `(q_before, step_length, delta_E)` triples, one per
    collision -- shared raw material for both
    `sample_equilibrium_charge_fractions` and
    `sample_equilibrium_stopping_power` below, so the two are read off the
    same random walk rather than two independently-noisy runs."""
    q = q0
    history = []
    for _ in range(n_collisions):
        p = Particle(kind="proton", position=np.zeros(3), direction=np.array([0.0, 0.0, 1.0]),
                     energy=E, q=q)
        s = sampler.sample_step(E, q, engine.medium.N, engine.xs, engine.rng)
        channel = sampler.select_channel(E, q, engine.xs, engine.rng)
        e_before = p.energy
        engine._interact(p, channel)
        history.append((q, s, e_before - p.energy))
        q = p.q
    return history


def sample_equilibrium_charge_fractions(
    engine: TrackEngine, E: float, q0: int, n_collisions: int,
) -> tuple[float, float]:
    """Empirical (f_H+, f_H0) at E: path-length-weighted fraction of the
    stationary-mode collision chain (`sample_equilibrium_charge_chain`)
    spent in each charge state -- the Monte Carlo counterpart to
    `equilibrium_charge_fractions`' closed form, for validating the actual
    sampling machinery (`sampler.select_channel` +
    `TrackEngine._interact`'s q-flipping) against the analytic Eqs.
    III.3/III.4 prediction, rather than just re-deriving it.

    Only practical to converge within a reasonable `n_collisions` budget at
    energies where capture/loss aren't too rare (see `tests/test_stopping_power.py`
    for why the corresponding tests only probe 100 keV/u, not 500/2000)."""
    history = sample_equilibrium_charge_chain(engine, E, q0, n_collisions)
    path_by_q = {1: 0.0, 0: 0.0}
    for q, s, _delta_E in history:
        path_by_q[q] += s
    total = path_by_q[1] + path_by_q[0]
    return path_by_q[1] / total, path_by_q[0] / total


def sample_equilibrium_stopping_power(
    engine: TrackEngine, E: float, q0: int, n_collisions: int,
) -> float:
    """Fully emergent S_eff(E) [eV/nm]: sum(delta_E)/sum(s) over the same
    stationary-mode collision chain as `sample_equilibrium_charge_fractions`
    -- a single self-contained Monte Carlo campaign whose charge mix
    emerges from the real capture/loss dynamics, rather than being blended
    in via the closed-form weights the way `emergent_equilibrium_stopping_power`
    (or the purely analytic `analytic_equilibrium_stopping_power`) does."""
    history = sample_equilibrium_charge_chain(engine, E, q0, n_collisions)
    total_s = sum(s for _q, s, _delta_E in history)
    total_dE = sum(delta_E for _q, _s, delta_E in history)
    return total_dE / total_s
