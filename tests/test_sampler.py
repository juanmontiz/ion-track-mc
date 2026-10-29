"""Free-path sampling reproduces
the analytic mean-free-path distribution, and channel selection reproduces
the sigma_c/sigma_tot ratios. Must pass before any further validation
(stopping power, dose, ...) is trusted.
"""
import numpy as np
from scipy import stats

from ion_track_mc import constants
from ion_track_mc.physics.channels import Channel
from ion_track_mc.transport import sampler

N_SAMPLES = 20000
E_PROBE_EV = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT  # 500 keV/u, q=1
SEED = 12345


def test_sample_step_ks_test_against_analytic_exponential(full_xs, water):
    rng = np.random.default_rng(SEED)
    Sigma_tot = water.N * full_xs.sigma_total(E_PROBE_EV, q=1)
    mean_free_path = 1.0 / Sigma_tot

    steps = np.array([
        sampler.sample_step(E_PROBE_EV, 1, water.N, full_xs, rng)
        for _ in range(N_SAMPLES)
    ])

    # H0: steps ~ Exponential(scale=mean_free_path).
    ks_stat, p_value = stats.kstest(steps, "expon", args=(0.0, mean_free_path))
    assert p_value > 0.01, f"KS test rejects exponential(lambda={mean_free_path}): p={p_value}"


def test_sample_step_empirical_mean_matches_analytic(full_xs, water):
    rng = np.random.default_rng(SEED)
    Sigma_tot = water.N * full_xs.sigma_total(E_PROBE_EV, q=1)
    mean_free_path = 1.0 / Sigma_tot

    steps = np.array([
        sampler.sample_step(E_PROBE_EV, 1, water.N, full_xs, rng)
        for _ in range(N_SAMPLES)
    ])

    # Standard error of the mean of an Exponential(lambda) sample of size N
    # is lambda/sqrt(N); require agreement within 5 sigma (should virtually
    # never flake at N=20000, ~0.9% expected precision).
    sem = mean_free_path / np.sqrt(N_SAMPLES)
    assert abs(steps.mean() - mean_free_path) < 5 * sem


def test_select_channel_frequencies_match_sigma_ratios(full_xs):
    rng = np.random.default_rng(SEED)
    chans = full_xs.channels(1)
    sig = np.array([full_xs.sigma(c, E_PROBE_EV, 1) for c in chans])
    p_expected = sig / sig.sum()

    counts = {c: 0 for c in chans}
    for _ in range(N_SAMPLES):
        counts[sampler.select_channel(E_PROBE_EV, 1, full_xs, rng)] += 1

    p_empirical = np.array([counts[c] for c in chans]) / N_SAMPLES
    # Multinomial SEM per category ~ sqrt(p(1-p)/N); allow 5 sigma slack.
    sem = np.sqrt(p_expected * (1 - p_expected) / N_SAMPLES)
    assert np.all(np.abs(p_empirical - p_expected) < 5 * sem + 1e-6)


def test_select_channel_only_returns_open_channels(full_xs):
    rng = np.random.default_rng(SEED)
    for _ in range(200):
        c = sampler.select_channel(E_PROBE_EV, 0, full_xs, rng)
        assert c in (Channel.IONIZATION, Channel.EXCITATION, Channel.LOSS, Channel.ELASTIC)
