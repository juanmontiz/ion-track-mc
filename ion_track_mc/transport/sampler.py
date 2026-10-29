"""The two random draws made at every collision.

1. How far does the projectile fly before its next collision?
   (`sample_step`, transport/README_transport.md §3)
2. Which process happens there? (`select_channel`, §4)

They are plain functions, separate from the TrackEngine loop, so each can be
tested on its own against its analytic distribution.
"""
import numpy as np

from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import CrossSectionSet


def sample_step(E: float, q: int, N: float, xs: CrossSectionSet, rng: np.random.Generator) -> float:
    """Distance [nm] to the next collision.

    The probability of flying a distance s without colliding is
    exp(-Sigma_tot * s), with Sigma_tot = N * sigma_tot(E, q) the
    macroscopic total cross section [nm^-1]. The distance to the first
    collision therefore has the density p(s) = Sigma_tot * exp(-Sigma_tot * s). Inverting its cumulative distribution with a uniform random
    number xi in (0, 1) gives s = -ln(xi) / Sigma_tot. The mean step is the
    mean free path 1 / Sigma_tot.
    """
    # Total collision probability per nm, summed over every open channel.
    Sigma_tot = N * xs.sigma_total(E, q)
    # Inverse-transform sampling of the exponential free-path distribution.
    return -np.log(rng.random()) / Sigma_tot


def select_channel(E: float, q: int, xs: CrossSectionSet, rng: np.random.Generator) -> Channel:
    """Choose which process happens at this collision.

    Given that a collision occurs, channel c is chosen with probability
    P(c) = sigma_c(E, q) / sigma_tot(E, q). This is done by discrete
    inverse-transform sampling: stack the probabilities end to end on
    [0, 1] and see which interval a uniform random number lands in.
    """
    chans = xs.channels(q)          # channels open for this charge state
    sig = np.array([xs.sigma(c, E, q) for c in chans])
    cdf = np.cumsum(sig) / sig.sum()  # cumulative probabilities, last = 1
    idx = int(np.searchsorted(cdf, rng.random()))
    return chans[idx]
