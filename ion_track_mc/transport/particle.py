"""Particle: everything the transport loop needs to know about one particle.

See transport/README_transport.md §1.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class Particle:
    """Complete state of a particle under transport.

    The projectile is a hydrogen ion whose charge state `q` changes during
    the track: electron capture turns H+ (q=1) into H0 (q=0), electron loss
    turns it back. Each charge state has its own cross sections, so `q`
    decides which collisions are possible at every step. This back-and-forth
    is the central physics point of the code.

    Electrons (`kind='electron'`) are the delta-electrons ejected by
    ionization; for them `q` is unused.
    """
    kind: str                # 'proton' or 'electron'
    position: np.ndarray     # [nm], shape (3,)
    direction: np.ndarray    # unit vector, shape (3,)
    energy: float            # [eV]
    q: int = 1               # projectile charge state (1 = H+, 0 = H0)
    alive: bool = True       # False only if its energy reaches 0; a proton stopped
                             # at the lowest tabulated energy (8 keV) stays True
    parent_id: int = -1      # id() of the particle that created this one (-1 for a primary)
