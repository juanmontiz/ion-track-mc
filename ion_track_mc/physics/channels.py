"""Enumeration of the physical interaction channels.

Every collision of the projectile with a water molecule is one of these five
processes. A small closed enum means the transport engine, the cross-section
tables and the energy-loss formulas all name a process the same way; adding a
new process starts here. See physics/README_physics.md §1.
"""
from enum import Enum


class Channel(str, Enum):
    IONIZATION = "ionization"   # ejects a delta-electron; projectile keeps q
    EXCITATION = "excitation"   # discrete local deposit; no secondary
    CAPTURE    = "capture"      # q: 1 -> 0  (H+ captures an electron -> H0)
    LOSS       = "loss"         # q: 0 -> 1  (H0 loses its electron -> H+)
    ELASTIC    = "elastic"      # scattering off the nuclei: deflects the
                                 # projectile and transfers a little energy
                                 # (physics/energy_loss.py::elastic_scatter).
                                 # It is counted in the free path and in the
                                 # deposits, but left out of the analytic
                                 # "electronic stopping power" by convention
                                 # (see analysis/stopping_power.py).
