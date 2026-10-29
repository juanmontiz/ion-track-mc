"""What happens to the delta-electrons that ionizations eject.

Electrons move very differently from protons: they are about 1836 times
lighter, so they scatter through large angles and follow their own cross
sections. Their treatment therefore lives in this separate plug-in, not in
the proton engine (secondaries/README_secondaries.md).

Three levels of fidelity are possible:
  (a) deposit the electron's whole energy where it was created;
  (b) spread it over an electron range;
  (c) follow the electron collision by collision, like the proton.

**This code implements (a) only.** Its consequence: delta-electron energy
never leaves the track core. Dose very close to the track axis is
overestimated, and dose far from it (carried there by energetic
delta-rays) is missing. Bragg curves and LET, which integrate over the
track, are much less affected.

`TrackEngine` only calls `ElectronTransport.transport(electron)`, so a (b)
or (c) model could replace (a) without changing the engine.
"""
from dataclasses import dataclass
from typing import Protocol

import numpy as np

from ion_track_mc.transport.particle import Particle


@dataclass
class ElectronDeposit:
    """One energy deposit made by an electron: `energy` [eV] at `position` [nm]."""
    position: np.ndarray
    energy: float


class ElectronTransport(Protocol):
    """Anything with this method can be used as the electron model: it takes
    one electron and returns the deposits it makes."""

    def transport(self, electron: Particle) -> list[ElectronDeposit]: ...


class LocalDepositElectronTransport:
    """Fidelity level (a): the electron deposits its full kinetic energy at
    its creation point; no actual transport is simulated.

    The `e_cut` parameter is accepted so that this class has the same
    constructor as a model that does transport electrons, but it is unused:
    here the entire energy is always deposited on the spot.
    """

    def __init__(self, e_cut: float = 10.0):
        self.e_cut = e_cut

    def transport(self, electron: Particle) -> list[ElectronDeposit]:
        deposit = ElectronDeposit(position=electron.position.copy(), energy=electron.energy)
        electron.energy = 0.0
        electron.alive = False
        return [deposit]
