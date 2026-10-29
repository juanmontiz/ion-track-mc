"""The Scorer interface: anything that wants to see each energy deposit.

`TrackEngine` calls `score(position, energy)` on every attached scorer, once
per deposit, while the track is being simulated (the "observer" pattern).
Several independent scorers (a depth-dose grid, a radial-dose profile, a
microdosimetric site, ...) can therefore watch the same run, and none of
them needs to know how the engine stores its own data.
See scoring/README_scoring.md §1.
"""
from abc import ABC, abstractmethod

import numpy as np


class Scorer(ABC):
    """Observer notified once per local energy-deposit event."""

    @abstractmethod
    def score(self, position: np.ndarray, energy: float) -> None:
        """Record one deposit of `energy` [eV] at `position` [nm]."""
        ...
