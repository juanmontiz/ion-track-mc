"""Scorer observer pattern. Required test: total
deposited energy conserved vs energy lost by primaries+secondaries.
"""
import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.scoring.scorer import Scorer
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

E_START_EV = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT


class RecordingScorer(Scorer):
    def __init__(self):
        self.events: list[tuple[np.ndarray, float]] = []

    def score(self, position: np.ndarray, energy: float) -> None:
        self.events.append((position.copy(), energy))

    def total(self) -> float:
        return sum(e for _, e in self.events)


def test_scorer_is_abstract():
    with pytest.raises(TypeError):
        Scorer()


def test_engine_with_no_scorers_still_works(full_xs, water, full_edep, full_eemean):
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=0)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    engine.transport(p)  # must not raise


def test_scorer_notified_matches_total_deposited(full_xs, water, full_edep, full_eemean):
    scorer = RecordingScorer()
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, scorers=[scorer], seed=3)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    e0 = p.energy
    spent = engine.transport(p)

    assert scorer.total() == pytest.approx(sum(d.energy for d in engine.deposits))
    assert scorer.total() == pytest.approx(e0 - spent.energy)


def test_multiple_scorers_all_notified(full_xs, water, full_edep, full_eemean):
    s1, s2 = RecordingScorer(), RecordingScorer()
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, scorers=[s1, s2], seed=3)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    engine.transport(p)
    assert len(s1.events) == len(s2.events)
    assert s1.total() == pytest.approx(s2.total())
