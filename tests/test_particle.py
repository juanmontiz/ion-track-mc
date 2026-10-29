import numpy as np

from ion_track_mc.transport.particle import Particle


def test_particle_q_is_mutable():
    p = Particle(kind="proton", position=np.zeros(3), direction=np.array([0.0, 0.0, 1.0]),
                 energy=1000.0, q=1)
    assert p.q == 1
    p.q = 0
    assert p.q == 0


def test_particle_defaults():
    p = Particle(kind="proton", position=np.zeros(3), direction=np.array([0.0, 0.0, 1.0]),
                 energy=1000.0)
    assert p.q == 1
    assert p.alive is True
    assert p.parent_id == -1
