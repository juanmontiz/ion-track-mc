"""Secondary-electron transport (fidelity level (a) only -- see
secondaries/electron_transport.py module docstring for why (b)/(c) aren't
possible yet).
"""
import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.physics.channels import Channel
from ion_track_mc.secondaries.electron_transport import LocalDepositElectronTransport
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

E_START_EV = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT


def test_local_deposit_transport_deposits_full_energy():
    transport = LocalDepositElectronTransport()
    e = Particle(kind="electron", position=np.array([1.0, 2.0, 3.0]),
                 direction=np.array([0.0, 0.0, 1.0]), energy=42.0)
    deposits = transport.transport(e)
    assert len(deposits) == 1
    assert deposits[0].energy == pytest.approx(42.0)
    assert np.array_equal(deposits[0].position, np.array([1.0, 2.0, 3.0]))
    assert e.energy == 0.0
    assert e.alive is False


def test_engine_uses_local_deposit_by_default(full_xs, water, full_edep, full_eemean):
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=0)
    assert isinstance(engine.electron_transport, LocalDepositElectronTransport)


def test_full_transport_drains_stack_and_deposits_electron_energy(full_xs, water, full_edep, full_eemean):
    """End-to-end: an ionization event's delta-electron gets pushed, then
    drained by transport() itself, ending up as its own Deposit entry."""
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=11)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    e0 = p.energy
    spent = engine.transport(p)

    assert engine.stack == []  # fully drained
    electron_deposits = [d for d in engine.deposits if d.channel is None]
    ionization_deposits = [d for d in engine.deposits if d.channel == Channel.IONIZATION]
    assert len(electron_deposits) == len(ionization_deposits)  # one electron per ionization event
    assert len(electron_deposits) > 0  # this seed/energy is known to ionize at least once

    for d in electron_deposits:
        assert d.kind == "electron"
        assert d.q is None
    proton_deposits = [d for d in engine.deposits if d.channel is not None]
    for d in proton_deposits:
        assert d.kind == "proton"
        assert d.q in (0, 1)

    total_deposited = sum(d.energy for d in engine.deposits)
    assert e0 - spent.energy == pytest.approx(total_deposited)


def test_injected_electron_transport_is_used(full_xs, water, full_edep, full_eemean):
    """A custom ElectronTransport is honored instead of the default."""
    class CountingTransport:
        def __init__(self):
            self.calls = 0

        def transport(self, electron):
            self.calls += 1
            electron.energy = 0.0
            electron.alive = False
            from ion_track_mc.secondaries.electron_transport import ElectronDeposit
            return [ElectronDeposit(position=electron.position.copy(), energy=0.0)]

    counting = CountingTransport()
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, electron_transport=counting, seed=11)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    engine.transport(p)
    assert counting.calls > 0
