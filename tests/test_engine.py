import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.physics import energy_loss
from ion_track_mc.physics.channels import Channel
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

E_START_EV = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT  # 500 keV/u


def _make_proton(energy=E_START_EV):
    return Particle(kind="proton", position=np.zeros(3),
                     direction=np.array([0.0, 0.0, 1.0]), energy=energy, q=1)


def _make_engine(water, full_xs, full_edep, full_eemean, **kwargs):
    return TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, **kwargs)


def test_interact_excitation_deposits_locally(full_xs, full_edep, full_eemean, water):
    engine = _make_engine(water, full_xs, full_edep, full_eemean, seed=0)
    p = _make_proton()
    e0 = p.energy
    engine._interact(p, Channel.EXCITATION)
    expected = energy_loss.mean_energy_lookup(e0, full_edep[(1, Channel.EXCITATION)])
    assert p.energy == pytest.approx(e0 - expected)
    assert len(engine.deposits) == 1
    assert engine.deposits[0].channel == Channel.EXCITATION
    assert engine.deposits[0].energy == pytest.approx(expected)
    assert engine.deposits[0].kind == "proton"
    assert engine.deposits[0].q == 1


def test_interact_ionization_deposits_binding_and_pushes_electron(full_xs, full_edep, full_eemean, water):
    engine = _make_engine(water, full_xs, full_edep, full_eemean, seed=0)
    p = _make_proton()
    e0 = p.energy
    engine._interact(p, Channel.IONIZATION)

    expected_edep, expected_W = energy_loss.ionization_energy(
        e0, full_edep[(1, Channel.IONIZATION)], full_eemean[1]
    )
    assert len(engine.deposits) == 1
    assert engine.deposits[0].channel == Channel.IONIZATION
    assert engine.deposits[0].energy == pytest.approx(expected_edep)
    assert engine.deposits[0].kind == "proton"
    assert engine.deposits[0].q == 1

    assert len(engine.stack) == 1
    electron = engine.stack[0]
    assert electron.kind == "electron"
    assert electron.energy == pytest.approx(expected_W)

    total_lost = e0 - p.energy
    assert total_lost == pytest.approx(expected_edep + expected_W)


def test_interact_ionization_W_is_deterministic_eemeantot(full_xs, full_edep, full_eemean, water):
    """W is now Eemeantot(E) exactly, not sampled -- repeated calls at the
    same E give the same W (unlike the old Rutherford model)."""
    engine = _make_engine(water, full_xs, full_edep, full_eemean, seed=3)
    p = _make_proton()
    _, expected_W = energy_loss.ionization_energy(
        p.energy, full_edep[(1, Channel.IONIZATION)], full_eemean[1]
    )
    for _ in range(20):
        p2 = _make_proton()
        engine._interact(p2, Channel.IONIZATION)
        W = engine.stack.pop().energy
        assert W == pytest.approx(expected_W)
        assert W > 0.0


def test_interact_capture_flips_q_and_deposits(full_xs, full_edep, full_eemean, water):
    engine = _make_engine(water, full_xs, full_edep, full_eemean, seed=0)
    p = _make_proton()
    p.q = 1
    e0 = p.energy
    engine._interact(p, Channel.CAPTURE)
    expected = energy_loss.capture_delta_E(e0, full_edep[(1, Channel.CAPTURE)])
    assert p.q == 0
    assert p.energy == pytest.approx(e0 - expected)
    assert len(engine.deposits) == 1
    assert engine.deposits[0].channel == Channel.CAPTURE
    assert engine.deposits[0].energy == pytest.approx(expected)
    assert engine.deposits[0].kind == "proton"
    assert engine.deposits[0].q == 1


def test_interact_loss_flips_q_and_deposits(full_xs, full_edep, full_eemean, water):
    engine = _make_engine(water, full_xs, full_edep, full_eemean, seed=0)
    p = _make_proton()
    p.q = 0
    e0 = p.energy
    engine._interact(p, Channel.LOSS)
    expected = energy_loss.mean_energy_lookup(e0, full_edep[(0, Channel.LOSS)])
    assert p.q == 1
    assert p.energy == pytest.approx(e0 - expected)
    assert expected == pytest.approx(constants.H_BINDING_ENERGY_EV)
    assert len(engine.deposits) == 1
    assert engine.deposits[0].channel == Channel.LOSS
    assert engine.deposits[0].kind == "proton"
    assert engine.deposits[0].q == 0


def test_interact_elastic_deposits_and_direction_stays_unit(full_xs, full_edep, full_eemean, water):
    """Elastic deposits real energy and deflects the
    trajectory (previously a zero-deposit, direction-unchanged no-op)."""
    engine = _make_engine(water, full_xs, full_edep, full_eemean, seed=0)
    p = _make_proton()
    e0 = p.energy
    engine._interact(p, Channel.ELASTIC)
    assert p.energy < e0
    assert len(engine.deposits) == 1
    assert engine.deposits[0].channel == Channel.ELASTIC
    assert engine.deposits[0].energy == pytest.approx(e0 - p.energy)
    assert np.linalg.norm(p.direction) == pytest.approx(1.0)
    assert engine.deposits[0].kind == "proton"
    assert engine.deposits[0].q == 1


def test_deposit_local_never_exceeds_available_energy(full_xs, full_edep, full_eemean, water):
    """_deposit_local's own clamp (min(delta_E, p.energy)), tested directly
    against the method rather than via _interact. With real tables, every
    tabulated E_min (~8 keV) is now far above any real Edepmeantot value
    (tens to low-hundreds of eV), so this scenario is no longer reachable
    through _interact/a real table lookup -- interp_log_log would raise
    ValueError first for any p.energy low enough to need clamping. This
    remains a defensive property of _deposit_local itself, tested in
    isolation."""
    engine = _make_engine(water, full_xs, full_edep, full_eemean, seed=0)
    p = _make_proton(energy=5.0)
    engine._deposit_local(p, 75.0, Channel.EXCITATION)  # delta_E > p.energy, on purpose
    assert p.energy == 0.0
    assert engine.deposits[0].energy == pytest.approx(5.0)
    assert p.alive is False


def test_transport_energy_accounting_is_exact(full_xs, full_edep, full_eemean, water):
    engine = _make_engine(water, full_xs, full_edep, full_eemean, e_cut=10.0, seed=42)
    p = _make_proton()
    e0 = p.energy
    spent = engine.transport(p)
    total_deposited = sum(d.energy for d in engine.deposits)
    assert e0 - spent.energy == pytest.approx(total_deposited)


def test_transport_energy_never_negative(full_xs, full_edep, full_eemean, water):
    engine = _make_engine(water, full_xs, full_edep, full_eemean, e_cut=10.0, seed=7)
    p = _make_proton()
    spent = engine.transport(p)
    assert spent.energy >= 0.0


def test_transport_terminates_at_or_below_e_cut_or_out_of_range(full_xs, full_edep, full_eemean, water):
    engine = _make_engine(water, full_xs, full_edep, full_eemean, e_cut=10.0, seed=1)
    p = _make_proton()
    spent = engine.transport(p)
    e_min, _ = full_xs.energy_range(spent.q)
    assert (spent.energy <= engine.e_cut) or (spent.energy < e_min) or not spent.alive


def test_transport_stays_within_tabulated_range_when_starting_out_of_range(full_xs, full_edep, full_eemean, water):
    e_max = full_xs.energy_range(1)[1]
    engine = _make_engine(water, full_xs, full_edep, full_eemean, e_cut=10.0, seed=0)
    p = _make_proton(energy=e_max * 10.0)  # starts above tabulated range
    spent = engine.transport(p)
    assert spent.energy == e_max * 10.0  # never moved: loop breaks immediately
