"""The event-by-event transport engine: follows one proton collision by
collision until it stops, and every delta-electron it produces.

One step of the loop (transport/README_transport.md §5):
  1. sample the distance s to the next collision (transport/sampler.py);
  2. move the proton s nm in a straight line (no energy is lost in flight:
     all loss happens at collisions);
  3. choose which process happens (ionization, excitation, capture, loss
     or elastic scattering), with probability sigma_c / sigma_tot;
  4. apply it: subtract the energy lost, record the local deposit, change
     the charge state (capture/loss), change direction (elastic), or eject
     a delta-electron (ionization);
  5. repeat until the proton's energy drops below the tabulated range or
     below the cut-off e_cut.

Energy lost per process comes from the mean-energy tables through
physics/energy_loss.py; the delta-electron direction from
physics/kinematics.py.

Delta-electrons are put on a stack and handled after the proton, by an
`ElectronTransport` plug-in. The default one deposits the electron's whole
energy where it was created (secondaries/README_secondaries.md): electrons
are not transported.

Every deposit is (a) stored in `TrackEngine.deposits` and (b) passed to each
attached `Scorer` (dose grids, microdosimetric sites, ...). The proton's own
(step length, energy lost) pairs go to `TrackEngine.collision_log`, which the
LET functions in scoring/radiobiology.py read; electron deposits are left out
of it, because LET describes the primary particle's own track.
"""
from dataclasses import dataclass

import numpy as np

from ion_track_mc.physics import energy_loss, kinematics
from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import CrossSectionSet
from ion_track_mc.physics.medium import Medium
from ion_track_mc.scoring.scorer import Scorer
from ion_track_mc.secondaries.electron_transport import ElectronTransport, LocalDepositElectronTransport

from . import sampler
from .particle import Particle


@dataclass
class Deposit:
    """One local energy deposit: `energy` [eV] left at `position` [nm].

    `channel` is the process that caused it, or None for a deposit made by
    a delta-electron. `kind` ('proton' or 'electron') and `q` (the proton's
    charge state at that moment; None for electrons) say which particle
    made it; the track plots use them to colour the points.
    """
    position: np.ndarray
    energy: float
    channel: Channel | None
    kind: str = "proton"
    q: int | None = None


class TrackEngine:
    """Event-by-event transport engine.

    Everything physical is handed in at construction: the medium, the cross
    sections, the mean-energy tables, the electron model and the scorers.
    The engine itself only runs the loop, so swapping any of them needs no
    change here. A primary and every secondary it spawns are finished
    before `transport()` returns.

    edep_tables[(q, channel)] and eemean_tables[q] are the Edepmeantot and
    Eemeantot tables, (E_eV, value_eV) pairs from load_mean_energy_table.
    e_cut [eV]: below this energy a particle is considered stopped
    (thermalised). seed fixes the random-number stream, so a run is exactly
    reproducible.
    """

    def __init__(
        self,
        medium: Medium,
        xs: CrossSectionSet,
        edep_tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]],
        eemean_tables: dict[int, tuple[np.ndarray, np.ndarray]],
        electron_transport: ElectronTransport | None = None,
        scorers: list[Scorer] | None = None,
        e_cut: float = 10.0,
        seed: int = 0,
    ):
        self.medium = medium
        self.xs = xs
        self.edep_tables = edep_tables
        self.eemean_tables = eemean_tables
        self.electron_transport = electron_transport or LocalDepositElectronTransport(e_cut=e_cut)
        self.scorers = list(scorers) if scorers else []   # notified of every deposit
        self.e_cut = e_cut            # [eV] local-deposit cutoff (thermalization)
        self.rng = np.random.default_rng(seed)
        self.deposits: list[Deposit] = []
        self.stack: list[Particle] = []   # particles waiting to be transported (LIFO)
        self.collision_log: list[tuple[float, float]] = []   # (step_length_nm, delta_E_eV), primary only

    def _notify_scorers(self, position: np.ndarray, energy: float) -> None:
        for scorer in self.scorers:
            scorer.score(position, energy)

    def transport(self, primary: Particle) -> Particle:
        """Transport `primary` and all secondaries it spawns (depth-first).

        Returns the (now spent) primary; by the time this returns, `self.stack`
        is empty and every secondary's energy has been fully accounted for in
        `self.deposits`.
        """
        self.collision_log = []   # one call = one track
        self.stack.append(primary)
        # Last in, first out: each delta-electron is handled right after
        # the proton has finished.
        while self.stack:
            p = self.stack.pop()
            if p.kind == "proton":
                self._transport_proton(p)
            elif p.kind == "electron":
                # The electron model returns the deposits it makes; record
                # and score them like any other.
                for d in self.electron_transport.transport(p):
                    self.deposits.append(Deposit(position=d.position, energy=d.energy, channel=None,
                                                  kind=p.kind, q=None))
                    self._notify_scorers(d.position, d.energy)
            else:
                raise ValueError(f"unknown particle kind: {p.kind}")
        return primary

    def _transport_proton(self, p: Particle) -> None:
        while p.alive and p.energy > self.e_cut:
            e_min, e_max = self.xs.energy_range(p.q)
            if not (e_min <= p.energy <= e_max):
                # Outside the tabulated energy range the cross sections are
                # unknown: stop the track rather than extrapolate. Any
                # energy left is not deposited.
                break
            # 1. How far to the next collision (random, exponential).
            s = sampler.sample_step(p.energy, p.q, self.medium.N, self.xs, self.rng)
            # 2. Fly there in a straight line, losing no energy on the way.
            p.position = p.position + s * p.direction
            # 3. Which process happens (random, weighted by cross section).
            channel = sampler.select_channel(p.energy, p.q, self.xs, self.rng)
            # 4. Apply it; the energy difference is this collision's loss.
            e_before = p.energy
            self._interact(p, channel)
            self.collision_log.append((s, e_before - p.energy))

    def _interact(self, p: Particle, channel: Channel) -> None:
        if channel == Channel.EXCITATION:
            # Molecule left in an excited state: energy deposited on the spot.
            delta_E = energy_loss.mean_energy_lookup(p.energy, self.edep_tables[(p.q, channel)])
            self._deposit_local(p, delta_E, channel)
        elif channel == Channel.IONIZATION:
            # An electron is knocked out: binding energy stays here, the
            # rest (W) leaves with the delta-electron.
            edep, W = energy_loss.ionization_energy(
                p.energy, self.edep_tables[(p.q, channel)], self.eemean_tables[p.q]
            )
            direction = kinematics.sample_ejected_electron_direction(
                p.energy, edep + W, W, p.direction, self.rng
            )
            self._deposit_local(p, edep, channel)
            self._eject_electron(p, W, direction)
        elif channel == Channel.CAPTURE:
            # H+ picks up an electron from the molecule and becomes H0.
            delta_E = energy_loss.capture_delta_E(p.energy, self.edep_tables[(p.q, channel)])
            self._deposit_local(p, delta_E, channel)
            p.q = 0
        elif channel == Channel.LOSS:
            # H0 is stripped of its electron and becomes H+ again.
            delta_E = energy_loss.mean_energy_lookup(p.energy, self.edep_tables[(p.q, channel)])
            self._deposit_local(p, delta_E, channel)
            p.q = 1
        elif channel == Channel.ELASTIC:
            # Scattering off the nuclei: new direction, small recoil energy.
            E_elastic, new_direction = energy_loss.elastic_scatter(
                p.energy, p.direction, self.medium.target_nuclear_charge,
                self.medium.target_nucleon_count, self.rng,
            )
            p.direction = new_direction
            self._deposit_local(p, E_elastic, channel)
        else:
            raise ValueError(f"unhandled channel: {channel}")

    def _deposit_local(self, p: Particle, delta_E: float, channel: Channel) -> None:
        delta_E = min(delta_E, p.energy)  # never deposit more than the particle carries
        # Energy conservation: what the projectile loses is what is deposited.
        p.energy -= delta_E
        self.deposits.append(Deposit(position=p.position.copy(), energy=delta_E, channel=channel,
                                      kind=p.kind, q=p.q))
        self._notify_scorers(p.position, delta_E)
        if p.energy <= 0.0:
            p.alive = False

    def _eject_electron(self, p: Particle, W: float, direction: np.ndarray) -> None:
        """Create a delta-electron of kinetic energy W and put it on the
        stack. W is taken off the projectile, so in total an ionization
        removes (local deposit + W) from it."""
        W = min(W, p.energy)
        p.energy -= W
        if p.energy <= 0.0:
            p.alive = False
        self.stack.append(Particle(kind="electron", position=p.position.copy(),
                                    direction=direction, energy=W, parent_id=id(p)))
