"""Run many protons as a pencil beam with a Gaussian energy spread.

`TrackEngine.transport()` follows ONE primary (and its delta-electrons) per
call. A beam is simply many such calls on the same engine: the scorers and
`engine.deposits` keep adding up across calls, so after the loop they hold
the whole beam's dose. `run_beam` is that loop plus energy sampling and LET
bookkeeping (analysis/README_analysis.md §1).

"Pencil beam" means every proton starts at the same point with the same
direction. A broad beam can be built afterwards by superposition (see
notebook 02).
"""
from dataclasses import dataclass

import numpy as np

from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import CrossSectionSet
from ion_track_mc.physics.medium import Medium
from ion_track_mc.scoring.radiobiology import (
    LET_segment_sums,
    dose_averaged_LET_from_sums,
    track_averaged_LET_from_sums,
)
from ion_track_mc.scoring.scorer import Scorer
from ion_track_mc.secondaries.electron_transport import ElectronTransport
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle


@dataclass
class BeamResult:
    """Summary of a `run_beam` call.

    `engine` holds the accumulated deposits/scorers from every history in
    the beam; the arrays below are per-particle (length `n_particles`).

    LET is accumulated as running sums over every particle's
    `TrackEngine.collision_log`, added up right after each `transport()`
    call (the next call overwrites the log). Memory use therefore does not
    grow with the number of collisions. Use `track_averaged_LET()` /
    `dose_averaged_LET()` for the beam-averaged values.

    The dose-averaged sum depends on the segment length LET is measured over
    (see `scoring.radiobiology.dose_averaged_LET`), so it is only filled in
    when `run_beam` was given `let_segment_length_nm`; otherwise it stays
    None and `dose_averaged_LET()` raises.
    """
    engine: TrackEngine
    E_start_eV: np.ndarray
    E_final_eV: np.ndarray
    q_final: np.ndarray
    n_rejected_out_of_range: int
    total_path_length_nm: float = 0.0
    total_delta_E_eV: float = 0.0
    let_segment_length_nm: float | None = None
    total_segment_E2_over_length_eV2_per_nm: float | None = None

    def track_averaged_LET(self) -> float:
        """Beam-averaged LET_T [eV/nm] over every particle transported."""
        return track_averaged_LET_from_sums(self.total_path_length_nm, self.total_delta_E_eV)

    def dose_averaged_LET(self) -> float:
        """Beam-averaged LET_D [eV/nm] over every particle transported, each
        track cut into its own segments of `let_segment_length_nm`."""
        if self.total_segment_E2_over_length_eV2_per_nm is None:
            raise ValueError(
                "dose-averaged LET depends on the segment length it is measured "
                "over: pass let_segment_length_nm to run_beam to get it"
            )
        return dose_averaged_LET_from_sums(self.total_segment_E2_over_length_eV2_per_nm,
                                           self.total_delta_E_eV)


def run_beam(
    medium: Medium,
    xs: CrossSectionSet,
    edep_tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]],
    eemean_tables: dict[int, tuple[np.ndarray, np.ndarray]],
    n_particles: int,
    E_beam_eV: float,
    E_sigma_eV: float = 0.0,
    q: int = 1,
    direction: np.ndarray | None = None,
    position: np.ndarray | None = None,
    scorers: list[Scorer] | None = None,
    electron_transport: ElectronTransport | None = None,
    e_cut: float = 10.0,
    seed: int = 0,
    max_resample_attempts: int = 100,
    let_segment_length_nm: float | None = None,
) -> BeamResult:
    """Transport `n_particles` primaries, energies ~ Normal(E_beam_eV,
    E_sigma_eV) (or exactly `E_beam_eV` if `E_sigma_eV <= 0`, i.e.
    monoenergetic), through `medium` using cross sections `xs`.

    All particles feed the SAME `TrackEngine`, so `scorers` and
    `engine.deposits` accumulate over the whole beam.

    Energies sampled outside the tabulated `xs.energy_range(q)` are
    rejected and resampled (up to `max_resample_attempts` tries each) rather
    than silently clipped or allowed to immediately terminate the track
    (the cross sections are unknown outside that range). Raises
    `ValueError` if a valid energy can't be found within the attempt budget
    (e.g. `E_sigma_eV` too large relative to the tabulated range).

    `let_segment_length_nm` enables the beam-averaged dose-averaged LET
    (each track segmented separately, see `BeamResult`); leave it None if
    only track-averaged LET is needed.

    Energy sampling uses its own RNG stream (spawned from `seed`, independent
    of the engine's internal stream) so beam-energy statistics don't
    correlate with the per-collision sampling inside `TrackEngine`.
    """
    direction = np.array([0.0, 0.0, 1.0]) if direction is None else np.asarray(direction, dtype=float)
    position = np.zeros(3) if position is None else np.asarray(position, dtype=float)

    # Two independent random streams from one seed: one for beam energies,
    # one for the collisions inside the engine.
    beam_seed, engine_seed = np.random.SeedSequence(seed).spawn(2)
    beam_rng = np.random.default_rng(beam_seed)
    engine = TrackEngine(medium, xs, edep_tables=edep_tables, eemean_tables=eemean_tables,
                          electron_transport=electron_transport,
                          scorers=scorers, e_cut=e_cut, seed=engine_seed)
    e_min, e_max = xs.energy_range(q)

    E_start = np.empty(n_particles)
    E_final = np.empty(n_particles)
    q_final = np.empty(n_particles, dtype=int)
    n_rejected = 0
    total_s = 0.0
    total_dE = 0.0
    total_E2_over_length = None if let_segment_length_nm is None else 0.0

    for i in range(n_particles):
        if E_sigma_eV <= 0.0:
            E0 = E_beam_eV          # monoenergetic beam
        else:
            # Draw from the Gaussian until the energy is inside the tables.
            for _ in range(max_resample_attempts):
                E0 = float(beam_rng.normal(E_beam_eV, E_sigma_eV))
                if e_min <= E0 <= e_max:
                    break
                n_rejected += 1
            else:
                raise ValueError(
                    f"could not sample an in-range energy after {max_resample_attempts} "
                    f"attempts (E_beam={E_beam_eV}, E_sigma={E_sigma_eV}, "
                    f"valid range=[{e_min}, {e_max}])"
                )
        p = Particle(kind="proton", position=position.copy(), direction=direction.copy(),
                     energy=E0, q=q)
        spent = engine.transport(p)   # one full history: proton + its electrons
        E_start[i] = E0
        E_final[i] = spent.energy
        q_final[i] = spent.q

        # Add this track to the beam's LET sums before the next track
        # overwrites engine.collision_log. Each track is segmented on its own.
        for s, dE in engine.collision_log:
            total_s += s
            total_dE += dE
        if let_segment_length_nm is not None and engine.collision_log:
            track_E2_over_length, _ = LET_segment_sums(engine.collision_log, let_segment_length_nm)
            total_E2_over_length += track_E2_over_length

    return BeamResult(engine=engine, E_start_eV=E_start, E_final_eV=E_final,
                       q_final=q_final, n_rejected_out_of_range=n_rejected,
                       total_path_length_nm=total_s, total_delta_E_eV=total_dE,
                       let_segment_length_nm=let_segment_length_nm,
                       total_segment_E2_over_length_eV2_per_nm=total_E2_over_length)
