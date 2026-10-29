"""Plotting: the standard figures of a track-structure study.

Three pieces (analysis/README_analysis.md §3):
  - style presets (`use_paper_style`, `use_presentation_style`) that set
    Matplotlib rcParams once, for whatever gets plotted next;
  - plot functions for what the code produces: cross sections, mean free
    path, Monte Carlo vs analytic stopping power (a self-consistency check,
    not a comparison with PSTAR -- see analysis/stopping_power.py), Bragg
    curve, radial dose, and a single track drawn collision by collision;
  - `save_figure`, which writes a figure as PNG and PDF.

There is no f(y) plot: f(y) distributions are not computed by this code.

Every function takes an optional `ax`; if omitted, a new figure+axes is
created. All return the `Axes` they drew on, for further customization or
`fig.savefig(...)` by the caller. **Exception:** `plot_stopping_power_self_consistency`
takes/returns a pair, `axes`/`(ax_log, ax_lin)`, not a single `ax` -- see its
own docstring.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ion_track_mc import constants
from ion_track_mc.analysis.stopping_power import analytic_stopping_power, emergent_stopping_power
from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import CrossSectionSet
from ion_track_mc.physics.medium import Medium
from ion_track_mc.scoring.dose import RadialDoseScorer, VoxelDoseScorer
from ion_track_mc.transport import sampler
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle


def use_paper_style() -> None:
    """Serif fonts, small figure, tight margins -- for print/manuscript figures."""
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.figsize": (3.5, 2.8),
        "figure.dpi": 300,
        "lines.linewidth": 1.2,
        "axes.grid": False,
        "savefig.bbox": "tight",
    })


def use_presentation_style() -> None:
    """Sans-serif, large fonts and lines, bigger figure -- for slides/screens."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 16,
        "axes.labelsize": 18,
        "axes.titlesize": 18,
        "legend.fontsize": 14,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "figure.figsize": (8, 6),
        "figure.dpi": 150,
        "lines.linewidth": 2.5,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "savefig.bbox": "tight",
    })


def save_figure(fig, out_dir, name: str, dpi: int = 300) -> list[Path]:
    """Save `fig` as `out_dir/name.png` (raster, `dpi`) and `out_dir/name.pdf`
    (vector), creating `out_dir` if needed. Returns the two paths.

    Used by the notebooks to write every figure into their `output_NN/`
    directory, so figures are available as files, not only inside the
    notebook.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [out_dir / f"{name}.png", out_dir / f"{name}.pdf"]
    for path in paths:
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return paths


def _new_ax(ax):
    if ax is None:
        _, ax = plt.subplots()
    return ax


def plot_cross_sections(xs: CrossSectionSet, q: int, ax=None, n_points: int = 200,
                         line_width: float | None = None):
    """sigma_c(E) [nm^2] per channel open for charge state q, log-log."""
    ax = _new_ax(ax)
    e_min, e_max = xs.energy_range(q)
    E_grid = np.geomspace(e_min, e_max, n_points)
    line_kwargs = {} if line_width is None else {"linewidth": line_width}
    for channel in xs.channels(q):
        sigma = [xs.sigma(channel, E, q) for E in E_grid]
        ax.plot(E_grid, sigma, label=channel.value, **line_kwargs)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Energy [eV]")
    ax.set_ylabel(r"$\sigma(E)$ [nm$^2$]")
    ax.set_title(f"Cross sections, q={q}")
    ax.legend()
    return ax


def plot_mean_free_path(xs: CrossSectionSet, medium: Medium, q: int, ax=None, n_points: int = 200,
                         line_width: float | None = None):
    """lambda(E) = 1 / (N * sigma_tot(E, q))  [nm], log-log."""
    ax = _new_ax(ax)
    e_min, e_max = xs.energy_range(q)
    E_grid = np.geomspace(e_min, e_max, n_points)
    lam = [1.0 / (medium.N * xs.sigma_total(E, q)) for E in E_grid]
    line_kwargs = {} if line_width is None else {"linewidth": line_width}
    ax.plot(E_grid, lam, **line_kwargs)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Energy [eV]")
    ax.set_ylabel(r"$\lambda(E)$ [nm]")
    ax.set_title(f"Mean free path, q={q}")
    return ax


def plot_stopping_power_self_consistency(
    medium: Medium, xs: CrossSectionSet,
    edep_tables: dict[tuple[int, Channel], tuple[np.ndarray, np.ndarray]],
    eemean_tables: dict[int, tuple[np.ndarray, np.ndarray]],
    q: int, axes=None,
    n_energy_points: int = 12, n_samples_per_point: int = 8000, seed: int = 0,
    line_width: float | None = None, marker_size: float | None = None,
):
    """S(E) Monte Carlo vs analytic: two panels, same x
    (energy, always log-scale given the multi-decade tabulated range) and
    same two series, differing only in y-scale -- log (left) to compare
    curves spanning orders of magnitude, linear (right) to see absolute
    shape/gaps that a log axis compresses away. Deviates from this module's
    usual single-`Axes`-return convention (see module docstring) because
    the log/linear comparison is the point of this particular plot.

    NOT a PSTAR/ICRU comparison -- both curves are computed from the same
    real (not synthetic) Edepmeantot/Eemeantot tables and cross sections;
    this checks the Monte Carlo sampling machinery against the closed-form
    analytic expectation of those same tables/rules, not against reference
    PSTAR/ICRU data. See analysis/stopping_power.py's module docstring for
    the full self-consistency framing, including the deliberate small
    elastic-channel divergence between the two curves.

    Returns `(ax_log, ax_lin)`. Callers overlaying more series (PSTAR,
    experimental data, ...) should plot onto both, e.g. `for ax in
    (ax_log, ax_lin): ax.plot(...)`, to keep both panels complete.
    """
    if axes is None:
        _fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    ax_log, ax_lin = axes
    e_min, e_max = xs.energy_range(q)
    E_grid = np.geomspace(e_min, e_max, n_energy_points)
    engine = TrackEngine(medium, xs, edep_tables=edep_tables, eemean_tables=eemean_tables, seed=seed)
    S_mc = [emergent_stopping_power(engine, E, q, n_samples_per_point) for E in E_grid]
    S_analytic = [analytic_stopping_power(medium, xs, E, q, edep_tables, eemean_tables) for E in E_grid]

    line_kwargs = {} if line_width is None else {"linewidth": line_width}
    marker_kwargs = {} if marker_size is None else {"markersize": marker_size}
    for ax in (ax_log, ax_lin):
        ax.plot(E_grid, S_analytic, "-", label="Analytic", **line_kwargs)
        ax.plot(E_grid, S_mc, "o", label="Monte Carlo (emergent)", **marker_kwargs)
        ax.set_xscale("log")
        ax.set_xlabel("Energy [eV]")
        ax.set_ylabel("S(E) [eV/nm]")
        ax.legend()
    ax_log.set_yscale("log")
    ax_lin.set_yscale("linear")
    ax_log.set_title("Emergent stopping power: MC vs analytic self-consistency\n"
                      "(not a PSTAR comparison) -- log scale")
    ax_lin.set_title("Emergent stopping power: MC vs analytic self-consistency\n"
                      "(not a PSTAR comparison) -- linear scale")
    return ax_log, ax_lin


def plot_bragg_curve(scorer: VoxelDoseScorer, ax=None, line_width: float | None = None):
    """Dose [Gy] vs depth for a scorer built with `make_bragg_curve_scorer`."""
    ax = _new_ax(ax)
    z_centers = 0.5 * (scorer.z_edges[:-1] + scorer.z_edges[1:])
    dose = scorer.dose_Gy()[0, 0, :]
    line_kwargs = {} if line_width is None else {"linewidth": line_width}
    ax.plot(z_centers, dose, **line_kwargs)
    ax.set_xlabel("Depth [nm]")
    ax.set_ylabel("Dose [Gy]")
    ax.set_title("Bragg curve")
    return ax


def plot_radial_dose(scorer: RadialDoseScorer, ax=None,
                      line_width: float | None = None, marker_size: float | None = None):
    """Dose [Gy] vs radial distance from the track axis, log-log."""
    ax = _new_ax(ax)
    r_centers = 0.5 * (scorer.r_edges[:-1] + scorer.r_edges[1:])
    dose = scorer.dose_Gy()
    kwargs = {}
    if line_width is not None:
        kwargs["linewidth"] = line_width
    if marker_size is not None:
        kwargs["markersize"] = marker_size
    ax.plot(r_centers, dose, "o-", **kwargs)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Radial distance [nm]")
    ax.set_ylabel("Dose [Gy]")
    ax.set_title("Radial dose")
    return ax


_CHARGE_STATE_COLORS = {1: "tab:blue", 0: "tab:orange"}
_CHARGE_STATE_LABELS = {1: "H1+", 0: "H0"}
_CHANNEL_COLORS = {
    Channel.EXCITATION: "tab:green",
    Channel.IONIZATION: "tab:red",
    Channel.CAPTURE: "tab:purple",
    Channel.LOSS: "tab:brown",
    Channel.ELASTIC: "tab:cyan",
}
_ELECTRON_COLOR = "black"


def plot_track(entry_particle: Particle, engine: TrackEngine, ax=None,
                marker_size: float = 12, electron_marker_size: float = 25,
                line_width: float | None = None, aspect: float | str | None = None,
                channel_marker_sizes: dict[Channel, float] | None = None,
                legend_fontsize: float | None = None):
    """One primary's full simulated track, projected onto the plane of
    depth z (along the beam) vs transverse position x.

    The proton path is drawn as at most two lines -- one per distinct
    charge state (H1+/H0) that actually occurs -- each a single Line2D
    built from NaN-separated runs of consecutive same-charge-state
    segments (so non-contiguous runs of the same charge state are not
    spuriously connected, and drawing stays O(1) artists regardless of
    collision count). A marker sits at every collision, face-colored by
    reaction channel. A separate marker (no line) sits at every secondary
    electron's deposit point -- electron transport is fidelity level (a)
    (secondaries/electron_transport.py): the electron deposits its full
    energy at its creation point instantly, so it has no trajectory to draw.

    `entry_particle` must be the Particle object as it was BEFORE being
    passed to `engine.transport()` -- transport() mutates position in
    place, so the caller must keep a copy for the track's starting point
    (vertex 0 of the polyline). Reads `engine.deposits`; supports one
    track's worth of deposits at a time -- there is no per-track id, so a
    multi-primary engine's accumulated deposits (e.g. from `run_beam`)
    would interleave multiple tracks with no way to separate them here.

    `marker_size`/`electron_marker_size` are scatter marker **areas** in
    points^2 (matplotlib `scatter`'s `s`) -- passed to `ax.scatter(s=...)`
    for, respectively, the collision markers and the electron deposit
    markers -- NOT the linear `markersize` convention used by this module's
    other line-plot functions, so a value that looks right there will look
    very different here. `line_width`, unlike the two marker sizes above,
    defaults to `None`: when `None`, the proton-path line uses whatever
    `linewidth` the active rcParams/style preset (`use_paper_style()`/
    `use_presentation_style()`) specifies, the same convention this
    module's other functions follow. `aspect` is passed straight through to
    `Axes.set_aspect` when not `None` -- accepts `'equal'`, `'auto'`, or a
    numeric ratio (y-data-units per x-data-unit); `None` (default) means
    `set_aspect` is never called, preserving matplotlib's normal
    auto-scaled axes.

    `channel_marker_sizes` overrides `marker_size` per reaction channel --
    e.g. `{Channel.IONIZATION: 40}` draws ionization markers at size 40
    while every other channel keeps using `marker_size`. Channels absent
    from the dict, and every channel when the dict itself is `None`
    (default), fall back to `marker_size`. Does not affect
    `electron_marker_size`. `legend_fontsize` sets the legend text size
    (passed to `Axes.legend(fontsize=...)`); `None` (default) leaves it at
    the rcParams/style-preset default, same convention as `line_width`.
    """
    ax = _new_ax(ax)
    proton_deposits = [d for d in engine.deposits if d.kind == "proton"]
    electron_deposits = [d for d in engine.deposits if d.kind == "electron"]

    vertices = [entry_particle.position] + [d.position for d in proton_deposits]

    # Group consecutive proton_deposits sharing the same q into runs, then
    # draw ONE Line2D per distinct q (not one per segment -- with tens of
    # thousands of collisions at high energy, one artist per segment makes
    # both drawing and `ax.legend(loc="best")` extremely slow). Each run of
    # deposit indices [a, ..., b] needs vertex indices range(a, b + 2): the
    # entry vertex just before the run plus every vertex inside it, so runs
    # connect end-to-end with no visible gap at a capture/loss event. Runs
    # of the same q that are NOT adjacent are joined with a (nan, nan) pair
    # so matplotlib breaks the line between them instead of drawing a
    # spurious connector across an intervening run of the other charge state.
    q_runs: dict[int, list[range]] = {}
    if proton_deposits:
        run_start = 0
        current_q = proton_deposits[0].q
        for i in range(1, len(proton_deposits)):
            if proton_deposits[i].q != current_q:
                q_runs.setdefault(current_q, []).append(range(run_start, i + 1))
                run_start = i
                current_q = proton_deposits[i].q
        q_runs.setdefault(current_q, []).append(range(run_start, len(proton_deposits) + 1))

    line_kwargs = {} if line_width is None else {"linewidth": line_width}
    for q, runs in q_runs.items():
        z_vals: list[float] = []
        x_vals: list[float] = []
        for run_index, vertex_indices in enumerate(runs):
            if run_index > 0:
                z_vals.append(np.nan)
                x_vals.append(np.nan)
            for vi in vertex_indices:
                z_vals.append(vertices[vi][2])
                x_vals.append(vertices[vi][0])
        ax.plot(z_vals, x_vals, "-", color=_CHARGE_STATE_COLORS[q], label=_CHARGE_STATE_LABELS[q],
                **line_kwargs)

    for channel in (Channel.EXCITATION, Channel.IONIZATION, Channel.CAPTURE, Channel.LOSS, Channel.ELASTIC):
        pts = [d.position for d in proton_deposits if d.channel == channel]
        if pts:
            size = marker_size if channel_marker_sizes is None else channel_marker_sizes.get(channel, marker_size)
            ax.scatter([p[2] for p in pts], [p[0] for p in pts],
                       color=_CHANNEL_COLORS[channel], label=channel.value, zorder=3, s=size)

    if electron_deposits:
        ax.scatter([d.position[2] for d in electron_deposits], [d.position[0] for d in electron_deposits],
                   color=_ELECTRON_COLOR, marker="*", s=electron_marker_size, label="e- (deposit)", zorder=4)

    ax.set_xlabel("Depth z [nm]")
    ax.set_ylabel("Transverse x [nm]")
    ax.set_title("Simulated track")
    if aspect is not None:
        ax.set_aspect(aspect)
    if ax.get_legend_handles_labels()[0]:
        legend_kwargs = {} if legend_fontsize is None else {"fontsize": legend_fontsize}
        ax.legend(**legend_kwargs)
    return ax
