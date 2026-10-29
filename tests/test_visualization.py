import matplotlib
matplotlib.use("Agg")  # headless, no display in CI/test environments

import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.analysis.visualization import (
    plot_bragg_curve,
    plot_cross_sections,
    plot_mean_free_path,
    plot_radial_dose,
    plot_stopping_power_self_consistency,
    plot_track,
    save_figure,
    use_paper_style,
    use_presentation_style,
)
from ion_track_mc.physics.channels import Channel
from ion_track_mc.scoring.dose import RadialDoseScorer, make_bragg_curve_scorer
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

E_START_EV = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT


def test_use_paper_style_sets_rcparams():
    use_paper_style()
    import matplotlib.pyplot as plt
    assert plt.rcParams["font.family"] == ["serif"]
    assert plt.rcParams["figure.figsize"] == [3.5, 2.8]


def test_use_presentation_style_sets_rcparams():
    use_presentation_style()
    import matplotlib.pyplot as plt
    assert plt.rcParams["font.family"] == ["sans-serif"]
    assert plt.rcParams["figure.figsize"] == [8.0, 6.0]


def test_plot_cross_sections_one_line_per_channel(full_xs):
    ax = plot_cross_sections(full_xs, q=1, n_points=20)
    assert len(ax.lines) == len(full_xs.channels(1))


def test_plot_cross_sections_line_width_override(full_xs):
    ax = plot_cross_sections(full_xs, q=1, n_points=20, line_width=3.0)
    for line in ax.lines:
        assert line.get_linewidth() == 3.0


def test_plot_cross_sections_line_width_default_unchanged(full_xs):
    ax_default = plot_cross_sections(full_xs, q=1, n_points=20)
    ax_explicit_none = plot_cross_sections(full_xs, q=1, n_points=20, line_width=None)
    assert ax_default.lines[0].get_linewidth() == ax_explicit_none.lines[0].get_linewidth()


def test_plot_mean_free_path_draws_one_line(full_xs, water):
    ax = plot_mean_free_path(full_xs, water, q=1, n_points=20)
    assert len(ax.lines) == 1
    assert np.all(np.array(ax.lines[0].get_ydata()) > 0.0)


def test_plot_mean_free_path_line_width_override(full_xs, water):
    ax = plot_mean_free_path(full_xs, water, q=1, n_points=20, line_width=3.0)
    assert ax.lines[0].get_linewidth() == 3.0


def test_plot_stopping_power_self_consistency_returns_two_panels(full_xs, water, full_edep, full_eemean):
    # Log-scale and linear-scale panels, same underlying two series each.
    ax_log, ax_lin = plot_stopping_power_self_consistency(
        water, full_xs, full_edep, full_eemean, q=1, n_energy_points=3, n_samples_per_point=200, seed=0
    )
    assert len(ax_log.lines) == 2
    assert len(ax_lin.lines) == 2


def test_plot_stopping_power_self_consistency_panel_yscales_differ(full_xs, water, full_edep, full_eemean):
    ax_log, ax_lin = plot_stopping_power_self_consistency(
        water, full_xs, full_edep, full_eemean, q=1, n_energy_points=3, n_samples_per_point=200, seed=0
    )
    assert ax_log.get_yscale() == "log"
    assert ax_lin.get_yscale() == "linear"


def test_plot_stopping_power_self_consistency_both_panels_share_figure(full_xs, water, full_edep, full_eemean):
    ax_log, ax_lin = plot_stopping_power_self_consistency(
        water, full_xs, full_edep, full_eemean, q=1, n_energy_points=3, n_samples_per_point=200, seed=0
    )
    assert ax_log.figure is ax_lin.figure


def test_plot_stopping_power_self_consistency_accepts_existing_axes(full_xs, water, full_edep, full_eemean):
    import matplotlib.pyplot as plt
    fig, (given_log, given_lin) = plt.subplots(1, 2)
    ax_log, ax_lin = plot_stopping_power_self_consistency(
        water, full_xs, full_edep, full_eemean, q=1, axes=(given_log, given_lin),
        n_energy_points=3, n_samples_per_point=200, seed=0,
    )
    assert ax_log is given_log
    assert ax_lin is given_lin
    assert len(given_log.lines) == 2
    assert len(given_lin.lines) == 2


def test_plot_stopping_power_self_consistency_style_override(full_xs, water, full_edep, full_eemean):
    ax_log, ax_lin = plot_stopping_power_self_consistency(
        water, full_xs, full_edep, full_eemean, q=1, n_energy_points=3, n_samples_per_point=200, seed=0,
        line_width=3.0, marker_size=15.0,
    )
    for ax in (ax_log, ax_lin):
        analytic_line = [l for l in ax.lines if l.get_label() == "Analytic"][0]
        mc_line = [l for l in ax.lines if l.get_label() == "Monte Carlo (emergent)"][0]
        assert analytic_line.get_linewidth() == 3.0
        assert mc_line.get_markersize() == 15.0


def test_plot_bragg_curve(full_xs, water, full_edep, full_eemean):
    z_edges = np.linspace(-1.0e9, 1.0e9, 11)
    scorer = make_bragg_curve_scorer(water, z_edges=z_edges, transverse_half_width_nm=1.0e9)
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, scorers=[scorer], seed=1)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    engine.transport(p)

    ax = plot_bragg_curve(scorer)
    assert len(ax.lines) == 1
    assert len(ax.lines[0].get_xdata()) == len(z_edges) - 1


def test_plot_bragg_curve_line_width_override(full_xs, water, full_edep, full_eemean):
    z_edges = np.linspace(-1.0e9, 1.0e9, 11)
    scorer = make_bragg_curve_scorer(water, z_edges=z_edges, transverse_half_width_nm=1.0e9)
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, scorers=[scorer], seed=1)
    p = Particle(kind="proton", position=np.zeros(3),
                 direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    engine.transport(p)

    ax = plot_bragg_curve(scorer, line_width=3.0)
    assert ax.lines[0].get_linewidth() == 3.0


def test_plot_radial_dose(water):
    scorer = RadialDoseScorer(water, r_edges=np.geomspace(1.0, 1000.0, 10), length_nm=100.0)
    scorer.score(np.array([5.0, 0.0, 50.0]), 100.0)
    ax = plot_radial_dose(scorer)
    assert len(ax.lines) == 1


def test_plot_radial_dose_style_override(water):
    scorer = RadialDoseScorer(water, r_edges=np.geomspace(1.0, 1000.0, 10), length_nm=100.0)
    scorer.score(np.array([5.0, 0.0, 50.0]), 100.0)
    ax = plot_radial_dose(scorer, line_width=3.0, marker_size=15.0)
    assert ax.lines[0].get_linewidth() == 3.0
    assert ax.lines[0].get_markersize() == 15.0


def _build_transported_track(full_xs, water, full_edep, full_eemean, seed=1):
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=seed)
    entry = Particle(kind="proton", position=np.zeros(3),
                      direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    entry_copy = Particle(kind=entry.kind, position=entry.position.copy(),
                           direction=entry.direction.copy(), energy=entry.energy, q=entry.q)
    engine.transport(entry)
    return entry_copy, engine


def test_plot_track_draws_charge_state_segments_and_reaction_markers(full_xs, water, full_edep, full_eemean):
    engine = TrackEngine(water, full_xs, edep_tables=full_edep, eemean_tables=full_eemean, seed=1)
    entry = Particle(kind="proton", position=np.zeros(3),
                      direction=np.array([0.0, 0.0, 1.0]), energy=E_START_EV, q=1)
    entry_copy = Particle(kind=entry.kind, position=entry.position.copy(),
                           direction=entry.direction.copy(), energy=entry.energy, q=entry.q)
    engine.transport(entry)  # mutates `entry` in place -- entry_copy keeps the pre-transport state

    ax = plot_track(entry_copy, engine)

    proton_deposits = [d for d in engine.deposits if d.kind == "proton"]
    electron_deposits = [d for d in engine.deposits if d.kind == "electron"]
    assert len(proton_deposits) > 0  # this seed/energy is known to produce collisions

    # Real physical invariant: every ionization
    # ejects exactly one electron. Assert this on the raw engine data,
    # independent of the plotting call, so the electron-marker-count check
    # below verifies the plot reflects a real invariant rather than just
    # being internally consistent with its own filtering.
    ionization_deposits = [d for d in engine.deposits if d.channel == Channel.IONIZATION]
    assert len(electron_deposits) == len(ionization_deposits)

    # One Line2D per distinct charge state actually present (not one per
    # segment -- plot_track groups consecutive same-q segments into a
    # single NaN-separated artist per q for scalability).
    distinct_q = {d.q for d in proton_deposits}
    assert len(ax.lines) == len(distinct_q)

    # No data lost in the regrouping: every proton deposit's plotted
    # position (z, x) must appear somewhere among the finite (non-NaN)
    # points of the drawn lines. This is the simplest direct check of the
    # "no data lost" invariant -- it doesn't re-derive run boundaries, it
    # just confirms every vertex the old per-segment implementation would
    # have drawn is still present in the new grouped-line data.
    plotted_points = set()
    for line in ax.lines:
        for z, x in zip(line.get_xdata(), line.get_ydata()):
            if np.isfinite(z) and np.isfinite(x):
                plotted_points.add((float(z), float(x)))
    for d in proton_deposits:
        assert (float(d.position[2]), float(d.position[0])) in plotted_points

    electron_collections = [c for c in ax.collections if c.get_label() == "e- (deposit)"]
    if electron_deposits:
        assert len(electron_collections) == 1
        assert len(electron_collections[0].get_offsets()) == len(electron_deposits)
    else:
        assert electron_collections == []


def test_plot_track_default_markers_smaller_than_old_values(full_xs, water, full_edep, full_eemean):
    entry_copy, engine = _build_transported_track(full_xs, water, full_edep, full_eemean)
    ax = plot_track(entry_copy, engine)

    channel_collections = [c for c in ax.collections if c.get_label() != "e- (deposit)"]
    assert channel_collections, "expected at least one reaction-channel marker collection"
    for c in channel_collections:
        assert all(s < 36 for s in c.get_sizes())  # old rcParams-default scatter size was 36

    electron_collections = [c for c in ax.collections if c.get_label() == "e- (deposit)"]
    if electron_collections:
        assert all(s < 80 for s in electron_collections[0].get_sizes())  # old hardcoded s=80


def test_plot_track_style_override(full_xs, water, full_edep, full_eemean):
    entry_copy, engine = _build_transported_track(full_xs, water, full_edep, full_eemean)
    ax = plot_track(entry_copy, engine, marker_size=40.0, electron_marker_size=90.0, line_width=4.0)

    channel_collections = [c for c in ax.collections if c.get_label() != "e- (deposit)"]
    for c in channel_collections:
        assert all(s == 40.0 for s in c.get_sizes())

    electron_collections = [c for c in ax.collections if c.get_label() == "e- (deposit)"]
    if electron_collections:
        assert all(s == 90.0 for s in electron_collections[0].get_sizes())

    for line in ax.lines:
        assert line.get_linewidth() == 4.0


def test_plot_track_aspect_default_is_auto(full_xs, water, full_edep, full_eemean):
    entry_copy, engine = _build_transported_track(full_xs, water, full_edep, full_eemean)
    ax = plot_track(entry_copy, engine)
    assert ax.get_aspect() == "auto"


def test_plot_track_aspect_equal(full_xs, water, full_edep, full_eemean):
    entry_copy, engine = _build_transported_track(full_xs, water, full_edep, full_eemean)
    ax = plot_track(entry_copy, engine, aspect="equal")
    assert ax.get_aspect() == 1.0  # matplotlib normalizes "equal" to the float 1.0


def test_plot_track_aspect_numeric_ratio(full_xs, water, full_edep, full_eemean):
    entry_copy, engine = _build_transported_track(full_xs, water, full_edep, full_eemean)
    ax = plot_track(entry_copy, engine, aspect=5.0)
    assert ax.get_aspect() == 5.0


def test_plot_track_channel_marker_sizes_override(full_xs, water, full_edep, full_eemean):
    entry_copy, engine = _build_transported_track(full_xs, water, full_edep, full_eemean)
    ax = plot_track(entry_copy, engine, marker_size=9.0,
                     channel_marker_sizes={Channel.IONIZATION: 60.0})

    ionization_collections = [c for c in ax.collections if c.get_label() == Channel.IONIZATION.value]
    assert ionization_collections, "expected ionization markers for this track"
    assert all(s == 60.0 for s in ionization_collections[0].get_sizes())

    # Channels not named in channel_marker_sizes fall back to marker_size.
    other_collections = [c for c in ax.collections
                          if c.get_label() not in (Channel.IONIZATION.value, "e- (deposit)")]
    assert other_collections, "expected at least one other channel present to check fallback"
    for c in other_collections:
        assert all(s == 9.0 for s in c.get_sizes())


def test_plot_track_channel_marker_sizes_default_matches_marker_size(full_xs, water, full_edep, full_eemean):
    entry_copy, engine = _build_transported_track(full_xs, water, full_edep, full_eemean)
    ax_default = plot_track(entry_copy, engine, marker_size=17.0)
    ax_explicit_none = plot_track(entry_copy, engine, marker_size=17.0, channel_marker_sizes=None)

    default_sizes = sorted(s for c in ax_default.collections if c.get_label() != "e- (deposit)"
                            for s in c.get_sizes())
    explicit_sizes = sorted(s for c in ax_explicit_none.collections if c.get_label() != "e- (deposit)"
                             for s in c.get_sizes())
    assert default_sizes == explicit_sizes == [17.0] * len(default_sizes)


def test_plot_track_legend_fontsize(full_xs, water, full_edep, full_eemean):
    entry_copy, engine = _build_transported_track(full_xs, water, full_edep, full_eemean)
    ax = plot_track(entry_copy, engine, legend_fontsize=6.0)
    legend = ax.get_legend()
    assert legend is not None
    for text in legend.get_texts():
        assert text.get_fontsize() == 6.0


def test_save_figure_writes_png_and_pdf_and_creates_dir(tmp_path):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    out_dir = tmp_path / "output_99"            # does not exist yet
    paths = save_figure(fig, out_dir, "01_test_plot")
    plt.close(fig)

    assert paths == [out_dir / "01_test_plot.png", out_dir / "01_test_plot.pdf"]
    assert paths[0].read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert paths[1].read_bytes()[:5] == b"%PDF-"
