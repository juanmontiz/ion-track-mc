# `ion_track_mc.analysis` — beams, stopping-power validation, plots

← [Back to the main README](../../README.md)

This package sits on top of transport and scoring. It does three things:

- runs **many primaries** as a beam;
- checks the Monte Carlo machinery against **closed-form stopping power**,
  including charge-state equilibrium;
- draws the **standard figures**.

It contains no new physics. Every energy-loss rule it uses is the engine's
own (`TrackEngine._interact`), so there is exactly one place where physics
lives.

The stopping-power formulas follow TILDA-V (Quinto et al. 2017, Eqs. 15–16).

| File | Role |
|---|---|
| `beam.py` | `run_beam` → `BeamResult`: pencil beam with a Gaussian energy spread |
| `stopping_power.py` | Emergent (MC) and analytic $S(E)$, fixed-$q$ and charge-equilibrium |
| `visualization.py` | Style presets and standard plots |

---

## 1. Beams (`beam.py`)

`TrackEngine.transport()` runs exactly one primary. `run_beam` loops it
$N_p$ times into **one shared engine**, so scorers and `engine.deposits`
accumulate over the whole beam.

- **Energy spread.** $E_0 \sim \mathcal N(E_{\text{beam}}, \sigma_E)$, or
  exactly $E_{\text{beam}}$ if $\sigma_E \le 0$. Draws outside the tabulated
  range are **rejected and redrawn**, never clipped. Clipping would pile
  probability up at the range edges. The draw gives up after
  `max_resample_attempts`. The number of rejections is reported in
  `n_rejected_out_of_range`.
- **Independent random streams.** `SeedSequence(seed).spawn(2)` gives one
  stream for beam energies and one for the engine. The energy statistics
  therefore never correlate with collision sampling, and the run is
  reproducible from a single `seed`.
- **Beam LET in $O(1)$ memory.** After each primary, its `collision_log` is
  folded into running sums: $\sum s$ and $\sum\Delta E$, plus that track's
  own segment sum $\sum_k E_k^2/\ell_k$ if `let_segment_length_nm` was
  given. Segments never span two tracks.
  `BeamResult.track_averaged_LET()` and `dose_averaged_LET()` reduce these
  sums with the formulas in [`scoring/README_scoring.md`](../scoring/README_scoring.md) §3.
  `dose_averaged_LET()` raises if the beam was run without
  `let_segment_length_nm`, because $\bar L_D$ is only defined at a chosen
  segment length.

`BeamResult` also holds per-primary arrays `E_start_eV`, `E_final_eV` (the
residual energy when the track stopped at the tabulated floor) and
`q_final`.

## 2. Stopping power (`stopping_power.py`)

### 2.1 Fixed charge state: emergent vs. analytic

**Analytic.** The mean energy lost per unit path is the collision rate per
unit length times the mean loss per collision:

$$
S_{\text{analytic}}^{(q)}(E) \;=\; N\sum_c \sigma_c^{(q)}(E)\,\langle\Delta E_c\rangle(E).
$$

`mean_delta_E` supplies $\langle\Delta E_c\rangle$ using the same rules as
the engine: ionization $\bar E_{\text{dep}} + \bar E_e$, excitation
$\bar E_{\text{dep}}$, capture $\bar E_{\text{dep}} + E/M_P - 13.6$, loss
$13.6$ (all in eV). For **elastic** it returns $0$, by policy (see §2.3).

**Emergent (Monte Carlo).** `emergent_stopping_power` draws $n$
independent single collisions at fixed $(E, q)$, each with a sampled step
$s_i$ and a sampled channel, and the real loss $\Delta E_i$ from
`TrackEngine._interact`. It then forms the ratio estimator

$$
S_{\text{MC}}^{(q)}(E) \;=\; \frac{\sum_i \Delta E_i}{\sum_i s_i}.
$$

**Why the two must agree.** By the law of large numbers the ratio tends to
$\mathbb E[\Delta E]/\mathbb E[s]$. Step length and channel are drawn
independently, with $\mathbb E[s] = 1/\Sigma_{\text{tot}}$ and
$\mathbb E[\Delta E] = \sum_c p_c\langle\Delta E_c\rangle$, where
$p_c = \sigma_c/\sigma_{\text{tot}}$. So

$$
\frac{\mathbb E[\Delta E]}{\mathbb E[s]}
\;=\; N\sigma_{\text{tot}}\sum_c\frac{\sigma_c}{\sigma_{\text{tot}}}\langle\Delta E_c\rangle
\;=\; N\sum_c\sigma_c\langle\Delta E_c\rangle \;=\; S_{\text{analytic}}.
$$

This is a **self-consistency check**. It tests the sampling machinery
(`sample_step`, `select_channel`, `_interact`) against the exact
expectation of the *same* tables. It is not a comparison with PSTAR or ICRU
data. `visualization.plot_stopping_power_self_consistency` plots it.

### 2.2 Charge-state equilibrium

A real proton beam cycles between $\mathrm{H}^+$ and $\mathrm{H}^0$ through
capture and loss. Model the charge state as a two-state process along the
path. $\mathrm{H}^+$ leaves its state at rate $N\sigma_C$ per unit length
(capture) and $\mathrm{H}^0$ at rate $N\sigma_L$ (loss). In steady state the
flux out of each state balances the flux in:

$$
f_{\mathrm{H^+}}\,N\sigma_C \;=\; f_{\mathrm{H^0}}\,N\sigma_L,
\qquad f_{\mathrm{H^+}} + f_{\mathrm{H^0}} = 1,
$$

$$
\boxed{\;f_{\mathrm{H^+}} = \frac{\sigma_L}{\sigma_L + \sigma_C},\qquad
f_{\mathrm{H^0}} = \frac{\sigma_C}{\sigma_L + \sigma_C}\;}
$$

These are Alcocer-Ávila 2021, Eqs. III.3/III.4, implemented in
`equilibrium_charge_fractions`. $f$ is the **fraction of path length**
spent in each state; $N$ cancels. The equilibrium stopping power blends
the two fixed-$q$ curves with these weights (Eq. III.8, the quantity
TILDA-V reports):

$$
S_{\text{eq}}(E) \;=\; f_{\mathrm{H^+}}(E)\,S^{(1)}(E) \;+\; f_{\mathrm{H^0}}(E)\,S^{(0)}(E).
$$

**This is the curve to compare against measured "proton" stopping power**,
not $S^{(1)}$ alone. Capture matters most at low energy; at high $E$,
$\sigma_C\to0$ and $f_{\mathrm{H^+}}\to1$.

There are three ways to compute $S_{\text{eq}}$, from cheapest and least
independent to most independent:

| Function | Charge mix | Stopping powers |
|---|---|---|
| `analytic_equilibrium_stopping_power` | closed form | closed form |
| `emergent_equilibrium_stopping_power` | closed form | two fixed-$q$ MC campaigns |
| `sample_equilibrium_stopping_power` | **emerges** from real capture/loss sampling | same single MC chain |

The last one uses `sample_equilibrium_charge_chain`. It follows one
projectile through $n$ collisions in TILDA-V's **stationary mode**: $E$ is
held fixed and only $q$ evolves, flipping on real CAPTURE/LOSS events from
`_interact`. Its empirical fractions are path-length weighted,

$$
\hat f_{\mathrm{H^+}} = \frac{\sum_{i:\,q_i=1} s_i}{\sum_i s_i},
$$

which is the right weighting: an $\mathrm{H}^0$ has different mean free
paths, so counting collisions instead would bias $\hat f$. The ratio
$\sum\Delta E_i/\sum s_i$ over the same chain equals
$\sum_q \hat f_q\,\hat S^{(q)}$ identically, so it is the
fully emergent $S_{\text{eq}}$. At $100\,\mathrm{keV/u}$ it matches the
closed form to about $0.002$ absolute in $f$ and about $0.004\%$ relative in
$S$. Tests probe only $100\,\mathrm{keV/u}$, because at higher energy
capture becomes too rare to converge in a test-sized budget.

### 2.3 Elastic: a deliberate, permanent divergence

Elastic collisions **do** deposit recoil energy in the engine, so $S_{\text{MC}}$
includes them. `mean_delta_E` returns $0$ for elastic, so
$S_{\text{analytic}}$ excludes them. This matches TILDA-V's Eq. 16, which
sums only the electronic processes $j, k = 1..3$, never elastic. The gap
is about $0.008\%$–$0.08\%$ of $S$ over $8\,\mathrm{keV}$–$2\,\mathrm{MeV}$,
well inside test tolerance. It is policy, not a bug to drive to zero.

## 3. Plots (`visualization.py`)

**Style presets** set Matplotlib `rcParams` globally for whatever is
plotted next: `use_paper_style()` (serif, small, for print) and
`use_presentation_style()` (sans-serif, large, for slides).

| Function | Plots |
|---|---|
| `plot_cross_sections(xs, q)` | $\sigma_c^{(q)}(E)$ per open channel, log-log |
| `plot_mean_free_path(xs, medium, q)` | $\lambda^{(q)}(E) = 1/(N\sigma_{\text{tot}}^{(q)}(E))$, log-log |
| `plot_stopping_power_self_consistency(...)` | $S_{\text{MC}}$ vs. $S_{\text{analytic}}$ (§2.1); **two panels** `(ax_log, ax_lin)` |
| `plot_bragg_curve(scorer)` | $D(z)$ from `make_bragg_curve_scorer` |
| `plot_radial_dose(scorer)` | $D(r)$, log-log |
| `plot_track(entry_particle, engine)` | one primary's track, $z$ vs. $x$, coloured by charge state, markers by channel, electron deposit points |

Each function accepts an optional `ax` and returns the `Axes` it drew on.
`plot_stopping_power_self_consistency` is the exception: it takes and
returns a pair. `plot_track` needs a copy of the primary **before**
transport, because `transport()` mutates its position. It also supports
only one track per engine, since deposits carry no track id.

**Saving figures.** `save_figure(fig, out_dir, name)` writes
`out_dir/name.png` (300 dpi) and `out_dir/name.pdf`, creating `out_dir`
if needed. Every notebook saves each of its figures this way into its own
folder, `notebooks/output_01/`, `output_02/` or `output_03/`. File names
are numbered in notebook order, for example `02_pencil_beam_kernel.pdf`.
These folders are regenerated on every run and are git-ignored.

PSTAR and experimental overlays are kept in the notebooks (`notebooks/01`,
`02`), not here. `plot_f_y` does not exist yet, because there is no $f(y)$
to plot.

---

## Fidelity level

| Piece | Status |
|---|---|
| `run_beam` | **Implemented**: pencil beam, single entry point and direction, Gaussian $E$ |
| Fixed-$q$ MC vs. analytic $S(E)$ | **Implemented**, agrees (self-consistency) |
| Equilibrium fractions and $S_{\text{eq}}$ (closed form + two MC routes) | **Implemented**, cross-validated at 100 keV/u |
| Comparison with PSTAR / experiment | in notebook 01 |
| Bragg curve / LET vs. reference data | not done |
| Multi-primary track plotting | not supported (one track per plot) |

## Complexity

- `run_beam`: $N_p$ × (cost of one track); see
  [`transport/README_transport.md`](../transport/README_transport.md). Memory
  is dominated by `engine.deposits`, which grows with the total number of
  deposits.
- `emergent_stopping_power`: $O(n)$ single collisions, each
  $O(K\log n_{\text{grid}})$. `plot_stopping_power_self_consistency` costs
  `n_energy_points × n_samples_per_point` collisions (default 12 × 8000).
- `analytic_*`: $O(K)$ table lookups per energy, effectively instant.

## Tests

`tests/test_beam.py` (monoenergetic and Gaussian beams, rejection, energy
conservation, shared scorer, reproducibility, beam LET);
`test_stopping_power.py` (MC vs. analytic, elastic-excluded policy,
fractions summing to 1 and matching the cross-section ratio, all three
$S_{\text{eq}}$ routes); `test_visualization.py`.
