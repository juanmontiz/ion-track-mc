# `ion_track_mc.physics` — cross sections, medium, energy loss, kinematics

← [Back to the main README](../../README.md)

This package holds the **physics data and per-collision physics rules**:
what can happen to a projectile at energy $E$ in charge state $q$, how
likely each outcome is, how much energy it costs, and which way things fly
afterwards. It does **not** move particles. Transport (sampling step
lengths, choosing a channel, looping over collisions) lives in
[`transport/`](../transport/README_transport.md). This package only answers
questions like "what is $\sigma_c(E)$?" or "how much energy does this
capture event remove?".

The formulas follow the reference code TILDA-V (Quinto et al. 2017).
Raw data file format: [`data/README_data_format.md`](../data/README_data_format.md).

| File | Role |
|---|---|
| `channels.py` | `Channel` enum: the closed set of interaction channels |
| `cross_section.py` | `CrossSectionSet` interface, table loaders, log-log interpolation |
| `medium.py` | `Medium` dataclass (density, molar mass, target identity), `WATER` |
| `energy_loss.py` | Per-channel projectile energy loss; elastic scattering |
| `kinematics.py` | Direction of the ejected δ-electron in ionization |

---

## 1. Interaction channels (`channels.py`)

The projectile is a hydrogen ion whose **charge state $q$ changes during
transport** ($q=1$: $\mathrm{H}^+$, $q=0$: $\mathrm{H}^0$).
Each charge state has its own set of open channels:

| Channel | Open for | Effect on projectile | Secondary |
|---|---|---|---|
| `IONIZATION` | $q=0,1$ | loses $E_{\text{ion}}$, keeps $q$ | one δ-electron |
| `EXCITATION` | $q=0,1$ | loses $E_{\text{exc}}$, keeps $q$ | none |
| `CAPTURE` | $q=1$ only | loses $E_{\text{cap}}$, $q: 1\to 0$ | none |
| `LOSS` | $q=0$ only | loses $E_{\text{loss}}$, $q: 0\to 1$ | none (see §4) |
| `ELASTIC` | $q=0,1$ | loses $E_{\text{elastic}}$, **deflected** | none |

Which channels are "open" for a given $q$ is not hard-coded anywhere: it is
whatever $(q, \text{channel})$ tables were loaded into the
`CrossSectionSet` (§2). Adding a channel is a one-line change to the enum
plus a branch in `TrackEngine._interact`.

## 2. Cross sections (`cross_section.py`)

### 2.1 The interface

`CrossSectionSet` is an abstract base class, the contract the transport
engine depends on:

| Method | Returns |
|---|---|
| `channels(q)` | list of open channels for charge state $q$ |
| `sigma(channel, E, q)` | $\sigma_c^{(q)}(E)$ in $\mathrm{nm}^2$ |
| `sigma_total(E, q)` | $\sigma_{\text{tot}}^{(q)}(E) = \sum_c \sigma_c^{(q)}(E)$ |
| `energy_range(q)` | $(E_{\min}, E_{\max})$, the validity range for $q$ |

To swap in a different medium or a different theory, write a new
`CrossSectionSet` and inject it into `TrackEngine`. The engine does not
change.

### 2.2 Tabulated implementation and log-log interpolation

`TabulatedCrossSectionSet` takes one $(E_i, \sigma_i)$ table per
$(q, \text{channel})$ and interpolates **linearly in log-log space**. For
$E_i \le E \le E_{i+1}$:

$$
\ln\sigma(E) \;=\; \ln\sigma_i \;+\;
\frac{\ln E - \ln E_i}{\ln E_{i+1} - \ln E_i}\,
\bigl(\ln\sigma_{i+1} - \ln\sigma_i\bigr).
$$

This is exact for any piece of the cross section that behaves as a power
law $\sigma \propto E^{a}$ between grid points. Over a decade of energy,
cross sections are close to power laws, so the interpolation is accurate on
a coarse grid. Linear-linear interpolation is not, on a 40-point grid
spanning $8\,\mathrm{keV}$–$300\,\mathrm{MeV}$. The same rule, exposed as
the standalone function `interp_log_log`, is reused for the mean-energy
tables in §4.

**No extrapolation.** Any query with $E$ outside the tabulated
$[E_{\min}, E_{\max}]$ raises `ValueError`. For one
charge state, `energy_range(q)` is the **intersection** of all its channels'
ranges,

$$
E_{\min}^{(q)} = \max_c E_{\min,c}, \qquad E_{\max}^{(q)} = \min_c E_{\max,c},
$$

so every open channel is valid at every energy inside it. All current
proton tables cover $8\,\mathrm{keV}$–$300\,\mathrm{MeV}$ with 40 points.
Log-log interpolation also needs **strictly positive** values: a table
containing a zero (for example the loss channel's all-zero `Eemeantot` file)
must never be queried, and `interp_log_log` raises if it is.

### 2.3 Loading and units

`load_cross_section_table` and `load_mean_energy_table` parse the
two-column `.dat` files. Each converts to internal units **once, at load
time**, using the factor tables in `constants.py` (see
[`README_constants.md`](../README_constants.md)). Current data: energies in
$\mathrm{eV}$, cross sections in $10^{-16}\,\mathrm{cm}^2$, which is
$10^{-2}\,\mathrm{nm}^2$. Energies must be strictly increasing; otherwise
loading raises.

## 3. Medium (`medium.py`)

Water is treated as **one effective molecular target**, with molecular
cross sections, not decomposed into H and O atoms. A
`Medium` supplies the **molecular number density** used for the
macroscopic cross section $\Sigma = N\sigma$:

$$
N \;=\; \frac{\rho\,N_A}{M_{\text{mol}}}
\qquad\Longrightarrow\qquad
N_{\mathrm{H_2O}} = \frac{1\,\mathrm{g\,cm^{-3}} \times 6.022\times10^{23}\,\mathrm{mol^{-1}}}{18.015\,\mathrm{g\,mol^{-1}}}
\approx 3.343\times10^{22}\,\mathrm{cm^{-3}} = 33.43\,\mathrm{nm^{-3}}.
$$

`Medium.N` returns this value in $\mathrm{nm}^{-3}$.

Two optional fields feed elastic scattering only (§4.5):
`target_nuclear_charge` $Z_T$ and `target_nucleon_count` $A_T$. For water,
$Z_T = 10$ ($1+1+8$) and $A_T = 18$. They default to `None`, so a custom
medium that forgets them fails loudly at the first elastic collision rather
than silently.

`check_medium_matches_filename` cross-checks a data file's `{Target}` token
(for example `H2O` in `H1+_+_H2O_ionization_TCS_CDW-EIS.dat`) against the
medium's `target_label`. This guards against wiring a PMMA table to water.

## 4. Projectile energy loss per channel (`energy_loss.py`)

Every loss is **deterministic for a given $E$** except elastic. The
deterministic losses come from two kinds of tables tabulated per
(projectile, channel), both log-log interpolated (§2.2):

- $\bar E_{\text{dep}}(E)$ (`Edepmeantot`): mean energy deposited at the
  collision site;
- $\bar E_{e}(E)$ (`Eemeantot`): mean kinetic energy carried away by the
  emitted electron.

Notation: $M_P = m_p/m_e = 1836.15267$ is the projectile mass in electron
masses, and $I_H = 13.6\,\mathrm{eV}$ is the H ground-state binding energy.
The formulas follow TILDA-V (Quinto et al. 2017).

### 4.1 Ionization ($q = 0, 1$)

$$
E_{\text{ion}}(E) \;=\; \underbrace{\bar E_{\text{dep}}(E)}_{\text{local deposit}}
\;+\; \underbrace{\bar E_{e}(E)}_{W,\ \text{δ-electron energy}}.
$$

`ionization_energy` returns the pair $(\bar E_{\text{dep}}, W)$ separately.
The engine deposits the first locally and gives the second to a new
δ-electron as its kinetic energy $W$.

### 4.2 Excitation ($q = 0, 1$)

$$
E_{\text{exc}}(E) \;=\; \bar E_{\text{dep}}(E).
$$

### 4.3 Electron capture ($q = 1 \to 0$)

$$
E_{\text{cap}}(E) \;=\; \bar E_{\text{dep}}(E) \;+\; \frac{E}{M_P} \;-\; I_H .
$$

The middle term is the kinetic energy needed to bring the captured electron
up to the projectile's velocity: an electron moving at the speed of a
projectile with energy $E$ has kinetic energy $E\,m_e/m_p = E/M_P$. The
$-I_H$ term is the binding energy released when that electron settles into
the H ground state. Checked over the full real data grid, $E_{\text{cap}}$
is **positive everywhere**: $5.36\,\mathrm{eV}$ at $8\,\mathrm{keV}$, rising
to $164\,\mathrm{keV}$ at $300\,\mathrm{MeV}$. There is therefore no
"negative deposit" case to handle (test:
`test_capture_delta_E_positive_across_full_real_grid`).

### 4.4 Electron loss ($q = 0 \to 1$)

$$
E_{\text{loss}} \;=\; I_H \;=\; 13.6\,\mathrm{eV}.
$$

This is a **lower bound**. The H atom must at least pay its electron's
binding energy. What the electron does afterwards is unknown in this model,
so it is not tracked. In code it is still a table lookup
(`mean_energy_lookup` on a file that is flat at $13.6\,\mathrm{eV}$), so no
special case is needed.

### 4.5 Elastic scattering ($q = 0, 1$): the only stochastic loss

`elastic_scatter` runs the full chain: momentum, screening cutoff angle,
sampled scattering angle, energy transfer. Relativistic projectile
momentum (in $\mathrm{eV}/c$):

$$
p c \;=\; \sqrt{\bigl(E + M_P\,m_ec^2\bigr)^2 - \bigl(M_P\,m_ec^2\bigr)^2},
\qquad m_ec^2 = 510\,999.0645\,\mathrm{eV}.
$$

Screening cutoff angle (Thomas–Fermi-like screening of the target nucleus by
its electrons):

$$
\theta_{\text{cut}} \;=\; \alpha\,Z_T^{1/3}\,\frac{m_ec^2}{pc},
\qquad \alpha = 1/137 .
$$

The scattering angle is sampled from a uniform $\eta\sim U(0,1)$ by

$$
\theta \;=\; 2\arcsin\Bigl\{\eta\bigl[1 - \sin^{-2}(\theta_{\text{cut}}/2)\bigr]
+ \sin^{-2}(\theta_{\text{cut}}/2)\Bigr\}^{-1/2}.
$$

**Why this formula?** Write $u = \sin^{-2}(\theta/2)$. The formula says $u$
is uniform between $u(\theta_{\text{cut}})$ (at $\eta=0$) and $u(\pi)=1$ (at
$\eta=1$). Since
$\mathrm{d}u = -\sin^{-3}(\theta/2)\cos(\theta/2)\,\mathrm{d}\theta$ and
$\sin\theta = 2\sin(\theta/2)\cos(\theta/2)$,

$$
\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega} \propto \frac{1}{\sin^4(\theta/2)}
\quad\Longrightarrow\quad
\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega}\,2\pi\sin\theta\,\mathrm{d}\theta
\;\propto\; \frac{\sin\theta}{\sin^4(\theta/2)}\,\mathrm{d}\theta
\;=\; -2\,\mathrm{d}u ,
$$

So a uniform $u$ is exactly **inverse-transform sampling of the Rutherford
distribution**, truncated below $\theta_{\text{cut}}$. Consequently
$\theta \in [\theta_{\text{cut}}, \pi]$, never $0$: every elastic event
deflects the projectile by at least $\theta_{\text{cut}}$.

The energy transferred to the recoiling target (two-body elastic
kinematics) is

$$
E_{\text{elastic}} \;=\; 4E\,\frac{M_P\,M_T}{(M_P + M_T)^2}\,\sin^2(\theta/2),
\qquad M_T = A_T\,M_P \;\;(=18\,M_P\ \text{for water}).
$$

The new direction is the incident one rotated by polar angle $\theta$ and
azimuth $\phi \sim U(0, 2\pi)$, using the same frame rotation as §5.

*Approximations to be aware of.* $M_T = A_T M_P$ ignores nuclear binding
and the neutron–proton mass difference. The single $\theta$ is used both in
the energy-transfer formula (where, strictly, the centre-of-mass angle
belongs) and as the lab-frame deflection. For a projectile 18× lighter than
its target, the two angles differ only at order $M_P/M_T \approx 5\%$.

**Policy:** elastic energy *is* deposited and does deflect the track, but
it is **excluded from the electronic stopping power** reported by
`analysis.stopping_power.analytic_stopping_power`, matching TILDA-V's Eq. 16
(see [`analysis/README_analysis.md`](../analysis/README_analysis.md)).

## 5. δ-electron emission direction (`kinematics.py`)

The δ-electron's energy is fixed, $W = \bar E_e(E)$ (§4.1). Its polar
angle relative to the projectile direction comes from momentum
conservation in a **binary encounter**, $\vec p_1 = \vec p_s + \vec p_e$
(incident = scattered projectile + electron). With non-relativistic
$p^2 = 2mE$, in units where $m_e = 1$:

$$
k_1 = 2M_PE, \qquad k_s = 2M_P\,(E - E_{\text{ion}}), \qquad k_e = 2W,
$$

where each $k$ is a squared momentum. The law of cosines on the momentum
triangle, $k_s = k_1 + k_e - 2\sqrt{k_1k_e}\cos\theta_e$, gives

$$
\cos\theta_e \;=\; \frac{k_1 + k_e - k_s}{2\sqrt{k_1\,k_e}} .
$$

**Domain failure at low $E$ and the isotropic fallback.** At low
projectile energy this argument exceeds 1 (for $\mathrm{H}^+$: at
$8$, $10$ and $20\,\mathrm{keV}$; for $\mathrm{H}^0$: up to $30\,\mathrm{keV}$).
This is a real triangle-inequality violation, not a bug: a slow, heavy
projectile must lose momentum $\Delta p \approx E_{\text{ion}}\sqrt{M_P/2E}$,
which grows as $E$ drops, and that exceeds the electron's momentum
$p_e = \sqrt{k_e}$. Closing the momentum balance then needs a third body
(the recoiling target ion), which a two-body formula cannot represent.
**Decision:** whenever $|\cos\theta_e| > 1$, sample isotropically,
$\cos\theta_e \sim U(-1, 1)$. No energy threshold is hard-coded; the check
runs per event.

**Frame rotation.** `_orthonormal_basis(w)` builds a right-handed
orthonormal frame $\{\hat u, \hat v, \hat w\}$ around the incident direction
$\hat w$, and the new direction is

$$
\hat d \;=\; \sin\theta\cos\phi\,\hat u + \sin\theta\sin\phi\,\hat v + \cos\theta\,\hat w,
\qquad \phi\sim U(0,2\pi).
$$

The helper's reference axis is $\hat x$ unless $|w_x| \ge 0.9$, in which
case it switches to $\hat y$. This avoids a near-degenerate cross product
when $\hat w$ is nearly parallel to the reference. Elastic scattering
(§4.5) reuses the same helper.

---

## Fidelity level

| Piece | Status |
|---|---|
| Cross sections (CDW-EIS ionization/capture, Miller–Green excitation/loss, classical elastic) | **Implemented**, real tabulated data, 8 keV–300 MeV |
| Log-log interpolation, no extrapolation | **Implemented** |
| Water as one effective molecular target | **Model choice**; composite media (tissue, PMMA) not supported |
| Ionization / excitation / capture losses from `Edepmeantot`/`Eemeantot` | **Implemented** |
| Loss channel $E_{\text{loss}} = 13.6\,\mathrm{eV}$ | **Model choice**, a lower bound |
| Elastic angle sampling + recoil energy | **Implemented**; deflection magnitude not yet validated against range-straggling data |
| δ-electron direction (binary encounter + isotropic fallback) | **Implemented** |
| δ-electron energy $W$ | deterministic mean $\bar E_e(E)$, not sampled from a distribution |

## Complexity and performance

Let $n$ be the number of grid points per table ($n = 40$ now) and $K$ the
number of open channels for one $q$ ($K = 4$ now).

- `TabulatedCrossSectionSet.sigma`: $O(\log n)$ per query (binary search
  inside `np.interp`). $\ln E_i$ and $\ln\sigma_i$ are precomputed at
  construction.
- `sigma_total`: $K$ calls to `sigma`, plus one call to `channels`.
- `channels(q)`: $O(K)$, a copy of a sorted list built once per $q$ at
  construction. It is called several times per collision (tens of
  thousands of times per track), so rebuilding it on every call was a
  measurable cost before it was cached.
- `interp_log_log` (used for mean-energy tables): recomputes the logs of the
  whole table on each call, so $O(n)$. That cost is negligible at $n = 40$.
- `_orthonormal_basis`: runs once per δ-electron and per elastic event.
  Its two cross products are written out component by component, because
  numpy's generic `np.cross` has high fixed overhead for 3-vectors. With
  `np.cross` this was the single largest hot spot in a 500 keV/u track.
  A test checks that the hand-written version matches `np.cross`.

## Tests

`tests/test_cross_section.py`, `test_medium.py`, `test_energy_loss.py`,
`test_kinematics.py`.
