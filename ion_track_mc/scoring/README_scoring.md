# `ion_track_mc.scoring` — turning deposits into dose, LET, and microdosimetry

← [Back to the main README](../../README.md)

The transport engine produces a stream of **local energy deposits**: an
energy $\varepsilon_i$ [eV] at a position $\vec r_i$ [nm]. This package turns
that stream into the quantities a radiobiologist reports:

- **dose** on a voxel grid or in radial shells;
- **LET**, track-averaged and dose-averaged;
- **single-event lineal energy $y$ and specific energy $z$** in a spherical
  microdosimetric site.

Where the deposits come from: [`transport/README_transport.md`](../transport/README_transport.md) §5.4.

| File | Role |
|---|---|
| `scorer.py` | `Scorer` abstract base class, the observer interface |
| `dose.py` | `VoxelDoseScorer`, `RadialDoseScorer`, `make_bragg_curve_scorer` |
| `radiobiology.py` | LET functions, $y$, $z$, `SphericalSiteScorer` |

---

## 1. The observer pattern (`scorer.py`)

A `Scorer` has one method:

```python
def score(self, position: np.ndarray, energy: float) -> None
```

`TrackEngine` calls it on **every injected scorer, once per local deposit**.
This covers ion deposits and electron deposits alike. Scorers never read
each other or the engine's internal lists, so any number of them can watch
the same run independently: a Bragg-curve grid, a radial-dose profile and a
microdosimetric site, all at once. Scorers **accumulate across
`transport()` calls**, so a beam is scored by injecting the scorers once and
transporting many primaries.

To add a new observable, subclass `Scorer`, implement `score`, and inject
it. Neither the engine nor the other scorers change.

## 2. Dose (`dose.py`)

Absorbed dose in a scoring region $V$ of density $\rho$:

$$
D \;=\; \frac{\sum_{i\in V}\varepsilon_i}{m_V} \;=\; \frac{\sum_{i\in V}\varepsilon_i}{\rho\,V}
\qquad [\mathrm{Gy} = \mathrm{J\,kg^{-1}}].
$$

The unit conversion is
$\varepsilon[\mathrm{J}] = \varepsilon[\mathrm{eV}]\times1.602176634\times10^{-19}$,
$V[\mathrm{cm^3}] = V[\mathrm{nm^3}]\times10^{-21}$, and
$m[\mathrm{kg}] = \rho[\mathrm{g\,cm^{-3}}]\,V[\mathrm{cm^3}]\times10^{-3}$.

At track-structure scale, single tracks give very large doses in
nanometre-sized regions. That is real physics, not a bug: a few eV in a
$\mathrm{nm}^3$ of water is of order $10^5$–$10^6\,\mathrm{Gy}$.

### 2.1 Cartesian voxels: `VoxelDoseScorer`

The caller gives explicit bin edges per axis (`x_edges`, `y_edges`,
`z_edges`, in nm), so voxels may be anisotropic. Each deposit is binned by
binary search on each axis. Voxel $(i,j,k)$ has volume
$V_{ijk} = \Delta x_i\,\Delta y_j\,\Delta z_k$.

**Deposits outside the grid are silently dropped.** To check whether the
grid captured everything, compare `total_deposited_eV()` with the run's
known total energy loss.

**Bragg curve.** `make_bragg_curve_scorer(medium, z_edges, w)` builds a
grid with a single transverse bin $[-w, w]\times[-w, w]$ and many depth
bins. The result `dose_Gy()[0, 0, :]` against `z_edges` is the depth–dose
curve $D(z)$.

### 2.2 Radial shells: `RadialDoseScorer`

This scorer bins deposits by **perpendicular distance from a track axis**
(an infinite line through `axis_point` along unit vector $\hat a$):

$$
r \;=\; \bigl\|\,(\vec r_i - \vec r_0) - [(\vec r_i - \vec r_0)\cdot\hat a]\,\hat a\,\bigr\| .
$$

Shell $k$ is the annular cylinder $r_k \le r < r_{k+1}$ of length $L$:

$$
V_k \;=\; \pi\bigl(r_{k+1}^2 - r_k^2\bigr)\,L,
\qquad D(r) \text{ from } D_k = \frac{\sum_{i\in k}\varepsilon_i}{\rho\,V_k}.
$$

$L$ (`length_nm`) is used **only to normalize volume**. The scorer does not
filter deposits by their position along the axis, so choose $L$ to match
the axial extent the deposits actually cover. With the local-deposit electron
model (δ-electrons deposited at their point of creation), the radial profile reflects
where δ-electrons are *created*, not where they would have travelled. The
real $D(r)$ tail out to the δ-ray range is therefore missing.

## 3. Linear energy transfer (`radiobiology.py`)

LET uses the primary's own collisions only: the `(s_i, ΔE_i)` pairs in
`TrackEngine.collision_log`. Here $\Delta E_i$ is the projectile's **total**
loss in collision $i$, including the δ-electron energy, so these are
**unrestricted** LETs. Define the local LET of collision $i$ as
$L_i = \Delta E_i / s_i$.

**Track-averaged LET**, the path-length-weighted mean of $L_i$:

$$
\bar L_T \;=\; \frac{\sum_i s_i\,L_i}{\sum_i s_i} \;=\; \frac{\sum_i \Delta E_i}{\sum_i s_i}.
$$

**Dose-averaged LET** is the energy-weighted mean of the local LET. It
is measured over **track segments**, not single collisions (see the next
box for why). Cut the track into segments of length $\Delta$
(`segment_length_nm`, a required argument). Collision $i$ sits at path
position $x_i=\sum_{j\le i}s_j$ and belongs to segment $k=\lfloor x_i/\Delta\rfloor$.
With $E_k=\sum_{i\in k}\Delta E_i$ and $L_k=E_k/\ell_k$,

$$
\bar L_D(\Delta) \;=\; \frac{\sum_k E_k\,L_k}{\sum_k E_k}
\;=\; \frac{\sum_k E_k^2/\ell_k}{\sum_k E_k} .
$$

Every segment has $\ell_k=\Delta$, with two exceptions. The leftover piece
at the end of the track, shorter than $\Delta$, is **merged into the last
full segment**, because a tiny last segment would bring back the
$1/\ell$ blow-up. A track shorter than $\Delta$ is one segment of its own
length. Segments with no collisions have $E_k=0$ and add nothing.

$\bar L_D \ge \bar L_T$ always holds; it is a weighted Cauchy–Schwarz
inequality (tested).

> **Why segments: the single-collision version does not converge.** One
> could try $L_i=\Delta E_i/s_i$ for each collision. In
> event-by-event transport the step $s_i$ is exponentially distributed
> (see the transport README, §3), and for an exponential variable
> $\mathbb{E}[1/s] = \int_0^\infty \Sigma e^{-\Sigma s}\,s^{-1}\,\mathrm{d}s$
> **diverges** logarithmically at $s\to0$. The sum $\sum_i\Delta E_i^2/s_i$
> was dominated by rare, very short steps and kept drifting upwards, with
> jumps, as collisions were added. On a synthetic track, $10^5$ and $10^6$
> collisions disagreed by 13%. With $\Delta=100\,\mathrm{nm}$ they agree
> within 2%, and a test enforces this. A fixed $\ell_k$ removes the random
> $1/s$.

**$\Delta$ sets the scale, so it has no default.** $\bar L_D(\Delta)$
depends on $\Delta$ by construction. As $\Delta$ grows to the whole track,
it tends to $\bar L_T$. As $\Delta$ shrinks to the mean free path (a few
nm), it goes back to the divergent single-collision value. Choose
$\Delta$ much larger than the mean free path and much smaller than the
distance over which LET changes. For example, $100\,\mathrm{nm}$ suits
tracks a few µm long; the notebooks use this value. Always report
$\bar L_D$ together with its $\Delta$.

**Units.** $1\,\mathrm{eV/nm} = 1\,\mathrm{keV/\mu m}$ exactly, so values in
eV/nm can be read directly as keV/µm, the conventional LET unit.

**Running sums.** The `*_from_sums` variants take
$\left(\sum s_i,\ \sum\Delta E_i\right)$ for $\bar L_T$ and
$\left(\sum_k E_k^2/\ell_k,\ \sum\Delta E_i\right)$ for $\bar L_D$.
`LET_segment_sums(collision_log, Δ)` computes the segment sum for one
track. A beam (`analysis.beam.run_beam`) can therefore accumulate LET in
$O(1)$ memory: each track is segmented on its own and its sums are added
to the running totals.

## 4. Single-event microdosimetry (`radiobiology.py`)

Microdosimetry looks at the energy $\varepsilon$ deposited in **one site**
(here a sphere of diameter $d$, modelling a cell nucleus subunit) by **one
event**.

**Mean chord length.** By Cauchy's theorem, a convex body crossed by
isotropic uniform random lines has mean chord $\bar\ell = 4V/S$. For a
sphere:

$$
\bar\ell \;=\; \frac{4\cdot\frac{4}{3}\pi R^3}{4\pi R^2} \;=\; \frac{4R}{3} \;=\; \frac{2}{3}\,d .
$$

**Lineal energy** (ICRU), in eV/nm, which equals keV/µm:

$$
y \;=\; \frac{\varepsilon}{\bar\ell}.
$$

**Specific energy** (ICRU), in Gy:

$$
z \;=\; \frac{\varepsilon}{m}, \qquad m = \rho\,\tfrac{4}{3}\pi R^3 .
$$

`SphericalSiteScorer(center, d)` sums every deposit with
$\|\vec r_i - \vec c\| \le R$ into `epsilon_eV`. It exposes `y()`, `z(medium)`
and `reset()`. One instance (or one `reset()` cycle) corresponds to one
event.

**Not implemented:** the distributions $f(y)$ and $d(y) = y f(y)/\bar y_F$.
Building them needs a Monte Carlo campaign with many primaries, randomized
impact parameters and an event definition; this toolkit does not include
it. `SphericalSiteScorer` is the per-event building block such a campaign
would reuse. MKM / RBE is not implemented either.

---

## Fidelity level

| Piece | Status |
|---|---|
| Observer interface | **Implemented** |
| Voxel / radial dose, Bragg-curve scorer | **Implemented**, textbook formulas |
| $\bar L_T$ | **Implemented**, unrestricted, primary only |
| $\bar L_D(\Delta)$ | **Implemented**, over fixed-length track segments; $\Delta$ required |
| $\bar\ell$, $y$, $z$ for one spherical site | **Implemented**, ICRU definitions |
| $f(y)$, $d(y)$, MKM | **Not implemented** |
| Radial dose tail | missing: δ-electrons are not transported |

## Complexity

- `VoxelDoseScorer.score`: three binary searches, $O(\log n_x + \log n_y + \log n_z)$
  per deposit. Memory $O(n_xn_yn_z)$.
- `RadialDoseScorer.score`: one projection plus a binary search,
  $O(\log n_r)$ per deposit.
- `SphericalSiteScorer.score`: $O(1)$ per deposit.
- Every scorer runs on every deposit, so the total scoring cost is
  (number of deposits) × (number of scorers).

## Tests

`tests/test_scorer.py`, `test_dose.py` (binning, manual-dose checks,
full-run energy conservation), `test_radiobiology.py` (known-value LET,
$\bar L_D \ge \bar L_T$, chord length, $y$, $z$, site scorer).
