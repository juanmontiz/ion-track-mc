# `ion_track_mc.transport` — particle state, sampling, the collision loop

← [Back to the main README](../../README.md)

This package **moves particles**. It samples how far a projectile travels
before its next collision, picks which interaction happens there, applies
the physics rules from [`physics/`](../physics/README_physics.md), and keeps
following the primary and every secondary it creates until all of their
energy is accounted for.

The sampling formulas are those of TILDA-V (Quinto et al. 2017, Eqs. 1–3).

| File | Role |
|---|---|
| `particle.py` | `Particle`: complete mutable state of one particle |
| `sampler.py` | `sample_step` (free path) and `select_channel` (which interaction) |
| `engine.py` | `TrackEngine`: the event-by-event loop, secondary stack, deposits, scorer notification |

---

## 1. Particle state (`particle.py`)

A `Particle` carries `kind` (`"proton"` or `"electron"`), `position`
$\vec r$ [nm], unit `direction` $\hat d$, kinetic `energy` $E$ [eV], charge
state `q`, `alive`, and `parent_id`.

**`q` is mutable on purpose.** Capture sets $q: 1\to0$ and loss sets
$q: 0\to1$. Charge state is part of the particle's evolving state, not a
fixed property of its species. This is the central physics point of the
code: the cross sections, and therefore the mean free
path, change every time $q$ flips.

## 2. Transport model: event-by-event, no continuous loss

The engine uses **event-by-event** (detailed) transport, not condensed
history. Every collision is simulated individually.
Between collisions the projectile flies **in a straight line with no energy
loss**:

$$
\vec r \;\leftarrow\; \vec r + s\,\hat d .
$$

All energy loss is discrete, at collision points. All direction change
comes from elastic collisions (and, for the δ-electron, from its emission
angle). This is the TILDA-V approach, and it is what gives the code
nanometre-scale track structure.

## 3. Sampling the distance to the next collision (`sampler.sample_step`)

For charge state $q$ at energy $E$ in a medium of molecular number density
$N$, the **macroscopic total cross section** (inverse mean free path) is

$$
\Sigma_{\text{tot}}^{(q)}(E) \;=\; N\sum_c \sigma_c^{(q)}(E)
\;=\; \frac{1}{\lambda^{(q)}(E)} .
$$

The sum runs over every open channel, **elastic included**: elastic
collisions are real collisions and shorten the mean free path.

Collisions along a straight path form a Poisson process with rate
$\Sigma_{\text{tot}}$ per unit length. The probability of travelling $s$
without a collision is therefore $e^{-\Sigma_{\text{tot}} s}$, and the step
length has cumulative distribution

$$
F(s) \;=\; 1 - e^{-\Sigma_{\text{tot}}\,s}.
$$

**Inverse-transform sampling.** Set $F(s) = 1 - \xi$ with
$\xi\sim U(0,1)$. This is valid because $1-\xi$ is also uniform. Solving
for $s$:

$$
\boxed{\; s \;=\; -\frac{\ln\xi}{\Sigma_{\text{tot}}^{(q)}(E)} \;=\; -\lambda^{(q)}(E)\,\ln\xi \;}
$$

This is the TILDA-V paper's Eq. 2. $E$ and $q$ do not change between
collisions, so $\Sigma_{\text{tot}}$ is constant along the step, which is
exactly what makes this closed form valid.

## 4. Choosing the interaction (`sampler.select_channel`)

Given that a collision happened, channel $c$ is chosen with probability

$$
p_c \;=\; \frac{\sigma_c^{(q)}(E)}{\sum_{c'}\sigma_{c'}^{(q)}(E)} .
$$

This uses **discrete inverse-transform sampling**. Build the cumulative
distribution $P_c = \sum_{c'\le c} p_{c'}$, draw $\xi'\sim U(0,1)$, and pick
the first $c$ with $\xi' < P_c$:

$$
P_{c-1} \;\le\; \xi' \;<\; P_c .
$$

In code: `np.cumsum` followed by `np.searchsorted`. The channel order is
fixed (alphabetical by enum value, from `CrossSectionSet.channels`), which
makes runs with the same seed reproducible. The order does not affect the
probabilities. This is TILDA-V's Eq. 3.

`sample_step` and `select_channel` are kept as **separate, standalone
functions**. Each one maps one-to-one onto one formula above and is tested
on its own (`tests/test_sampler.py`). Fusing them would save repeated
cross-section lookups. They are kept separate anyway, so the
code reads formula by formula; see "Complexity" below.

## 5. The collision loop (`engine.py`)

### 5.1 One primary per call, secondaries depth-first

`TrackEngine.transport(primary)` follows **one** primary and its **entire**
cascade before returning. Secondaries go on a LIFO stack (a plain Python
list, `self.stack`), so the cascade is explored depth-first:

```
push primary
while stack not empty:
    p = pop()
    if p is a proton:   run the ion loop (§5.2)
    if p is an electron: hand it to the injected ElectronTransport
```

Electrons are not transported here. They are delegated to an injected
`ElectronTransport` plug-in (see
[`secondaries/README_secondaries.md`](../secondaries/README_secondaries.md)),
because electron physics is different enough to live in its own module.

To simulate a beam, call `transport()` once per primary. Deposits and
scorers **accumulate across calls**; `analysis.beam.run_beam` is the thin
wrapper that does this loop for you.

### 5.2 Ion loop (`_transport_proton`)

```
while p.alive and E > e_cut:
    if E outside xs.energy_range(q): stop          # no extrapolation
    s       = sample_step(E, q, N, xs, rng)          # §3
    r      += s * d                                  # straight-line flight
    channel = select_channel(E, q, xs, rng)          # §4
    _interact(p, channel)                            # apply physics, may change E, q, d
    collision_log.append((s, E_before - E))          # for LET (scoring/)
```

The projectile stops when its energy drops below `e_cut` (default
$10\,\mathrm{eV}$), or below the tabulated validity floor (currently
$8\,\mathrm{keV}$), whichever it reaches first. With the current data the
validity floor is always reached first. The energy remaining at that point
is **not** deposited. It is visible as the returned particle's `energy`
(`BeamResult.E_final_eV` in `analysis/`).

### 5.3 Applying one interaction (`_interact`)

| Channel | Energy removed from projectile $\Delta E$ | Other effects |
|---|---|---|
| excitation | $\bar E_{\text{dep}}(E)$, deposited locally | — |
| ionization | $\bar E_{\text{dep}}(E)$ deposited locally **plus** $W = \bar E_e(E)$ given to a new δ-electron | electron pushed on stack, direction from `kinematics` |
| capture | $\bar E_{\text{dep}}(E) + E/M_P - 13.6\,\mathrm{eV}$, deposited locally | $q \leftarrow 0$ |
| loss | $13.6\,\mathrm{eV}$, deposited locally | $q \leftarrow 1$ |
| elastic | $E_{\text{elastic}}(\theta)$, deposited locally | $\hat d \leftarrow$ rotated by $(\theta,\phi)$ |

Formulas and derivations: [`physics/README_physics.md`](../physics/README_physics.md) §4–5.

**Energy conservation guard.** Every deposit is clamped,
$\Delta E \leftarrow \min(\Delta E, E)$, so a particle never loses more
energy than it carries. If $E$ reaches $0$, the particle is marked dead.

### 5.4 Outputs

Each local deposit is recorded in three ways:

1. A `Deposit(position, energy, channel, kind, q)` appended to
   `engine.deposits`. This is internal bookkeeping, used by `plot_track`.
2. A call to `scorer.score(position, energy)` on **every injected
   `Scorer`**. This is the observer pattern; see
   [`scoring/README_scoring.md`](../scoring/README_scoring.md).
3. For the **primary ion only**, a `(step_length, delta_E)` pair in
   `engine.collision_log`, reset at the start of each `transport()` call.
   This is the input for LET. Electron deposits are deliberately left out,
   because LET describes the primary's own track.

Note that `delta_E` in `collision_log` is the projectile's **total** energy
loss in that collision, including the δ-electron's $W$. LET computed from
it is therefore the *unrestricted* LET.

### 5.5 Dependency injection

`TrackEngine` is constructed with a `Medium`, a `CrossSectionSet`, the
mean-energy tables, an optional `ElectronTransport`, and optional
`Scorer`s. It only knows these through their interfaces. A new medium,
cross-section model, electron model, or scorer is injected; the engine is
never edited.

### 5.6 Randomness

There is one `np.random.default_rng(seed)` per engine, and all sampling
inside the engine draws from it. With the same seed and the same inputs,
the track is identical.

---

## Fidelity level

| Piece | Status |
|---|---|
| Event-by-event transport, straight-line steps, no continuous loss | **Implemented** |
| Free-path sampling, channel selection | **Implemented**; match TILDA-V Eqs. 2–3 ✅ |
| Dynamic charge state $q$ | **Implemented** |
| Projectile angular deflection | via elastic channel only; not yet validated against range-straggling data |
| Tracking cutoff | tabulated range floor (8 keV) in practice; residual energy not deposited |
| Electron transport | delegated to `secondaries/`; local-deposit model |

## Complexity and performance

The cost per primary is proportional to the number of collisions, which
grows with initial energy, since a faster projectile needs more collisions
before it reaches the floor. Per ion collision, with $K$ open channels:

- `energy_range(q)` (once), `sample_step` and `select_channel` each call
  `channels(q)`. That is cheap, because the list is cached per $q$.
- `sample_step` and `select_channel` each evaluate all $K$ cross sections at
  the **same** $E$, so there are $2K$ `sigma` calls where $K$ would suffice.
  This duplication is the price of keeping the two formulas as separate,
  readable functions. It is accepted on purpose: in this codebase,
  simplicity beats performance.
- Each ionization or elastic event builds an orthonormal frame
  (`kinematics._orthonormal_basis`, with the cross products written out by
  hand).

Measured on a single 500 keV/u track (seed 0, 24 410 deposits): 3.80 s
before the two readability-neutral fixes (hand-written cross product and
cached `channels(q)`), 2.39 s after, with an identical seeded
result. Most of the remaining time is cross-section lookups and the
Python-level collision loop itself. Further speed-ups were deliberately
not made, because they would cost readability.

Memory: `engine.deposits` grows by one entry per deposit and is **never
cleared** between `transport()` calls. For long beam runs, rely on scorers
instead of keeping every deposit.

## Tests

`tests/test_particle.py`; `test_sampler.py` (Kolmogorov–Smirnov test of
sampled step lengths against the exponential law, empirical mean vs.
$\lambda(E)$, channel frequencies vs. $p_c$); `test_engine.py` (each
`_interact` branch, deposit clamping, exact energy accounting over a full
track).
