# `ion_track_mc.secondaries` — δ-electron transport plug-in

← [Back to the main README](../../README.md)

Every ionization creates a **δ-electron** with kinetic energy
$W = \bar E_e(E)$ and a sampled direction (see
[`physics/README_physics.md`](../physics/README_physics.md) §4.1, §5).
`TrackEngine` pushes it on the secondary stack. When the electron is popped,
the engine does **not** transport it itself. It hands the electron to an
injected `ElectronTransport` object, and this package defines that
interface and its current implementation.

| File | Role |
|---|---|
| `electron_transport.py` | `ElectronTransport` protocol, `ElectronDeposit`, `LocalDepositElectronTransport` |

---

## 1. Why a separate plug-in?

Electron transport is physically different from ion transport. Electrons
scatter elastically at large angles, change direction many times along
their path, and need their own cross sections. Folding them into the ion
engine's `CrossSectionSet` would mix two kinds of physics in one loop.
Keeping them behind an interface means a better electron model only
touches this package.

## 2. The interface

```python
class ElectronTransport(Protocol):
    def transport(self, electron: Particle) -> list[ElectronDeposit]: ...
```

The method receives one electron and returns **every** deposit it makes,
each as `ElectronDeposit(position, energy)`. The engine records each one
(`Deposit(kind="electron", channel=None, q=None)`) and notifies all
scorers. The contract is energy conservation:

$$
\sum_{j}\varepsilon_j \;=\; W
$$

for each electron. All its energy must end up in the returned deposits,
unless an implementation documents an escape or cutoff term.

## 3. Current implementation: fidelity level (a), local deposit

`LocalDepositElectronTransport` deposits the electron's **entire** kinetic
energy at its creation point, in a single deposit:

$$
\varepsilon = W, \qquad \vec r_{\text{deposit}} = \vec r_{\text{ionization}} .
$$

No transport is simulated. The `e_cut` argument exists only to match the
signature of future implementations and is ignored.

The three fidelity levels considered:

| Level | What it does | Needs | Status |
|---|---|---|---|
| (a) | deposit $W$ locally at creation point | nothing | **implemented, default** |
| (b) | deposit along a straight path of length given by a range–energy relation (CSDA) | electron range data | not implemented |
| (c) | full event-by-event electron track structure | electron elastic + inelastic cross sections in water | not implemented |

Some electron tables do sit in `data/cross_sections/` (the `e-_+_H2O_*`
files). They are **not wired** into any transport model.

### What (a) gets right and wrong

- **Right:** total energy. Integral quantities such as the depth–dose
  (Bragg) curve at scales much larger than δ-electron ranges, the stopping
  power, and $\bar L_T$ are unaffected, because the energy is deposited
  either way.
- **Wrong:** where the energy goes, at nanometre scale. The radial dose
  profile $D(r)$ loses its δ-ray tail. Microdosimetric $\varepsilon$ in
  small sites is biased: energy that should leak out of, or into, a site is
  kept at the ionization point. Every deposit an electron makes is a single
  point.

### Why (a)

Moving to (b) or (c) would make every track **slower**: each δ-electron
would run its own collision loop. The code is meant to run on students'
laptops, so (a) is kept. Its consequence is that δ-electron energy never
leaves the track core: dose very close to the track axis is overestimated
and dose far from it is missing, while track-integrated quantities (Bragg
curve, LET) are much less affected. Electron cross sections for water
(`e-_+_H2O_*` tables) are shipped but not used.

## Adding a new electron model

1. Write a class with `transport(self, electron) -> list[ElectronDeposit]`.
   No inheritance is needed; it is a `typing.Protocol`.
2. Make sure its deposits conserve $W$.
3. Inject it with `TrackEngine(..., electron_transport=MyModel())`, or the
   same argument of `run_beam`.

## Complexity

(a) is $O(1)$ per electron: one deposit. (c) would cost
$O(\text{collisions per electron})$, which is why this code does not implement it.

## Tests

`tests/test_electron_transport.py`: full-energy local deposit, default
selection in the engine, the stack draining with electron energy accounted
for, and use of an injected custom model.
