# `ion_track_mc.constants` — unit policy and physical constants

← [Back to the main README](../README.md)

`constants.py` is a single module, not a package. It fixes the **internal
units** for the whole code and holds every physical constant and
raw-data unit-conversion factor. Nothing else in the codebase hard-codes a
unit conversion.

Raw data file format: [`data/README_data_format.md`](data/README_data_format.md).

---

## 1. Internal unit policy

| Quantity | Internal unit | Why |
|---|---|---|
| energy $E$ | $\mathrm{eV}$ | natural scale of single collisions |
| length $s$, $\vec r$ | $\mathrm{nm}$ | track structure is resolved at nm scale |
| cross section $\sigma$ | $\mathrm{nm}^2$ | so that $N\sigma$ is directly in $\mathrm{nm}^{-1}$ |
| number density $N$ | $\mathrm{nm}^{-3}$ | same reason |

With these choices the core transport formula needs no conversion factors:

$$
\Sigma_{\text{tot}}\;[\mathrm{nm}^{-1}] \;=\; N\;[\mathrm{nm}^{-3}]\;\times\;\sigma_{\text{tot}}\;[\mathrm{nm}^2],
\qquad s = -\ln\xi/\Sigma_{\text{tot}}\;[\mathrm{nm}].
$$

Stopping power and LET then come out in $\mathrm{eV/nm}$, which equals
$\mathrm{keV/\mu m}$ exactly. Conversion to SI happens only at the
dose-scoring boundary, through `EV_TO_JOULE` and `CM3_TO_NM3`.

## 2. Convert once, at load time

Raw tables arrive in whatever units their source used. The loaders in
`physics/cross_section.py` convert **once, when a table is read**, by
looking up a unit *name*:

```python
ENERGY_UNIT_TO_EV        = {"eV": 1.0, "keV": 1e3, "keV/u": 1e3 * PROTON_NUCLEON_COUNT}
CROSS_SECTION_UNIT_TO_NM2 = {"cm2": 1e14, "1e-16cm2": 1e-2}
DEFAULT_ENERGY_UNIT        = "eV"
DEFAULT_CROSS_SECTION_UNIT = "1e-16cm2"
```

To support a data set tabulated differently, pass (or default to) a
different unit name. The loader code never changes. The relevant identities
are $1\,\mathrm{cm} = 10^7\,\mathrm{nm}$, so $1\,\mathrm{cm^2} = 10^{14}\,\mathrm{nm^2}$
and $10^{-16}\,\mathrm{cm^2} = 10^{-2}\,\mathrm{nm^2}$; likewise
$1\,\mathrm{cm^3} = 10^{21}\,\mathrm{nm^3}$.

Mean-energy tables (`Edepmeantot`, `Eemeantot`) have an energy in their
value column too, so both columns use `ENERGY_UNIT_TO_EV`.

## 3. Energy per nucleon: keV/u means divide by $A$, not by the atomic mass

Accelerator and radiobiology convention quotes ion energies **per
nucleon**:

$$
E\;[\mathrm{keV/u}] \;=\; \frac{E_{\text{kin}}\;[\mathrm{keV}]}{A},
$$

where $A$ is the integer **mass number** (nucleon count) of the
projectile. For hydrogen, $A = 1$, which is `PROTON_NUCLEON_COUNT`. So
$1\,\mathrm{keV/u} = 10^3\,\mathrm{eV}$ of kinetic energy, exactly.

It is **not** the isotope's atomic mass in unified mass units:
$m(^1\mathrm{H}) = 1.0078250319\,\mathrm{u}$. That number differs from 1
because $1\,\mathrm{u}$ is defined as $1/12$ of a $^{12}\mathrm{C}$ atom
(whose nucleons are bound, so lighter per nucleon), and because the
$^1\mathrm{H}$ atom includes its electron. It is only *numerically close*
to $A$: using the atomic mass as the conversion factor would be a unit
error of about $0.78\%$. The code therefore uses
`PROTON_NUCLEON_COUNT = 1` and the unit string "keV/u".

**Electrons have no nucleon count**, so `e-` tables are in plain
$\mathrm{eV}$, never $\mathrm{keV/u}$.

## 4. Physical constants

| Name | Value | Used for |
|---|---|---|
| `EV_PER_HARTREE` | $27.211386245988\,\mathrm{eV}$ | atomic-unit energies → eV (defined, not currently used) |
| `AVOGADRO` | $6.02214076\times10^{23}\,\mathrm{mol^{-1}}$ | $N = \rho N_A/M_{\text{mol}}$ |
| `EV_TO_JOULE` | $1.602176634\times10^{-19}\,\mathrm{J/eV}$ (exact, SI 2019) | dose in Gy |
| `PROTON_TO_ELECTRON_MASS_RATIO` | $M_P = m_p/m_e = 1836.15267343$ | capture drag $E/M_P$, kinematics, elastic |
| `ELECTRON_TO_PROTON_MASS_RATIO` | $1/M_P$ | same, written as $m_e/m_p$ |
| `H_BINDING_ENERGY_EV` | $I_H = 13.6\,\mathrm{eV}$ | capture and loss energy balance |
| `ELECTRON_REST_ENERGY_EV` | $m_ec^2 = 510\,999.0645\,\mathrm{eV}$ | elastic: projectile momentum |
| `FINE_STRUCTURE_CONSTANT` | $\alpha = 1/137$ | elastic: screening angle $\theta_{\text{cut}}$ |
| `N_WATER_NM3` | $33.4\,\mathrm{nm^{-3}}$ | reference value only; `Medium.N` computes $33.43$ from density and molar mass |

Note that $\alpha$ is the rounded $1/137$, as in the source formula, not
CODATA $1/137.036$. The difference is $0.03\%$ in $\theta_{\text{cut}}$.

Where each constant enters the physics:
[`physics/README_physics.md`](physics/README_physics.md) §3–5.

## Tests

`tests/test_constants.py`.
