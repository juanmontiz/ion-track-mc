# Cross-section and mean-energy data files

← [Back to the main README](../../README.md)

The physics of every collision comes from tables in
`ion_track_mc/data/cross_sections/`. This page describes how those files are
named, what is inside them, how they are read, and which ones the code
actually uses.

There are two kinds of table:

- **cross sections** $\sigma(E)$, which say how *likely* a process is at
  projectile energy $E$;
- **mean energies**, which say how much energy a process *costs* on average:
  - `Edepmeantot(E)`: mean energy left in the target molecule at the
    collision site;
  - `Eemeantot(E)`: mean kinetic energy of the electron the collision
    emits.

How the code uses them is explained in
[`physics/README_physics.md`](../physics/README_physics.md) (§2 for cross
sections, §4 for energy loss).

## File naming

```
{Projectile}_+_{Target}_{process}_{datatype}_{model}.dat
```

- **`Projectile`**: the projectile species *and* its charge state, written
  explicitly:
  - `H1+`: a proton, H⁺ ($q = 1$);
  - `H0`: a neutral hydrogen atom, H⁰ ($q = 0$);
  - `e-`: an electron.

  The charge state is part of the name because each charge state has its
  own set of possible processes.
- **`Target`**: the medium's molecule, here always `H2O`. The code checks
  that it matches the `Medium` in use (`physics/medium.py`,
  `check_medium_matches_filename`).
- **`process`**: one of `ionization`, `excitation`, `capture`, `loss`,
  `elastic`, the same names as the `Channel` enum in `physics/channels.py`.
  - Capture ($q: 1 \to 0$) files are named after the projectile *before*
    the capture, so they start with `H1+`.
  - Loss ($q: 0 \to 1$) files are likewise named after the projectile before
    the loss, so they start with `H0`.
- **`datatype`**: what the value column holds:
  - `TCS`: total cross section;
  - `Edepmeantot`: mean energy deposited at the collision site;
  - `Eemeantot`: mean energy of the emitted electron.
- **`model`**: the theory the values come from:
  - `CDW-EIS`: Continuum-Distorted-Wave Eikonal-Initial-State;
  - `Miller-Green`: a semi-empirical model;
  - `Classical`: classical screened-Coulomb scattering;
  - `PW`, `Olivero`, `BEB`: electron-impact models;
  - `constant`: a flat value, see [Loss energetics](#loss-energetics).

Example: `H1+_+_H2O_ionization_TCS_CDW-EIS.dat` is the total ionization
cross section of water by protons, computed with CDW-EIS.

The H⁺/H⁰ cross sections were computed for water **vapour**; the code uses
them for liquid water, as TILDA-V does.

## File structure

- Plain text, two whitespace-separated columns (tabs or spaces).
- Column 1 is the **collision energy**. Column 2 is the **value**.
- Blank lines and lines starting with `#` are skipped.
- An optional one-line text header is detected automatically. If the first
  data line does not parse as two numbers, it is skipped.
- Energies must be **strictly increasing**. The loader checks this, because
  interpolation needs it.

## Units and conversion on load

Inside the code, energies are in eV and cross sections in nm²
(see [`README_constants.md`](../README_constants.md)). Tables are converted
**once, when they are loaded**, so nothing downstream ever handles raw units.

| Column | Unit in the shipped files | Conversion |
|---|---|---|
| collision energy | eV | ×1 |
| `TCS` value | $10^{-16}\,\mathrm{cm^2}$ | ×$10^{-2}$ → nm² |
| `Edepmeantot`, `Eemeantot` value | eV | ×1 |

The loaders `load_cross_section_table` and `load_mean_energy_table`
(`physics/cross_section.py`) take the units as arguments. They default to
`constants.DEFAULT_ENERGY_UNIT` (`"eV"`) and
`constants.DEFAULT_CROSS_SECTION_UNIT` (`"1e-16cm2"`). Their conversion
factors are listed in `constants.ENERGY_UNIT_TO_EV` and
`constants.CROSS_SECTION_UNIT_TO_NM2`.

A `"keV/u"` energy unit (energy per nucleon) also exists, for proton tables
tabulated that way. For hydrogen it is simply keV, since $A = 1$. It must
never be used for electron tables: an electron has no nucleons.

The proton and electron grids look very different only because the physics
is different. The H⁺/H⁰ grids run from 8 keV to 300 MeV (beam energies);
the e⁻ grids run from 7.4 eV (13.1 eV for ionization) to 8.7 MeV (2.6 MeV
for elastic scattering): electron collision energies.
Both are in eV.

## Validity range: no extrapolation

A table is valid from its first to its last energy. Asking for a value
outside that range **raises `ValueError`**; the code never extrapolates. In
a track, this is why a proton stops being followed when its energy drops
below the lowest tabulated energy (8 keV). Its remaining energy is not
deposited.

## Inventory

Where it applies, the energy grid has 40 points from 8 keV to 300 MeV for
H⁺/H⁰ and 101 points for e⁻.

| File | Process | Projectile | Datatype | Used by the code? |
|---|---|---|---|---|
| `H1+_+_H2O_ionization_TCS_CDW-EIS.dat` | ionization | H⁺ | TCS | yes |
| `H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat` | ionization | H⁺ | Edepmeantot | yes |
| `H1+_+_H2O_ionization_Eemeantot_CDW-EIS.dat` | ionization | H⁺ | Eemeantot | yes |
| `H1+_+_H2O_excitation_TCS_Miller-Green.dat` | excitation | H⁺ | TCS | yes |
| `H1+_+_H2O_excitation_Edepmeantot_Miller-Green.dat` | excitation | H⁺ | Edepmeantot | yes |
| `H1+_+_H2O_capture_TCS_CDW-EIS.dat` | capture (σ₁₀) | H⁺ | TCS | yes |
| `H1+_+_H2O_capture_Edepmeantot_CDW-EIS.dat` | capture | H⁺ | Edepmeantot | yes |
| `H1+_+_H2O_elastic_TCS_Classical.dat` | elastic | H⁺ | TCS | yes |
| `H0_+_H2O_ionization_TCS_CDW-EIS.dat` | ionization | H⁰ | TCS | yes |
| `H0_+_H2O_ionization_Edepmeantot_CDW-EIS.dat` | ionization | H⁰ | Edepmeantot | yes |
| `H0_+_H2O_ionization_Eemeantot_CDW-EIS.dat` | ionization | H⁰ | Eemeantot | yes |
| `H0_+_H2O_excitation_TCS_Miller-Green.dat` | excitation | H⁰ | TCS | yes |
| `H0_+_H2O_excitation_Edepmeantot_Miller-Green.dat` | excitation | H⁰ | Edepmeantot | yes |
| `H0_+_H2O_loss_TCS_Miller-Green.dat` | loss (σ₀₁) | H⁰ | TCS | yes |
| `H0_+_H2O_loss_Edepmeantot_constant.dat` | loss | H⁰ | Edepmeantot | yes (flat 13.6 eV) |
| `H0_+_H2O_loss_Eemeantot_constant.dat` | loss | H⁰ | Eemeantot | no (flat 0 eV) |
| `H0_+_H2O_elastic_TCS_Classical.dat` | elastic | H⁰ | TCS | yes |
| `e-_+_H2O_elastic_TCS_PW.dat` | elastic | e⁻ | TCS | no |
| `e-_+_H2O_excitation_TCS_Olivero.dat` | excitation | e⁻ | TCS | no |
| `e-_+_H2O_excitation_Edepmeantot_Olivero.dat` | excitation | e⁻ | Edepmeantot | no |
| `e-_+_H2O_ionization_TCS_BEB.dat` | ionization | e⁻ | TCS | no |
| `e-_+_H2O_ionization_Edepmeantot_BEB.dat` | ionization | e⁻ | Edepmeantot | no |
| `e-_+_H2O_ionization_Eemeantot_BEB.dat` | ionization | e⁻ | Eemeantot | no |

**Which processes each charge state has.** H⁺ undergoes ionization,
excitation, capture and elastic scattering; H⁰ undergoes ionization,
excitation, loss and elastic scattering. Capture has no `Eemeantot` table
because no electron is emitted: it is captured.

**Elastic scattering** has only a cross section. Its energy transfer and
deflection angle are computed from collision kinematics
(`physics/energy_loss.py::elastic_scatter`; physics README §4.5), not from
tables.

**Electron (e⁻) tables are shipped but not used.** In this code a
δ-electron deposits all its energy where it is created
([`secondaries/README_secondaries.md`](../secondaries/README_secondaries.md)),
so electron cross sections are never needed. They are kept for a future
electron-transport model. One detail matters for that future work: the
excitation `Edepmeantot` grid (13.1 eV–8.7 MeV) is the ionization grid, not
the excitation cross-section grid (7.4 eV–8.7 MeV).

### Loss energetics

When H⁰ loses its electron, the code charges it the hydrogen ground-state
binding energy, 13.6 eV, at every energy. This is a lower bound: the atom
must at least pay that to free the electron, and what the freed electron
does afterwards is not tracked. To keep a single loader, this value is
still stored as an ordinary two-column table (`..._loss_Edepmeantot_constant.dat`,
flat 13.6 eV on the loss cross-section grid), together with a flat 0 eV
`Eemeantot` table.

## Using different data

To use other tables, replace the files, keeping the same naming and
two-column format; no code changes are needed. If the new files use
different units, pass `energy_unit=` / `cross_section_unit=` to the loaders,
or add the unit to the conversion dictionaries in `constants.py`.
