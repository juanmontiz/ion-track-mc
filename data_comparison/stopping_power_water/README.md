# Proton stopping power in water: data for comparison

Published stopping powers of protons in water, used in
`notebooks/01_full_pipeline_walkthrough.ipynb` (§6) to compare with the
stopping power that emerges from the simulation.

**Format.** Plain text, whitespace-separated columns:

| Column | Quantity | Unit |
|---|---|---|
| 1 | proton kinetic energy | keV (= keV/u for protons) |
| 2 | stopping power | keV/µm (= eV/nm) |
| 3 | absolute uncertainty (measurements only) | keV/µm |

## Measurements

| File | Phase | Energy range (keV) | Reference |
|---|---|---|---|
| `Exp_SP_Phillips_1953.dat` | vapour | 10.6–73.4 | J. A. Phillips, «The energy loss of low energy protons in some gases», *Phys. Rev.* **90**, 532 (1953). doi:10.1103/PhysRev.90.532 |
| `Exp_SP_Reynolds_1953.dat` | vapour | 40–550 | H. K. Reynolds, D. N. F. Dunbar, W. A. Wenzel and W. Whaling, «The stopping cross section of gases for protons, 30–600 keV», *Phys. Rev.* **92**, 742 (1953). doi:10.1103/PhysRev.92.742 |
| `Exp_SP_Mitterschiffthaler_1990.dat` | vapour | 26.4–360.2 | C. Mitterschiffthaler and P. Bauer, «Stopping cross section of water vapor for hydrogen ions», *Nucl. Instrum. Methods B* **48**, 58 (1990). doi:10.1016/0168-583X(90)90073-4 |
| `Exp_SP_Baek_2006.dat` | vapour | 1.1–90.8 | W. Y. Baek, B. Grosswendt and G. Willems, «Ionization ranges of protons in water vapour in the energy range 1–100 keV», *Radiat. Prot. Dosim.* **122**, 32 (2006). doi:10.1093/rpd/ncl514 |
| `Exp_SP_Siiskonen_2011.dat` | liquid | 4758–15243 | T. Siiskonen *et al.*, «Energy loss measurement of protons in liquid water», *Phys. Med. Biol.* **56**, 2367 (2011). doi:10.1088/0031-9155/56/8/003 |

The vapour-phase values are given in keV/µm, that is, scaled to a density
of 1 g/cm³.

The Baek *et al.* values are stopping powers derived from measured
ionization ranges with semi-empirical detour factors (see the abstract).

*Note (2026-10-29).* In the development copy of this code the Reynolds and
Baek files had each other's names. They were swapped back after checking
the energy ranges stated in the abstracts (Reynolds: 30–600 keV; Baek:
1–100 keV).

## Model (not a measurement)

| File | Energy range (keV) | Reference |
|---|---|---|
| `SP_2017_Quinto_EJPD.dat` | 10–100000 | M. A. Quinto *et al.*, «Monte Carlo simulation of proton track structure in biological matter», *Eur. Phys. J. D* **71**, 130 (2017). doi:10.1140/epjd/e2017-70709-6. The TILDA-V model whose cross sections `ion_track_mc` uses. |
