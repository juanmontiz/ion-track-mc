"""Stopping power of protons in water: model vs NIST PSTAR vs measurements (README figure).

Model: charge-equilibrium expectation of the code's collision rules
(analytic_equilibrium_stopping_power), which the Monte Carlo reproduces.
Measurements: data_comparison/stopping_power_water/ (references in its README).
Run from the repository root:  python docs/stopping_power_plot.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ion_track_mc.analysis.stopping_power import analytic_equilibrium_stopping_power
from ion_track_mc.physics.medium import WATER

from _tables import PSTAR_KEV, PSTAR_MEV_CM2_G, edep, eemean, xs

HERE = Path(__file__).parent
DATA = HERE.parent / "data_comparison" / "stopping_power_water"
MEASUREMENTS = {   # label: (file, marker)
    "Phillips 1953 (vapour)": ("Exp_SP_Phillips_1953.dat", "v"),
    "Reynolds et al. 1953 (vapour)": ("Exp_SP_Reynolds_1953.dat", "^"),
    "Mitterschiffthaler & Bauer 1990 (vapour)": ("Exp_SP_Mitterschiffthaler_1990.dat", "D"),
    "Baek et al. 2006 (vapour)": ("Exp_SP_Baek_2006.dat", "s"),
    "Siiskonen et al. 2011 (liquid)": ("Exp_SP_Siiskonen_2011.dat", "o"),
}

E_keV = np.geomspace(8, 1e5, 160)
model = [analytic_equilibrium_stopping_power(WATER, xs, E * 1e3, edep, eemean) for E in E_keV]

fig, ax = plt.subplots(figsize=(8, 4.8))
for label, (fname, marker) in MEASUREMENTS.items():
    E, S, err = np.loadtxt(DATA / fname, unpack=True)
    ax.errorbar(E, S, yerr=err, fmt=marker, ms=4, lw=0.8, capsize=0, alpha=0.85, label=label)
ax.plot(PSTAR_KEV, PSTAR_MEV_CM2_G * WATER.density_g_cm3 * 0.1, "k--", lw=1.4,
        label="NIST PSTAR (liquid water)")
ax.plot(E_keV, model, color="tab:blue", lw=2.6, label="ion-track-mc (this code)", zorder=5)
ax.set_xscale("log")
ax.set_xlim(1, 1e5)
ax.set_ylim(0, 115)
ax.set_xlabel("Proton energy [keV]")
ax.set_ylabel("Stopping power [keV/µm]")
ax.set_title("Proton stopping power in water")
ax.grid(which="major", alpha=0.3)
ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
fig.tight_layout()
fig.savefig(HERE / "stopping_power.png", dpi=130)
