"""Stopping power of the model vs NIST PSTAR at a few energies (README table).

The model value is the charge-equilibrium expectation of the code's own
collision rules (analytic_equilibrium_stopping_power); the Monte Carlo
reproduces it (tests/test_stopping_power.py).

Run from the repository root:  python docs/stopping_power_table.py
"""
from ion_track_mc.analysis.stopping_power import analytic_equilibrium_stopping_power
from ion_track_mc.physics.medium import WATER

from _tables import PSTAR_KEV, PSTAR_MEV_CM2_G, edep, eemean, xs

ENERGIES_KEV = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 100000]

print("| Energy (keV) | Model (keV/µm) | PSTAR (keV/µm) | Model / PSTAR |")
print("|---:|---:|---:|---:|")
for E_keV, s in zip(PSTAR_KEV, PSTAR_MEV_CM2_G):
    if E_keV not in ENERGIES_KEV:
        continue
    model = analytic_equilibrium_stopping_power(WATER, xs, E_keV * 1e3, edep, eemean)
    pstar = s * WATER.density_g_cm3 * 0.1
    print(f"| {E_keV:g} | {model:.1f} | {pstar:.1f} | {model / pstar:.2f} |")
