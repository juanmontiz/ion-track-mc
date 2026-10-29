"""Bragg curve for the README: depth-dose curve of 50 protons of 1 MeV in water.

Run from the repository root:  python docs/make_bragg_figure.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ion_track_mc.physics.medium import WATER
from ion_track_mc.scoring.dose import make_bragg_curve_scorer
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

from _tables import edep, eemean, xs

N_PROTONS = 50
scorer = make_bragg_curve_scorer(WATER, np.linspace(0, 30e3, 121), transverse_half_width_nm=1e6)
engine = TrackEngine(WATER, xs, edep_tables=edep, eemean_tables=eemean, scorers=[scorer], seed=1)
for _ in range(N_PROTONS):
    engine.transport(Particle(kind="proton", position=np.zeros(3),
                              direction=np.array([0.0, 0.0, 1.0]), energy=1e6, q=1))
    engine.deposits.clear()

z = 0.5 * (scorer.z_edges[1:] + scorer.z_edges[:-1]) / 1e3
dose = scorer.dose_Gy()[0, 0, :]
fig, ax = plt.subplots(figsize=(7, 3.6))
ax.step(z, dose / dose.max(), where="mid", color="tab:blue", lw=1.8)
ax.fill_between(z, dose / dose.max(), step="mid", color="tab:blue", alpha=0.12)
ax.set_xlim(0, 30)
ax.set_ylim(0, 1.08)
ax.set_xlabel("Depth [µm]")
ax.set_ylabel("Dose [relative]")
ax.set_title(f"Bragg curve: {N_PROTONS} protons of 1 MeV in water")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(Path(__file__).with_name("bragg_curve.png"), dpi=130)
