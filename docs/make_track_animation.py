"""Animated GIF for the README: a proton track built collision by collision.

The first 25 nm of a 15 keV proton track in water (same track as the README
figure, seed 1). Collisions appear in the order the code samples them, coloured
by process with the colours of plot_track; the path is blue while the projectile
is H+ and orange while it is H0.
Run from the repository root:  python docs/make_track_animation.py
"""
import copy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D

from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.medium import WATER
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

from _tables import edep, eemean, xs

Z_MAX = 25.0   # nm shown
COLORS = {Channel.EXCITATION: "tab:green", Channel.IONIZATION: "tab:red",
          Channel.CAPTURE: "tab:purple", Channel.LOSS: "tab:brown", Channel.ELASTIC: "tab:cyan"}
CHARGE = {1: "tab:blue", 0: "tab:orange"}

engine = TrackEngine(WATER, xs, edep_tables=edep, eemean_tables=eemean, seed=1)
p = Particle(kind="proton", position=np.zeros(3), direction=np.array([0.0, 0.0, 1.0]),
             energy=15e3, q=1)
start = copy.deepcopy(p)
engine.transport(p)
dep = [d for d in engine.deposits if d.kind == "proton" and d.position[2] <= Z_MAX]
pos = np.array([start.position] + [d.position for d in dep])
# A deposit records the charge state the projectile had when it arrived,
# so it is also the state on the segment that ends there.
q_seg = [d.q for d in dep]
# Energy after each collision (collision_log has one entry per primary collision, in order).
E_after = 15.0 - np.cumsum([c[1] for c in engine.collision_log])[:len(dep)] / 1e3

fig, ax = plt.subplots(figsize=(9, 2.9), dpi=100)
ax.set_xlim(0, Z_MAX)
ax.set_ylim(-0.25, 0.5)
ax.set_xlabel("Depth z [nm]")
ax.set_ylabel("x [nm]")
ax.set_title("A 15 keV proton in water, one collision at a time")
handles = [Line2D([], [], color=c, lw=2, label=l) for l, c in (("H⁺", CHARGE[1]), ("H⁰", CHARGE[0]))]
handles += [Line2D([], [], ls="", marker="o", color=c, label=ch.value) for ch, c in COLORS.items()]
ax.legend(handles=handles, loc="upper left", fontsize=8.5, ncol=7, frameon=False, columnspacing=1.2, handletextpad=0.4)
counter = ax.text(0.99, 0.04, "", transform=ax.transAxes, ha="right", fontsize=9)
fig.tight_layout()

N = len(dep)


def frame(k):
    """Draw segment k (start → collision k) and the marker of collision k."""
    if k < N:
        a, b = pos[k], pos[k + 1]
        ax.plot([a[2], b[2]], [a[0], b[0]], color=CHARGE[q_seg[k]], lw=1.6, zorder=1)
        ax.scatter(b[2], b[0], s=36, color=COLORS[dep[k].channel], zorder=3, ec="white", lw=0.4)
        counter.set_text(f"collision {k + 1}   E = {E_after[k]:.2f} keV")
    return []


OUT = Path(__file__).with_name("track_animation.gif")
FuncAnimation(fig, frame, frames=N, interval=110, blit=False).save(OUT, writer=PillowWriter(fps=9))

# Hold the finished track for 2.5 s before the loop restarts.
from PIL import Image, ImageSequence  # noqa: E402
gif = Image.open(OUT)
frames = [f.copy() for f in ImageSequence.Iterator(gif)]
durations = [110] * (len(frames) - 1) + [2500]
frames[0].save(OUT, save_all=True, append_images=frames[1:], duration=durations, loop=0)
print(N, "collisions,", len(frames), "frames")
