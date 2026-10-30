"""Figures for the tutorial slides (docs/tutorial/tutorial.md).

Every figure is computed with ion_track_mc itself, using the same tables as the
README figures (docs/_tables.py). Course palette, large fonts, SVG output.
Run from the repository root:
  .venv/bin/python docs/tutorial/img/make_figures.py        # English labels, <name>.svg
  .venv/bin/python docs/tutorial/img/make_figures.py es     # Spanish labels, <name>_es.svg
"""
import copy
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "docs"))

from _tables import PSTAR_KEV, PSTAR_MEV_CM2_G, edep, eemean, xs  # noqa: E402

from ion_track_mc.analysis.beam import run_beam  # noqa: E402
from ion_track_mc.analysis.stopping_power import (  # noqa: E402
    analytic_equilibrium_stopping_power, equilibrium_charge_fractions)
from ion_track_mc.physics.channels import Channel  # noqa: E402
from ion_track_mc.physics.medium import WATER  # noqa: E402
from ion_track_mc.scoring.dose import RadialDoseScorer, make_bragg_curve_scorer  # noqa: E402
from ion_track_mc.transport import sampler  # noqa: E402
from ion_track_mc.transport.engine import TrackEngine  # noqa: E402
from ion_track_mc.transport.particle import Particle  # noqa: E402

NAVY, NAVY2, AMBER, CYAN, RED, MUTED = "#0b1d3a", "#16305c", "#f5a524", "#1597b5", "#c8364b", "#5b6475"
GREEN, VIOLET, BROWN = "#2e8b57", "#7b4fa0", "#8c5a2b"
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 22, "mathtext.fontset": "cm",
                     "axes.linewidth": 1.2, "xtick.major.width": 1.2, "ytick.major.width": 1.2,
                     "xtick.labelsize": 20, "ytick.labelsize": 20})
COLOR = {Channel.IONIZATION: RED, Channel.EXCITATION: GREEN, Channel.CAPTURE: VIOLET,
         Channel.LOSS: BROWN, Channel.ELASTIC: CYAN}
LANG = sys.argv[1] if len(sys.argv) > 1 else "en"


def T(en, es):
    """The label in the language of this run."""
    return es if LANG == "es" else en


def num(x, fmt):
    """A number as text, with a decimal comma in Spanish."""
    t = format(x, fmt)
    return t.replace(".", ",") if LANG == "es" else t


def decimal_comma(axis):
    """Decimal comma on a linear axis (Spanish only)."""
    if LANG == "es":
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",")))


NAME = {Channel.IONIZATION: T("ionization", "ionización"), Channel.EXCITATION: T("excitation", "excitación"),
        Channel.CAPTURE: T("capture", "captura"), Channel.LOSS: T("loss", "pérdida"),
        Channel.ELASTIC: T("elastic", "elástica")}


def save(fig, name):
    if LANG != "en":
        name = name.replace(".svg", f"_{LANG}.svg")
    fig.savefig(HERE / name, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print("wrote", name)


# ------------------------------------------------ cross sections, H+ and H0
fig, axs = plt.subplots(1, 2, figsize=(15, 4.6), sharey=True)
for ax, q, title in ((axs[0], 1, "H$^+$"), (axs[1], 0, "H$^0$")):
    e_min, e_max = xs.energy_range(q)
    E = np.geomspace(max(e_min, 8e3), min(e_max, 1e7), 300)
    for ch in [c for c in COLOR if c in xs.channels(q)]:
        ax.plot(E / 1e3, [xs.sigma(ch, e, q) for e in E], color=COLOR[ch], lw=3, label=NAME[ch])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(T("Projectile energy (keV)", "Energía del proyectil (keV)"))
    ax.set_title(title, color=NAVY)
axs[0].set_ylabel(T("Cross section $\\sigma$ (nm$^2$)", "Sección eficaz $\\sigma$ (nm$^2$)"))
axs[0].set_ylim(1e-6, 1)
fig.tight_layout(rect=(0, 0, 0.84, 1))
fig.legend(handles=[Line2D([], [], color=COLOR[c], lw=3, label=NAME[c]) for c in COLOR],
           frameon=False, fontsize=21, loc="center left", bbox_to_anchor=(0.84, 0.5))
save(fig, "cross_sections.svg")

# ------------------------------------------------ sampled steps vs exponential law
rng = np.random.default_rng(0)
E0 = 1e6                                             # 1 MeV, H+
lam = 1.0 / (WATER.N * xs.sigma_total(E0, 1))        # mean free path, nm
s = np.array([sampler.sample_step(E0, 1, WATER.N, xs, rng) for _ in range(20000)])
fig, ax = plt.subplots(figsize=(8.5, 6))
ax.hist(s, bins=60, range=(0, 6 * lam), density=True, color=CYAN, alpha=0.55, label=T("20 000 sampled steps", "20 000 pasos muestreados"))
x = np.linspace(0, 6 * lam, 200)
ax.plot(x, np.exp(-x / lam) / lam, color=NAVY, lw=3, label="$p(s)=e^{-s/\\lambda}/\\lambda$")
ax.axvline(lam, color=RED, lw=2, ls="--")
ax.text(lam * 1.08, 0.75 / lam, f"$\\lambda$ = {num(lam, '.2f')} nm", color=RED, fontsize=21)
ax.set_xlabel(T("Step length $s$ (nm)", "Longitud del paso $s$ (nm)"))
ax.set_ylabel(T("Probability density (nm$^{-1}$)", "Densidad de probabilidad (nm$^{-1}$)"))
decimal_comma(ax.yaxis)
ax.set_xlim(0, 6 * lam)
ax.legend(frameon=False, fontsize=19)
ax.text(0.97, 0.45, "H$^+$, 1 MeV", transform=ax.transAxes, ha="right", color=MUTED, fontsize=20)
fig.tight_layout()
save(fig, "step_distribution.svg")
print(f"  mean free path at 1 MeV: {lam:.3f} nm; sample mean {s.mean():.3f} nm")
for e in (10e3, 100e3, 1e6, 10e6):
    print(f"  lambda({e/1e3:g} keV) = {1/(WATER.N*xs.sigma_total(e,1)):.3f} nm")

# ------------------------------------------------ equilibrium charge fractions
E = np.geomspace(8e3, 1e6, 200)
f = np.array([equilibrium_charge_fractions(xs, e) for e in E])
fig, ax = plt.subplots(figsize=(8.5, 6))
ax.plot(E / 1e3, f[:, 0], color=NAVY, lw=3)
ax.plot(E / 1e3, f[:, 1], color=AMBER, lw=3)
ax.set_xscale("log")
ax.set_ylim(0, 1.02)
ax.set_xlabel(T("Projectile energy (keV)", "Energía del proyectil (keV)"))
ax.set_ylabel(T("Equilibrium fraction", "Fracción de equilibrio"))
decimal_comma(ax.yaxis)
i = np.argmin(abs(E - 2e5))
ax.text(200, f[i, 0] - 0.06, "H$^+$", color=NAVY, fontsize=26, va="top")
ax.text(200, f[i, 1] + 0.06, "H$^0$", color=AMBER, fontsize=26, va="bottom")
fig.tight_layout()
save(fig, "charge_fractions.svg")
for e in (10e3, 30e3, 100e3, 300e3):
    print(f"  fractions at {e/1e3:g} keV: H+ {equilibrium_charge_fractions(xs, e)[0]:.2f}")

# ------------------------------------------------ one track, 15 keV, first 25 nm
engine = TrackEngine(WATER, xs, edep_tables=edep, eemean_tables=eemean, seed=1)
p = Particle(kind="proton", position=np.zeros(3), direction=np.array([0.0, 0.0, 1.0]), energy=15e3, q=1)
start = copy.deepcopy(p)
engine.transport(p)
ZMAX = 25.0
dep = [d for d in engine.deposits if d.kind == "proton" and d.position[2] <= ZMAX]
pos = np.array([start.position] + [d.position for d in dep])
fig, ax = plt.subplots(figsize=(15, 4.4))
for k, d in enumerate(dep):
    ax.plot(pos[k:k + 2, 2], pos[k:k + 2, 0], color=NAVY if d.q == 1 else AMBER,
            lw=2.4, zorder=1)
for ch in COLOR:
    pts = np.array([d.position for d in dep if d.channel == ch]).reshape(-1, 3)
    ax.scatter(pts[:, 2], pts[:, 0], s=120, color=COLOR[ch], edgecolor="white", lw=0.8, zorder=3)
ax.set_xlim(0, ZMAX)
ax.set_xlabel(T("Depth $z$ (nm)", "Profundidad $z$ (nm)"))
decimal_comma(ax.yaxis)
ax.set_ylabel("$x$ (nm)")
h = [Line2D([], [], color=NAVY, lw=3, label="H$^+$"), Line2D([], [], color=AMBER, lw=3, label="H$^0$")]
h += [Line2D([], [], ls="", marker="o", ms=12, color=COLOR[c], label=NAME[c]) for c in COLOR]
ax.legend(handles=h, ncol=7, frameon=False, fontsize=21, loc="lower center", bbox_to_anchor=(0.5, 1.0),
          columnspacing=1.2, handletextpad=0.3)
fig.tight_layout()
save(fig, "track_15keV.svg")
print(f"  15 keV track: {len(engine.collision_log)} collisions, first {ZMAX} nm: {len(dep)}")

# ------------------------------------------------ stopping power: model vs PSTAR
Ek = np.geomspace(8, 1e5, 160)
model = [analytic_equilibrium_stopping_power(WATER, xs, e * 1e3, edep, eemean) for e in Ek]
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(PSTAR_KEV, PSTAR_MEV_CM2_G * WATER.density_g_cm3 * 0.1, color=MUTED, lw=2.6, ls="--")
ax.plot(Ek, model, color=RED, lw=3.2)
ax.set_xscale("log")
ax.set_xlim(8, 1e5)
ax.set_ylim(0, 115)
ax.set_xlabel(T("Proton energy (keV)", "Energía del protón (keV)"))
ax.set_ylabel(T("Stopping power (keV/µm)", "Poder de frenado (keV/µm)"))
ax.legend(handles=[Line2D([], [], color=RED, lw=3.2, label="ion-track-mc"),
                   Line2D([], [], color=MUTED, lw=2.6, ls="--", label=T("NIST PSTAR (liquid water)", "NIST PSTAR (agua líquida)"))],
          frameon=False, fontsize=21, loc="upper right")
fig.tight_layout()
save(fig, "stopping_power.svg")

# ------------------------------------------------ Bragg curve + radial dose, 1 MeV beam
N_P = 30
z_edges = np.linspace(0, 26e3, 53)                     # nm (range ≈ 23 µm)
bragg = make_bragg_curve_scorer(WATER, z_edges=z_edges, transverse_half_width_nm=1e6)
r_edges = np.geomspace(0.05, 200, 31)
radial = RadialDoseScorer(WATER, r_edges=r_edges, length_nm=23e3)
res = run_beam(WATER, xs, edep, eemean, n_particles=N_P, E_beam_eV=1e6,
               scorers=[bragg, radial], seed=0, let_segment_length_nm=1000.0)
print(f"  beam: LET_T = {res.track_averaged_LET():.2f} eV/nm, LET_D = {res.dose_averaged_LET():.2f} eV/nm")

zc = 0.5 * (z_edges[1:] + z_edges[:-1]) / 1e3
D = bragg.dose_Gy()[0, 0, :]
fig, ax = plt.subplots(figsize=(10, 6))
ax.step(zc, D / D.max(), where="mid", color=NAVY, lw=3)
ax.set_xlabel(T("Depth (µm)", "Profundidad (µm)"))
ax.set_ylabel(T("Dose (relative)", "Dosis (relativa)"))
decimal_comma(ax.yaxis)
ax.set_xlim(0, 26)
ax.set_ylim(0, 1.08)
ax.text(0.03, 0.95, T(f"{N_P} protons, 1 MeV", f"{N_P} protones, 1 MeV"), transform=ax.transAxes, color=MUTED, fontsize=21, va="top")
fig.tight_layout()
save(fig, "bragg_1MeV.svg")
print(f"  Bragg peak at {zc[np.argmax(D)]:.1f} µm")

rc = np.sqrt(r_edges[1:] * r_edges[:-1])
Dr = radial.dose_Gy() / N_P
fig, ax = plt.subplots(figsize=(10, 6))
m = Dr > 0
ax.plot(rc[m], Dr[m], "o-", color=CYAN, lw=2.6, ms=8)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(T("Distance from the track axis $r$ (nm)", "Distancia al eje de la traza $r$ (nm)"))
ax.set_ylabel(T("Dose per proton (Gy)", "Dosis por protón (Gy)"))
ax.text(0.97, 0.95, T("1 MeV protons, whole track", "protones de 1 MeV, traza completa"), transform=ax.transAxes, ha="right", va="top",
        color=MUTED, fontsize=21)
fig.tight_layout()
save(fig, "radial_dose_1MeV.svg")
