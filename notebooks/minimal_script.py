import numpy as np
from pathlib import Path
from ion_track_mc import constants
from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import (
    TabulatedCrossSectionSet,
    load_cross_section_table,
    load_mean_energy_table,
)
from ion_track_mc.physics.medium import WATER
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle
from ion_track_mc.scoring.dose import make_bragg_curve_scorer
from ion_track_mc.analysis.visualization import plot_bragg_curve
import matplotlib.pyplot as plt
import ion_track_mc

# 1. load cross sections (all 8 charge/channel tables, incl. elastic)
DATA_DIR = Path(ion_track_mc.__file__).parent / "data" / "cross_sections"
FILES = {
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_TCS_CDW-EIS.dat",
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_TCS_Miller-Green.dat",
    (1, Channel.CAPTURE):    "H1+_+_H2O_capture_TCS_CDW-EIS.dat",
    (1, Channel.ELASTIC):    "H1+_+_H2O_elastic_TCS_Classical.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_TCS_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_TCS_Miller-Green.dat",
    (0, Channel.LOSS):       "H0_+_H2O_loss_TCS_Miller-Green.dat",
    (0, Channel.ELASTIC):    "H0_+_H2O_elastic_TCS_Classical.dat",
}
xs = TabulatedCrossSectionSet({k: load_cross_section_table(DATA_DIR / f) for k, f in FILES.items()})

EDEP_FILES = {
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_Edepmeantot_Miller-Green.dat",
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat",
    (1, Channel.CAPTURE):    "H1+_+_H2O_capture_Edepmeantot_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_Edepmeantot_Miller-Green.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_Edepmeantot_CDW-EIS.dat",
    (0, Channel.LOSS):       "H0_+_H2O_loss_Edepmeantot_constant.dat",
}
EEMEAN_FILES = {
    1: "H1+_+_H2O_ionization_Eemeantot_CDW-EIS.dat",
    0: "H0_+_H2O_ionization_Eemeantot_CDW-EIS.dat",
}
edep_tables = {k: load_mean_energy_table(DATA_DIR / f) for k, f in EDEP_FILES.items()}
eemean_tables = {k: load_mean_energy_table(DATA_DIR / f) for k, f in EEMEAN_FILES.items()}

# 2. pick a scorer (Bragg curve here)
scorer = make_bragg_curve_scorer(WATER, z_edges=np.linspace(-1e9, 1e9, 41), transverse_half_width_nm=1e9)

# 3. run one primary
engine = TrackEngine(WATER, xs, edep_tables=edep_tables, eemean_tables=eemean_tables, scorers=[scorer], seed=0)
p = Particle(kind="proton", position=np.zeros(3), direction=np.array([0.,0.,1.]),
             energy=1000e3 * constants.PROTON_NUCLEON_COUNT, q=1)  # 1 MeV/u
engine.transport(p)

# 4. plot
plot_bragg_curve(scorer)
plt.savefig("bragg_curve.png")
