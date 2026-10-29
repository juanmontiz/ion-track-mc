"""Water cross-section and mean-energy tables, loaded once for the docs/ scripts.

Exposes `xs`, `edep`, `eemean` (the three inputs of TrackEngine) and
`PSTAR_KEV`/`PSTAR_MEV_CM2_G`: NIST PSTAR electronic stopping power of liquid
water (material 276, default energy grid), as quoted in notebook 01.
"""
from pathlib import Path

import numpy as np

import ion_track_mc
from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import (
    TabulatedCrossSectionSet, load_cross_section_table, load_mean_energy_table)

DATA = Path(ion_track_mc.__file__).parent / "data" / "cross_sections"
xs = TabulatedCrossSectionSet({k: load_cross_section_table(DATA / f) for k, f in {
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_TCS_CDW-EIS.dat",
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_TCS_Miller-Green.dat",
    (1, Channel.CAPTURE): "H1+_+_H2O_capture_TCS_CDW-EIS.dat",
    (1, Channel.ELASTIC): "H1+_+_H2O_elastic_TCS_Classical.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_TCS_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_TCS_Miller-Green.dat",
    (0, Channel.LOSS): "H0_+_H2O_loss_TCS_Miller-Green.dat",
    (0, Channel.ELASTIC): "H0_+_H2O_elastic_TCS_Classical.dat"}.items()})
edep = {k: load_mean_energy_table(DATA / f) for k, f in {
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_Edepmeantot_Miller-Green.dat",
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat",
    (1, Channel.CAPTURE): "H1+_+_H2O_capture_Edepmeantot_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_Edepmeantot_Miller-Green.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_Edepmeantot_CDW-EIS.dat",
    (0, Channel.LOSS): "H0_+_H2O_loss_Edepmeantot_constant.dat"}.items()}
eemean = {k: load_mean_energy_table(DATA / f) for k, f in {
    1: "H1+_+_H2O_ionization_Eemeantot_CDW-EIS.dat",
    0: "H0_+_H2O_ionization_Eemeantot_CDW-EIS.dat"}.items()}

PSTAR_KEV = np.array([8, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000,
                      50000, 100000, 200000, 300000], float)
PSTAR_MEV_CM2_G = np.array([378.2, 422.9, 567.3, 774.0, 814.5, 660.4, 412.8, 260.6,
                            158.5, 79.06, 45.64, 26.05, 12.44, 7.286, 4.491, 3.519])
