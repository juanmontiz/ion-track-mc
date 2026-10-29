from pathlib import Path

import pytest

from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import (
    TabulatedCrossSectionSet,
    load_cross_section_table,
    load_mean_energy_table,
)
from ion_track_mc.physics.medium import WATER

DATA_DIR = Path(__file__).resolve().parents[1] / "ion_track_mc" / "data" / "cross_sections"

FULL_TABLE_FILES = {
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_TCS_CDW-EIS.dat",
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_TCS_Miller-Green.dat",
    (1, Channel.CAPTURE): "H1+_+_H2O_capture_TCS_CDW-EIS.dat",
    (1, Channel.ELASTIC): "H1+_+_H2O_elastic_TCS_Classical.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_TCS_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_TCS_Miller-Green.dat",
    (0, Channel.LOSS): "H0_+_H2O_loss_TCS_Miller-Green.dat",
    (0, Channel.ELASTIC): "H0_+_H2O_elastic_TCS_Classical.dat",
}

FULL_EDEP_TABLE_FILES = {
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_Edepmeantot_Miller-Green.dat",
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat",
    (1, Channel.CAPTURE): "H1+_+_H2O_capture_Edepmeantot_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_Edepmeantot_Miller-Green.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_Edepmeantot_CDW-EIS.dat",
    (0, Channel.LOSS): "H0_+_H2O_loss_Edepmeantot_constant.dat",
}

FULL_EEMEAN_TABLE_FILES = {
    1: "H1+_+_H2O_ionization_Eemeantot_CDW-EIS.dat",
    0: "H0_+_H2O_ionization_Eemeantot_CDW-EIS.dat",
}


@pytest.fixture
def full_xs() -> TabulatedCrossSectionSet:
    tables = {
        key: load_cross_section_table(DATA_DIR / fname)
        for key, fname in FULL_TABLE_FILES.items()
    }
    return TabulatedCrossSectionSet(tables)


@pytest.fixture
def full_edep():
    return {
        key: load_mean_energy_table(DATA_DIR / fname)
        for key, fname in FULL_EDEP_TABLE_FILES.items()
    }


@pytest.fixture
def full_eemean():
    return {
        q: load_mean_energy_table(DATA_DIR / fname)
        for q, fname in FULL_EEMEAN_TABLE_FILES.items()
    }


@pytest.fixture
def water():
    return WATER
