import pytest

from ion_track_mc import constants
from ion_track_mc.physics.medium import (
    WATER,
    Medium,
    check_medium_matches_filename,
    target_label_from_filename,
)


def test_water_number_density_matches_design_doc_value():
    # N_water ~ 33.4 nm^-3 at 1 g/cm^3.
    assert WATER.N == pytest.approx(constants.N_WATER_NM3, rel=1e-2)


def test_N_scales_with_density():
    dense = Medium(name="dense_water", density_g_cm3=2.0,
                    molar_mass_g_mol=18.015, target_label="H2O")
    assert dense.N == pytest.approx(2.0 * WATER.N)


def test_N_scales_inversely_with_molar_mass():
    heavy = Medium(name="heavy", density_g_cm3=1.0,
                    molar_mass_g_mol=2.0 * 18.015, target_label="H2O")
    assert heavy.N == pytest.approx(0.5 * WATER.N)


def test_target_label_from_filename():
    assert target_label_from_filename(
        "H1+_+_H2O_ionization_TCS_CDW-EIS.dat"
    ) == "H2O"
    assert target_label_from_filename(
        "H0_+_H2O_loss_TCS_Miller-Green.dat"
    ) == "H2O"


def test_target_label_from_filename_rejects_bad_convention():
    with pytest.raises(ValueError):
        target_label_from_filename("not_a_valid_filename.dat")


def test_check_medium_matches_filename_ok():
    check_medium_matches_filename(WATER, "H1+_+_H2O_ionization_TCS_CDW-EIS.dat")


def test_check_medium_matches_filename_mismatch_raises():
    pmma = Medium(name="pmma", density_g_cm3=1.18, molar_mass_g_mol=100.12,
                  target_label="PMMA")
    with pytest.raises(ValueError):
        check_medium_matches_filename(pmma, "H1+_+_H2O_ionization_TCS_CDW-EIS.dat")


def test_water_target_nuclear_charge_and_nucleon_count():
    assert WATER.target_nuclear_charge == pytest.approx(10.0)
    assert WATER.target_nucleon_count == pytest.approx(18.0)


def test_medium_without_target_fields_defaults_to_none():
    m = Medium(name="dense_water", density_g_cm3=2.0, molar_mass_g_mol=18.015, target_label="H2O")
    assert m.target_nuclear_charge is None
    assert m.target_nucleon_count is None
