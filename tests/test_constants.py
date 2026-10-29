"""Locks the pytest harness and sanity-checks the unit-policy constants."""
import pytest

from ion_track_mc import constants


def test_harness_works():
    assert 1 + 1 == 2


def test_cm2_to_nm2_conversion():
    assert constants.CM2_TO_NM2 == 1.0e14


def test_water_number_density_is_positive():
    assert constants.N_WATER_NM3 > 0


def test_h_binding_energy_ev():
    assert constants.H_BINDING_ENERGY_EV == pytest.approx(13.6)


def test_electron_rest_energy_ev():
    assert constants.ELECTRON_REST_ENERGY_EV == pytest.approx(510999.0645)


def test_fine_structure_constant():
    assert constants.FINE_STRUCTURE_CONSTANT == pytest.approx(1.0 / 137.0)


def test_proton_to_electron_mass_ratio_is_reciprocal_of_existing_constant():
    assert constants.PROTON_TO_ELECTRON_MASS_RATIO == pytest.approx(
        1.0 / constants.ELECTRON_TO_PROTON_MASS_RATIO
    )


def test_proton_nucleon_count_is_one():
    # ^1H has exactly one nucleon (the proton itself) -- this is an integer
    # mass number A, distinct from the isotope's physical mass in amu
    # (1.0078250319 u, which includes the proton-electron mass difference
    # and is NOT what "energy per nucleon" divides by).
    assert constants.PROTON_NUCLEON_COUNT == 1


def test_kev_per_u_divides_by_nucleon_count_not_atomic_mass():
    # "keV/u" is the standard accelerator-physics "energy per nucleon"
    # convention: total kinetic energy = (keV/u value) * A, where A is the
    # integer nucleon count -- NOT the isotope's actual mass in amu. For the
    # ^1H projectile (A=1) this must come out to exactly 1000 eV per 1
    # keV/u, not 1000 * 1.0078250319 eV (the old, incorrect "keV/amu"
    # reading that divided by atomic mass instead of nucleon count).
    assert constants.ENERGY_UNIT_TO_EV["keV/u"] == pytest.approx(1.0e3)


def test_kev_per_amu_key_removed():
    # Renamed to "keV/u" (see above) -- the old key name conflated "per
    # nucleon" with "per atomic mass unit", which are only numerically
    # close (not equal) for a real isotope.
    assert "keV/amu" not in constants.ENERGY_UNIT_TO_EV
