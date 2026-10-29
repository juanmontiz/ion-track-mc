"""Physical constants and unit converters.

Every quantity inside the code uses one fixed set of units, so no function
ever has to guess what unit a number is in.

UNIT POLICY (internal):
    energy    -> eV
    length    -> nm
    cross sec -> nm^2   (raw tables come in whatever RAW_CROSS_SECTION_UNIT is;
                          convert on load)
    number density N -> nm^-3

Raw data tables (cross sections, mean deposited/emitted energies) are
provided in whatever units the source (theory/experiment) used; conversion
to internal units happens once, at load time, inside the CrossSectionSet
loader (`physics/cross_section.py`). Which raw unit a given table is in is
never hardcoded in the loader itself -- it looks up the conversion factor in
the dicts below, keyed by a unit name string. To support a differently-
tabulated data set, add/point at a different unit name; no loader code
changes needed.
"""

EV_PER_HARTREE = 27.211386245988   # a.u. of energy -> eV
CM2_TO_NM2     = 1.0e14            # 1 cm^2 = 1e14 nm^2
CM3_TO_NM3     = 1.0e21            # 1 cm^3 = 1e21 nm^3
NM_PER_CM      = 1.0e7             # 1 cm = 1e7 nm
AVOGADRO       = 6.02214076e23     # mol^-1 (CODATA)

# m_e / M_p (CODATA proton-electron mass ratio 1836.15267343), used for the
# non-relativistic binary-encounter kinematic ceiling W_max = 4*(m/M)*E.
ELECTRON_TO_PROTON_MASS_RATIO = 1.0 / 1836.15267343

EV_TO_JOULE = 1.602176634e-19      # exact, SI 2019 redefinition (elementary charge)

# Nucleon count (integer mass number A) of the projectile, always ^1H here.
# This is what "energy per nucleon" units (keV/u, MeV/u -- the standard
# accelerator-physics convention) divide by. It is NOT the isotope's actual
# physical mass in atomic mass units (1.0078250319 u for ^1H, which includes
# the proton-electron mass difference) -- those are only numerically close,
# not equal, and conflating them was a real (if small, ~0.78%) unit bug.
PROTON_NUCLEON_COUNT = 1

# Water number density at 1 g/cm^3 (molecules / nm^3):
#   (1 g/cm^3 / 18.015 g/mol) * N_A  ->  ~3.34e22 cm^-3  ->  33.4 nm^-3
N_WATER_NM3 = 33.4

# --- Raw data-table units --------------------------------------------------
#
# Every table file under ion_track_mc/data/cross_sections/ has a collision-
# energy column and a value column (cross section, or a mean deposited/
# emitted energy). Both are tabulated in whatever unit the source used, not
# necessarily the internal unit. The loader (`physics/cross_section.py`)
# converts once, on load, by looking up the raw unit's name in the dicts
# below -- so switching a table's unit is a one-line change to the unit name
# passed to the loader (or to the DEFAULT_* below), never a loader edit.
#
# Energy-column units -> eV.
ENERGY_UNIT_TO_EV = {
    "eV": 1.0,
    "keV": 1.0e3,
    "keV/u": 1.0e3 * PROTON_NUCLEON_COUNT,  # energy per nucleon; projectile is always H (A=1) here
}

# Cross-section-column units -> nm^2.
CROSS_SECTION_UNIT_TO_NM2 = {
    "cm2": CM2_TO_NM2,                     # 1 cm^2      = 1e14 nm^2
    "1e-16cm2": CM2_TO_NM2 * 1.0e-16,      # 1e-16 cm^2  = 1e-2 nm^2
}

# The units actually used by the current data drop under
# data/cross_sections/: collision energy in eV, cross sections in 1e-16 cm^2
# (i.e. values are already "in units of 1e-16 cm^2" -- multiply the raw
# number by CROSS_SECTION_UNIT_TO_NM2["1e-16cm2"]). Mean deposited/emitted
# energy columns (Edepmeantot, Eemeantot) are energies, so they use
# ENERGY_UNIT_TO_EV too, not a separate table.
DEFAULT_ENERGY_UNIT = "eV"
DEFAULT_CROSS_SECTION_UNIT = "1e-16cm2"

# --- Projectile energy-loss physics ----------------------------------------
#
# Constants used by the per-collision energy-loss and emission-angle formulas
# in physics/energy_loss.py and physics/kinematics.py (derivations in
# physics/README_physics.md §4-§5).
H_BINDING_ENERGY_EV = 13.6             # H ground-state binding energy (capture/loss Q-value term)
ELECTRON_REST_ENERGY_EV = 510999.0645  # m_e * c^2, used by the elastic-scattering momentum formula
FINE_STRUCTURE_CONSTANT = 1.0 / 137.0  # alpha, used by the elastic-scattering screening-angle formula
PROTON_TO_ELECTRON_MASS_RATIO = 1836.15267343  # M_P; reciprocal of ELECTRON_TO_PROTON_MASS_RATIO
# above, kept as its own named constant since the energy-loss formulas use
# M_P directly (not 1/m_e_to_M_p) at every call site.
