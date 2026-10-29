"""How much energy the projectile loses in each kind of collision.

Two tabulated quantities, both functions of the projectile energy E, drive
all of it (physics/README_physics.md §4):
  - Edepmeantot(E): mean energy left in the water molecule at the collision
    site (excitation energy, binding energy, ...);
  - Eemeantot(E): mean kinetic energy W of the delta-electron ejected by an
    ionization.

Ionization, excitation, capture and loss use these mean values, so their
energy loss is deterministic at a given E. Elastic scattering is the one
stochastic loss: its energy transfer follows from a randomly sampled
scattering angle.

`table`/`edep_table`/`eemean_table` parameters throughout are always
(E_eV, value_eV) pairs as returned by
physics.cross_section.load_mean_energy_table.
"""
import numpy as np

from ion_track_mc import constants
from ion_track_mc.physics import kinematics
from ion_track_mc.physics.cross_section import interp_log_log


def mean_energy_lookup(E: float, table: tuple[np.ndarray, np.ndarray]) -> float:
    """Edepmeantot/Eemeantot table lookup, log-log interpolated.

    Used directly for excitation (E_exc = Edepmeantot(E)) and loss
    (E_loss = Edepmeantot(E); that table is a flat 13.6 eV, the binding
    energy of the electron H0 gives up).
    """
    E_table, value_table = table
    return interp_log_log(E, E_table, value_table)


def capture_delta_E(E: float, edep_table: tuple[np.ndarray, np.ndarray]) -> float:
    """Energy lost by H+ when it captures an electron and becomes H0:

        E_cap(E) = Edepmeantot(E) + E*(m_e/M_p) - 13.6 eV.

    Three terms (physics/README_physics.md §4.3): the energy left in the
    ionised molecule; the kinetic energy the captured electron needs to
    travel along with the proton (it moves at the same speed, so it carries
    E*m_e/M_p); minus the 13.6 eV released when it binds into the H ground
    state. The result is positive over the whole tabulated energy range
    (checked in tests/test_energy_loss.py).
    """
    edep = mean_energy_lookup(E, edep_table)
    return edep + E * constants.ELECTRON_TO_PROTON_MASS_RATIO - constants.H_BINDING_ENERGY_EV


def ionization_energy(
    E: float,
    edep_table: tuple[np.ndarray, np.ndarray],
    eemean_table: tuple[np.ndarray, np.ndarray],
) -> tuple[float, float]:
    """Returns (deposit_local, W) = (Edepmeantot(E), Eemeantot(E)).

    E_ion(E) = Edepmeantot(E) + Eemeantot(E) is the total energy the
    projectile loses in one ionization: deposit_local (binding energy,
    left at the collision site) plus W (kinetic energy handed to the
    ejected delta-electron, which carries it away). W is the tabulated mean,
    not sampled from a distribution.
    """
    deposit_local = mean_energy_lookup(E, edep_table)
    W = mean_energy_lookup(E, eemean_table)
    return deposit_local, W


def elastic_scatter(
    E: float,
    incident_direction: np.ndarray,
    target_nuclear_charge: float | None,
    target_nucleon_count: float | None,
    rng: np.random.Generator,
) -> tuple[float, np.ndarray]:
    """Elastic scattering off the target nuclei: deflection and energy loss.

    Chain of steps (physics/README_physics.md §4.5): projectile momentum ->
    screening cut-off angle -> random scattering angle -> energy transferred
    to the recoiling target. Returns (E_elastic [eV], new_direction).

    The angle is drawn from the Rutherford distribution, proportional to
    1/sin^4(theta/2), cut off below theta_cut where the target's electrons
    screen the nuclear charge. theta therefore ranges from theta_cut
    (eta=0) to pi (eta=1) and is never exactly 0. Small angles dominate, so
    most elastic collisions barely deflect the proton.

    Raises `ValueError` if `target_nuclear_charge`/`target_nucleon_count` are
    `None` (e.g. a custom `Medium` that didn't set them -- see
    `physics/medium.py`'s `Medium` docstring) rather than letting the
    `** (1.0 / 3.0)` below fail with a confusing `TypeError` deep into a beam
    run.
    """
    if target_nuclear_charge is None or target_nucleon_count is None:
        raise ValueError(
            "elastic_scatter requires Medium.target_nuclear_charge and "
            "target_nucleon_count to be set (WATER sets both; a custom "
            "Medium must set them explicitly before elastic scattering will "
            "work)."
        )
    M_P = constants.PROTON_TO_ELECTRON_MASS_RATIO
    m_e_c2 = constants.ELECTRON_REST_ENERGY_EV
    alpha = constants.FINE_STRUCTURE_CONSTANT

    # Relativistic momentum p*c [eV] of the proton (rest energy M_P*m_e*c^2).
    rest_energy = M_P * m_e_c2
    p = np.sqrt((E + rest_energy) ** 2 - rest_energy ** 2)
    # Screening angle: below it the atomic electrons hide the nuclear charge
    # (Thomas-Fermi screening radius ~ a0 * Z^(-1/3)).
    theta_cut = alpha * target_nuclear_charge ** (1.0 / 3.0) * (m_e_c2 / p)

    # Inverse-transform sampling of the Rutherford distribution truncated at
    # theta_cut: in the variable 1/sin^2(theta/2) the cumulative distribution
    # is linear, so a uniform eta maps straight onto it.
    sin2inv = 1.0 / np.sin(theta_cut / 2.0) ** 2
    eta = rng.random()
    bracket = eta * (1.0 - sin2inv) + sin2inv
    theta = 2.0 * np.arcsin(bracket ** -0.5)

    # Two-body elastic kinematics: energy given to a target of mass M_T
    # (binding of the nucleons ignored, M_T = A_T * M_P).
    M_T = target_nucleon_count * M_P
    E_elastic = 4.0 * E * (M_P * M_T) / (M_P + M_T) ** 2 * np.sin(theta / 2.0) ** 2

    # New direction: polar angle theta about the old direction, azimuth phi
    # uniform, expressed in a basis whose w axis is the old direction.
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    phi = rng.uniform(0.0, 2.0 * np.pi)
    u, v, w = kinematics._orthonormal_basis(incident_direction)
    new_direction = sin_theta * np.cos(phi) * u + sin_theta * np.sin(phi) * v + cos_theta * w

    return float(E_elastic), new_direction
