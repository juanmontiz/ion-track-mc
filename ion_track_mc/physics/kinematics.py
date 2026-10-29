"""Direction of the delta-electron ejected in an ionization.

The electron's energy W comes from physics/energy_loss.py::ionization_energy;
this module picks the direction it flies off in, from momentum conservation
(physics/README_physics.md §5). It also provides the small vector helper
used to rotate any (theta, phi) angle pair into the frame of a moving
particle.
"""
import numpy as np

from ion_track_mc import constants


def _orthonormal_basis(w: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a right-handed {u, v, w} orthonormal basis with unit vector `w`.

    With w along the particle's current direction, a new direction at polar
    angle theta and azimuth phi is
    sin(theta)cos(phi) u + sin(theta)sin(phi) v + cos(theta) w.
    """
    # u = np.cross(reference, w) and v = np.cross(w, u), written out by hand --
    # numpy's generic cross() has per-call overhead unsuited to fixed
    # 3-vectors, and this runs on every ionization and elastic event.
    wx, wy, wz = w
    # Pick a reference axis not nearly parallel to w, so the cross product
    # with it is well conditioned.
    if abs(wx) < 0.9:
        u = np.array([0.0, -wz, wy])    # reference = x-hat
    else:
        u = np.array([wz, 0.0, -wx])    # reference = y-hat
    u = u / np.linalg.norm(u)
    ux, uy, uz = u
    v = np.array([wy * uz - wz * uy, wz * ux - wx * uz, wx * uy - wy * ux])
    return u, v, w


def sample_ejected_electron_direction(
    E: float,
    E_ion: float,
    W: float,
    incident_direction: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Direction of the delta-electron ejected in an ionization event.

    Momentum conservation: the incoming proton momentum p1 equals the
    scattered proton momentum ps plus the electron momentum pe (the target
    molecule's recoil is neglected). These three vectors form a triangle;
    the law of cosines then fixes the electron's polar angle:

        cos(theta_e) = (k1 + ke - ks) / (2 sqrt(k1 ke)),

    with k = p^2 in atomic units: k1 = 2 M_P E, ks = 2 M_P (E - E_ion),
    ke = 2 W.

    This two-body picture breaks down at low collision energy, where the
    formula gives |cos(theta_e)| > 1 (no such triangle exists). In that
    case the electron is emitted isotropically, cos(theta) ~ U(-1, 1). The
    check is made per collision; there is no fixed energy cut-off.

    The azimuth phi ~ U(0, 2 pi) is always random (momentum conservation
    only fixes the polar angle).
    """
    M_P = constants.PROTON_TO_ELECTRON_MASS_RATIO
    k1 = 2.0 * M_P * E
    ks = 2.0 * M_P * (E - E_ion)
    ke = 2.0 * W
    arg = (k1 + ke - ks) / (2.0 * np.sqrt(k1 * ke))

    if abs(arg) <= 1.0:
        cos_theta = float(arg)                      # momentum triangle closes
    else:
        cos_theta = float(rng.uniform(-1.0, 1.0))   # isotropic fallback

    sin_theta = float(np.sqrt(max(0.0, 1.0 - cos_theta ** 2)))
    phi = rng.uniform(0.0, 2.0 * np.pi)
    u, v, w = _orthonormal_basis(incident_direction)
    return sin_theta * np.cos(phi) * u + sin_theta * np.sin(phi) * v + cos_theta * w
