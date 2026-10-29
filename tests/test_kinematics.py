"""Electron-emission-angle formula: deterministic where the classical
binary-encounter formula is valid, isotropic fallback where it isn't (the
cosine argument can exceed 1 at low collision energy).
"""
import numpy as np
import pytest

from ion_track_mc import constants
from ion_track_mc.physics.kinematics import sample_ejected_electron_direction

DIRECTION = np.array([0.0, 0.0, 1.0])
SEED = 99


def test_direction_is_unit_vector():
    rng = np.random.default_rng(SEED)
    E = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    for _ in range(200):
        direction = sample_ejected_electron_direction(E, 50.0, 20.0, DIRECTION, rng)
        assert np.linalg.norm(direction) == pytest.approx(1.0)


def test_deterministic_angle_matches_formula_when_in_domain():
    E = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT  # far from the low-E domain-violation regime
    E_ion = 50.0
    W = 20.0
    M_P = constants.PROTON_TO_ELECTRON_MASS_RATIO
    k1 = 2.0 * M_P * E
    ks = 2.0 * M_P * (E - E_ion)
    ke = 2.0 * W
    expected_cos_theta = (k1 + ke - ks) / (2.0 * np.sqrt(k1 * ke))
    assert abs(expected_cos_theta) <= 1.0  # sanity: this test case must be in-domain

    rng = np.random.default_rng(1)
    got = sample_ejected_electron_direction(E, E_ion, W, DIRECTION, rng)
    assert float(np.dot(got, DIRECTION)) == pytest.approx(expected_cos_theta, abs=1e-9)


def test_isotropic_fallback_when_out_of_domain():
    """At E=8000 eV (H1+), the deterministic formula is known out-of-domain
    (cosine argument ~ 2.05). Check the result is
    isotropic (cos(theta) uniform on [-1, 1]) instead of raising."""
    E = 8000.0
    edep = 14.25434242818631070  # H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat at E=8000 eV
    eemean = 5.11618732            # H1+_+_H2O_ionization_Eemeantot_CDW-EIS.dat at E=8000 eV
    E_ion = edep + eemean
    W = eemean

    rng = np.random.default_rng(2)
    cos_thetas = np.array([
        float(np.dot(sample_ejected_electron_direction(E, E_ion, W, DIRECTION, rng), DIRECTION))
        for _ in range(5000)
    ])
    assert cos_thetas.min() >= -1.0 and cos_thetas.max() <= 1.0
    # isotropic in cos(theta) -> mean ~ 0, std ~ 1/sqrt(3)
    assert abs(cos_thetas.mean()) < 0.05
    assert abs(cos_thetas.std() - 1.0 / np.sqrt(3)) < 0.05


def test_azimuthal_angle_is_uniform():
    """Rotate each sampled direction back into the (x, y) plane relative to
    DIRECTION and check the azimuthal angle is uniform on [0, 2*pi) --
    exercised via the in-domain deterministic branch, whose polar angle is
    fixed, so only phi varies."""
    E = 500.0 * 1.0e3 * constants.PROTON_NUCLEON_COUNT
    rng = np.random.default_rng(3)
    xs, ys = [], []
    for _ in range(2000):
        d = sample_ejected_electron_direction(E, 50.0, 20.0, DIRECTION, rng)
        xs.append(d[0])
        ys.append(d[1])
    phis = np.arctan2(ys, xs)
    # a uniform angle has zero mean resultant vector length
    resultant = np.hypot(np.mean(np.cos(phis)), np.mean(np.sin(phis)))
    assert resultant < 0.05


def test_orthonormal_basis_matches_np_cross_reference():
    """`_orthonormal_basis` writes its cross products out by hand (speed);
    check it gives the same frame as the np.cross version it replaced, for
    random directions and for both reference-axis branches (|w_x| < 0.9 and
    >= 0.9)."""
    from ion_track_mc.physics.kinematics import _orthonormal_basis

    def reference(w):
        ref = np.array([1.0, 0.0, 0.0]) if abs(w[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        u = np.cross(ref, w)
        u = u / np.linalg.norm(u)
        return u, np.cross(w, u), w

    rng = np.random.default_rng(4)
    ws = rng.normal(size=(500, 3))
    ws /= np.linalg.norm(ws, axis=1)[:, None]
    ws = np.vstack([ws, [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [-0.95, 0.0, np.sqrt(1 - 0.95**2)]]])
    for w in ws:
        for got, expected in zip(_orthonormal_basis(w), reference(w)):
            np.testing.assert_allclose(got, expected, rtol=0, atol=1e-15)
