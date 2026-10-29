"""Monte Carlo de fotones en una placa de agua, COMPLETO (versión básica del TP4-A).

Solo hay que ejecutarlo: `simular_placa(E0, espesor, n_fotones)` devuelve la
fracción transmitida en haz estrecho y en haz ancho. El coeficiente Compton
del agua se calcula con la fórmula de Klein–Nishina integrada; el efecto
fotoeléctrico se desprecia (válido por encima de ~100 keV en agua).

Unidades: energías en MeV, longitudes en cm.
"""
import numpy as np

ME_C2 = 0.51099895   # MeV


def interp_loglog(x_tabla, y_tabla, x):
    return np.exp(np.interp(np.log(x), np.log(x_tabla), np.log(y_tabla)))


# ---------------------------------------------------------------------------
# 1. Física de la dispersión Compton
# ---------------------------------------------------------------------------
def energia_dispersada(E, cos_t):
    """Fórmula de Compton: energía del fotón dispersado [MeV]."""
    return E / (1.0 + (E / ME_C2) * (1.0 - cos_t))


def klein_nishina_rel(E, cos_t):
    """(dσ/dΩ) / r0², relativa a su valor en θ = 0 (que vale 1).

    """
    r = energia_dispersada(E, cos_t) / E
    return 0.5 * r**2 * (r + 1/r - (1 - cos_t**2))


def muestrear_cos_compton(E, rng):
    """cos θ con la distribución de Klein–Nishina, por rechazo: se propone
    cos θ uniforme en [−1, 1] y se acepta con probabilidad klein_nishina_rel."""
    while True:
        c = rng.uniform(-1, 1)
        if rng.random() < klein_nishina_rel(E, c):
            return c


# ---------------------------------------------------------------------------
# 2. Geometría
# ---------------------------------------------------------------------------
def rotar(u, cos_t, phi):
    """Nueva dirección a un ángulo polar θ (cos_t) y azimutal phi respecto de u."""
    sin_t = np.sqrt(max(0.0, 1.0 - cos_t**2))
    ux, uy, uz = u
    if abs(uz) > 0.99999:
        return np.array([sin_t * np.cos(phi), sin_t * np.sin(phi), np.sign(uz) * cos_t])
    d = np.sqrt(1.0 - uz**2)
    return np.array([
        ux * cos_t + sin_t * (ux * uz * np.cos(phi) - uy * np.sin(phi)) / d,
        uy * cos_t + sin_t * (uy * uz * np.cos(phi) + ux * np.sin(phi)) / d,
        uz * cos_t - sin_t * np.cos(phi) * d,
    ])


# ---------------------------------------------------------------------------
# 3. Transporte en la placa
# ---------------------------------------------------------------------------
def transportar_foton(E0, espesor, mu_compton, mu_foto, rng):
    """Sigue un fotón que entra en z = 0 según +z en una placa [0, espesor].

    `mu_compton(E)` y `mu_foto(E)` son funciones que devuelven 1/cm.
    Devuelve (resultado, n_interacciones, E_final, E_depositada), con
    resultado en {'transmitido', 'retrodispersado', 'absorbido'}.
    """
    z, u, E = 0.0, np.array([0.0, 0.0, 1.0]), E0
    n, E_dep = 0, 0.0
    while True:
        mu_c, mu_f = mu_compton(E), mu_foto(E)
        mu_tot = mu_c + mu_f
        # distancia hasta la próxima interacción: s = −ln ξ / μ
        s = -np.log(rng.random()) / mu_tot
        z += s * u[2]
        if z >= espesor:
            return "transmitido", n, E, E_dep
        if z < 0.0:
            return "retrodispersado", n, E, E_dep
        n += 1
        # fotoeléctrico con probabilidad mu_f/mu_tot
        if rng.random() < mu_f / mu_tot:
            E_dep += E                      # fotoeléctrico: todo se deposita localmente
            return "absorbido", n, 0.0, E_dep
        # Compton
        cos_t = muestrear_cos_compton(E, rng)
        E_nueva = energia_dispersada(E, cos_t)
        E_dep += E - E_nueva                # energía cedida al electrón Compton (local)
        E = E_nueva
        u = rotar(u, cos_t, 2 * np.pi * rng.random())


def simular_placa(E0, espesor, n_fotones, mu_compton, mu_foto, seed=0):
    """Transmisión en geometría estrecha y ancha para n_fotones."""
    rng = np.random.default_rng(seed)
    estrecho = ancho = 0
    for _ in range(n_fotones):
        res, n, E, _ = transportar_foton(E0, espesor, mu_compton, mu_foto, rng)
        if res == "transmitido":
            ancho += 1
            if n == 0:
                estrecho += 1
    return estrecho / n_fotones, ancho / n_fotones


R0, NE = 2.8179403262e-13, 3.3428e23      # radio clásico del electrón (cm), electrones/cm³ en agua


def mu_compton(E):
    """Coeficiente de atenuación Compton del agua [1/cm] (Klein–Nishina integrada)."""
    k = E / ME_C2
    s = 2 * np.pi * R0**2 * ((1 + k) / k**2 * (2 * (1 + k) / (1 + 2 * k) - np.log(1 + 2 * k) / k)
                             + np.log(1 + 2 * k) / (2 * k) - (1 + 3 * k) / (1 + 2 * k) ** 2)
    return NE * s


def mu_foto(E):
    return 0.0


_simular_placa_general = simular_placa


def simular_placa(E0, espesor, n_fotones, seed=0):
    """(fracción transmitida en haz estrecho, en haz ancho) para n_fotones de E0 [MeV]."""
    return _simular_placa_general(E0, espesor, n_fotones, mu_compton, mu_foto, seed=seed)
