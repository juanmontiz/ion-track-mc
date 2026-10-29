"""TP4-A: esqueleto del Monte Carlo de fotones en una placa de agua.

Completar las funciones marcadas con TODO. El resto (geometría, rotación de
la dirección, conteo) ya está escrito.

Unidades: energías en MeV, longitudes en cm, coeficientes de atenuación
lineales en 1/cm.
"""
import numpy as np

ME_C2 = 0.51099895   # MeV


# ---------------------------------------------------------------------------
# 1. Datos de atenuación (XCOM)
# ---------------------------------------------------------------------------
def cargar_xcom(archivo):
    """Lee una tabla exportada de XCOM para el agua.

    TODO: adaptar a las columnas del archivo descargado. Se necesitan:
    E [MeV], mu/rho incoherente (Compton) y mu/rho fotoeléctrico [cm²/g].
    Devuelve tres arrays (E, mu_compton, mu_foto) en 1/cm (agua: rho = 1 g/cm³).
    """
    raise NotImplementedError


def interp_loglog(x_tabla, y_tabla, x):
    return np.exp(np.interp(np.log(x), np.log(x_tabla), np.log(y_tabla)))


# ---------------------------------------------------------------------------
# 2. Física de la dispersión Compton
# ---------------------------------------------------------------------------
def energia_dispersada(E, cos_t):
    """Fórmula de Compton: energía del fotón dispersado [MeV]."""
    return E / (1.0 + (E / ME_C2) * (1.0 - cos_t))


def klein_nishina_rel(E, cos_t):
    """(dσ/dΩ) / r0², relativa a su valor en θ = 0 (que vale 1).

    TODO: implementar
        (1/2) (E'/E)² (E'/E + E/E' − sin²θ)
    """
    raise NotImplementedError


def muestrear_cos_compton(E, rng):
    """TODO: muestrear cos θ con la distribución de Klein–Nishina por rechazo:
    proponer cos θ uniforme en [−1, 1] y aceptarlo con probabilidad
    klein_nishina_rel(E, cos θ) (que es ≤ 1)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# 3. Geometría (ya escrito)
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
# 4. Transporte en la placa
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
        # TODO: muestrear la distancia s hasta la próxima interacción.
        s = ...
        z += s * u[2]
        if z >= espesor:
            return "transmitido", n, E, E_dep
        if z < 0.0:
            return "retrodispersado", n, E, E_dep
        n += 1
        # TODO: elegir el proceso con probabilidad mu_f/mu_tot (fotoeléctrico).
        if ...:
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
