"""Utilidades comunes para los TP del curso de Radiobiología.

Envuelve la carga de datos y el arranque del código Monte Carlo `ion_track_mc`
para que cada TP no repita ~40 líneas de preparación. Requiere tener
`ion_track_mc` instalado (ver 00_Guia_general.md).

Unidades: energías en eV, longitudes en nm (las del toolkit).
"""
import copy
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import ion_track_mc
from ion_track_mc.physics.channels import Channel
from ion_track_mc.physics.cross_section import (
    TabulatedCrossSectionSet,
    load_cross_section_table,
    load_mean_energy_table,
)
from ion_track_mc.physics.medium import WATER
from ion_track_mc.transport.engine import TrackEngine
from ion_track_mc.transport.particle import Particle

DIR_DATOS = Path(ion_track_mc.__file__).parent / "data" / "cross_sections"

_TCS = {
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_TCS_CDW-EIS.dat",
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_TCS_Miller-Green.dat",
    (1, Channel.CAPTURE):    "H1+_+_H2O_capture_TCS_CDW-EIS.dat",
    (1, Channel.ELASTIC):    "H1+_+_H2O_elastic_TCS_Classical.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_TCS_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_TCS_Miller-Green.dat",
    (0, Channel.LOSS):       "H0_+_H2O_loss_TCS_Miller-Green.dat",
    (0, Channel.ELASTIC):    "H0_+_H2O_elastic_TCS_Classical.dat",
}
_EDEP = {
    (1, Channel.EXCITATION): "H1+_+_H2O_excitation_Edepmeantot_Miller-Green.dat",
    (1, Channel.IONIZATION): "H1+_+_H2O_ionization_Edepmeantot_CDW-EIS.dat",
    (1, Channel.CAPTURE):    "H1+_+_H2O_capture_Edepmeantot_CDW-EIS.dat",
    (0, Channel.EXCITATION): "H0_+_H2O_excitation_Edepmeantot_Miller-Green.dat",
    (0, Channel.IONIZATION): "H0_+_H2O_ionization_Edepmeantot_CDW-EIS.dat",
    (0, Channel.LOSS):       "H0_+_H2O_loss_Edepmeantot_constant.dat",
}
_EEMEAN = {
    1: "H1+_+_H2O_ionization_Eemeantot_CDW-EIS.dat",
    0: "H0_+_H2O_ionization_Eemeantot_CDW-EIS.dat",
}


def cargar_agua():
    """Secciones eficaces y tablas de energía media de H+/H0 en agua.

    Devuelve un dict con claves 'xs', 'edep', 'eemean' y 'medio'.
    """
    xs = TabulatedCrossSectionSet(
        {k: load_cross_section_table(DIR_DATOS / f) for k, f in _TCS.items()})
    edep = {k: load_mean_energy_table(DIR_DATOS / f) for k, f in _EDEP.items()}
    eemean = {k: load_mean_energy_table(DIR_DATOS / f) for k, f in _EEMEAN.items()}
    return {"xs": xs, "edep": edep, "eemean": eemean, "medio": WATER}


def nuevo_motor(datos, seed=0, scorers=None):
    """Un TrackEngine de agua listo para transportar protones."""
    return TrackEngine(datos["medio"], datos["xs"], edep_tables=datos["edep"],
                       eemean_tables=datos["eemean"], scorers=scorers, seed=seed)


def proton(E_keV, posicion=(0.0, 0.0, 0.0), direccion=(0.0, 0.0, 1.0), q=1):
    """Un protón (q=1) o hidrógeno neutro (q=0) de energía E_keV."""
    return Particle(kind="proton", position=np.array(posicion, dtype=float),
                    direction=np.array(direccion, dtype=float), energy=E_keV * 1e3, q=q)


def simular_proton(datos, E_keV, seed=0, scorers=None):
    """Transporta un protón que entra en el origen según +z.

    Devuelve (motor, entrada, final): `entrada` es una copia de la partícula
    antes del transporte (la pide `plot_track`) y `final` la partícula al
    terminar (su `position` es el punto donde se detuvo).
    """
    motor = nuevo_motor(datos, seed=seed, scorers=scorers)
    p = proton(E_keV)
    entrada = copy.deepcopy(p)
    final = motor.transport(p)
    return motor, entrada, final


def ionizaciones(motor):
    """Posiciones [nm] de las ionizaciones del ion, array (n, 3)."""
    pos = [d.position for d in motor.deposits if d.channel == Channel.IONIZATION]
    return np.array(pos).reshape(-1, 3)


def ionizaciones_en_blancos(posiciones, centros, diametro_nm=3.4):
    """Número de ionizaciones dentro de cada blanco esférico.

    Modelo de «blanco del tamaño del ADN» (Goodhead): una esfera de 3,4 nm
    contiene unos 10 pares de bases. `centros` es un array (m, 3) [nm];
    devuelve un array de m enteros.
    """
    if len(posiciones) == 0:
        return np.zeros(len(centros), dtype=int)
    arbol = cKDTree(posiciones)
    return arbol.query_ball_point(centros, r=diametro_nm / 2, return_length=True)


def centros_cerca_del_eje(n, z_min, z_max, r_max, rng):
    """n puntos al azar, uniformes en volumen, en un cilindro de radio r_max
    alrededor del eje z (la dirección del haz), entre z_min y z_max [nm]."""
    r = r_max * np.sqrt(rng.random(n))
    phi = 2 * np.pi * rng.random(n)
    z = rng.uniform(z_min, z_max, n)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def tabla_electrones(nombre):
    """Lee una tabla e-_+_H2O_* del toolkit: columnas (E [eV], valor).

    `nombre` es el resto del nombre de archivo, p. ej. 'ionization_TCS_BEB'.
    Secciones eficaces en 1e-16 cm^2; energías medias en eV.
    """
    return np.loadtxt(DIR_DATOS / f"e-_+_H2O_{nombre}.dat")


# ---------------------------------------------------------------------------
# Funciones de los notebooks guiados (TP3-B y TP5-B)
# ---------------------------------------------------------------------------

def resumen_traza(motor, E_keV, ventana_nm=100.0):
    """Datos de una traza para los notebooks guiados.

    Usa el orden de las colisiones del ion (`motor.collision_log` y sus
    depósitos, uno por colisión). Devuelve un dict con: camino recorrido [nm],
    posiciones de las ionizaciones, número de ionizaciones en los primeros y en
    los últimos `ventana_nm` de camino, distancia media entre ionizaciones
    consecutivas en cada ventana [nm] y energías [keV] a las que el protón
    capturó un electrón (pasó a H⁰).
    """
    dep = [d for d in motor.deposits if d.kind == "proton"]
    pasos = np.array([c[0] for c in motor.collision_log])
    perdidas = np.array([c[1] for c in motor.collision_log])
    s = np.cumsum(pasos)                                   # camino al llegar a cada colisión
    E_antes = E_keV - np.concatenate([[0.0], np.cumsum(perdidas)[:-1]]) / 1e3
    es_ion = np.array([d.channel == Channel.IONIZATION for d in dep])
    es_cap = np.array([d.channel == Channel.CAPTURE for d in dep])
    L = s[-1] if len(s) else 0.0

    def ventana(m):
        si = s[es_ion & m]
        return len(si), (np.diff(si).mean() if len(si) > 1 else np.nan)

    n_ini, d_ini = ventana(s <= ventana_nm)
    n_fin, d_fin = ventana(s >= L - ventana_nm)
    return {"camino_nm": L,
            "profundidad_nm": dep[-1].position[2] if dep else 0.0,
            "posiciones": np.array([d.position for d in dep]).reshape(-1, 3),
            "es_ionizacion": es_ion,
            "carga": np.array([d.q for d in dep]),
            "n_ionizaciones": int(es_ion.sum()),
            "n_inicio": n_ini, "dist_inicio_nm": d_ini,
            "n_final": n_fin, "dist_final_nm": d_fin,
            "E_captura_keV": E_antes[es_cap]}


def graficar_traza(res, titulo="", ax=None):
    """Proyección (z, x) de la traza: ionizaciones en rojo, otras colisiones
    en gris; en naranja las colisiones hechas como H⁰. Rótulos en español."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 3.5))
    p = res["posiciones"]
    h0 = res["carga"] == 0
    ax.plot(p[:, 2], p[:, 0], color="0.75", lw=0.6, zorder=1)
    otras = ~res["es_ionizacion"] & ~h0
    ax.scatter(p[otras, 2], p[otras, 0], s=4, color="0.55", label="otras colisiones (H⁺)")
    ax.scatter(p[res["es_ionizacion"] & ~h0, 2], p[res["es_ionizacion"] & ~h0, 0], s=6,
               color="#c8364b", label="ionizaciones (H⁺)")
    if h0.any():
        ax.scatter(p[h0, 2], p[h0, 0], s=10, color="#f5a524", label="colisiones como H⁰")
    ax.set_xlabel("Profundidad z (nm)")
    ax.set_ylabel("x (nm)")
    ax.set_title(titulo)
    ax.legend(loc="best", fontsize=9, markerscale=2)
    return ax


def curva_bragg(datos, E_keV, n, z_max_nm, n_bins=150, seed=0):
    """Simula n protones y devuelve (bordes [nm], energía depositada por
    intervalo de profundidad [eV], profundidades finales [nm])."""
    bordes = np.linspace(0.0, z_max_nm, n_bins + 1)
    dosis = np.zeros(n_bins)
    finales = []
    motor = nuevo_motor(datos, seed=seed)
    for _ in range(n):
        p = proton(E_keV)
        motor.transport(p)
        pos = np.array([d.position[2] for d in motor.deposits])
        en = np.array([d.energy for d in motor.deposits])
        dosis += np.histogram(pos, bins=bordes, weights=en)[0]
        finales.append(p.position[2])
        motor.deposits.clear()
    return bordes, dosis, np.array(finales)


def bragg_bortfeld(z_cm, R0_cm, sigma_E_rel=0.01, p=1.77, alfa=0.0022, beta=0.012, gamma=0.6):
    """Curva de Bragg analítica de Bortfeld (Med. Phys. 24, 2024, 1997) en agua.

    Dosis por unidad de fluencia [MeV g⁻¹ cm²] en las profundidades z_cm, para
    un haz de alcance R0_cm con dispersión de energía relativa sigma_E_rel.
    Sin la cola de fluencia (epsilon = 0). Lejos del pico (zeta > 10) se usa la
    forma sin dispersión de alcance de la misma referencia.
    """
    from scipy.special import gamma as Gamma, pbdv
    z = np.atleast_1d(np.asarray(z_cm, float))
    E0 = (R0_cm / alfa) ** (1 / p)
    sigma = np.sqrt((0.012 * R0_cm ** 0.935) ** 2
                    + (sigma_E_rel * E0) ** 2 * alfa ** 2 * p ** 2 * E0 ** (2 * p - 2))
    zeta = (R0_cm - z) / sigma
    D = np.zeros_like(z)
    cerca = zeta <= 10
    zc = zeta[cerca]
    pref = (np.exp(-zc ** 2 / 4) * sigma ** (1 / p) * Gamma(1 / p)
            / (np.sqrt(2 * np.pi) * p * alfa ** (1 / p) * (1 + beta * R0_cm)))
    D[cerca] = pref * (pbdv(-1 / p, -zc)[0] / sigma
                       + (beta / p + gamma * beta) * pbdv(-1 / p - 1, -zc)[0])
    r = R0_cm - z[~cerca]
    D[~cerca] = ((r ** (1 / p - 1) + (beta + gamma * beta * p) * r ** (1 / p))
                 / (p * alfa ** (1 / p) * (1 + beta * R0_cm)))
    return D


def segmento_traza(datos, E_keV, longitud_nm=100.0, seed=0):
    """Primeros `longitud_nm` de la traza de un protón que entra según +z.

    Para no simular la traza completa (10 MeV ≈ 10 s), el motor se detiene
    cuando el protón perdió 30 keV, mucho más de lo que pierde en 100 nm a
    cualquier energía del rango del código. Devuelve (posiciones de las
    ionizaciones con z < longitud_nm [nm], TEL del segmento [keV/µm]).
    """
    E_eV = E_keV * 1e3
    motor = TrackEngine(datos["medio"], datos["xs"], edep_tables=datos["edep"],
                        eemean_tables=datos["eemean"], e_cut=max(10.0, E_eV - 30e3), seed=seed)
    motor.transport(proton(E_keV))
    s = np.cumsum([c[0] for c in motor.collision_log])
    dE = np.array([c[1] for c in motor.collision_log])
    TEL = dE[s <= longitud_nm].sum() / longitud_nm          # eV/nm = keV/µm
    ion = ionizaciones(motor)
    return ion[ion[:, 2] < longitud_nm], TEL


def ajuste_lq(D, S):
    """Ajuste LQ por cuadrados mínimos de -ln S = α D + β D² (sin término
    independiente). Devuelve (α [Gy⁻¹], β [Gy⁻²])."""
    D, S = np.asarray(D, float), np.asarray(S, float)
    A = np.column_stack([D, D ** 2])
    (alfa, beta), *_ = np.linalg.lstsq(A, -np.log(S), rcond=None)
    return alfa, beta


def dosis_para_supervivencia(alfa, beta, S=0.1):
    """Dosis [Gy] que da la supervivencia S en el modelo LQ."""
    L = -np.log(S)
    if abs(beta) < 1e-12:
        return L / alfa
    return (-alfa + np.sqrt(alfa ** 2 + 4 * beta * L)) / (2 * beta)


# ---------------------------------------------------------------------------
# Rótulos en español para los gráficos de ion_track_mc.analysis.visualization
# ---------------------------------------------------------------------------

_TRADUCCION = {
    "Energy [eV]": "Energía [eV]",
    "Depth [nm]": "Profundidad [nm]",
    "Dose [Gy]": "Dosis [Gy]",
    "Bragg curve": "Curva de Bragg",
    "Radial distance [nm]": "Distancia radial [nm]",
    "Radial dose": "Dosis radial",
    "Depth z [nm]": "Profundidad z [nm]",
    "Transverse x [nm]": "Transversal x [nm]",
    "Simulated track": "Traza simulada",
    "Analytic": "Analítico",
    "Monte Carlo (emergent)": "Monte Carlo (emergente)",
    "e- (deposit)": "e⁻ (depósito)",
    "H1+": "H⁺",
    "H0": "H⁰",
    "ionization": "ionización",
    "excitation": "excitación",
    "capture": "captura",
    "loss": "pérdida",
    "elastic": "elástica",
}


def _traducir(texto):
    import re
    if texto in _TRADUCCION:
        return _TRADUCCION[texto]
    m = re.fullmatch(r"(Cross sections|Mean free path), q=(\d)", texto)
    if m:
        return {"Cross sections": "Secciones eficaces", "Mean free path": "Recorrido libre medio"}[m.group(1)] + f", q = {m.group(2)}"
    return texto


def en_espanol(fig=None):
    """Traduce al español los rótulos de los ejes, los títulos y las leyendas de
    los gráficos de ion_track_mc (plot_track, plot_bragg_curve, plot_cross_sections,
    plot_mean_free_path, plot_radial_dose). Llamarla después de graficar:

        plot_track(entrada, motor); en_espanol(); plt.show()

    Sin argumento actúa sobre la figura actual.
    """
    import matplotlib.pyplot as plt
    fig = fig or plt.gcf()
    for ax in fig.get_axes():
        ax.set_xlabel(_traducir(ax.get_xlabel()))
        ax.set_ylabel(_traducir(ax.get_ylabel()))
        ax.set_title(_traducir(ax.get_title()))
        ley = ax.get_legend()
        if ley is not None:
            for t in ley.get_texts():
                t.set_text(_traducir(t.get_text()))
    return fig
