# Utilidades del curso de Radiobiología

Archivos de los trabajos prácticos del curso *Interacción radiación–materia. Aplicación a la radiobiología* (UNR):

| Archivo | Para qué |
|---|---|
| `tp_utils.py` | Carga de las tablas y arranque de `ion_track_mc`, y funciones de los notebooks guiados |
| `mc_fotones.py` | Monte Carlo de fotones completo (TP4-A, versión básica) |
| `esqueleto_tp4a.py` | Esqueleto del Monte Carlo de fotones para completar (TP4-A, versiones avanzada y media) |
| `verificar_instalacion.py` | Comprueba que la instalación sirve para todos los TP |
| `TP5-B_supervivencia_V79_avanzada.xlsx` | Curvas de supervivencia de células V79: rayos X de 240 kVp y protones de 10, 18 y 28 keV/µm, más α de deuterones y ³He (TP5-B, versión avanzada) |
| `TP5-B_supervivencia_V79_media.xlsx` | Rayos X y protones de 28 keV/µm (TP5-B, versión media) |

Datos de las planillas: puntos digitalizados de la figura 3a y valores de la tabla 1 de M. Folkard *et al.*, «Inactivation of V79 cells by low-energy protons, deuterons and helium-3 ions», *Int. J. Radiat. Biol.* **69**, 729–738 (1996), doi:[10.1080/095530096145472](https://doi.org/10.1080/095530096145472). Cada planilla trae una hoja «Fuente» con los detalles.

Las copias de trabajo están en el material del curso: el código en `Trabajos practicos/codigo/`, que el generador de los notebooks copia acá; las planillas en `Trabajos practicos/datos/supervivencia_V79/`, que copia `generar_planillas.py`. No editar estos archivos.
