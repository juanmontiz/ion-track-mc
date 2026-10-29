"""Verificación de la instalación para los trabajos prácticos del curso de Radiobiología.

Uso:  python verificar_instalacion.py
(en la carpeta donde están tp_utils.py, mc_fotones.py y esqueleto_tp4a.py).
Imprime «TODO EN ORDEN» si se puede hacer cualquier TP; si no, dice qué falta.
Tarda unos segundos.
"""
import importlib
import sys
import time

problemas = []
print(f"Python {sys.version.split()[0]}")
if sys.version_info < (3, 10):
    problemas.append("hace falta Python 3.10 o más nuevo")

for nombre in ("numpy", "scipy", "matplotlib", "ion_track_mc"):
    try:
        modulo = importlib.import_module(nombre)
        print(f"  {nombre:18s} {getattr(modulo, '__version__', 'instalado')}")
    except ImportError:
        problemas.append(f"falta el paquete {nombre}")
        print(f"  {nombre:18s} NO ENCONTRADO")

for archivo in ("tp_utils", "mc_fotones", "esqueleto_tp4a"):
    try:
        importlib.import_module(archivo)
        print(f"  {archivo + '.py':18s} encontrado")
    except ImportError:
        problemas.append(f"falta el archivo {archivo}.py en esta carpeta")
        print(f"  {archivo + '.py':18s} NO ENCONTRADO")

if not problemas:
    import tp_utils as tp
    t0 = time.time()
    datos = tp.cargar_agua()
    _, _, final = tp.simular_proton(datos, 50, seed=0)
    z = final.position[2] / 1e3
    print(f"  prueba: un protón de 50 keV se detuvo a {z:.2f} µm ({time.time() - t0:.1f} s)")
    if not 0.4 < z < 0.7:
        problemas.append("la simulación de prueba dio un resultado inesperado")

print()
if problemas:
    print("HAY PROBLEMAS:")
    for p in problemas:
        print("  -", p)
    print("Ver «Problemas frecuentes» en las instrucciones de instalación.")
else:
    print("TODO EN ORDEN: la instalación sirve para todos los TP.")
