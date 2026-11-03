"""
Gráficas de una corrida guardada.

    python plot.py                        # resultados del material (output_dir)
    python plot.py resultados/onsager_FM  # cualquier otra carpeta de resultados
"""
import sys

import parametros
from isingmc.flujo import config_desde_modulo
from isingmc.graficas import graficar_histograma, graficar_observables

if __name__ == "__main__":
    cfg = config_desde_modulo(parametros)
    carpeta = sys.argv[1] if len(sys.argv) > 1 else cfg["output_dir"]
    if carpeta == cfg["output_dir"]:
        graficar_histograma(cfg, carpeta)
    graficar_observables(carpeta, cfg["k_B"], cfg["n_bloques"])
