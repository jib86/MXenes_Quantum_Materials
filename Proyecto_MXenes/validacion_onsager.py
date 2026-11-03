"""
================================================================
 validacion_onsager.py — Prueba del código contra la solución exacta
================================================================
Red cuadrada, un átomo por celda, J = ±1 meV a primeros vecinos.
Onsager:  k_B T_c / |J| = 2 / ln(1 + sqrt 2) = 2.269185...
          -> con |J| = 1 meV,  T_c = 26.333 K.
Cumulante en el cruce (red cuadrada, PBC): U* ≈ 0.6107.

Se corre el MISMO flujo que para el material (run.ejecutar) con:
  * J = +1 (ferromagnético: temperatura de Curie)
  * J = -1 (antiferromagnético: temperatura de Néel; valida que el
    estado fundamental y el parámetro de orden con eps funcionan)
Criterio de aprobación: |T_c - T_c^exacto| / T_c^exacto < 0.5 %.

    python validacion_onsager.py
"""
import json
import os
import numpy as np

from run import ejecutar
import parametros

TC_EXACTO_ADIM = 2.0 / np.log(1.0 + np.sqrt(2.0))
TOL_REL = 0.005


def config_cuadrada(J, salida):
    return dict(
        k_B=parametros.k_B,
        cell=np.diag([1.0, 1.0, 20.0]), posiciones=np.array([[0.0, 0.0, 0.0]]),
        modo_J="anillos", anillos=[{"d": 1.0, "J": J}], tol=0.05, enlaces=[],
        factor_J=1.0, epsilon_base=None, periodo_epsilon=None, MAE_meV=None,
        L_list=[16, 24, 32, 48], aspecto=1,
        T_min=None, T_max=None,
        N_temps_grueso=12, N_temps_fino=12, ventana_fina=0.04, fraccion_grueso=0.2,
        algoritmo="hibrido", N_term=2000, N_skip=1, N_meas=30000, n_bloques=20,
        exponente_fss=1.0, seed=2026, n_procesos=min(4, os.cpu_count() or 1),
        output_dir=salida, forzar=False)


def main():
    Tc_exacto = TC_EXACTO_ADIM * 1.0 / parametros.k_B
    informe = {}
    for J, nombre in ((+1.0, "FM"), (-1.0, "AFM")):
        print(f"\n######## Validación Onsager {nombre} (J = {J:+.0f} meV) ########")
        r = ejecutar(config_cuadrada(J, f"./resultados_onsager_{nombre}"))
        if not r.get("exito"):
            informe[nombre] = dict(aprobado=False, mensaje=r.get("mensaje"))
            print(f"\n{nombre}: NO APROBADO ({r.get('mensaje')})")
            continue
        err = abs(r["Tc_K"] - Tc_exacto) / Tc_exacto
        ok = err < TOL_REL
        informe[nombre] = dict(Tc_K=float(r["Tc_K"]), sTc_K=float(r["sTc_K"]),
                               Tc_alt_K=float(r["Tc_alt_K"]), sTc_alt_K=float(r["sTc_alt_K"]),
                               Tc_exacto_K=float(Tc_exacto), U_estrella=float(r["U_estrella"]),
                               error_rel=float(err), aprobado=bool(ok))
        print(f"\n{nombre}: T_c = {r['Tc_K']:.3f} ± {r['sTc_K']:.3f} K   exacto = "
              f"{Tc_exacto:.3f} K   error = {100*err:.2f} %  ->  "
              + ("APROBADO" if ok else "NO APROBADO"))
        print(f"{nombre}: U* = {r['U_estrella']:.4f}  (esperado ~ 0.6107)")
    with open("resultados_onsager.json", "w") as f:
        json.dump(informe, f, indent=2)
    return informe


if __name__ == "__main__":
    main()
