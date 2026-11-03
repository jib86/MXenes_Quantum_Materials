"""
================================================================
 plot.py — Gráficas a partir de lo que guardó run.py
================================================================
    python plot.py                      # usa output_dir de parametros.py
    python plot.py resultados_onsager_FM

Genera en la carpeta de resultados:
  observables.png    U_L, chi, C y <|M|> vs T, con barras de error (pasada fina)
  extrapolacion.png  T_x(L) vs L_eff^-x y la recta hasta L -> infinito
  histograma_distancias.png  (solo para el material de parametros.py)
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import parametros
from analisis import tabla_observables
from red import distancias_celda


def cargar(carpeta, etiqueta="fino"):
    datos = {}
    for f in glob.glob(os.path.join(carpeta, f"{etiqueta}_L*.npz")):
        d = np.load(f)
        datos[int(d["L1"])] = {k: d[k] for k in d.files}
    return dict(sorted(datos.items()))


def graficar_observables(carpeta, k_B=parametros.k_B, n_bloques=20):
    datos = cargar(carpeta)
    if not datos:
        print(f"No hay datos 'fino_L*.npz' en {carpeta}")
        return
    resumen = None
    ruta = os.path.join(carpeta, "resumen.json")
    if os.path.exists(ruta):
        with open(ruta) as f:
            resumen = json.load(f)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    paneles = [("U", r"$U_L$"), ("chi", r"$\chi$ (1/meV)"),
               ("C", r"$C/k_B$ por sitio"), ("Mabs", r"$\langle|M_{ord}|\rangle$")]
    for L, d in datos.items():
        tab = tabla_observables(d["T"], d["E"], d["M"], int(d["N"]), k_B, n_bloques)
        o = np.argsort(d["T"])
        for a, (clave, nombre) in zip(ax.ravel(), paneles):
            v, e = tab[clave]
            a.errorbar(d["T"][o], v[o], yerr=e[o], marker="o", ms=3, capsize=2,
                       label=f"L = {L}")
            a.set_ylabel(nombre)
    if resumen and resumen.get("exito"):
        for a in ax.ravel():
            a.axvline(resumen["Tc_K"], color="k", ls="--", lw=1)
        ax[0, 0].set_title(f"{resumen['nombre_Tc']} = {resumen['Tc_K']:.2f} ± "
                           f"{resumen['sTc_K']:.2f} K", fontsize=10)
    for a in ax[1]:
        a.set_xlabel("T (K)")
    ax[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(carpeta, "observables.png"), dpi=150)
    plt.close(fig)

    if resumen and resumen.get("exito") and len(resumen["cruces"]) >= 2:
        cr = resumen["cruces"]
        x_exp = resumen["exponente_fss"]
        Le = np.array([c["L_eff"] for c in cr])
        Tx = np.array([c["T"] for c in cr])
        sT = np.array([c["sT"] for c in cr])
        x = Le ** (-x_exp)
        fig, a = plt.subplots(figsize=(5.5, 4))
        a.errorbar(x, Tx, yerr=sT, fmt="o", capsize=3, label="cruces")
        c1 = np.polyfit(x, Tx, 1, w=1 / sT)
        xx = np.linspace(0, x.max() * 1.1, 50)
        a.plot(xx, np.polyval(c1, xx), "k--", lw=1)
        a.errorbar([0], [resumen["Tc_K"]], yerr=[resumen["sTc_K"]], fmt="s", color="C3",
                   capsize=3, label=r"$T_c(\infty)$")
        a.set_xlabel(rf"$L_{{eff}}^{{-{x_exp:g}}}$")
        a.set_ylabel(r"$T_\times$ (K)")
        a.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(carpeta, "extrapolacion.png"), dpi=150)
        plt.close(fig)
    print(f"Gráficas guardadas en {carpeta}/")


def graficar_histograma(carpeta, d_max=7.0):
    ds = distancias_celda(parametros.cell, parametros.posiciones, d_max)
    fig, a = plt.subplots(figsize=(7, 4))
    a.hist(ds, bins=300)
    if parametros.modo_J == "anillos":
        for an in parametros.anillos:
            for d in np.atleast_1d(an["d"]):
                a.axvspan(d - parametros.tol, d + parametros.tol, color="C1", alpha=0.3)
    a.set_xlabel("Distancia (Å)")
    a.set_ylabel("Conteo (pares k -> k'+R)")
    a.set_title("Distancias entre átomos magnéticos; en naranja las ventanas de cada J",
                fontsize=9)
    fig.tight_layout()
    os.makedirs(carpeta, exist_ok=True)
    fig.savefig(os.path.join(carpeta, "histograma_distancias.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    carpeta = sys.argv[1] if len(sys.argv) > 1 else parametros.output_dir
    if carpeta == parametros.output_dir:
        graficar_histograma(carpeta)
    graficar_observables(carpeta)
