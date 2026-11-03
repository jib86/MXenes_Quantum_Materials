"""
================================================================
 run.py — Flujo completo: de la red y los J hasta T_c en Kelvin
================================================================

    python run.py                 # usa parametros.py

Pasos:
  0. Lista de enlaces + verificaciones (coordinación por subred, L mínimo).
  1. Estado fundamental exacto -> patrón eps, tipo de orden, T_c de campo medio.
  2. Pasada GRUESA en T (poca estadística) para ubicar T_c.
  3. Pasada FINA alrededor de T_c con toda la estadística.
  4. Observables con errores (jackknife), tau_int, cruces de Binder por pares
     de L, extrapolación a L -> infinito.
  5. Diagnósticos: U*, |M| a baja T, tau_int, criterio de anisotropía.
Guarda las series en output_dir/ y un resumen en output_dir/resumen.json.
"""

import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from red import (enlaces_por_anillos, enlaces_explicitos, resumen_enlaces,
                 construir_libreta, L_minimo)
from estado_fundamental import (buscar_estado_fundamental, construir_epsilon,
                                resumen_estado_fundamental, fraccion_satisfechos)
from montecarlo import barrido_temperaturas
from analisis import (tabla_observables, cruces_binder, extrapolar_Tc,
                      refinar_cruce, U_ESTRELLA_ISING_2D)

CLAVES = ["k_B", "cell", "posiciones", "modo_J", "anillos", "tol", "enlaces",
          "factor_J", "epsilon_base", "periodo_epsilon", "MAE_meV", "L_list",
          "aspecto", "T_min", "T_max", "N_temps_grueso", "N_temps_fino",
          "ventana_fina", "fraccion_grueso", "algoritmo", "N_term", "N_skip",
          "N_meas", "n_bloques", "exponente_fss", "seed", "n_procesos",
          "output_dir", "forzar"]


def config_desde_modulo(mod):
    return {c: getattr(mod, c, None) for c in CLAVES}


# ============================ PREPARACIÓN ===================================

def preparar(cfg, verbose=True):
    pos = np.asarray(cfg["posiciones"], float)
    nat = len(pos)
    if cfg["modo_J"] == "anillos":
        enl, avisos = enlaces_por_anillos(cfg["cell"], pos, cfg["anillos"],
                                          cfg["tol"], cfg["factor_J"])
    elif cfg["modo_J"] == "enlaces":
        enl, avisos = enlaces_explicitos(cfg["cell"], pos, cfg["enlaces"], cfg["factor_J"])
    else:
        raise ValueError("modo_J debe ser 'anillos' o 'enlaces'")

    gs = buscar_estado_fundamental(enl, nat, cfg["k_B"],
                                   epsilon_usuario=cfg["epsilon_base"],
                                   periodo_usuario=cfg["periodo_epsilon"])
    sat = fraccion_satisfechos(enl, gs)

    # Tamaños: mínimos y múltiplos del periodo magnético
    m1, m2 = gs["periodo"]
    Lmin = L_minimo(enl)
    tamanos = []
    for L in cfg["L_list"]:
        L1, L2 = int(L), int(round(cfg["aspecto"] * L))
        if L1 < Lmin or L2 < Lmin:
            raise ValueError(f"L={L} es menor que el mínimo {Lmin}.")
        if L1 % m1 or L2 % m2:
            raise ValueError(f"L={L}: con un orden de periodo {m1}x{m2}, L1={L1} y "
                             f"L2={L2} deben ser múltiplos de {m1} y {m2}.")
        tamanos.append((L1, L2))
    if len(tamanos) < 3:
        avisos.append("Con menos de 3 tamaños no hay extrapolación con error confiable.")

    if verbose:
        print("=" * 72)
        print(resumen_enlaces(enl, nat))
        for a in avisos:
            print("  AVISO:", a)
        print()
        print(resumen_estado_fundamental(gs))
        print("  enlaces satisfechos por anillo:",
              {f"anillo {a+1}": round(f, 3) for a, f in sat.items()})
        print("=" * 72)
    return dict(enl=enl, gs=gs, nat=nat, tamanos=tamanos, avisos=avisos + gs["avisos"])


# ============================ SIMULACIÓN ====================================

def _correr_un_L(args):
    (L1, L2), nat, enl, gs, temps, k_B, n_term, n_skip, n_meas, alg, seed = args
    vec, Jv, nvec = construir_libreta(L1, L2, nat, enl)
    eps = construir_epsilon(L1, L2, nat, gs)
    t0 = time.time()
    E, M, _ = barrido_temperaturas(vec, Jv, nvec, eps, temps, k_B,
                                   n_term, n_skip, n_meas, alg, seed=seed)
    return (L1, L2), E, M, time.time() - t0


def pasada(cfg, prep, temps, n_meas, etiqueta, verbose=True):
    """Simula todos los L en las temperaturas dadas (de mayor a menor)."""
    temps = np.sort(np.asarray(temps, float))[::-1]
    trabajos = [((L1, L2), prep["nat"], prep["enl"], prep["gs"], temps, cfg["k_B"],
                 cfg["N_term"], cfg["N_skip"], n_meas, cfg["algoritmo"],
                 cfg["seed"] + 1000 * i + (0 if etiqueta == "grueso" else 7))
                for i, (L1, L2) in enumerate(prep["tamanos"])]
    if cfg["n_procesos"] > 1:
        with ProcessPoolExecutor(cfg["n_procesos"]) as ex:
            res = list(ex.map(_correr_un_L, trabajos))
    else:
        res = [_correr_un_L(t) for t in trabajos]
    datos = {}
    os.makedirs(cfg["output_dir"], exist_ok=True)
    for (L1, L2), E, M, dt in res:
        N = L1 * L2 * prep["nat"]
        datos[L1] = dict(T=temps, E=E, M=M, N=N, L2=L2)
        np.savez_compressed(os.path.join(cfg["output_dir"], f"{etiqueta}_L{L1}.npz"),
                            T=temps, E=E, M=M, N=N, L1=L1, L2=L2)
        if verbose:
            print(f"  [{etiqueta}] L={L1}x{L2} (N={N}) listo en {dt:.1f} s")
    return datos


# ============================ ANÁLISIS ======================================

def analizar(cfg, datos, reponderar=True):
    """
    Observables con errores para cada L y cruces de Binder por pares
    consecutivos. Los cruces robustos se refinan por reponderación de
    histogramas (curva continua, error jackknife).
    """
    k_B = cfg["k_B"]
    Ls = sorted(datos)
    tablas = {L: tabla_observables(datos[L]["T"], datos[L]["E"], datos[L]["M"],
                                   datos[L]["N"], k_B, cfg["n_bloques"]) for L in Ls}
    pares = []
    for L1, L2 in zip(Ls[:-1], Ls[1:]):
        T = datos[L1]["T"]
        cr = cruces_binder(T, *tablas[L1]["U"], *tablas[L2]["U"], rng=cfg["seed"])
        if reponderar and cr and cr[0]["robusto"]:
            Ts = np.sort(T)
            j = np.searchsorted(Ts, cr[0]["T"])
            a, b = Ts[max(0, j - 2)], Ts[min(len(Ts) - 1, j + 1)]
            ref = refinar_cruce(datos[L1], datos[L2], k_B, a, b, cfg["n_bloques"])
            if ref is not None:
                cr[0].update(T=ref[0], sT=ref[1], U=ref[2], sU=ref[3], reponderado=True)
        pares.append(dict(L1=L1, L2=L2, cruces=cr))
    return tablas, pares


def estimar_Tc_rapido(datos, tablas, pares):
    """Para la pasada gruesa: cruces robustos o, si no hay, el pico de chi."""
    robustos = [p["cruces"][0]["T"] for p in pares if p["cruces"] and p["cruces"][0]["robusto"]]
    if robustos:
        return float(np.median(robustos)), "cruces de Binder"
    Lmax = max(datos)
    T = datos[Lmax]["T"]
    i = int(np.nanargmax(tablas[Lmax]["chi"][0]))
    if i in (0, len(T) - 1):
        return None, "pico de chi en el borde del rango"
    return float(T[i]), "pico de la susceptibilidad (sin cruces robustos)"


# ============================ FLUJO PRINCIPAL ===============================

def ejecutar(cfg, verbose=True):
    t_inicio = time.time()
    prep = preparar(cfg, verbose)
    gs = prep["gs"]
    avisos = list(prep["avisos"])
    problema = gs["n_clases"] > 1 or any("NULO" in a for a in gs["avisos"])
    if problema and not cfg.get("forzar"):
        msg = ("El estado fundamental es degenerado o tiene subredes desacopladas: con "
               "estos J no hay un único orden que medir y el T_c del cumulante no "
               "tendría sentido. Revise los J (o ponga forzar = True en parametros.py "
               "para simular de todos modos).")
        print("\nDETENIDO:", msg)
        return dict(exito=False, mensaje=msg, avisos=avisos)

    # --- Rango de temperaturas ---
    T_min = cfg["T_min"] if cfg["T_min"] is not None else 0.25 * gs["T_MF"]
    T_max = cfg["T_max"] if cfg["T_max"] is not None else 1.05 * gs["T_MF"]
    if verbose:
        print(f"\nPasada gruesa: {cfg['N_temps_grueso']} temperaturas entre "
              f"{T_min:.1f} y {T_max:.1f} K")
    n_g = max(2000, int(cfg["fraccion_grueso"] * cfg["N_meas"]))
    datos_g = pasada(cfg, prep, np.linspace(T_min, T_max, cfg["N_temps_grueso"]),
                     n_g, "grueso", verbose)
    tab_g, pares_g = analizar(cfg, datos_g)
    Tc_est, metodo = estimar_Tc_rapido(datos_g, tab_g, pares_g)
    if Tc_est is None:
        msg = (f"No se pudo ubicar T_c en [{T_min:.1f}, {T_max:.1f}] K ({metodo}). "
               "Amplíe T_min/T_max o revise el estado fundamental.")
        print("\nERROR:", msg)
        return dict(exito=False, mensaje=msg, avisos=avisos)
    if verbose:
        print(f"  T_c estimada (gruesa) ~ {Tc_est:.2f} K  [{metodo}]")

    # --- Pasadas intermedias: acercarse a T_c con poca estadística ---
    # Cada pasada reduce la ventana alrededor de los cruces hasta llegar a
    # T_c*(1 ± ventana_fina). Luego una pasada FINAL con toda la estadística.
    v = cfg["ventana_fina"]
    nT = cfg["N_temps_fino"]
    dT = (T_max - T_min) / max(1, cfg["N_temps_grueso"] - 1)
    semi = 1.5 * dT
    for intento in range(4):
        if semi <= v * Tc_est * 1.01:
            break
        temps = np.linspace(max(Tc_est - semi, 0.05 * Tc_est), Tc_est + semi, nT)
        if verbose:
            print(f"\nPasada intermedia: {temps[0]:.2f} - {temps[-1]:.2f} K")
        d_i = pasada(cfg, prep, temps, n_g, "intermedio", verbose)
        t_i, p_i = analizar(cfg, d_i, reponderar=False)
        rob = [p["cruces"][0]["T"] for p in p_i if p["cruces"] and p["cruces"][0]["robusto"]]
        paso = temps[1] - temps[0]
        if rob:
            Tc_est = float(np.median(rob))
            semi = max(v * Tc_est, 1.5 * (max(rob) - min(rob)), 2 * paso)
        else:
            nuevo, _ = estimar_Tc_rapido(d_i, t_i, p_i)
            Tc_est = nuevo if nuevo is not None else Tc_est
            semi *= 1.5
        if verbose:
            print(f"  T_c estimada ~ {Tc_est:.3f} K")

    # --- Pasada final (toda la estadística), con re-centrado si hace falta ---
    semi = max(semi, v * Tc_est)
    for intento in range(3):
        temps_f = np.linspace(max(Tc_est - semi, 0.05 * Tc_est), Tc_est + semi, nT)
        if verbose:
            print(f"\nPasada final: {nT} temperaturas entre {temps_f[0]:.2f} y "
                  f"{temps_f[-1]:.2f} K, N_meas = {cfg['N_meas']}")
        datos = pasada(cfg, prep, temps_f, cfg["N_meas"], "fino", verbose)
        tablas, pares = analizar(cfg, datos)
        rob = [p["cruces"][0]["T"] for p in pares
               if p["cruces"] and p["cruces"][0]["robusto"]]
        interior = all(temps_f[1] <= t <= temps_f[-2] for t in rob)
        if (len(rob) == len(pares) and interior) or intento == 2:
            break
        if rob:
            Tc_est = float(np.median(rob))
            if interior:
                semi *= 1.5
        else:
            nuevo, _ = estimar_Tc_rapido(datos, tablas, pares)
            Tc_est = nuevo if nuevo is not None else Tc_est
            semi *= 1.5
        if verbose:
            print(f"  -> faltan cruces robustos o están en el borde: re-centrando en "
                  f"{Tc_est:.2f} K (semiancho {semi:.2f} K)")

    # --- Cruces y extrapolación ---
    usados = []
    if verbose:
        print("\nCruces del cumulante de Binder:")
    for p in pares:
        cr = p["cruces"]
        if not cr:
            if verbose:
                print(f"  L={p['L1']} / {p['L2']}: sin cruce en la ventana fina")
            continue
        c = cr[0]
        etiqueta = "robusto" if c["robusto"] else "NO robusto (descartado)"
        if verbose:
            print(f"  L={p['L1']} / {p['L2']}: T_x = {c['T']:.3f} ± {c['sT']:.3f} K, "
                  f"U* = {c['U']:.4f} ± {c['sU']:.4f}  [{etiqueta}]"
                  + (f"  (+{len(cr)-1} cruces más)" if len(cr) > 1 else ""))
        if c["robusto"] and np.isfinite(c["sT"]):
            usados.append((np.sqrt(p["L1"] * p["L2"]), c["T"], c["sT"], c["U"], c["sU"]))

    if not usados:
        msg = "Ningún par de tamaños dio un cruce robusto en la pasada fina."
        print("\nERROR:", msg)
        return dict(exito=False, mensaje=msg, avisos=avisos)

    Leff, Tx, sTx, Ux, sUx = map(np.array, zip(*usados))
    x_alt = 3.0 if cfg["exponente_fss"] == 1.0 else 1.0
    ext = extrapolar_Tc(Leff, Tx, sTx, cfg["exponente_fss"])
    ext_alt = extrapolar_Tc(Leff, Tx, sTx, x_alt)
    sistematico = abs(ext["Tc"] - ext_alt["Tc"])

    # --- Diagnósticos ---
    Lmax = max(datos)
    i_frio = int(np.argmin(datos[Lmax]["T"]))
    M_frio = tablas[Lmax]["Mabs"][0][i_frio]
    if M_frio < 0.5:
        avisos.append(f"<|M_ord|> = {M_frio:.2f} a la T más baja (L={Lmax}): el sistema no "
                      "se ordena según eps. Revise el estado fundamental.")
    if not all(np.all(tablas[L]["bloques_ok"]) for L in tablas):
        avisos.append("En algunas temperaturas los bloques del jackknife son cortos "
                      "frente a tau_int: aumente N_meas o N_skip.")
    tau_max = max(float(np.max(np.maximum(tablas[L]["tau_E"], tablas[L]["tau_M2"]))) for L in tablas)
    U_med = float(np.average(Ux, weights=1 / sUx ** 2))
    if verbose:
        print(f"\nExtrapolación ({ext['n']} pares, x = {ext['exponente']:g}): "
              f"T_c(inf) = {ext['Tc']:.3f} ± {ext['sTc']:.3f} K   (chi2_red = {ext['chi2_red']:.2f})")
        print(f"Alternativa (x = {x_alt:g}):            T_c(inf) = {ext_alt['Tc']:.3f} ± "
              f"{ext_alt['sTc']:.3f} K")
        print(f"Diferencia entre ambas formas (sistemático): {sistematico:.3f} K")
        print(f"U* medio = {U_med:.4f}  (Ising 2D, red cuadrada isótropa con PBC: "
              f"{U_ESTRELLA_ISING_2D}; con otra geometría de celda o acoplamientos "
              "anisótropos U* cambia, así que es una guía, no un criterio)")
        print(f"tau_int máximo = {tau_max:.1f} medidas  (N_skip = {cfg['N_skip']})")

    if cfg["MAE_meV"] is not None:
        r = cfg["MAE_meV"] / (cfg["k_B"] * ext["Tc"])
        txt = (f"MAE / (k_B T_c) = {r:.1f}: "
               + ("Ising es razonable (>= 10)." if r >= 10 else
                  "MENOR que 10: Ising sobreestima T_c; use Heisenberg con anisotropía."))
        if r < 10:
            avisos.append(txt)
        if verbose:
            print(txt)
    else:
        avisos.append("MAE_meV no dado: no se verificó si el modelo de Ising está justificado.")

    resumen = dict(
        exito=True, tipo_orden=gs["tipo"], nombre_Tc=gs["nombre_Tc"],
        Tc_K=ext["Tc"], sTc_K=ext["sTc"], sistematico_K=sistematico,
        exponente_fss=ext["exponente"], Tc_alt_K=ext_alt["Tc"], sTc_alt_K=ext_alt["sTc"],
        chi2_red=ext["chi2_red"], T_MF_K=gs["T_MF"], U_estrella=U_med,
        cruces=[dict(L_eff=float(a), T=float(b), sT=float(c), U=float(d), sU=float(e))
                for a, b, c, d, e in usados],
        patron_eps=gs["patron"].tolist(), periodo=list(gs["periodo"]),
        tau_int_max=tau_max, avisos=avisos, minutos=(time.time() - t_inicio) / 60)
    with open(os.path.join(cfg["output_dir"], "resumen.json"), "w") as f:
        json.dump(resumen, f, indent=2, ensure_ascii=False)

    if verbose:
        print("\n" + "=" * 72)
        print(f"RESULTADO: {gs['nombre_Tc']} = {ext['Tc']:.2f} ± {ext['sTc']:.2f} (estad.) "
              f"± {sistematico:.2f} (sist.) K")
        for a in avisos:
            print("  AVISO:", a)
        print("=" * 72)
    return resumen


def main():
    import parametros
    ejecutar(config_desde_modulo(parametros))


if __name__ == "__main__":
    main()
