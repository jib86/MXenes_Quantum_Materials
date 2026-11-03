"""
================================================================
 estado_fundamental.py — Encontrar el orden magnético a T = 0
================================================================

Por qué hace falta
------------------
El parámetro de orden M_ord = (1/N) sum_i eps_i S_i solo mide algo útil si
eps es el estado fundamental REAL del Hamiltoniano con los J que se usan.
Si eps se escribe a mano y no es el estado fundamental (o si el estado
fundamental es degenerado), el cumulante de Binder no mide la transición.

Qué hace
--------
Enumera EXACTAMENTE todas las configuraciones periódicas con celdas
magnéticas pequeñas (1x1, 2x1, 1x2, 2x2, 3x1, ... celdas unidad) y se queda
con la de menor energía por sitio. Luego:
  * cuenta cuántos estados fundamentales distintos hay (descontando
    traslaciones y la inversión global S -> -S); si hay más de uno, avisa;
  * calcula el campo de canje sobre cada subred en ese estado; si alguna
    subred siente campo NULO, avisa (subred desacoplada: no puede ordenar);
  * clasifica el orden (FM / AFM / ferri) para decir si la temperatura de
    transición es de Curie o de Néel;
  * da una estimación de campo medio de T_c a partir del máximo autovalor
    de la matriz de acoplamientos.

Limitación: un estado fundamental con periodo mayor que los enumerados no
se detecta (el programa lo dice). Se puede ampliar con `max_espines`.
"""

import numpy as np
from numba import njit

from red import matriz_acoplamiento


@njit(cache=True)
def _enumerar_gray(A, tol_rel, max_guardar):
    """
    Recorre las 2^(N-1) configuraciones (espín 0 fijo en +1) en código Gray:
    entre una configuración y la siguiente cambia un solo espín, así que la
    energía se actualiza en O(N). Devuelve la energía mínima total y hasta
    `max_guardar` configuraciones (como máscaras de bits) que la alcanzan.
    """
    N = A.shape[0]
    S = np.ones(N)
    h = A @ S
    E = -0.5 * np.dot(S, h)
    escala = np.sum(np.abs(A)) + 1e-300
    tol = tol_rel * escala
    E_min = E
    guard = np.zeros(max_guardar, dtype=np.int64)
    n_min = 1
    mascara = np.int64(0)
    guard[0] = mascara
    for it in range(1, 2 ** (N - 1)):
        # bit que cambia = número de ceros a la derecha de it, +1 (espín 0 fijo)
        t = it
        p = 0
        while (t & 1) == 0:
            t >>= 1
            p += 1
        p += 1
        sp = S[p]
        dE = 2.0 * sp * (h[p] - A[p, p] * sp)
        E += dE
        for j in range(N):
            h[j] -= 2.0 * A[j, p] * sp
        S[p] = -sp
        mascara ^= (np.int64(1) << p)
        if E < E_min - tol:
            E_min = E
            n_min = 1
            guard[0] = mascara
        elif abs(E - E_min) <= tol:
            if n_min < max_guardar:
                guard[n_min] = mascara
            n_min += 1
    return E_min, guard[:min(n_min, max_guardar)], n_min


def _mascara_a_espines(m, N):
    return np.array([-1 if (m >> i) & 1 else 1 for i in range(N)], dtype=np.int64)


def _canonica(S, m1, m2, nat):
    """Representante único de S bajo traslaciones de celda y S -> -S."""
    T = S.reshape(m1, m2, nat)
    mejor = None
    for t1 in range(m1):
        for t2 in range(m2):
            R = np.roll(np.roll(T, t1, axis=0), t2, axis=1).ravel()
            for signo in (1, -1):
                c = tuple(signo * R)
                if mejor is None or c < mejor:
                    mejor = c
    return mejor


def periodos_candidatos(nat, max_espines, m_max=4):
    per = [(m1, m2) for m1 in range(1, m_max + 1) for m2 in range(1, m_max + 1)
           if m1 * m2 * nat <= max_espines]
    return sorted(per, key=lambda p: (p[0] * p[1], p))


def buscar_estado_fundamental(enl, nat, k_B, max_espines=24, epsilon_usuario=None,
                              periodo_usuario=None):
    """
    Devuelve un diccionario con:
      periodo (m1, m2), patron (m1*m2*nat signos), e_sitio (meV/sitio),
      n_clases (estados fundamentales distintos), tipo ("FM"/"AFM"/"ferri"),
      nombre_Tc, T_MF (K), avisos (lista de textos), campo_subred.
    """
    avisos = []
    resultados = []
    for (m1, m2) in periodos_candidatos(nat, max_espines):
        A = matriz_acoplamiento(m1, m2, nat, enl)
        N = A.shape[0]
        E_min, mascaras, n_min = _enumerar_gray(A, 1e-10, 256)
        resultados.append(dict(periodo=(m1, m2), A=A, e=E_min / N,
                               mascaras=mascaras, n_min=n_min, N=N))

    escala = np.max(np.abs(enl["J"]))
    e_glob = min(r["e"] for r in resultados)
    tol_e = 1e-8 * escala
    elegido = next(r for r in resultados if r["e"] <= e_glob + tol_e)
    (m1, m2), N = elegido["periodo"], elegido["N"]

    # Clases distintas de estados fundamentales. Se cuentan en TODAS las
    # celdas magnéticas que alcanzan la energía mínima (p. ej. en la red
    # triangular AFM las franjas en distintas orientaciones aparecen en
    # celdas distintas) y se toma el máximo.
    n_clases = 1
    for r in resultados:
        if r["e"] <= e_glob + tol_e:
            a1, a2 = r["periodo"]
            cl = {_canonica(_mascara_a_espines(m, r["N"]), a1, a2, nat) for m in r["mascaras"]}
            n_clases = max(n_clases, len(cl))
            if r["n_min"] > len(r["mascaras"]):
                avisos.append("Más de 256 configuraciones degeneradas: degeneración masiva.")
    clases = {_canonica(_mascara_a_espines(m, N), m1, m2, nat) for m in elegido["mascaras"]}
    patron = np.array(sorted(clases)[-1], dtype=np.int64)  # determinista
    if patron[0] < 0:
        patron = -patron

    if n_clases > 1:
        avisos.append(
            f"El estado fundamental es DEGENERADO: {n_clases} configuraciones "
            "distintas (sin contar traslaciones ni S->-S) tienen la misma energía. "
            "Puede haber frustración o subredes desacopladas; el parámetro de "
            "orden con un solo patrón eps no es confiable.")

    # Campo de canje sobre cada sitio en el estado fundamental
    A = elegido["A"]
    # (la diagonal de A son enlaces con imágenes de la misma subred en otras
    #  celdas: en una configuración periódica tienen el mismo espín, cuentan)
    h = A @ patron
    campo = (h * patron).reshape(m1 * m2, nat)   # >0: el sitio está "contento"
    campo_subred = campo.mean(axis=0)
    for k in range(nat):
        if np.all(np.abs(campo[:, k]) < 1e-6 * escala):
            avisos.append(
                f"La subred {k} siente campo de canje NULO en el estado "
                "fundamental: está desacoplada del orden y no tiene T_c propia.")

    # Tipo de orden
    s = patron.sum()
    if abs(s) == N:
        tipo, nombre = "FM", "temperatura de Curie"
    elif s == 0:
        tipo, nombre = "AFM", "temperatura de Néel"
    else:
        tipo, nombre = "ferri", "temperatura de Curie (ferrimagnético)"

    # Campo medio: máximo autovalor de J(q) en los q de la celda elegida
    lam = np.linalg.eigvalsh(A).max()
    T_MF = lam / k_B

    # Si el usuario dio un eps, comparar con el encontrado
    if epsilon_usuario is not None:
        eu = np.asarray(epsilon_usuario, dtype=np.int64)
        per_u = (1, 1) if len(eu) == nat else periodo_usuario
        if per_u is None or per_u[0] * per_u[1] * nat != len(eu):
            raise ValueError("epsilon_base no cuadra con nat; indique periodo_epsilon=(m1, m2).")
        Au = matriz_acoplamiento(per_u[0], per_u[1], nat, enl)
        e_u = -0.5 * eu @ Au @ eu / len(eu)
        if e_u > e_glob + tol_e:
            avisos.append(
                f"El epsilon_base dado tiene energía {e_u:.6f} meV/sitio, MAYOR que "
                f"la del estado fundamental encontrado ({e_glob:.6f}). No es el "
                "estado fundamental: se usará el encontrado.")
        else:
            patron, (m1, m2) = eu, tuple(per_u)

    avisos_lim = (f"(Se enumeraron celdas magnéticas de hasta {max_espines} espines; "
                  f"un orden con periodo mayor no se detectaría.)")
    return dict(periodo=(m1, m2), patron=patron, e_sitio=e_glob,
                n_clases=n_clases, tipo=tipo, nombre_Tc=nombre, T_MF=T_MF,
                campo_subred=campo_subred, avisos=avisos, nota=avisos_lim,
                energias={r["periodo"]: r["e"] for r in resultados})


def construir_epsilon(L1, L2, nat, gs):
    """Repite el patrón del estado fundamental sobre la supercelda L1 x L2."""
    m1, m2 = gs["periodo"]
    if L1 % m1 or L2 % m2:
        raise ValueError(
            f"El estado fundamental tiene periodo {m1}x{m2} celdas: L1 y L2 deben "
            f"ser múltiplos de {m1} y {m2} (si no, el orden queda frustrado por "
            "las condiciones periódicas).")
    P = gs["patron"].reshape(m1, m2, nat)
    eps = np.empty(L1 * L2 * nat, dtype=np.int64)
    for n1 in range(L1):
        for n2 in range(L2):
            eps[(n1 * L2 + n2) * nat:(n1 * L2 + n2 + 1) * nat] = P[n1 % m1, n2 % m2]
    return eps


def fraccion_satisfechos(enl, gs):
    """Para cada anillo: fracción de enlaces satisfechos (J S_i S_j > 0)."""
    m1, m2 = gs["periodo"]
    nat = int(max(enl["k"].max(), enl["kp"].max())) + 1
    P = gs["patron"].reshape(m1, m2, nat)
    out = {}
    for a in np.unique(enl["anillo"]):
        sel = np.where(enl["anillo"] == a)[0]
        sat = []
        for c1 in range(m1):
            for c2 in range(m2):
                for e in sel:
                    si = P[c1, c2, enl["k"][e]]
                    sj = P[(c1 + enl["n1"][e]) % m1, (c2 + enl["n2"][e]) % m2, enl["kp"][e]]
                    sat.append(enl["J"][e] * si * sj > 0)
        out[int(a)] = float(np.mean(sat))
    return out


def resumen_estado_fundamental(gs):
    l = [f"Estado fundamental: periodo {gs['periodo'][0]}x{gs['periodo'][1]} celdas, "
         f"orden {gs['tipo']}  ->  la transición es una {gs['nombre_Tc']}.",
         f"  patrón eps = {gs['patron'].tolist()}",
         f"  energía = {gs['e_sitio']:.6f} meV/sitio,  clases degeneradas = {gs['n_clases']}",
         f"  campo de canje medio por subred (meV) = {np.round(gs['campo_subred'], 4).tolist()}",
         f"  T_c de campo medio ~ {gs['T_MF']:.1f} K (sobreestima; en 2D el real es menor)",
         "  " + gs["nota"]]
    for a in gs["avisos"]:
        l.append("  AVISO: " + a)
    return "\n".join(l)
