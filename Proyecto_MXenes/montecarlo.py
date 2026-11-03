"""
================================================================
 montecarlo.py — El motor: Metrópolis, Wolff y medición de series
================================================================

Unidades: J en meV, T en K, k_B en meV/K  ->  beta = 1/(k_B T) en 1/meV.
Hamiltoniano: H = - sum_<ij> J_ij S_i S_j   (ver red.py).

"Un paso Monte Carlo" (la unidad de tiempo de todo el programa) depende del
algoritmo elegido:
  * "metropolis": 1 barrido = N intentos de volteo de un espín.
  * "wolff"     : un número FIJO de clusters de Wolff que en promedio
                  voltea ~N espines (se calibra en la termalización).
  * "hibrido"   : 1 barrido de Metrópolis + esos clusters de Wolff.
                  Es el recomendado: Wolff elimina el "critical slowing down"
                  cerca de T_c y Metrópolis ayuda donde hay acoplamientos
                  frustrados o débiles que los clusters no mueven bien.

Wolff con acoplamientos de signo arbitrario: desde un sitio i del cluster se
agrega el vecino j SOLO si el enlace está satisfecho (J_ij S_i S_j > 0), con
probabilidad p = 1 - exp(-2 beta |J_ij|). Esto cumple balance detallado para
cualquier signo de J (representación de Fortuin–Kasteleyn); con frustración
fuerte pierde eficiencia, por eso existe el modo híbrido.

La función central es `simular_temperatura`: termaliza, mide y devuelve las
SERIES TEMPORALES completas E(t) y M(t). Guardar las series (y no solo los
promedios) es lo que permite calcular el tiempo de autocorrelación y los
errores jackknife en analisis.py.
"""

import numpy as np
from numba import njit

ALGORITMOS = {"metropolis": 0, "wolff": 1, "hibrido": 2}


@njit(cache=True)
def sembrar(seed):
    """El generador de numba es independiente del de numpy: hay que sembrarlo aparte."""
    np.random.seed(seed)


@njit(cache=True)
def campo_local(i, spins, vec, Jv, nvec):
    """h_i = sum_j J_ij S_j  (meV)."""
    h = 0.0
    for c in range(nvec[i]):
        h += Jv[i, c] * spins[vec[i, c]]
    return h


@njit(cache=True)
def barrido_metropolis(spins, vec, Jv, nvec, beta):
    N = spins.size
    for _ in range(N):
        i = np.random.randint(N)
        dE = 2.0 * spins[i] * campo_local(i, spins, vec, Jv, nvec)
        if dE <= 0.0 or np.random.random() < np.exp(-beta * dE):
            spins[i] = -spins[i]


@njit(cache=True)
def cluster_wolff(spins, vec, Jv, nvec, beta, pila, en_cluster):
    """Construye y voltea UN cluster de Wolff. Devuelve su tamaño."""
    N = spins.size
    semilla = np.random.randint(N)
    en_cluster[semilla] = True
    pila[0] = semilla
    tope = 1
    tam = 0
    # 'pila' guarda también a todos los miembros: [0, tam) ya procesados
    while tam < tope:
        i = pila[tam]
        tam += 1
        si = spins[i]
        for c in range(nvec[i]):
            j = vec[i, c]
            if en_cluster[j]:
                continue
            JSS = Jv[i, c] * si * spins[j]
            if JSS > 0.0 and np.random.random() < 1.0 - np.exp(-2.0 * beta * JSS):
                en_cluster[j] = True
                pila[tope] = j
                tope += 1
    for m in range(tam):
        i = pila[m]
        spins[i] = -spins[i]
        en_cluster[i] = False
    return tam


@njit(cache=True)
def paso_mc(spins, vec, Jv, nvec, beta, alg, pila, en_cluster, n_cl):
    """
    Un paso MC. Con Wolff se voltean EXACTAMENTE n_cl clusters (número fijo).
    OJO: no se debe usar "clusters hasta voltear N espines" como regla de
    parada durante la medición: ese tiempo de parada depende del estado y
    sesga los promedios (lo comprobamos: ~6 % en E contra Onsager).
    Devuelve el total de espines volteados por Wolff.
    """
    if alg == 0 or alg == 2:
        barrido_metropolis(spins, vec, Jv, nvec, beta)
    volteados = 0
    if alg == 1 or alg == 2:
        for _ in range(n_cl):
            volteados += cluster_wolff(spins, vec, Jv, nvec, beta, pila, en_cluster)
    return volteados


@njit(cache=True)
def energia_total(spins, vec, Jv, nvec):
    """E = -(1/2) sum_i S_i h_i   (el 1/2 corrige que cada par aparece dos veces)."""
    s = 0.0
    for i in range(spins.size):
        s += spins[i] * campo_local(i, spins, vec, Jv, nvec)
    return -0.5 * s


@njit(cache=True)
def orden(spins, eps):
    """M_ord = (1/N) sum_i eps_i S_i."""
    s = 0.0
    for i in range(spins.size):
        s += eps[i] * spins[i]
    return s / spins.size


@njit(cache=True)
def simular_temperatura(spins, vec, Jv, nvec, eps, beta, n_term, n_skip, n_meas, alg):
    """
    Termaliza n_term pasos y luego toma n_meas medidas separadas por n_skip
    pasos. Modifica 'spins' en el lugar (para el annealing a la siguiente T).
    Devuelve (E[n_meas] en meV, M[n_meas]).
    """
    N = spins.size
    pila = np.empty(N, dtype=np.int64)
    en_cluster = np.zeros(N, dtype=np.bool_)
    # Termalización. Con Wolff, el número de clusters por paso se ajusta para
    # voltear ~N espines por paso usando el tamaño MEDIO ACUMULADO de los
    # clusters (se reinicia a mitad de la termalización para olvidar el
    # arranque). No usar solo el último paso: en la fase ordenada los tamaños
    # son bimodales (gotas pequeñas o casi todo el sistema) y un paso con una
    # gota diminuta dispararía miles de clusters gigantes en el siguiente.
    # Durante la medición n_cl queda FIJO.
    n_cl = 1
    suma, cuenta = 0.0, 0
    for t in range(n_term):
        if t == n_term // 2:
            suma, cuenta = 0.0, 0
        v = paso_mc(spins, vec, Jv, nvec, beta, alg, pila, en_cluster, n_cl)
        if alg != 0:
            suma += v
            cuenta += n_cl
            n_cl = min(N, max(1, int(round(N / (suma / cuenta)))))
    E = np.empty(n_meas)
    M = np.empty(n_meas)
    for t in range(n_meas):
        for _ in range(n_skip):
            paso_mc(spins, vec, Jv, nvec, beta, alg, pila, en_cluster, n_cl)
        E[t] = energia_total(spins, vec, Jv, nvec)
        M[t] = orden(spins, eps)
    return E, M


def barrido_temperaturas(vec, Jv, nvec, eps, temperaturas, k_B, n_term, n_skip,
                         n_meas, algoritmo="hibrido", seed=None, spins0=None):
    """
    Recorre las temperaturas EN EL ORDEN DADO reutilizando la configuración
    (annealing: conviene ir de mayor a menor T). Devuelve
        E[nT, n_meas], M[nT, n_meas], spins_final
    """
    alg = ALGORITMOS[algoritmo]
    N = len(nvec)
    if seed is not None:
        np.random.seed(seed)
        sembrar(seed)
    spins = (np.random.choice(np.array([-1, 1]), size=N).astype(np.int64)
             if spins0 is None else spins0.copy())
    eps = np.asarray(eps, dtype=np.int64)
    nT = len(temperaturas)
    E = np.empty((nT, n_meas))
    M = np.empty((nT, n_meas))
    for it, T in enumerate(temperaturas):
        beta = 1.0 / (k_B * T)
        E[it], M[it] = simular_temperatura(spins, vec, Jv, nvec, eps, beta,
                                           n_term, n_skip, n_meas, alg)
    return E, M, spins
