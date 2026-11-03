"""
================================================================
 analisis.py — De las series temporales a T_c con barras de error
================================================================

1. tau_int (Madras–Sokal): cuántas medidas consecutivas están correlacionadas.
2. Jackknife por bloques: errores de cantidades NO lineales (U, chi, C).
3. Cruces del cumulante de Binder entre pares de tamaños, con error
   (bootstrap paramétrico sobre las barras de error de U).
4. Extrapolación al límite termodinámico con ajuste ponderado.

Definiciones (por sitio, N espines, beta = 1/(k_B T) en 1/meV):
    U_L  = 1 - <M^4> / (3 <M^2>^2)
    chi  = N beta (<M^2> - <|M|>^2)            [1/meV]
    C/kB = beta^2 (<E^2> - <E>^2) / N            [adimensional]
"""

import numpy as np

U_ESTRELLA_ISING_2D = 0.6107   # red cuadrada, PBC, celda cuadrada (Kamieniarz–Blöte)


# ============================ AUTOCORRELACIÓN ================================

def funcion_autocorrelacion(x):
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    return acf / acf[0] if acf[0] > 0 else np.zeros(n)


def tau_int(x, c=6.0):
    """
    Tiempo de autocorrelación integrado (en unidades de MEDIDAS) con la
    ventana automática de Madras–Sokal: se suma C(t) hasta la primera
    ventana W que cumple W >= c * tau(W).
    """
    x = np.asarray(x, float)
    if len(x) < 4 or np.var(x) == 0:
        return 0.5
    rho = funcion_autocorrelacion(x)
    tau_W = 0.5 + np.cumsum(rho[1:])
    W = np.arange(1, len(x))
    ok = W >= c * tau_W
    return float(tau_W[np.argmax(ok)]) if ok.any() else float(tau_W[-1])


# ============================ OBSERVABLES + JACKKNIFE =======================

def _derivadas(m, N, beta):
    """m = [<|M|>, <M^2>, <M^4>, <E>, <E^2>] -> cantidades físicas."""
    mabs, m2, m4, e, e2 = m
    U = 1.0 - m4 / (3.0 * m2 ** 2) if m2 > 0 else np.nan
    chi = N * beta * (m2 - mabs ** 2)
    C = beta ** 2 * (e2 - e ** 2) / N
    return np.array([U, chi, C, mabs, e / N])


NOMBRES = ["U", "chi", "C", "Mabs", "E"]


def observables(E, M, N, T, k_B, n_bloques=20):
    """
    Calcula U, chi, C, <|M|>, E/N con errores jackknife por bloques, y los
    tiempos de autocorrelación de E y de M^2.
    Devuelve un diccionario {nombre: (valor, error), ..., 'tau_E', 'tau_M2',
    'n_eff', 'bloques_ok'}.
    """
    beta = 1.0 / (k_B * T)
    E = np.asarray(E, float)
    M = np.asarray(M, float)
    n = len(E)
    K = int(n_bloques)
    b = n // K
    datos = np.stack([np.abs(M), M ** 2, M ** 4, E, E ** 2])[:, :K * b]
    bloques = datos.reshape(5, K, b).mean(axis=2)            # (5, K)
    total = bloques.mean(axis=1)
    loo = (K * total[:, None] - bloques) / (K - 1)            # deja uno fuera
    theta = np.array([_derivadas(loo[:, k], N, beta) for k in range(K)])
    valor = _derivadas(total, N, beta)
    err = np.sqrt((K - 1) / K * np.sum((theta - theta.mean(axis=0)) ** 2, axis=0))

    tE, tM = tau_int(E), tau_int(M ** 2)
    tmax = max(tE, tM)
    out = {nom: (valor[i], err[i]) for i, nom in enumerate(NOMBRES)}
    out.update(tau_E=tE, tau_M2=tM, n_eff=n / (2 * tmax), largo_bloque=b,
               bloques_ok=bool(b >= 10 * tmax))
    return out


def tabla_observables(T, E, M, N, k_B, n_bloques=20):
    """Aplica `observables` a cada temperatura; devuelve arreglos (valor, error)."""
    filas = [observables(E[i], M[i], N, T[i], k_B, n_bloques) for i in range(len(T))]
    tab = {nom: (np.array([f[nom][0] for f in filas]), np.array([f[nom][1] for f in filas]))
           for nom in NOMBRES}
    for clave in ("tau_E", "tau_M2", "n_eff", "bloques_ok"):
        tab[clave] = np.array([f[clave] for f in filas])
    return tab


# ============================ CRUCES DEL CUMULANTE ==========================

def _ajuste_local(T, U, s, i, ancho=3):
    """Polinomio de grado 2 ponderado alrededor del intervalo [T_i, T_i+1]."""
    lo, hi = max(0, i - ancho + 1), min(len(T), i + ancho + 1)
    w = 1.0 / np.maximum(s[lo:hi], 1e-12)
    grado = 2 if hi - lo >= 4 else 1
    T0 = 0.5 * (T[i] + T[i + 1])
    return np.polyfit(T[lo:hi] - T0, U[lo:hi], grado, w=w), T0


def _raiz_en(p1, p2, T0, a, b):
    d = np.polysub(p1, p2)
    r = np.roots(d) + T0 if len(d) > 1 else np.array([])
    r = r[np.isreal(r)].real
    r = r[(r >= a) & (r <= b)]
    return r[np.argmin(np.abs(r - 0.5 * (a + b)))] if len(r) else np.nan


def cruces_binder(T, U1, s1, U2, s2, n_boot=400, z_min=2.0, rng=None):
    """
    Cruces entre U_L1(T) y U_L2(T) con L2 > L1.

    Solo se aceptan cruces con la orientación física: por debajo de T_c el
    sistema grande está MÁS ordenado (U_L2 > U_L1) y por encima, menos.
    Un cruce es 'robusto' si a ambos lados la diferencia es significativa
    (|dU|/sigma > z_min); así se descartan cruces de ruido en las mesetas
    U ~ 2/3 (baja T) y U ~ 0 (alta T).

    Devuelve una lista de diccionarios {T, sT, U, sU, robusto}, ordenada con
    los robustos primero y, entre ellos, por significancia.
    """
    rng = np.random.default_rng(rng)
    o = np.argsort(T)
    T, U1, s1, U2, s2 = (np.asarray(x, float)[o] for x in (T, U1, s1, U2, s2))
    dU = U2 - U1
    sd = np.sqrt(s1 ** 2 + s2 ** 2) + 1e-12
    z = dU / sd
    salida = []
    for i in range(len(T) - 1):
        if not (dU[i] > 0 and dU[i + 1] <= 0):
            continue
        z_abajo = np.max(z[max(0, i - 2):i + 1])
        z_arriba = -np.min(z[i + 1:i + 4])
        robusto = bool(z_abajo > z_min and z_arriba > z_min)
        a, b = T[max(0, i - 1)], T[min(len(T) - 1, i + 2)]
        p1, T0 = _ajuste_local(T, U1, s1, i)
        p2, _ = _ajuste_local(T, U2, s2, i)
        Tx = _raiz_en(p1, p2, T0, a, b)
        if not np.isfinite(Tx):
            Tx = T[i] + (T[i + 1] - T[i]) * dU[i] / (dU[i] - dU[i + 1])
        Ux = np.polyval(p1, Tx - T0)
        bT, bU = [], []
        for _ in range(n_boot):
            V1 = U1 + s1 * rng.standard_normal(len(T))
            V2 = U2 + s2 * rng.standard_normal(len(T))
            q1, _ = _ajuste_local(T, V1, s1, i)
            q2, _ = _ajuste_local(T, V2, s2, i)
            r = _raiz_en(q1, q2, T0, a, b)
            if np.isfinite(r):
                bT.append(r)
                bU.append(np.polyval(q1, r - T0))
        sT = float(np.std(bT)) if len(bT) > 10 else np.nan
        sU = float(np.std(bU)) if len(bU) > 10 else np.nan
        salida.append(dict(T=float(Tx), sT=sT, U=float(Ux), sU=sU, robusto=robusto,
                           signif=float(min(z_abajo, z_arriba)),
                           frac_boot=len(bT) / n_boot))
    salida.sort(key=lambda c: (not c["robusto"], -c["signif"]))
    return salida


# ============================ EXTRAPOLACIÓN =================================

def ajuste_lineal_ponderado(x, y, s):
    """y = c0 + c1 x con pesos 1/s^2. Devuelve c, cov, chi2_red."""
    x, y, s = (np.asarray(v, float) for v in (x, y, s))
    A = np.vstack([np.ones_like(x), x]).T
    W = np.diag(1.0 / s ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    c = cov @ A.T @ W @ y
    dof = len(x) - 2
    chi2 = float(np.sum(((y - A @ c) / s) ** 2))
    chi2_red = chi2 / dof if dof > 0 else np.nan
    if dof > 0 and chi2_red > 1:        # errores subestimados: se inflan
        cov = cov * chi2_red
    return c, cov, chi2_red


def extrapolar_Tc(L_eff, T_cruce, s_cruce, exponente=1.0):
    """
    T_x(L) = T_c(inf) + a * L_eff^(-exponente).

    exponente = 1/nu = 1     : forma usada en el protocolo.
    exponente = 1/nu + omega : corrimiento teórico de los cruces del cumulante
                               (el término L^-1/nu se cancela en U_L = U_L'),
                               ~3 para Ising 2D.
    Con un solo par no se extrapola: se devuelve el cruce del par.
    """
    L_eff, T_cruce, s_cruce = (np.asarray(v, float) for v in (L_eff, T_cruce, s_cruce))
    if len(T_cruce) == 1:
        return dict(Tc=float(T_cruce[0]), sTc=float(s_cruce[0]), a=np.nan,
                    chi2_red=np.nan, exponente=exponente, n=1)
    x = L_eff ** (-exponente)
    c, cov, chi2r = ajuste_lineal_ponderado(x, T_cruce, s_cruce)
    return dict(Tc=float(c[0]), sTc=float(np.sqrt(cov[0, 0])), a=float(c[1]),
                chi2_red=chi2r, exponente=exponente, n=len(x))
