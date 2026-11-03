import numpy as np

from isingmc.analisis import ajuste_lineal_ponderado, cruces_binder, extrapolar_Tc, tau_int


def test_tau_ruido_blanco():
    x = np.random.default_rng(0).standard_normal(20000)
    assert abs(tau_int(x) - 0.5) < 0.1


def test_tau_AR1():
    rng = np.random.default_rng(1)
    rho, n = 0.8, 200000
    x = np.empty(n)
    x[0] = 0
    for t in range(1, n):
        x[t] = rho * x[t - 1] + rng.standard_normal()
    esperado = 0.5 * (1 + rho) / (1 - rho)      # 4.5
    assert abs(tau_int(x) - esperado) / esperado < 0.1


def test_extrapolacion_recta_exacta():
    L = np.array([10.0, 20.0, 40.0])
    r = extrapolar_Tc(L, 100 + 50 / L, np.full(3, 0.1), exponente=1.0)
    assert abs(r["Tc"] - 100) < 1e-8


def test_cruce_sintetico():
    T = np.linspace(0.8, 1.2, 21)
    U = lambda L: 0.6 - np.tanh(L * (T - 1.0)) * 0.3
    s = np.full_like(T, 1e-3)
    cr = cruces_binder(T, U(10), s, U(20), s, n_boot=50, rng=0)
    assert cr and cr[0]["robusto"] and abs(cr[0]["T"] - 1.0) < 0.01
