import numpy as np

from isingmc import K_B
from isingmc.analisis import observables
from isingmc.montecarlo import barrido_temperaturas, energia_total
from isingmc.red import construir_libreta
from redes_prueba import cuadrada

TC = 2.0 / np.log(1 + np.sqrt(2)) / K_B        # Onsager, |J| = 1 meV


def test_energia_FM_ordenado():
    vec, Jv, nvec = construir_libreta(6, 6, 1, cuadrada())
    assert np.isclose(energia_total(np.ones(36, dtype=np.int64), vec, Jv, nvec), -72.0)


def test_algoritmos_coinciden_y_reproducen_onsager_lejos_de_Tc():
    # A T = 0.8 Tc y L = 16 la energía casi no depende del tamaño: Onsager exacto
    from scipy.special import ellipk
    T = 0.8 * TC
    b = 1 / (K_B * T)
    k = 2 * np.sinh(2 * b) / np.cosh(2 * b) ** 2
    u_exacta = -1 / np.tanh(2 * b) * (1 + 2 / np.pi * (2 * np.tanh(2 * b) ** 2 - 1) * ellipk(k ** 2))
    vec, Jv, nvec = construir_libreta(16, 16, 1, cuadrada())
    eps = np.ones(256, dtype=np.int64)
    for alg in ("metropolis", "wolff", "hibrido"):
        E, M, _ = barrido_temperaturas(vec, Jv, nvec, eps, [T], K_B, 1000, 2, 4000, alg, seed=7)
        o = observables(E[0], M[0], 256, T, K_B)
        u, su = o["E"]
        assert abs(u - u_exacta) < 4 * su + 2e-3, (alg, u, su, u_exacta)
