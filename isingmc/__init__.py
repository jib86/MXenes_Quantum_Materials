"""
isingmc — Monte Carlo de Ising para materiales magnéticos 2D
=============================================================

Predice la temperatura de transición magnética (Curie o Néel) de una red 2D
a partir de sus acoplamientos de canje J, con el cumulante de Binder y
escalamiento de tamaño finito.

Módulos
-------
red                 enlaces (k, k', n1, n2, J), supercelda y libreta de vecinos
estado_fundamental  orden magnético a T = 0 por enumeración exacta
montecarlo          Metrópolis, Wolff e híbrido (numba); series E(t), M(t)
analisis            tau_int, jackknife, cruces de Binder, extrapolación
flujo               flujo completo: de la configuración a T_c ± error
graficas            figuras de los resultados

Convención: H = - sum_<ij> J_ij S_i S_j, S = ±1, cada par una vez, J en meV.
"""

__version__ = "1.0.0"

K_B = 0.0861733  # constante de Boltzmann, meV/K
