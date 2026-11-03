"""
================================================================
 parametros.py — ÚNICO archivo que hay que editar
================================================================
 Material (celda, posiciones, J) + parámetros de la simulación.
 Convención:  H = - sum_<ij> J_ij S_i S_j ,  S = ±1, cada par UNA vez,
              J > 0 ferromagnético, J en meV  (ver red.py).
"""
import numpy as np

# constante de Boltzmann
k_B = 0.0861733   # meV/K

# ===================== DATOS DEL MATERIAL =====================================
# MXene Nb–Ti 30 % DESNUDA (10 Ti en una supercelda 4x2), J de TB2J tomados del
# informe de Luis (13 oct 2026). OJO: el informe no trae la lista completa de
# enlaces ni las posiciones, así que modelo_desnuda.py RECONSTRUYE la red y los
# enlaces a partir de lo que el informe sí da (ver los supuestos allí). Control:
# reproduce J0 por capa y lambda_max = 49.2 meV del informe.
# Cuando se tenga el exchange.out de TB2J, reemplazar por los enlaces reales.
import modelo_desnuda
cell, posiciones, _enlaces_desnuda, info_modelo = modelo_desnuda.construir()

# ---------------- Acoplamientos -----------------------------------------------
# modo_J = "anillos": se da una distancia (o una LISTA de distancias) por J y
#                     el programa encuentra todos los pares a esa distancia.
# modo_J = "enlaces": lista explícita (k, k', n1, n2, J_meV): átomo k de la
#                     celda (0,0) con el átomo k' de la celda (n1, n2). Es el
#                     formato natural de TB2J. Se puede dar cada par en los
#                     dos sentidos o en uno solo (no mezclar).
modo_J = "enlaces"
enlaces = _enlaces_desnuda

anillos = [             # solo se usa con modo_J = "anillos"
    # {"d": 3.13, "J": 16.6},
    # {"d": [5.40, 5.45], "J": 0.96},
]
tol = 0.03        # Å, ventana de cada distancia (verificar con el histograma)

# Factor que multiplica todos los J para llevarlos a la convención de arriba.
#   1.0 : sus J ya usan H = -sum_<ij> J S_i S_j con S = ±1
#   2.0 : su fuente usa la suma doble H = -sum_{i≠j} J S_i S_j (cada par 2 veces)
#         -> ES EL CASO DE TB2J (espines unitarios, i ≠ j)
#  -1.0 : su fuente usa H = +sum J S S (J > 0 antiferro)
# (si su mapeo DFT usó |S| ≠ 1, multiplique además por |S|^2)
factor_J = 2.0

# Patrón del estado fundamental. None = se busca automáticamente (recomendado).
# Si se da a mano, el programa verifica que de verdad sea el de menor energía.
epsilon_base = None
periodo_epsilon = None   # (m1, m2) si epsilon_base cubre más de una celda

# Energía de anisotropía magnética por átomo magnético (meV), del DFT.
# Ising solo es cuantitativo si MAE >~ 10 k_B T_c. None = no se verifica.
# Desnuda: escenario (ii) 10 meV/celda = 1.0 meV/Ti; escenario (i) ~ 8-10 meV/Ti.
MAE_meV = 1.0

# ===================== PARÁMETROS DE LA SIMULACIÓN ===========================

# Tamaños de supercelda L x (aspecto*L). Para escalamiento de tamaño finito
# hacen falta >= 3 tamaños; si el estado fundamental tiene periodo m, los L
# deben ser múltiplos de m.
# Desnuda: celda 12.52 x 6.26 Å, con aspecto = 2 la supercelda es L x 2L celdas
# (casi isótropa): L = 6, 9, 12, 18 -> 7.5 a 22.5 nm, N = 720 a 6480 espines.
L_list  = [6, 9, 12, 18]
aspecto = 2

# Rango de temperaturas (K). None = automático a partir del campo medio.
T_min = None
T_max = None

# Barrido en dos pasadas: una gruesa para ubicar T_c y una fina alrededor.
N_temps_grueso = 14
N_temps_fino   = 14
ventana_fina   = 0.04     # la pasada fina cubre T_c_estimada * (1 ± ventana)
fraccion_grueso = 0.2     # estadística de la pasada gruesa relativa a la fina

# Algoritmo: "hibrido" (recomendado), "wolff" o "metropolis"
algoritmo = "hibrido"

# Pasos MC por temperatura: N_term + N_skip * N_meas
N_term = 2000      # termalización (pasos descartados)
N_skip = 1         # pasos entre medidas (verifique N_skip >= 2 tau_int en la salida)
N_meas = 20000     # medidas por temperatura (pasada fina)
n_bloques = 20     # bloques para el jackknife

# Extrapolación T_x(L) = T_c + a L_eff^-x.  x = 1 (protocolo) ó 3 (1/nu + omega)
exponente_fss = 1.0

seed = 42
n_procesos = 1     # > 1 corre los distintos L en paralelo
forzar = False     # True: simular aunque el estado fundamental sea degenerado
output_dir = "./resultados"
