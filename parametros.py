"""
================================================================
 parametros.py — ÚNICO archivo que hay que editar
================================================================
 Material (celda, posiciones, J) + parámetros de la simulación.
 Convención:  H = - sum_<ij> J_ij S_i S_j ,  S = ±1, cada par UNA vez,
              J > 0 ferromagnético, J en meV  (ver isingmc/red.py y docs/metodo.md).
"""
import numpy as np

# constante de Boltzmann
k_B = 0.0861733   # meV/K

# ===================== DATOS DEL MATERIAL =====================================
# MXene Nb–Ti 30 % DESNUDA: supercelda 4x2 con 10 Ti magnéticos (6 en la capa
# metálica inferior, 4 en la superior). J de TB2J tomados del informe del
# 13 oct 2026. El informe no trae la lista completa de enlaces ni las
# posiciones: esta lista es la RECONSTRUCCIÓN de materiales/nbti30_desnuda.py
# (ver allí los supuestos). Control: reproduce lambda_max = 49.2 meV del informe.

# ---------------- Celda y átomos magnéticos -----------------------------------
# Vectores de red (Å). Filas 1 y 2: plano del material. Fila 3: vacío.
cell = np.array([
    [  12.5200,    0.0000,    0.000],   # a1
    [  -3.1300,    5.4213,    0.000],   # a2
    [   0.0000,    0.0000,   24.120],   # c (vacío, no periódica)
])

# Posiciones (Å) de los átomos magnéticos. El número k de cada átomo es su
# lugar en esta lista (empieza en 0) y es el que se usa en la tabla de J.
posiciones = np.array([
    [   0.5217,    0.9036,   15.800],   # k = 0:  Ti1 (capa superior)
    [  -1.0433,    3.6142,   15.800],   # k = 1:  Ti2 (capa superior)
    [  -1.5650,    2.7107,    8.200],   # k = 2:  Ti3 (capa inferior)
    [   4.6950,    2.7107,    8.200],   # k = 3:  Ti4 (capa inferior)
    [   3.1300,    0.0000,    8.200],   # k = 4:  Ti5 (capa inferior)
    [   1.5650,    2.7107,    8.200],   # k = 5:  Ti6 (capa inferior)
    [   6.7817,    0.9036,   15.800],   # k = 6:  Ti7 (capa superior)
    [   5.2167,    3.6142,   15.800],   # k = 7:  Ti8 (capa superior)
    [   9.3900,    0.0000,    8.200],   # k = 8:  Ti9 (capa inferior)
    [   7.8250,    2.7107,    8.200],   # k = 9: Ti10 (capa inferior)
])

# ==============================================================================
#  ACOPLAMIENTOS DE CANJE J (meV)  —  AQUÍ SE ESCRIBEN LOS J DEL DFT
# ==============================================================================
# Cada línea es un enlace:   (k, k', n1, n2, J)
#   k, k'  : números de los dos átomos (su lugar en 'posiciones', desde 0)
#   n1, n2 : celda en la que está el átomo k' respecto del átomo k:
#            (0, 0) = misma celda; (1, 0) = celda vecina a lo largo de a1;
#            (0, -1) = celda anterior a lo largo de a2; etc.
#            Es la columna R del archivo exchange.out de TB2J.
#   J      : en meV, tal como lo da la fuente (TB2J). J > 0 ferromagnético.
#            La conversión de convención la hace factor_J, más abajo.
# Cada par se escribe UNA sola vez: el programa agrega el sentido inverso.
# Un enlace que no aparece en la lista tiene J = 0.
modo_J = "enlaces"
enlaces = [
    # --- Capa inferior, primeros vecinos (3.13 Å) ---
    (2, 4,  0,  1,   16.60),   #  Ti3 – Ti5  celda vecina (+0, +1)
    (2, 5,  0,  0,   16.60),   #  Ti3 – Ti6  misma celda
    (3, 8,  0,  1,   16.60),   #  Ti4 – Ti9  celda vecina (+0, +1)
    (3, 9,  0,  0,   16.60),   #  Ti4 – Ti10 misma celda
    (4, 5,  0,  0,   16.60),   #  Ti5 – Ti6  misma celda
    (8, 9,  0,  0,   16.60),   #  Ti9 – Ti10 misma celda
    (2, 8, -1,  0,    8.90),   #  Ti3 – Ti9  celda vecina (-1, +0)
    (2, 9, -1,  0,    8.90),   #  Ti3 – Ti10 celda vecina (-1, +0)
    (3, 4,  0,  0,    8.90),   #  Ti4 – Ti5  misma celda
    (3, 5,  0,  0,    8.90),   #  Ti4 – Ti6  misma celda
    (4, 5,  0, -1,    8.90),   #  Ti5 – Ti6  celda vecina (+0, -1)
    (8, 9,  0, -1,    8.90),   #  Ti9 – Ti10 celda vecina (+0, -1)
    # --- Capa inferior, terceros vecinos (6.26 Å) ---
    (2, 3, -1, -1,   -0.65),   #  Ti3 – Ti4  celda vecina (-1, -1)
    (2, 3, -1,  0,   -0.65),   #  Ti3 – Ti4  celda vecina (-1, +0)
    (2, 3,  0,  0,   -0.65),   #  Ti3 – Ti4  misma celda
    (2, 3,  0,  1,   -0.65),   #  Ti3 – Ti4  celda vecina (+0, +1)
    (4, 8, -1, -1,   -0.65),   #  Ti5 – Ti9  celda vecina (-1, -1)
    (4, 8, -1,  0,   -0.65),   #  Ti5 – Ti9  celda vecina (-1, +0)
    (4, 8,  0,  0,   -0.65),   #  Ti5 – Ti9  misma celda
    (4, 8,  0,  1,   -0.65),   #  Ti5 – Ti9  celda vecina (+0, +1)
    (5, 9, -1, -1,   -0.65),   #  Ti6 – Ti10 celda vecina (-1, -1)
    (5, 9, -1,  0,   -0.65),   #  Ti6 – Ti10 celda vecina (-1, +0)
    (5, 9,  0,  0,   -0.65),   #  Ti6 – Ti10 misma celda
    (5, 9,  0,  1,   -0.65),   #  Ti6 – Ti10 celda vecina (+0, +1)
    # --- Capa superior, primeros vecinos (3.13 Å) ---
    (0, 1,  0,  0,   14.00),   #  Ti1 – Ti2  misma celda
    (6, 7,  0,  0,   14.00),   #  Ti7 – Ti8  misma celda
    (0, 1,  0, -1,   11.80),   #  Ti1 – Ti2  celda vecina (+0, -1)
    (6, 7,  0, -1,   11.80),   #  Ti7 – Ti8  celda vecina (+0, -1)
    # --- Capa superior, segundos vecinos (5.42 Å) ---
    (0, 7, -1, -1,    0.96),   #  Ti1 – Ti8  celda vecina (-1, -1)
    (0, 7,  0,  0,    0.96),   #  Ti1 – Ti8  misma celda
    (1, 6, -1,  0,    0.96),   #  Ti2 – Ti7  celda vecina (-1, +0)
    (1, 6,  0,  1,    0.96),   #  Ti2 – Ti7  celda vecina (+0, +1)
    # --- Capa superior, terceros vecinos (6.26 Å) ---
    (0, 6, -1, -1,   -0.50),   #  Ti1 – Ti7  celda vecina (-1, -1)
    (0, 6, -1,  0,   -0.50),   #  Ti1 – Ti7  celda vecina (-1, +0)
    (0, 6,  0,  0,   -0.50),   #  Ti1 – Ti7  misma celda
    (0, 6,  0,  1,   -0.50),   #  Ti1 – Ti7  celda vecina (+0, +1)
    (1, 7, -1, -1,   -0.50),   #  Ti2 – Ti8  celda vecina (-1, -1)
    (1, 7, -1,  0,   -0.50),   #  Ti2 – Ti8  celda vecina (-1, +0)
    (1, 7,  0,  0,   -0.50),   #  Ti2 – Ti8  misma celda
    (1, 7,  0,  1,   -0.50),   #  Ti2 – Ti8  celda vecina (+0, +1)
    # --- Entre capas (7.7–13.4 Å) ---
    # (el informe no da estos J; valor uniforme que reproduce J0 = 31 meV en la capa superior)
    (0, 2,  0,  0,    0.88),   #  Ti1 – Ti3  misma celda
    (0, 3,  0,  0,    0.88),   #  Ti1 – Ti4  misma celda
    (0, 4,  0,  0,    0.88),   #  Ti1 – Ti5  misma celda
    (0, 5,  0,  0,    0.88),   #  Ti1 – Ti6  misma celda
    (0, 8,  0,  0,    0.88),   #  Ti1 – Ti9  misma celda
    (0, 9,  0,  0,    0.88),   #  Ti1 – Ti10 misma celda
    (1, 2,  0,  0,    0.88),   #  Ti2 – Ti3  misma celda
    (1, 3,  0,  0,    0.88),   #  Ti2 – Ti4  misma celda
    (1, 4,  0,  0,    0.88),   #  Ti2 – Ti5  misma celda
    (1, 5,  0,  0,    0.88),   #  Ti2 – Ti6  misma celda
    (1, 8,  0,  0,    0.88),   #  Ti2 – Ti9  misma celda
    (1, 9,  0,  0,    0.88),   #  Ti2 – Ti10 misma celda
    (2, 6,  0,  0,    0.88),   #  Ti3 – Ti7  misma celda
    (2, 7,  0,  0,    0.88),   #  Ti3 – Ti8  misma celda
    (3, 6,  0,  0,    0.88),   #  Ti4 – Ti7  misma celda
    (3, 7,  0,  0,    0.88),   #  Ti4 – Ti8  misma celda
    (4, 6,  0,  0,    0.88),   #  Ti5 – Ti7  misma celda
    (4, 7,  0,  0,    0.88),   #  Ti5 – Ti8  misma celda
    (5, 6,  0,  0,    0.88),   #  Ti6 – Ti7  misma celda
    (5, 7,  0,  0,    0.88),   #  Ti6 – Ti8  misma celda
    (6, 8,  0,  0,    0.88),   #  Ti7 – Ti9  misma celda
    (6, 9,  0,  0,    0.88),   #  Ti7 – Ti10 misma celda
    (7, 8,  0,  0,    0.88),   #  Ti8 – Ti9  misma celda
    (7, 9,  0,  0,    0.88),   #  Ti8 – Ti10 misma celda
]

# Alternativa más simple: modo_J = "anillos", un J por distancia; el programa
# busca todos los pares a esa distancia. Solo se usa si modo_J = "anillos".
anillos = [
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

# Rango de temperaturas (K). None = automático: de 0.25 a 1.05 veces lambda_max/k_B
# (cota superior de T_c; no se usa como estimación).
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
output_dir = "resultados/material"
