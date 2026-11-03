"""
================================================================
 modelo_desnuda.py — MXene Nb–Ti 30 % DESNUDA: modelo RECONSTRUIDO
================================================================

Fuente: informe "MXene Nb–Ti 30 %: magnetismo de la monocapa mixta y la
desnuda" (Luis, 13 oct 2026; TB2J 0.9.12 + Vampire). El informe da los J
resumidos por capa y distancia, pero NO la lista completa de enlaces ni las
posiciones de los Ti. Este módulo reconstruye un modelo consistente con todo
lo que el informe sí dice. CUANDO SE TENGA EL exchange.out DE TB2J, ESTE
MODELO SE REEMPLAZA POR LOS ENLACES REALES.

Lo que se toma TAL CUAL del informe
-----------------------------------
* Supercelda 4x2 de la red triangular de metales: a1 = 12.52 Å, a2 = 6.26 Å,
  120°; 10 Ti: 6 en la capa metálica inferior (z ≈ 8.2 Å), 4 en la superior
  (z ≈ 15.8 Å). Estado fundamental FM (confirmado por DFT, J(q) y 2^10
  configuraciones).
* Capa inferior, primeros vecinos: 12 enlaces FM entre +8.9 y +16.9 meV. Los
  más fuertes (+16.3 a +16.9 meV): Ti3–Ti6, Ti4–Ti10, Ti3–Ti5, Ti4–Ti9,
  Ti9–Ti10, Ti5–Ti6.
* Capa inferior, 6.26 Å (2a): 14 enlaces AFM entre −0.3 y −1.0 meV (la
  red reconstruida tiene 12 pares distintos a esa distancia: se les da el
  valor medio);
  5.4 Å: menores de 0.3 meV (se toman 0).
* Capa superior: Ti1–Ti2 y Ti7–Ti8, +11.8 a +14.0 meV (dos imágenes cada uno,
  3.12–3.14 Å); 5.4 Å: 4 FM de +0.96 meV; 6.26 Å: 8 AFM de −0.47 a −0.54 meV.
* J0 = sum_j J_ij = +31 meV en la capa superior y +49 a +53 en la inferior.
* Convención TB2J: E = −sum_{i≠j} J S_i·S_j (cada par DOS veces, espines
  unitarios)  ->  en nuestro código factor_J = 2.
* MAE: escenario (ii) 10 meV/celda = 1.0 meV/Ti; escenario (i) 78–106
  meV/celda = 7.8–10.6 meV/Ti.

Lo que se SUPONE (reconstrucción)
---------------------------------
1. Sitios de Ti en la red triangular ideal (a = 3.13 Å), sin relajación.
   Las vacantes de Ti en la capa inferior (ocupadas por Nb) están en sitios
   no vecinos: es la ÚNICA forma de obtener exactamente 12 enlaces de primeros
   vecinos con 6 Ti en 8 sitios, y reproduce la lista de pares del informe.
2. En la capa superior los 4 Ti forman dos columnas no vecinas a lo largo de
   a2 (única forma de tener solo 4 enlaces de primeros vecinos, 4 a 5.4 Å y
   8 a 6.26 Å, como dice el informe).
3. Los 6 primeros vecinos inferiores que no están en la lista de "más
   fuertes" reciben un mismo valor J_debil, y el canje entre capas (que el
   informe no da, salvo "sin AFM mayores de 0.5 meV") se representa como un
   acoplamiento uniforme J_perp entre cada Ti superior y cada Ti inferior.
   J_perp se fija para reproducir J0 = 31 meV en la capa superior. J_debil
   se fija para J0 = 51 meV en la inferior, pero sin bajar del extremo
   inferior del rango del informe (8.9 meV); con 8.9 el J0 inferior da
   51.9 meV, dentro de los 49–53 del informe.
4. Valores dentro de cada rango: fuertes 16.6, AFM inferior −0.65, AFM
   superior −0.50 meV; Ti1–Ti2 y Ti7–Ti8: 14.0 (R = 0) y 11.8 (imagen en a2).
Control independiente: el máximo autovalor de J(q=0) debe dar ~49.2 meV,
como en el informe (no se usó para ajustar).
"""

import numpy as np

A = 3.13                                     # Å, parámetro de la red triangular
a = np.array([A, 0.0, 0.0])
b = np.array([-A / 2, A * np.sqrt(3) / 2, 0.0])
cell = np.array([4 * a, 2 * b, [0.0, 0.0, 24.12]])
Z_INF, Z_SUP = 8.2, 15.8

# Ti -> (capa, i, j) en coordenadas de la red primitiva (i mod 4, j mod 2)
SITIOS = {
    "Ti1": ("sup", 0, 0), "Ti2": ("sup", 0, 1),
    "Ti7": ("sup", 2, 0), "Ti8": ("sup", 2, 1),
    "Ti3": ("inf", 0, 1), "Ti4": ("inf", 2, 1),
    "Ti5": ("inf", 1, 0), "Ti6": ("inf", 1, 1),
    "Ti9": ("inf", 3, 0), "Ti10": ("inf", 3, 1),
}
NOMBRES = [f"Ti{n}" for n in range(1, 11)]          # k = número de Ti - 1
DESPL_SUP = (a + b) / 3                              # registro de la capa superior (supuesto)

J_FUERTE = 16.6
J_AFM_INF = -0.65
J_AFM_SUP = -0.50
J_SUP_R0, J_SUP_A2 = 14.0, 11.8
J_2DO_SUP = 0.96
J0_INF, J0_SUP = 51.0, 31.0
J_DEBIL_MIN = 8.9
FUERTES = {frozenset(p) for p in [("Ti3", "Ti6"), ("Ti4", "Ti10"), ("Ti3", "Ti5"),
                                  ("Ti4", "Ti9"), ("Ti9", "Ti10"), ("Ti5", "Ti6")]}

V1 = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]       # primeros vecinos
V3 = [(2 * x, 2 * y) for (x, y) in V1]                          # 2a
V2 = [(2, 1), (1, 2), (1, -1), (-2, -1), (-1, -2), (-1, 1)]     # sqrt(3) a


def posiciones():
    pos = []
    for nom in NOMBRES:
        capa, i, j = SITIOS[nom]
        r = i * a + j * b
        if capa == "sup":
            r = r + DESPL_SUP
        pos.append(r + np.array([0, 0, Z_SUP if capa == "sup" else Z_INF]))
    return np.array(pos)


def _ocupante(capa, i, j):
    ii, jj = i % 4, j % 2
    for nom, (c, x, y) in SITIOS.items():
        if c == capa and x == ii and y == jj:
            return nom, (i - ii) // 4, (j - jj) // 2
    return None


def _vecinos(nom, vectores):
    capa, i, j = SITIOS[nom]
    for (dx, dy) in vectores:
        oc = _ocupante(capa, i + dx, j + dy)
        if oc and not (oc[0] == nom and oc[1] == 0 and oc[2] == 0):
            yield oc


def enlaces_dentro_de_capa(J_debil):
    """Lista (k, k', n1, n2, J) en los DOS sentidos, J en convención TB2J."""
    lista = []
    for nom in NOMBRES:
        capa = SITIOS[nom][0]
        k = NOMBRES.index(nom)
        for (otro, n1, n2) in _vecinos(nom, V1):
            if capa == "inf":
                fuerte = frozenset((nom, otro)) in FUERTES
                # en un par vertical (misma columna), la imagen R=0 es la fuerte
                misma_col = SITIOS[nom][1] == SITIOS[otro][1]
                if misma_col and (n1, n2) != (0, 0):
                    fuerte = False
                J = J_FUERTE if fuerte else J_debil
            else:
                J = J_SUP_R0 if (n1, n2) == (0, 0) else J_SUP_A2
            lista.append((k, NOMBRES.index(otro), n1, n2, J))
        for (otro, n1, n2) in _vecinos(nom, V2):
            if capa == "sup":
                lista.append((k, NOMBRES.index(otro), n1, n2, J_2DO_SUP))
        for (otro, n1, n2) in _vecinos(nom, V3):
            if otro == nom:
                continue                     # imagen de sí mismo a 6.26 Å: no listada, J = 0
            lista.append((k, NOMBRES.index(otro), n1, n2,
                          J_AFM_INF if capa == "inf" else J_AFM_SUP))
    return lista


def construir():
    """
    Devuelve cell, posiciones, enlaces (convención TB2J, ambos sentidos) y un
    diccionario con los parámetros ajustados y los controles.
    """
    sup = [n for n in NOMBRES if SITIOS[n][0] == "sup"]
    inf = [n for n in NOMBRES if SITIOS[n][0] == "inf"]
    # J0 superior sin el canje entre capas -> J_perp
    base = enlaces_dentro_de_capa(0.0)
    J0 = np.zeros(10)
    for (k, kp, n1, n2, J) in base:
        J0[k] += J
    J_perp = (J0_SUP - np.mean([J0[NOMBRES.index(n)] for n in sup])) / len(inf)
    # J0 inferior: cada Ti inferior tiene 2 enlaces débiles
    J0_inf_sin = np.mean([J0[NOMBRES.index(n)] for n in inf]) + len(sup) * J_perp
    J_debil = max((J0_INF - J0_inf_sin) / 2.0, J_DEBIL_MIN)
    lista = enlaces_dentro_de_capa(J_debil)
    for s in sup:
        for i_ in inf:
            ks, ki = NOMBRES.index(s), NOMBRES.index(i_)
            lista.append((ks, ki, 0, 0, J_perp))
            lista.append((ki, ks, 0, 0, J_perp))

    # Controles contra el informe
    J0 = np.zeros(10)
    A_q0 = np.zeros((10, 10))
    cuenta = {}
    pos = posiciones()
    for (k, kp, n1, n2, J) in lista:
        J0[k] += J
        A_q0[k, kp] += J
        d = np.linalg.norm(pos[kp] + n1 * cell[0] + n2 * cell[1] - pos[k])
        capas = tuple(sorted((SITIOS[NOMBRES[k]][0], SITIOS[NOMBRES[kp]][0])))
        clave = capas + ((round(d, 2) if capas[0] == capas[1] else "7.7-13.4"),)
        cuenta[clave] = cuenta.get(clave, 0) + 1
    info = dict(J_debil=J_debil, J_perp=J_perp,
                J0={n: float(J0[NOMBRES.index(n)]) for n in NOMBRES},
                lambda_max=float(np.linalg.eigvalsh(A_q0).max()),
                enlaces_por_capa_y_distancia={f"{c1}-{c2} {d} Å": v // 2
                                              for (c1, c2, d), v in sorted(cuenta.items(), key=str)})
    return cell, pos, lista, info


if __name__ == "__main__":
    cell_, pos_, lista_, info_ = construir()
    print(f"J_debil = {info_['J_debil']:.2f} meV (informe: entre 8.9 y 16.9)")
    print(f"J_perp  = {info_['J_perp']:.3f} meV por par superior-inferior")
    print("J0 por sitio (meV):", {k: round(v, 2) for k, v in info_["J0"].items()})
    print(f"lambda_max J(q=0) = {info_['lambda_max']:.2f} meV (informe: 49.2)")
    print("Enlaces únicos por capa y distancia:")
    for k, v in info_["enlaces_por_capa_y_distancia"].items():
        print(f"   {k}: {v}")
