"""
================================================================
 red.py — La red magnética: enlaces, supercelda y libreta de vecinos
================================================================

IDEA GENERAL
------------
Toda la información del material se resume en una LISTA DE ENLACES de la
celda unidad. Cada enlace es

        (k, k', n1, n2, J)

y significa: "el átomo k de la celda (0,0) interactúa con el átomo k' de la
celda (n1, n2) con un acoplamiento J (meV)". Esta es exactamente la forma en
que TB2J y los ajustes DFT (energy mapping) entregan los acoplamientos, y
tiene dos ventajas sobre "una distancia por J":

  1. Enlaces a la misma distancia pero cristalográficamente distintos pueden
     tener J distintos (p. ej. 3.133 y 3.146 Å en el Ti6).
  2. Al conocer (n1, n2) explícitamente, la supercelda periódica se arma con
     aritmética modular: no hace falta buscar vecinos por distancia en la
     supercelda, y se puede verificar que no hay imágenes duplicadas.

La lista se puede construir de dos maneras (ver parametros.py):
  - modo "anillos": se dan distancias y J por anillo, y el programa busca
    todos los pares (k, k', n1, n2) cuya distancia cae en cada anillo.
  - modo "enlaces": se da la lista explícita (k, k', n1, n2, J).

CONVENCIÓN DEL HAMILTONIANO (la que usa TODO el programa)
--------------------------------------------------------
        H = - sum_{<i,j>}  J_ij  S_i S_j ,      S_i = ±1

  * cada PAR se cuenta UNA sola vez,
  * J > 0  ferromagnético,  J < 0  antiferromagnético,
  * J en meV.

Si sus J vienen con otra convención, conviértalos con `factor_J` en
parametros.py (por ejemplo: factor 2 si su fuente usa la suma doble
sum_{i≠j}, o el cuadrado del espín si su mapeo usó |S| ≠ 1).

ÍNDICES DE SITIO EN LA SUPERCELDA L1 x L2
-----------------------------------------
        sitio(n1, n2, k) = (n1 * L2 + n2) * nat + k

que coincide con el orden de ase.Atoms.repeat((L1, L2, 1)).
"""

import numpy as np


# ============================ ENLACES DE LA CELDA UNIDAD ===================

def anchos_celda(cell):
    """
    Anchos perpendiculares de la celda 2D (Å): distancia entre las dos caras
    paralelas a b (ancho en la dirección de a) y entre las dos caras paralelas
    a a (ancho en la dirección de b). Sirven para saber cuántas celdas vecinas
    hay que recorrer para encontrar todos los pares hasta un cierto corte.
    """
    a, b = np.asarray(cell[0][:2], float), np.asarray(cell[1][:2], float)
    area = abs(a[0] * b[1] - a[1] * b[0])
    return area / np.linalg.norm(b), area / np.linalg.norm(a)


def _normalizar_anillos(anillos):
    """Acepta d como número o lista de números; devuelve lista de (ds, J)."""
    salida = []
    for an in anillos:
        ds = np.atleast_1d(np.asarray(an["d"], float))
        salida.append((ds, float(an["J"])))
    return salida


def enlaces_por_anillos(cell, posiciones, anillos, tol, factor_J=1.0,
                        avisar_no_asignados=True):
    """
    Construye la lista de enlaces buscando, para cada par de átomos (k, k')
    y cada traslación (n1, n2), si su distancia cae en algún anillo.

    anillos : lista de diccionarios {"d": distancia o lista de distancias (Å),
                                     "J": acoplamiento (meV)}
    tol     : tolerancia en Å para asignar una distancia a un anillo.

    Devuelve un diccionario con los arreglos
        k, kp, n1, n2  (enteros),  J (meV, ya multiplicado por factor_J),
        d (Å), anillo (índice del anillo)
    y la lista de avisos encontrados.
    """
    cell = np.asarray(cell, float)
    pos = np.asarray(posiciones, float)
    nat = len(pos)
    anillos_n = _normalizar_anillos(anillos)
    avisos = []

    # --- Las ventanas de distancia no deben solaparse ---
    todas = sorted((d, ia) for ia, (ds, _) in enumerate(anillos_n) for d in ds)
    for (d1, a1), (d2, a2) in zip(todas[:-1], todas[1:]):
        if a1 != a2 and d2 - d1 < 2 * tol:
            raise ValueError(
                f"Las ventanas de los anillos {a1+1} (d={d1}) y {a2+1} (d={d2}) "
                f"se solapan con tol={tol} Å. Reduzca 'tol'.")

    d_max = max(d for d, _ in todas) + tol
    w1, w2 = anchos_celda(cell)
    m1, m2 = int(np.ceil(d_max / w1)) + 1, int(np.ceil(d_max / w2)) + 1

    lista = []
    no_asignadas = []
    for k in range(nat):
        for kp in range(nat):
            for n1 in range(-m1, m1 + 1):
                for n2 in range(-m2, m2 + 1):
                    if k == kp and n1 == 0 and n2 == 0:
                        continue
                    r = pos[kp] + n1 * cell[0] + n2 * cell[1] - pos[k]
                    d = np.linalg.norm(r)
                    if d > d_max:
                        continue
                    asignado = False
                    for ia, (ds, J) in enumerate(anillos_n):
                        if np.any(np.abs(ds - d) < tol):
                            lista.append((k, kp, n1, n2, J * factor_J, d, ia))
                            asignado = True
                            break
                    if not asignado:
                        no_asignadas.append(d)

    if avisar_no_asignados and no_asignadas:
        u = np.unique(np.round(no_asignadas, 3))
        avisos.append(
            "Hay pares por debajo del corte que NO cayeron en ningún anillo "
            f"(distancias en Å: {u.tolist()}). Si deberían tener J, agregue "
            "esas distancias a un anillo o aumente 'tol'.")

    if not lista:
        raise ValueError("No se encontró ningún enlace. Revise d y tol.")
    enl = _a_arreglos(lista)
    return enl, avisos


def enlaces_explicitos(cell, posiciones, lista, factor_J=1.0):
    """
    Lista explícita de enlaces: cada elemento es (k, k', n1, n2, J_meV).

    Se aceptan dos formas de entrada:
      * cada par aparece en LOS DOS sentidos (k->k',R) y (k'->k,-R), como en
        la salida de TB2J; los dos J deben coincidir.
      * cada par aparece UNA sola vez; el programa agrega el sentido inverso.
    Mezclar las dos formas es un error (daría doble conteo parcial).
    """
    cell = np.asarray(cell, float)
    pos = np.asarray(posiciones, float)
    dic = {}
    for (k, kp, n1, n2, J) in lista:
        clave = (int(k), int(kp), int(n1), int(n2))
        if clave in dic:
            raise ValueError(f"Enlace repetido en la lista: {clave}")
        dic[clave] = float(J)

    con_inverso = sum(1 for (k, kp, n1, n2) in dic if (kp, k, -n1, -n2) in dic)
    if con_inverso == len(dic):
        for (k, kp, n1, n2), J in dic.items():
            Jinv = dic[(kp, k, -n1, -n2)]
            if not np.isclose(J, Jinv, rtol=1e-6, atol=1e-9):
                raise ValueError(f"J asimétrico: {(k,kp,n1,n2)}={J}, inverso={Jinv}")
    elif con_inverso == 0:
        for (k, kp, n1, n2), J in list(dic.items()):
            dic[(kp, k, -n1, -n2)] = J
    else:
        raise ValueError(
            "La lista de enlaces mezcla pares con y sin su sentido inverso. "
            "Dé todos los pares en los dos sentidos, o todos en uno solo.")

    filas = []
    for (k, kp, n1, n2), J in dic.items():
        if k == kp and n1 == 0 and n2 == 0:
            raise ValueError("Un átomo no puede acoplarse consigo mismo en (0,0).")
        d = np.linalg.norm(pos[kp] + n1 * cell[0] + n2 * cell[1] - pos[k])
        filas.append((k, kp, n1, n2, J * factor_J, d, -1))
    enl = _a_arreglos(filas)
    # Agrupar en "anillos" por valor de J (solo para los reportes)
    _, enl["anillo"] = np.unique(np.round(enl["J"], 9), return_inverse=True)
    return enl, []


def _a_arreglos(filas):
    filas = sorted(filas, key=lambda f: (f[6], f[0], f[1], f[2], f[3]))
    a = np.array(filas, dtype=object)
    return {
        "k":      a[:, 0].astype(np.int64),
        "kp":     a[:, 1].astype(np.int64),
        "n1":     a[:, 2].astype(np.int64),
        "n2":     a[:, 3].astype(np.int64),
        "J":      a[:, 4].astype(np.float64),
        "d":      a[:, 5].astype(np.float64),
        "anillo": a[:, 6].astype(np.int64),
    }


# ============================ VERIFICACIONES ================================

def coordinacion_por_subred(enl, nat):
    """
    Tabla z[k, anillo]: cuántos vecinos de cada anillo tiene el átomo k.
    A diferencia del promedio, muestra si los sitios NO son equivalentes
    (en el Ti6 original, los Ti 0 y 3 no tenían ningún vecino J2).
    """
    n_an = enl["anillo"].max() + 1
    z = np.zeros((nat, n_an), dtype=int)
    for k, a in zip(enl["k"], enl["anillo"]):
        z[k, a] += 1
    return z


def L_minimo(enl):
    """
    Tamaño mínimo de la supercelda para que dos traslaciones distintas de un
    enlace nunca caigan en la misma imagen (si no, aparecen vecinos
    duplicados, como pasaba con L=1). Condición: L > 2 * max|n|.
    """
    return 2 * int(max(np.abs(enl["n1"]).max(), np.abs(enl["n2"]).max())) + 1


def validar_tamano(L1, L2, enl):
    Lmin = L_minimo(enl)
    if L1 < Lmin or L2 < Lmin:
        raise ValueError(
            f"Supercelda {L1}x{L2} demasiado pequeña: con estos enlaces el "
            f"mínimo es {Lmin}x{Lmin} para no duplicar vecinos por imágenes "
            f"periódicas. Para escalamiento de tamaño finito use L >= 12.")


# ============================ SUPERCELDA ====================================

def indice(n1, n2, k, L2, nat):
    return (n1 * L2 + n2) * nat + k


def construir_libreta(L1, L2, nat, enl, permitir_plegado=False):
    """
    Libreta de vecinos de la supercelda L1 x L2 con condiciones periódicas.

    Devuelve
        vec  (N, z_max) int64  : índices de los vecinos de cada sitio
        Jv   (N, z_max) float64: J (meV) de cada uno de esos vecinos
        nvec (N,)       int64  : cuántos vecinos tiene realmente cada sitio
    Las filas se rellenan con J = 0 hasta z_max, pero los bucles del Monte
    Carlo solo recorren los primeros nvec[i] elementos.

    permitir_plegado=True solo se usa para enumerar estados fundamentales en
    celdas magnéticas pequeñas (ahí sí es correcto que una traslación caiga
    sobre el mismo sitio: es la energía de una configuración periódica).
    """
    if not permitir_plegado:
        validar_tamano(L1, L2, enl)
    N = L1 * L2 * nat
    z_por_k = np.bincount(enl["k"], minlength=nat)
    z_max = int(z_por_k.max())
    vec = np.zeros((N, z_max), dtype=np.int64)
    Jv = np.zeros((N, z_max), dtype=np.float64)
    nvec = np.zeros(N, dtype=np.int64)
    for n1 in range(L1):
        for n2 in range(L2):
            for e in range(len(enl["k"])):
                k = enl["k"][e]
                i = indice(n1, n2, k, L2, nat)
                j = indice((n1 + enl["n1"][e]) % L1, (n2 + enl["n2"][e]) % L2,
                           enl["kp"][e], L2, nat)
                c = nvec[i]
                vec[i, c] = j
                Jv[i, c] = enl["J"][e]
                nvec[i] += 1
    if not permitir_plegado:
        # Comprobación de seguridad: ningún vecino repetido ni autointeracción
        for i in range(N):
            v = vec[i, :nvec[i]]
            if len(np.unique(v)) != len(v) or np.any(v == i):
                raise RuntimeError(f"Vecinos duplicados en el sitio {i}.")
    return vec, Jv, nvec


def matriz_acoplamiento(L1, L2, nat, enl):
    """Matriz densa J_ij (meV) de una supercelda pequeña, con plegado."""
    vec, Jv, nvec = construir_libreta(L1, L2, nat, enl, permitir_plegado=True)
    N = L1 * L2 * nat
    A = np.zeros((N, N))
    for i in range(N):
        for c in range(nvec[i]):
            A[i, vec[i, c]] += Jv[i, c]
    return A


def posiciones_supercelda(L1, L2, cell, posiciones):
    cell = np.asarray(cell, float)
    pos = np.asarray(posiciones, float)
    out = []
    for n1 in range(L1):
        for n2 in range(L2):
            for k in range(len(pos)):
                out.append(pos[k] + n1 * cell[0] + n2 * cell[1])
    return np.array(out)


def distancias_celda(cell, posiciones, d_max):
    """Todas las distancias k -> k'+R menores que d_max (para el histograma)."""
    cell = np.asarray(cell, float)
    pos = np.asarray(posiciones, float)
    w1, w2 = anchos_celda(cell)
    m1, m2 = int(np.ceil(d_max / w1)) + 1, int(np.ceil(d_max / w2)) + 1
    ds = []
    for k in range(len(pos)):
        for kp in range(len(pos)):
            for n1 in range(-m1, m1 + 1):
                for n2 in range(-m2, m2 + 1):
                    if k == kp and n1 == 0 and n2 == 0:
                        continue
                    d = np.linalg.norm(pos[kp] + n1 * cell[0] + n2 * cell[1] - pos[k])
                    if d <= d_max:
                        ds.append(d)
    return np.array(ds)


def resumen_enlaces(enl, nat):
    """Texto con la tabla de coordinación por subred y los J por anillo."""
    z = coordinacion_por_subred(enl, nat)
    n_an = z.shape[1]
    lineas = ["Anillos (J en meV, ya con factor_J):"]
    for a in range(n_an):
        m = enl["anillo"] == a
        ds = np.unique(np.round(enl["d"][m], 3))
        Js = np.unique(np.round(enl["J"][m], 6))
        lineas.append(f"  anillo {a+1}: J = {Js.tolist()}  d = {ds.tolist()} Å")
    lineas.append("Coordinación por subred z[k, anillo]:")
    lineas.append("   k  " + "  ".join(f"z{a+1}" for a in range(n_an)))
    for k in range(nat):
        lineas.append(f"  {k:2d}  " + "  ".join(f"{z[k,a]:2d}" for a in range(n_an)))
    if len({tuple(f) for f in z}) > 1:
        lineas.append("  AVISO: los sitios no son equivalentes (coordinación distinta).")
    return "\n".join(lineas)
