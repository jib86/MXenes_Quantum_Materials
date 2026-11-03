import numpy as np

from isingmc import K_B
from isingmc.estado_fundamental import buscar_estado_fundamental, construir_epsilon
from isingmc.red import enlaces_explicitos
from redes_prueba import CUADRADA, cuadrada, triangular


def test_cuadrada_FM():
    gs = buscar_estado_fundamental(cuadrada(1.0), 1, K_B)
    assert gs["tipo"] == "FM" and gs["n_clases"] == 1
    assert np.isclose(gs["e_sitio"], -2.0)


def test_cuadrada_AFM_tablero():
    gs = buscar_estado_fundamental(cuadrada(-1.0), 1, K_B)
    assert gs["tipo"] == "AFM" and gs["periodo"] == (2, 2) and gs["n_clases"] == 1
    eps = construir_epsilon(4, 4, 1, gs).reshape(4, 4)
    assert np.all(eps * np.roll(eps, 1, axis=0) == -1)


def test_triangular_AFM_degenerado():
    gs = buscar_estado_fundamental(triangular(-1.0), 1, K_B)
    assert gs["n_clases"] > 1 and any("DEGENERADO" in a for a in gs["avisos"])


def test_subred_desacoplada():
    # Subred 0: red cuadrada AFM (tablero). Subred 1, en el centro de cada
    # plaqueta, solo se acopla a sus 4 esquinas, cuya suma es cero en el
    # tablero: siente campo nulo (el defecto del código original).
    pos = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]])
    lista = [(0, 0, 1, 0, -1.0), (0, 0, 0, 1, -1.0),
             (1, 0, 0, 0, 0.3), (1, 0, 1, 0, 0.3), (1, 0, 0, 1, 0.3), (1, 0, 1, 1, 0.3)]
    enl, _ = enlaces_explicitos(CUADRADA, pos, lista)
    gs = buscar_estado_fundamental(enl, 2, K_B)
    assert any("NULO" in a for a in gs["avisos"])
