import numpy as np
import pytest

from isingmc.red import (construir_libreta, coordinacion_por_subred, enlaces_explicitos,
                         enlaces_por_anillos, L_minimo, validar_tamano)
from redes_prueba import CUADRADA, ORIGEN, cuadrada


def test_coordinacion_cuadrada():
    enl = cuadrada()
    assert coordinacion_por_subred(enl, 1).tolist() == [[4]]
    assert L_minimo(enl) == 3


def test_supercelda_sin_duplicados():
    vec, Jv, nvec = construir_libreta(8, 8, 1, cuadrada())
    assert vec.shape == (64, 4) and np.all(nvec == 4)
    assert np.allclose(Jv, 1.0)


def test_rechaza_L_pequeno():
    with pytest.raises(ValueError):
        validar_tamano(2, 2, cuadrada())


def test_ventanas_solapadas():
    with pytest.raises(ValueError):
        enlaces_por_anillos(CUADRADA, ORIGEN, [{"d": 1.0, "J": 1}, {"d": 1.05, "J": 2}], 0.05)


def test_enlaces_explicitos_un_sentido_se_completa():
    lista = [(0, 0, 1, 0, 1.0), (0, 0, 0, 1, 1.0)]
    enl, _ = enlaces_explicitos(CUADRADA, ORIGEN, lista)
    assert len(enl["J"]) == 4


def test_enlaces_explicitos_mezcla_es_error():
    lista = [(0, 0, 1, 0, 1.0), (0, 0, -1, 0, 1.0), (0, 0, 0, 1, 1.0)]
    with pytest.raises(ValueError):
        enlaces_explicitos(CUADRADA, ORIGEN, lista)


def test_factor_J():
    enl, _ = enlaces_por_anillos(CUADRADA, ORIGEN, [{"d": 1.0, "J": 1.5}], 0.05, factor_J=2.0)
    assert np.allclose(enl["J"], 3.0)
