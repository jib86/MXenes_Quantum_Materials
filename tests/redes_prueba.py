"""Redes pequeñas con solución conocida para las pruebas."""
import numpy as np

from isingmc.red import enlaces_por_anillos

CUADRADA = np.diag([1.0, 1.0, 20.0])
TRIANGULAR = np.array([[1.0, 0.0, 0.0], [0.5, np.sqrt(3) / 2, 0.0], [0.0, 0.0, 20.0]])
ORIGEN = np.array([[0.0, 0.0, 0.0]])


def cuadrada(J=1.0):
    enl, _ = enlaces_por_anillos(CUADRADA, ORIGEN, [{"d": 1.0, "J": J}], 0.05)
    return enl


def triangular(J=-1.0):
    enl, _ = enlaces_por_anillos(TRIANGULAR, ORIGEN, [{"d": 1.0, "J": J}], 0.05)
    return enl
