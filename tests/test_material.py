import numpy as np

from materiales import nbti30_desnuda


def test_desnuda_reproduce_informe():
    _, _, _, info = nbti30_desnuda.construir()
    assert abs(info["lambda_max"] - 49.2) < 0.2           # informe TB2J: 49.2 meV
    assert np.isclose(info["J0"]["Ti1"], 31.0)              # capa superior
    assert 49 <= info["J0"]["Ti3"] <= 53                    # capa inferior
    c = info["enlaces_por_capa_y_distancia"]
    assert c["inf-inf 3.13 Å"] == 12 and c["sup-sup 3.13 Å"] == 4
