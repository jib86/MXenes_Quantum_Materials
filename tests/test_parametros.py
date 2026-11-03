from isingmc.flujo import config_desde_modulo, preparar
from isingmc.red import enlaces_explicitos


def test_parametros_se_cargan_y_la_red_es_valida():
    import parametros
    cfg = config_desde_modulo(parametros)
    if cfg["modo_J"] == "enlaces":
        enlaces_explicitos(cfg["cell"], cfg["posiciones"], cfg["enlaces"], cfg["factor_J"])
    prep = preparar(cfg, verbose=False)
    assert prep["nat"] == len(cfg["posiciones"])
