"""
Simula el material definido en parametros.py y reporta su T_c.

    python run.py
"""
import parametros
from isingmc.flujo import config_desde_modulo, ejecutar

if __name__ == "__main__":
    ejecutar(config_desde_modulo(parametros))
