# Monte Carlo de Ising para materiales magnéticos 2D

Cálculo de la temperatura de transición magnética (**Curie o Néel**) de una
monocapa a partir de sus acoplamientos de canje *J* (por ejemplo, de TB2J), con
el cumulante de Binder y escalamiento de tamaño finito. Aplicado a MXenes
Nb–Ti.

```
red + J  →  estado fundamental  →  Monte Carlo (varios L)  →  cruce de Binder  →  L → ∞  →  T_c ± error
```

**Validado** contra la solución exacta de Onsager (red cuadrada FM y AFM):
error de 0.02–0.37 %, U* = 0.612–0.615. Ver [docs/resultados.md](docs/resultados.md).

**Documentos en PDF:**
[Guía del código](docs/Guia_del_codigo.pdf) (qué hace, cómo funciona y cómo correrlo paso a paso) ·
[Informe de correcciones](docs/Informe_de_correcciones.pdf) (qué estaba mal en el código original y qué se corrigió).

---

## Instalación

```bash
pip install -r requirements.txt      # numpy, scipy, numba, matplotlib
pip install pytest                   # opcional, para las pruebas
```

## Uso

```bash
python validacion_onsager.py   # 1. comprobar el código contra Onsager (~20 min)
python run.py                  # 2. simular el material definido en parametros.py
python plot.py                 # 3. gráficas en resultados/material/
pytest                         # pruebas unitarias (segundos)
```

Lo **único** que se edita es [`parametros.py`](parametros.py): celda, posiciones
de los átomos magnéticos, acoplamientos y parámetros de la simulación. El
cuaderno [`notebooks/MonteCarlo.ipynb`](notebooks/MonteCarlo.ipynb) hace el
mismo flujo paso a paso.

Salida típica de `run.py`:

```
Estado fundamental: periodo 1x1 celdas, orden FM  ->  la transición es una temperatura de Curie.
...
Cruces del cumulante de Binder:
  L=12 / 18: T_x = 590.975 ± 1.748 K, U* = 0.5983 ± 0.0155  [robusto]
Extrapolación (3 pares, x = 1): T_c(inf) = 591.724 ± 3.955 K
RESULTADO: temperatura de Curie = 591.72 ± 3.95 (estad.) ± 0.77 (sist.) K
  AVISO: MAE / (k_B T_c) = 0.0196: MENOR que 10: Ising sobreestima T_c; use Heisenberg con anisotropía.
```

## Dónde se ponen los J (valores de DFT)

En [`parametros.py`](parametros.py), en la sección marcada
**«ACOPLAMIENTOS DE CANJE J (meV) — AQUÍ SE ESCRIBEN LOS J DEL DFT»**. Cada línea
de la lista `enlaces` es un par de átomos y su J, con un comentario que dice qué
par es:

```python
enlaces = [
    # --- Capa inferior, primeros vecinos (3.13 Å) ---
    (2, 5,  0,  0,   16.60),   #  Ti3 – Ti6  misma celda
    (2, 4,  0,  1,   16.60),   #  Ti3 – Ti5  celda vecina (+0, +1)
    ...
]
```

- **k, k′**: los dos átomos, por su lugar en la lista `posiciones` (desde 0).
- **n1, n2**: celda del segundo átomo respecto del primero (columna `R` de TB2J).
- **J**: en meV, tal como lo da la fuente; `factor_J` hace la conversión (2 para TB2J).

Cada par se escribe una vez; un enlace que no aparece vale 0. Hoy la lista tiene
los 64 enlaces de la MXene desnuda, reconstruidos del informe TB2J
(ver [`materiales/nbti30_desnuda.py`](materiales/nbti30_desnuda.py)).

**Convención:** H = −Σ<sub>⟨ij⟩</sub> J<sub>ij</sub> S<sub>i</sub> S<sub>j</sub>,
S = ±1, cada par contado una vez, J > 0 ferromagnético, J en meV.
Para J de TB2J use `factor_J = 2`.

## Estructura

```
├── parametros.py            configuración del material y de la simulación (único archivo a editar)
├── run.py                   corre el flujo completo
├── plot.py                  gráficas de una corrida guardada
├── validacion_onsager.py    prueba contra la solución exacta (FM y AFM)
├── isingmc/                 paquete
│   ├── red.py                 enlaces, coordinación por subred, supercelda verificada
│   ├── estado_fundamental.py  orden a T = 0 por enumeración exacta; Curie o Néel
│   ├── montecarlo.py          Metrópolis, Wolff e híbrido (numba)
│   ├── analisis.py            τ_int, jackknife, cruces de Binder con reponderación, extrapolación
│   ├── flujo.py               barrido grueso → intermedio → final; T_c ± error
│   └── graficas.py            U_L, χ, C, ⟨|M|⟩, extrapolación, histograma de distancias
├── materiales/              modelos de materiales
│   └── nbti30_desnuda.py      MXene Nb–Ti 30 % desnuda (J de TB2J)
├── notebooks/               cuaderno de trabajo y rendimiento
├── tests/                   pruebas unitarias (pytest)
├── docs/
│   ├── Guia_del_codigo.pdf    guía en lenguaje llano: qué hace y cómo correrlo
│   ├── Informe_de_correcciones.pdf  qué estaba mal y qué se corrigió
│   ├── metodo.md              método, convenciones y verificaciones
│   ├── resultados.md          validación y resultados del material
│   └── figuras/
└── resultados/              salida de las simulaciones (no se versiona)
```

## Qué verifica el programa

| Etapa | Verificación |
|---|---|
| Red | ventanas de distancia solapadas, pares sin asignar, coordinación por subred, L mínimo, vecinos duplicados |
| Estado fundamental | degeneración y subredes desacopladas (**se detiene**), L múltiplo del periodo magnético |
| Monte Carlo | τ<sub>int</sub> frente al largo de los bloques del jackknife |
| Cruces de Binder | orientación física y significancia frente al ruido |
| Modelo | MAE ≳ 10 k<sub>B</sub>T<sub>c</sub> para que Ising sea cuantitativo |

## Alcance

- Es un modelo de **Ising** (espines ±1). Es cuantitativo con anisotropía de
  eje fácil fuerte. Con anisotropía débil, como en las MXenes Nb–Ti
  (MAE / k<sub>B</sub>T<sub>c</sub> ≈ 0.02–0.2), sobreestima T<sub>c</sub> y da una cota
  superior; ahí corresponde Heisenberg con anisotropía.
- El valor universal U* ≈ 0.611 es el de la red cuadrada isótropa; con otra
  geometría de celda U* cambia y sirve solo como guía.

Detalles del método y referencias en [docs/metodo.md](docs/metodo.md).

## Licencia

Distribuido bajo la licencia MIT. Ver [LICENSE](LICENSE).
