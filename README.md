# Monte Carlo de Ising para MXenes 2D

Predice la temperatura de transición magnética (Curie o Néel) de un material 2D
a partir de su red y de los acoplamientos de canje *J* obtenidos por DFT.

**La idea en una línea:** red + *J* → estado fundamental → Monte Carlo en varios
tamaños *L* → cruce del cumulante de Binder → extrapolación a *L* → ∞ → *T*c en K.

## Uso

```bash
pip install -r requirements.txt
python validacion_onsager.py   # 1) probar el código contra la solución exacta
python run.py                  # 2) simular el material definido en parametros.py
python plot.py                 # 3) gráficas en resultados/
```

Lo **único** que se edita es `Proyecto_MXenes/parametros.py`: celda, posiciones de
los átomos magnéticos, los *J* (por anillos de distancia o como lista explícita
de enlaces tipo TB2J) y los parámetros de la simulación.

## Convención

H = − Σ<sub>⟨ij⟩</sub> *J*<sub>ij</sub> S<sub>i</sub> S<sub>j</sub>, con S = ±1,
cada par contado **una** vez, *J* > 0 ferromagnético, *J* en meV.
Si sus *J* vienen con otra convención use `factor_J` (p. ej. 2 si la fuente
suma sobre i ≠ j).

## Estructura (`Proyecto_MXenes/`)

| Archivo | Qué hace |
|---|---|
| `parametros.py` | Datos del material y de la simulación (único archivo a editar) |
| `red.py` | Lista de enlaces (k, k′, n1, n2, J), coordinación por subred, supercelda y libreta de vecinos verificadas |
| `estado_fundamental.py` | Enumeración exacta del orden a T = 0: patrón ε, degeneración, subredes desacopladas, FM/AFM/ferri, campo medio |
| `montecarlo.py` | Metrópolis, Wolff e híbrido (numba); series temporales E(t), M(t) |
| `analisis.py` | τ<sub>int</sub> (Madras–Sokal), errores jackknife, cruces de Binder con error, extrapolación ponderada |
| `run.py` | Flujo completo: pasada gruesa → pasada fina (con re-centrado) → *T*c ± error |
| `validacion_onsager.py` | Red cuadrada FM y AFM: *T*c = 2.269 \|J\|/k<sub>B</sub>, U* ≈ 0.611 |
| `plot.py` | U<sub>L</sub>, χ, C, ⟨\|M\|⟩ con barras de error; extrapolación; histograma de distancias |
| `MonteCarlo.ipynb` | Cuaderno que usa los módulos (sin código duplicado) |
| `Benchmark.ipynb` | Python puro vs Numba |

## Verificaciones que hace el programa

- Ventanas de distancia que se solapan o pares sin asignar dentro del corte.
- Tamaño mínimo de supercelda (sin vecinos duplicados por imágenes periódicas).
- Estado fundamental degenerado o subredes con campo de canje nulo: **se detiene**
  (no hay un orden único que medir).
- *L* múltiplo del periodo magnético.
- Bloques del jackknife largos frente a τ<sub>int</sub>.
- Cruces de Binder con la orientación física y significativos frente al ruido.
- Criterio de validez de Ising: MAE ≳ 10 k<sub>B</sub>*T*c (si se da `MAE_meV`).

## Limitaciones

- Modelo de Ising (espines ±1): solo cuantitativo con anisotropía de eje fácil
  fuerte. Si no, Ising sobreestima *T*c y hace falta Heisenberg con anisotropía.
- El valor universal U* ≈ 0.611 es para la red cuadrada isótropa; con otra
  geometría de celda U* cambia y sirve solo como guía.
