# Método

Resumen del camino: **red + J → estado fundamental → Monte Carlo en varios L →
cruce del cumulante de Binder → extrapolación a L → ∞ → T_c**. Cada etapa tiene
una verificación que detiene el programa o emite un aviso si algo no cuadra.

## 1. Modelo y convención

$$H = -\sum_{\langle ij\rangle} J_{ij}\, S_i S_j, \qquad S_i = \pm 1$$

- Cada par se cuenta **una** vez; J > 0 ferromagnético; J en meV, T en K.
- `factor_J` convierte otras convenciones: **2 para TB2J** (E = −Σ_{i≠j} J S·S,
  espines unitarios), −1 si la fuente usa H = +Σ J S S, y |S|² si el mapeo DFT
  usó espines no unitarios.
- Validez: Ising describe un material con eje fácil fuerte. Criterio práctico,
  que el programa verifica si se da `MAE_meV`: **MAE ≳ 10 k_B T_c**. Si no se
  cumple, Ising sobreestima T_c.

## 2. Red y enlaces (`isingmc/red.py`)

Los acoplamientos se describen como enlaces de la celda unidad
**(k, k′, n1, n2, J)**: el átomo k de la celda (0, 0) con el átomo k′ de la
celda (n1, n2). Es el formato de TB2J y permite J distintos para enlaces a la
misma distancia. Se pueden dar explícitamente o por anillos de distancia.

Verificaciones: ventanas de distancia que se solapan, pares sin asignar dentro
del corte, coordinación **por subred** (el promedio oculta sitios no
equivalentes), tamaño mínimo de supercelda (L > 2 max|n|, si no hay vecinos
duplicados por imágenes periódicas) y ausencia de vecinos repetidos.

## 3. Estado fundamental (`isingmc/estado_fundamental.py`)

El parámetro de orden es M = (1/N) Σ ε_i S_i, con ε el estado fundamental. En
vez de suponerlo, se enumeran **exactamente** todas las configuraciones
periódicas en celdas magnéticas de hasta 24 espines (código Gray en numba).
El programa:

- elige la de menor energía por sitio y clasifica el orden (FM → Curie,
  AFM → Néel, ferri → Curie);
- **se detiene** si el estado fundamental es degenerado o si una subred siente
  campo de canje nulo (subred desacoplada: no hay un orden único que medir);
- calcula λ_max de J(q) solo como **cota superior** para fijar el rango del
  barrido en T (no se usa como estimación de T_c).

## 4. Monte Carlo (`isingmc/montecarlo.py`)

- **Metrópolis**: un barrido = N intentos de volteo.
- **Wolff**: se agrega un vecino al cluster solo si el enlace está satisfecho
  (J S_i S_j > 0), con p = 1 − exp(−2β|J|); válido para J de cualquier signo.
  El número de clusters por paso se calibra en la termalización (tamaño medio
  acumulado) y queda **fijo** durante la medición: una regla de parada que
  dependa del estado sesga los promedios.
- **Híbrido** (recomendado): un barrido de Metrópolis + los clusters de Wolff.

Se recorren las temperaturas de mayor a menor (annealing) y se guardan las
series completas E(t) y M(t).

## 5. Análisis (`isingmc/analisis.py`)

- **τ_int** de E y M² con la ventana automática de Madras–Sokal.
- **Errores jackknife** por bloques de U, χ, C, ⟨|M|⟩ y E, con aviso si los
  bloques son cortos frente a τ_int.
- **Cumulante de Binder** U_L = 1 − ⟨M⁴⟩ / (3⟨M²⟩²). Los cruces entre tamaños
  consecutivos se aceptan solo con la orientación física (el sistema grande
  más ordenado debajo de T_c) y si la diferencia es significativa a ambos
  lados; luego se ubican con **reponderación de histogramas**
  (Ferrenberg–Swendsen) y error jackknife.
- **Extrapolación** ponderada T_×(L_eff) = T_c + a L_eff^(−x), con
  L_eff = √(L L′). Se reportan x = 1 (protocolo) y x = 3 (1/ν + ω); su
  diferencia es el error sistemático.

## 6. Barrido en temperatura (`isingmc/flujo.py`)

1. Pasada gruesa con poca estadística.
2. Pasadas intermedias que estrechan la ventana alrededor de los cruces.
3. Pasada final con toda la estadística en T_c (1 ± `ventana_fina`), con
   re-centrado automático si los cruces caen en el borde.

## Referencias

- L. Onsager, Phys. Rev. 65, 117 (1944).
- K. Binder, Z. Phys. B 43, 119 (1981).
- U. Wolff, Phys. Rev. Lett. 62, 361 (1989).
- A. M. Ferrenberg y R. H. Swendsen, Phys. Rev. Lett. 61, 2635 (1988).
- N. Madras y A. D. Sokal, J. Stat. Phys. 50, 109 (1988).
- G. Kamieniarz y H. W. J. Blöte, J. Phys. A 26, 201 (1993): U* de Ising 2D.
- X. He et al., Comput. Phys. Commun. 264, 107938 (2021): TB2J.
