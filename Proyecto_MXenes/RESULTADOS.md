# Resultados registrados

## Validación contra Onsager (red cuadrada, |J| = 1 meV, exacto 26.333 K)

| Caso | T_c (x = 1) | T_c (x = 3) | U* | Error (x = 1) |
|---|---|---|---|---|
| FM  | 26.236 ± 0.043 K | 26.298 ± 0.014 K | 0.612 | 0.37 % |
| AFM | 26.404 ± 0.051 K | 26.343 ± 0.020 K | 0.615 | 0.27 % |

Criterio (< 0.5 %): aprobado en los dos. La forma x = 3 (1/ν + ω) queda más
cerca del exacto en ambos casos (0.13 % y 0.04 %).

Repetición tras cambiar la calibración de Wolff (commit "Wolff: calibración..."):
FM 26.337 ± 0.084 K (0.02 %), AFM 26.235 ± 0.090 K (0.37 %); con x = 3,
26.312 ± 0.037 y 26.284 ± 0.037 K. Aprobado.

## MXene Nb–Ti 30 % desnuda — modelo de Ising con J reconstruidos (2 nov 2026)

- J: `modelo_desnuda.py`, reconstruidos del informe TB2J de Luis (13 oct 2026);
  `factor_J = 2` (convención TB2J i ≠ j). Control: λ_max J(q=0) = 49.19 meV
  (informe: 49.2).
- Estado fundamental encontrado por el código: FM, no degenerado.
- L = 6, 9, 12, 18 (supercelda L x 2L), híbrido, N_meas = 10 000.
- **T_c (Ising) = 592 ± 4 (estad.) ± 1 (sist.) K**; campo medio Ising 1142 K.
- Vampire (Heisenberg + anisotropía, mismo informe): 170–213 K.
- La diferencia es la esperada: MAE / (k_B T_c) ≈ 0.02–0.2, muy lejos del
  criterio de Ising (≳ 10). Con esta anisotropía débil el modelo de Ising no
  es cuantitativo para este material; da una cota superior.
- tau_int llegó a ~54 medidas en L = 18 cerca de T_c: para producción usar
  N_meas >= 50 000.
