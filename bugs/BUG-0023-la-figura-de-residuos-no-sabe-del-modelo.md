---
id: BUG-0023
title: La figura de residuos de un modelo no sabe del modelo — la Q resta todos los parámetros, el primer residuo lleva la fecha del primer dato y los retardos no son los de fug C
status: fixed
severity: high
component: plots
found_in: 0.1.14
fixed_in: 0.1.16
reported: 2026-09-25
reporter: David — revisión de las figuras de art (BUG-0165 de art)
tags:
  - plots
  - ljung-box
  - degrees-of-freedom
  - dates
references:
  - art/bugs/BUG-0166 (el mismo df en el texto de art; su parte en fue quedó sin registrar)
  - art/bugs/BUG-0165 (dos renderizadores de la misma diagnosis)
  - art/bugs/BUG-0172, BUG-0185 (el desfase de ifadf)
  - fug-1.14-proto/src/diagnose.c, default_lags()
  - tests/test_bug_0023_lo_que_los_residuos_de_un_modelo_necesitan.py
---

## Summary

`plot_model_diagnostics(model)` dibuja los residuos de un modelo con su ACF/PACF
y la Q de Ljung-Box. Tres cosas que no están en los residuos, sino en el modelo,
salían mal:

1. **Los grados de libertad de la Q.** El rótulo era `Q(lags − npar)` con
   `npar = model._result.npar`, que cuenta **todos** los parámetros: armónicos,
   intervenciones y media incluidos. La corrección de Ljung-Box sólo resta los
   parámetros ARMA **estimados**.
2. **La fecha del primer residuo.** `plot_residuals_ts` fechaba el primer
   residuo en `model.series.start`, pero la diferenciación consume
   d + D·s + las raíces de `ifadf`: el eje quedaba corrido hacia atrás.
3. **El número de retardos por defecto.** `_default_lags` no era la regla de
   fug C: le faltaba la rama anual larga (n > 200 → 45) y el tope n − 2. El
   `.out` (`report.py`) sí la tenía, escrita a mano en otro sitio.

## Impact

Alto: la etiqueta de la Q es un **número publicado** que se lee para decidir si
un modelo está terminado. Medido sobre el run 3 de SF_MEG (IPC de España):

| modelo | parámetros | ARMA estimados | rótulo | df correcto |
|---|---|---|---|---|
| A_m00 (10 armónicos, AR(1) fijo en 0) | 11 | 0 | **Q(28)** | 39 |
| A_m02 (10 armónicos, MA(1) libre) | 12 | 1 | **Q(27)** | 38 |

Lo ve cualquier usuario de `fue` (`Model.plot_residuals`) y, a través de
`plot_diagnosis`, las herramientas de art que dibujan con fue: `record_version`,
`build_model`, `full_report` y `save_diagnosis_report`.

## Reproduction

```python
import fue
from fue.plots import plot_model_diagnostics
m = fue.load("tests/data/bug_0023/ES_CPI_A_m00.pre")[1]; m.fit()
fig, _ = plot_model_diagnostics(m)
fig.axes[1].get_xlabel()          # 'Q(28) = 75.1'  — df correcto 39
```

## Root cause

**Una regla con tres copias, y un conocimiento del modelo que no vivía en
fue.** El contador correcto de ARMA libres existía (`report._count_nparma`, el
del `.out`), pero `plots.py` usaba `_result.npar`. La regla de retardos estaba
dos veces en fue y una en art, distintas. Y la cuenta del desfase de la
diferenciación vivía en art (`desfase_observaciones`), así que fue no podía
fechar sus propios residuos.

La parte de fue de art/BUG-0166 se anotó en aquel informe («Y `fue/plots.py`
igual — pero eso es un cambio en fue, no en art») y **nunca se registró aquí**.

## Fix

Las cuatro piezas viven ahora en `fue.diagnostics`, una vez:

* `default_lags(nobs, freq)` — la regla de fug C, exacta (diagnose.c).
* `free_arma_count(model)` — ARMA regulares y estacionales **libres** más los
  factores de frecuencia fija; `report._count_nparma` delega en ella.
* `differencing_offset(model)` — d + D·s + raíces de `ifadf` (la cuenta de
  art/BUG-0172, que baja al dueño del modelo).
* `residuals_start(model)` — la fecha del primer residuo.

`plots.py` las usa: el rótulo es `Q(df)` con df = retardos − ARMA estimados, el
mismo convenio que pyfug; con df ≤ 0 dice `Q = … (df ≤ 0: no test)` en vez de
un número que se lee como contraste. `report.py` usa `default_lags`.

## Validation

`tests/test_bug_0023_lo_que_los_residuos_de_un_modelo_necesitan.py`: la regla de
fug C caso por caso; el recuento de ARMA en A_m00 (0), A_m02 (1) y A_m06 (2, con
un MA(2) de frecuencia fija); el desfase contra la longitud real de los
residuos (1, 3 y 5); los rótulos Q(39) y Q(38); la fecha del primer residuo en
B_m11; y la guarda de df ≤ 0. Los cuatro tests de la figura fallan con el
`plots.py` anterior.
