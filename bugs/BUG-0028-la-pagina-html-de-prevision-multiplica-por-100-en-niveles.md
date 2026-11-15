---
id: BUG-0028
title: La página HTML de previsión multiplica por 100 las variaciones sea cual sea λ — en niveles el cambio sale cien veces mayor bajo una columna «(%)»
status: open
severity: medium
component: report_forecast
found_in: 0.1.21
fixed_in:
reported: 2026-11-15
reporter: David / Claude — revisión de incoherencias entre componentes antes de mtram
tags:
  - forecast
  - report
  - scale
  - lambda-1
references:
  - src/fue/report_forecast.py:_table_data (variaciones, ERR y level_std·100 de la tabla)
  - src/fue/report_forecast.py:_make_charts_svg_legacy (la tasa anual del panel superior, ERR y σ del inferior)
  - src/fue/forecast.py:562-565 (los campos diff1, seasonal_diff y sus d.t. de ForecastResult)
  - BUG-0027 (el mismo defecto en la tabla del .out, arreglado en 0.1.21)
  - atsw-gui engines/fuf BUG-0001 (el del C)
  - drtran src/drtran/report.py (adapta la previsión conjunta a esta página: hereda el defecto)
---

## Summary

`write_forecast_report` (la página HTML de fuf, y la de mtram, que la usa a
través de `drtran.report`) escribe las variaciones de periodo y anuales, sus
d.t., ERR y la d.t. del nivel como `100·Δ(refactor·BoxCox(y))/refactor`. Con
λ = 0 es un tanto por ciento. **Con λ = 1 es el cambio del nivel por cien**,
bajo una columna rotulada «Monthly (%)». BUG-0027 arregló la tabla del `.out`
(vs = 100/refactor con λ = 0, 1/refactor si no); la página y su gráfico siguen
en la regla vieja.

Serie sintética en torno a 100 en niveles, λ = 1, refactor = 100:

```
última observación: 98.1090 -> 100.5195, cambio +2.4104
columna «Monthly (%)» de la página: +241.04
primera previsión: nivel 100.12, su d.t. 87.87, variación −39.57 (d.t. 87.87)
```

Y las primeras filas observadas escriben como «variación anual» el nivel
transformado (k ≤ freq), el defecto que BUG-0027 quitó del `.out`.

## Impact

Medio: la página es lo que el analista lee y comparte; en niveles todas las
variaciones y d.t. salen cien veces mayores y la etiqueta dice «%». El `.out`
y la página de la misma previsión ya no dicen lo mismo desde 0.1.21.

## Reproduction

`bugs/BUG-0028-repro/repro.py`: AR(1) sintético en torno a 100, λ = 1,
refactor = 100; las filas que pinta la página (`_table_data`, la entrada de la
plantilla). Sale con 1 mientras la variación de la última observada sea cien
veces el cambio del nivel.

## Root cause

Los campos de `ForecastResult` (`diff1`, `seasonal_diff` y sus d.t.) se definen
como «cambios en porcentaje» —100·Δ/refactor— para todo λ, y la página los pinta
tal cual y calcula las filas observadas con la misma regla.

## Fix

La escala de BUG-0027 (`(100 si λ = 0, si no 1)/refactor`) en la página y el
gráfico, sin cambiar el significado de los campos de `ForecastResult` (los usan
drtran y pyfug): o se aplica al pintar, o `ForecastResult` gana un campo que
diga en qué unidades vienen. A decidir al arreglarlo; drtran `report.py` se
revisa en la misma tanda.

## Validation

Pendiente.
