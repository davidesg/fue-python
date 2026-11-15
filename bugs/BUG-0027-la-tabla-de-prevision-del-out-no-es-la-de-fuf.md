---
id: BUG-0027
title: La tabla de previsión del .out no es la de fuf — en niveles las variaciones y ERR salen ×100, los campos van pegados y las series anuales llenan columnas que no existen
status: fixed
severity: high
component: report
found_in: 0.1.20
fixed_in: 0.1.21 (unreleased)
reported: 2026-11-15
reporter: David / Claude — auditoría de la suite antes de mtram
tags:
  - forecast
  - report
  - fuf
  - legado-c
references:
  - src/fue/report.py (_fuf_section_forecast_table)
  - src/fue/forecast.py (forecast: los vectores del C, para la tabla)
  - atsw-gui engines/fuf/BUGS.md (BUG-0001, BUG-0002, BUG-0003, arreglados el 15-nov-2026)
  - atsw-gui engines/fuf/src/usfo.c (forecast_table_acii)
  - atsw-gui docs/TODO.md ("fue-python, la tabla de fuf (0.1.21), lo primero")
  - tests/test_bug_0027_tabla_de_fuf.py, tests/data/tablas_de_fuf/
---

## Summary

`write_fuf_out` escribe la tabla de previsión de fuf (`report._fuf_section_forecast_table`).
fuf arregló la suya el 15-nov-2026 (sus BUG-0001, 0002 y 0003) y el porte
seguía en la de antes, con sus defectos y alguno propio:

- **la escala:** multiplicaba por 100 las variaciones sea cual sea λ. En
  logaritmos es lo que se quiere (porcentaje); **en niveles (λ = 1) las deja
  cien veces más grandes**. ERR salía sin escalar (el residuo por refactor) y
  la d.t. del nivel como √v₁, sin la escala;
- **los campos pegados:** anchos fijos sin separador (`{:9.4f}{:7.4f}…`): un
  número que no cabe se junta con el de al lado;
- **la variación anual de las primeras observadas:** con k ≤ freq escribía el
  nivel transformado (`100·BoxCox(z_k)/refactor`) donde no hay dato un año
  antes; la del primer periodo, 0;
- **las series anuales:** llenaba las columnas anuales, que con un periodo por
  año repetirían las del periodo.

En `CPI_USA_model_AR` (λ = 1, refactor = 100), diciembre de 2017:

```
fuf : 12/2017    246.5240        -   -0.1450        -    5.0920        -    0.1899
fue : 12/2017 246.5240     - -14.5000      -  509.2000     -  18.9852
```

## Impact

Alto en niveles: la variación, sus d.t., la d.t. del nivel y ERR salen cien
veces más grandes, y nada lo dice. En logaritmos con refactor = 100 (el caso
habitual) los números ya eran los de fuf; cambian el formato y la columna
anual de las primeras filas. La tabla la leen a máquina el gestor de
previsiones y `lib/outfcst`, que parten por espacios.

`ForecastResult` no cambia: sus campos siguen siendo los de siempre (cambios
en porcentaje de la serie transformada).

## Reproduction

`bugs/BUG-0027-repro/repro.py`: dos casos, cada uno con el `.inp` de fuf y la
tabla que fuf escribe (en niveles, `CPI_USA_model_AR`; anual, `MINK_m01`).
Sale con 1 mientras la tabla no sea la de fuf; con 0 cuando lo es.

## Root cause

El porte copió la tabla de fuf-1.08 (`forecast_table_acii` antes del
15-nov-2026), que ya tenía los tres defectos, y añadió el de la columna anual.

## Fix

`_fuf_section_forecast_table` escribe ahora lo que escribe `forecast_table_acii`:
vs = 100/refactor con λ = 0 y 1/refactor si no; la cabecera de 76 columnas de
drtran, drvarma y drvec y su borde de cierre; un espacio garantizado entre
campos; «-» donde no hay dato; la fecha «1912» con un periodo por año. Los
números de las previsiones salen de los vectores del propio cálculo (f₁, f₂,
f₃, v₁, v₂, v₃), que `forecast` deja en el resultado: reconstruirlos de los
campos de `ForecastResult` movía la última cifra.

Con una Box-Cox negativa fuf escribe además una tabla propia (el nivel con
bandas de dos d.t.). No se porta: fue no la tuvo nunca.

## Validation

`tests/test_bug_0027_tabla_de_fuf.py` con 17 casos (`tests/data/tablas_de_fuf/`):
los 14 de `engines/fuf/tests/golden`, los dos de su BUG-0001 con el factor
bien puesto y el visón anual.

- Byte a byte en Linux x86_64 con glibc y el motor C: 11 casos. Los otros
  seis son ficheros del corpus con factor de reescalado 1 y parámetros
  estimados con 100, cuyas previsiones explotan: el porte y el C se apartan
  en la última cifra (4e-5 relativo como mucho) y en el signo de un cero.
- En todas partes y en los 17: las mismas filas y campos, los números a
  1e-4.

`tests/real_cases/forecast_S.2.out` regenerado con el fuf nuevo: solo cambia
la tabla.
