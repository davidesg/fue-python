---
id: BUG-0026
title: El Jarque-Bera del .out y de describe() trunca n/6 a entero — el puerto copió el C antiguo (nobs / 6, división entera)
status: fixed
severity: low
component: report
found_in: 0.1.16
fixed_in: 0.1.17 (unreleased)
reported: 2026-10-06
reporter: David / Claude — art/BUG-0214
tags:
  - jarque-bera
  - legado-c
references:
  - src/fue/report.py (_section_residual_stats)
  - src/fue/series.py (TimeSeries.describe)
  - atsw-gui engines/fue/src/diagnose.c (JarqueBera, ya con nobs / 6.0)
  - art-python/bugs/BUG-0214-el-jb-de-la-diagnosis-no-es-el-del-out.md
---

## Summary

`report.py` y `series.py` calculan `jb = (n // 6) * (S² + K²/4)`: copian el
`JarqueBera` del C antiguo (fue-1.14 `diagnose.c`), donde `nobs / 6` es una
división entera. El estadístico es n/6 · (S² + K²/4). Con n = 215 sale 35 en
lugar de 35,83: IPC_US AR(2) da 359,190 en el `.out` y 367,742 en la diagnosis
de art, que usa n/6.

## Impact

Bajo, y sólo cuando 6 no divide a n. Pero son dos cifras para el mismo
contraste, y el alumno no sabe cuál citar (art/BUG-0214).

## Root cause

Fidelidad a un defecto del C. El C de atsw-gui (`engines/fue`, `fuf`, `fug`)
ya lo corrigió a `nobs / 6.0`; scipy, pyfug y art usan n/6.

## Fix

`n / 6.0` en `report.py` y `series.py`. El `.out` y `describe()` dan el mismo JB
que scipy, pyfug, art y el C actual. Los `.out` escritos antes conservan el
valor truncado.

## Validation

`tests/test_bug_0026_jarque_bera_n_entre_6.py`: `describe()` de 215 datos da
n/6 · (S² + K²/4) y no ⌊n/6⌋; la fórmula es la de scipy. `test_api.py` cambia
su fórmula esperada.
