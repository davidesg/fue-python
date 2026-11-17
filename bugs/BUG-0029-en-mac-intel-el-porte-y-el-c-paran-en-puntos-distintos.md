---
id: BUG-0029
title: En macOS Intel el porte (raxopt) y el motor C paran en puntos distintos sobre un fichero fuera de la banda del optimizador
status: open
severity: low
component: qnewtopt
found_in: 0.1.22
fixed_in:
reported: 2026-11-17
reporter: David / Claude — preflight de atsw 2.1.3 con las suites (macos-15-intel)
tags:
  - optimizer
  - platform
  - homologation
  - scale
references:
  - tests/test_qnewtopt.py (test_raxopt_sigma2_matches_c[R.4], test_cov_raxopt_matches_c)
  - fue.scale (la banda [0.1, 1000]); atsw-gui docs/ESTUDIO-escala.md
  - art-python, run 44154374179 de preflight-atsw.yml, job «Suite fue — macos-15-intel · Python 3.12»
  - bugs/BUG-0029-repro/repro.py
---

## Resumen

`tests/test_qnewtopt.py` compara los dos motores de fue **en la misma
máquina**. En el preflight de atsw 2.1.3 (17-nov-2026, fue 0.1.22 de PyPI,
macos-15-intel, Python 3.12), dos comparaciones fallan allí y en ningún otro
sitio (Linux, macOS arm64, Windows):

```
test_raxopt_sigma2_matches_c[R.4]: raxopt sigma2 0.00009672 vs C 0.00009453
test_cov_raxopt_matches_c (GDP R.2): los errores típicos, rtol 1e-5
```

Hasta atsw 2.1.1 el Mac Intel no tenía rueda binaria, así que corría el porte
contra sí mismo y estas pruebas se saltaban: el fallo no es nuevo, es nuevo
verlo.

En Linux pasan, pero cerca del borde: con el motor C cargado, R.4 da logL
211.214577 en C y 211.214804 en el porte (el C para antes), y los errores
típicos de R.2 difieren en 3·10⁻⁶ relativos (la tolerancia es 10⁻⁵).

R.4 (IPC trimestral, Coint) tiene λ = 1 y refactor = 1: σ² ≈ 10⁻⁴, una
dispersión de ≈ 0.01, por debajo de la banda [0.1, 1000] de `fue.scale`. Es
donde el punto de parada pierde precisión.

## Impacto

Bajo. Los dos motores llegan a la misma verosimilitud salvo decimales; lo que
cambia es dónde certifica el optimizador, y sólo en un fichero fuera de la
banda. Pero son motores que se dicen homologados, y en una plataforma no lo
están dentro de su propia tolerancia.

## Reproducción

`bugs/BUG-0029-repro/repro.py`, con la rueda binaria de fue instalada (si el C
no carga, no juzga y sale 2). Sale 1 si los motores discrepan; 0 si coinciden
en macOS x86_64; 2 si coinciden en otra plataforma. Aquí (Linux x86_64) sale 2:
hay que correrlo en un Mac Intel.

## Por estudiar

Si es la libm o el FMA del x86_64 de Apple; si basta reescalar R.4 a la banda
(refactor 100) para que coincidan, que confirmaría la causa; o si las dos
pruebas deben llevar la tolerancia de fuera de Linux.
