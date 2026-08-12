---
id: BUG-0012
title: On Box-Jenkins Series A the ARMA(1,1) stops on the AR boundary, 6.86 in log-likelihood below what Mauricio's own C reaches from the same start
status: open
severity: medium
component: estimation
found_in: 0.1.9
fixed_in:
reported: 2026-08-12
reporter: al montar el banco de Box-Jenkins de DRVUS 1.2.01
tags:
  - optimizer
  - boundary
  - benchmark
references:
  - drvus-source/1.2.01/drvus/src/Box_y_Jenkings/SeriesA/a1.inp (el caso)
  - drvus-source/1.2.01/drvus/src/Box_y_Jenkings/SeriesA/a1.out (el recorrido de DRVUS)
  - tests/test_box_jenkins_series.py
  - BUG-0005 (óptimo espurio según el build; misma familia)
---

## Summary

Sobre la **Serie A de Box-Jenkins**, ARMA(1,1) en el nivel (`a1.inp`, n=197,
semillas φ=0.87 θ=0.48), fue se detiene en

    φ = 0.999978   θ = 0.699565   logL = -57.603862

y el C original de Mauricio, **desde las mismas semillas y sobre los mismos
datos**, alcanza

    φ = 0.908683   θ = 0.575839   logL = -50.745092

**6.86 de verosimilitud mejor**, y además es el valor publicado: Box y Jenkins
dan φ≈0.92, θ≈0.58 para esta especificación.

Las otras ocho comparaciones del banco coinciden a **1e-11 o mejor**, así que no
es un problema general del motor: es este caso.

## Impact

Medio. El punto donde fue se detiene no es un óptimo local cualquiera: es la
**frontera del AR**, φ→1, donde la ARMA(1,1) degenera en la IMA(1,1) —que es
justamente `a2.inp`, y ahí fue clava el resultado de DRVUS a 1.5e-11—. Así que
el modelo que devuelve es interpretable, pero es el equivocado, y lo devuelve sin
avisar.

Importa más de lo que su severidad sugiere porque **la frontera es donde esta
suite trabaja**: el MEG, el DCD y los contrastes de no invertibilidad viven ahí.

## Reproduction

```python
import fue
ts, m = fue.load(".../Box_y_Jenkings/SeriesA/a1.inp")   # requiere BUG-0011 arreglado
m.fit()
m.ar[0][0], m.ma[0][0], m.loglik      # 0.999978, 0.699565, -57.603862
```

El `.out` de DRVUS conserva el recorrido y es lo más informativo del caso:

```
iter  0   φ=0.870000  θ=0.480000  |g|=14.684
iter 10   φ=0.999550  θ=0.526354  |g|= 0.666
iter 30   φ=0.999979  θ=0.700636  |g|= 0.047     <- donde fue se queda
iter 50   φ=0.995054  θ=0.715333  |g|= 0.006
iter 64   φ=0.908683  θ=0.575839  |g|= 0.000     <- DRVUS escapa
```

DRVUS **pasa por el punto de fue** hacia la iteración 30 y sale de él.

## Lo descartado

- **No es el optimizador.** `drvus/src/qnewtop.c` y `csrc/internal/qnewtopt.c`
  difieren en 24 líneas y ninguna es numérica: cabecera GPL, `#include`, y
  `printf` → `if (outputv) fprintf(outputv, …)`.
- **No es el núcleo de la verosimilitud.** `elfvarma.c` es verbatim salvo la
  sustitución de autovalores por GSL (`docs/PROVENANCE.md` §2).

- **No es el porte a Python.** Los DOS motores de fue se detienen en el mismo
  sitio:

  | motor | φ | θ | logL |
  |---|---|---|---|
  | C compilado (`_fue_engine.abi3.so`, el que usa `.fit()`) | 0.999978 | 0.699565 | **−57.603862** |
  | `raxopt` en Python puro (`cast_us.estimate_py`) | — | — | **−57.649946** |
  | **DRVUS, el C original** | 0.908683 | 0.575839 | **−50.745092** |

  Las dos implementaciones de fue coinciden entre sí a 0.05 y difieren de DRVUS
  en ~6.9. Que el porte en Python —una reimplementación independiente del mismo
  algoritmo— caiga donde cae el C de fue dice que el comportamiento viene de lo
  que fue hace, no de un fallo de una de las dos rutas.

Quedan dos candidatos, y no se han medido:

1. **El presupuesto de iteraciones.** DRVUS necesitó 64 y las suyas se ven en el
   `.out`; fue no expone control de iteraciones — `fit()` no acepta argumentos —
   así que la hipótesis no se pudo contrastar desde la API. **Es lo primero que
   hay que mirar**, porque el recorrido de DRVUS muestra que el punto de fue es
   por donde él pasa hacia la iteración 30.
2. **`nlatools.c`, que sí se reescribió** (764 líneas contra 1355, alrededor de
   GSL). Es el único módulo del núcleo que no es de Mauricio, y el álgebra lineal
   cerca de la frontera es exactamente donde una reescritura puede notarse. Pero
   el porte en Python **no usa** `nlatools`, y cae en el mismo sitio — lo que
   debilita esta hipótesis y refuerza la primera.

## Fix

Primero medir cuál de los dos. Lo más barato es **exponer el límite de
iteraciones en `fit()`** —que hace falta de todas formas— y ver si con 100 fue
escapa como DRVUS. Si escapa, el arreglo es el presupuesto y el caso queda como
prueba de regresión. Si no, el sospechoso es `nlatools`.

Y con independencia de eso: **fue no debería devolver un óptimo de frontera sin
decirlo.** Es la misma carencia que BUG-0005 —`converged=True, ifault=0` sobre un
ajuste absurdo— y el principio que la arquitectura de la suite ya declara: el
motor lleva sus reservas como dato.

## Validation

`tests/test_box_jenkins_series.py` lleva el caso marcado como divergencia
conocida y las otras ocho como igualdad a 1e-9. Si `a1` empieza a coincidir, el
test lo dirá.
