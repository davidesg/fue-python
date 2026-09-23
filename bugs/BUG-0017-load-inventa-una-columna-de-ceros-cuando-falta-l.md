---
id: BUG-0017
title: load() inventa una columna de ceros cuando falta la de un determinista no estandar, en vez de fallar: convierte un fichero mutilado en un modelo con un regresor nulo
status: fixed
severity: high
component: inp-reader
found_in: fue-python
fixed_in: 0.1.16
reported: 2026-09-16
reporter: David (estudio del contrato de ficheros, bateria de conformidad)
tags:
  - inp
  - deterministas
  - contrato
references: []
---

## Summary

En la sección de datos del `.inp`, la primera columna es la serie y detrás va
una columna por cada determinista NO ESTÁNDAR. `load()` lee la columna k-ésima
así (`src/fue/inp.py:502`):

```python
float(toks[col_idx]) if col_idx < len(toks) else 0.0
```

Si la columna no está, **el valor es cero** y la lectura sigue. Un fichero al
que le falten las columnas de los regresores externos no da ningún error: da
un modelo en el que cada regresor externo es idénticamente nulo, con su ω
declarada y estimable.

Ese modelo estima, converge y diagnostica. El efecto del regresor sale cero
porque el regresor ES cero, y nada en el camino lo dice.

## Impact

Alto por lo silencioso. El lector es el único del formato en Python, así que
todo lo que hay por encima —art, drtran-python, el MCP— hereda la invención.

El motor en C no hace esto: `fue` 1.14 exige exactamente `1 + nstdet` reales
por observación y aborta si faltan («observation NN of MM expected»). El
validador `inpcheck.c:348-357` es igual de estricto. Python es aquí más
permisivo que el motor, y la permisividad no produce un aviso: produce datos.

Tiene además un productor conocido de esos ficheros:
`art/pipeline.py::_write_inp` no escribe las columnas (BUG-0187 del registro
de art). Los dos defectos por separado son un fichero mutilado y un lector
generoso; juntos son un regresor que desaparece sin rastro.

## Reproduction

Partiendo de un `.inp` con un determinista no estándar con datos
—`R.2_5.inp` del corpus de conformidad, 68 observaciones, regresor en la
posición 3— y quitándole la segunda columna:

```sh
awk '/Time series/{f=1} f&&NF>1{print $1; next} {print}' R.2_5.inp > mutilado.inp
```

```python
import fue, numpy as np
ts, m = fue.load("mutilado.inp")          # no lanza
print(np.asarray(m.interventions[3].data)[:4])
# [0. 0. 0. 0.]
```

El mismo fichero en el motor:

```
$ fue mutilado
Error in the input file mutilado.inp, line ...: observation 21 of 68 expected
```

## Root cause

`src/fue/inp.py:492-506`. El `else 0.0` está puesto para tolerar ficheros
antiguos sin columnas; pero la tolerancia no distingue «este fichero no
declara deterministas no estándar» de «este fichero declara tres y no trae
ninguna columna». El primero es legítimo; el segundo es un fichero roto.

## Fix

Que la generosidad sea condicional a que NO se hayan declarado columnas:

```python
n_cols_esperadas = 1 + len(nonstd)
if len(toks) < n_cols_esperadas:
    raise ValueError(
        f"{path}: observación {i+1}: se esperaban {n_cols_esperadas} columnas "
        f"(la serie y {len(nonstd)} determinista(s) no estándar) y hay "
        f"{len(toks)}"
    )
```

Es el mismo mensaje que da el motor, y en el mismo sitio. Un fichero que el
motor no puede leer no debería leerlo el puerto: esa es la promesa del puerto.

## Validation

La batería de conformidad (`atws/conformidad/bateria.sh`) usa `fue.load()`
como juez, así que este bug la ciega: un fichero al que le falten las columnas
«pasa» con ceros a los dos lados. Con el arreglo, la batería tiene que empezar
a acusar los ficheros que hoy da por buenos.

## Resolución (0.1.16)

Confirmado. Si el fichero declara `k` deterministas no estándar, cada
observación tiene que traer `1 + k` columnas; si no, `ValueError` con el número
de observación, las columnas esperadas y las que hay. Sin deterministas no
estándar la tolerancia de siempre no cambia.

Validación: `tests/test_bug_0017_0022_contrato_ficheros.py`.
