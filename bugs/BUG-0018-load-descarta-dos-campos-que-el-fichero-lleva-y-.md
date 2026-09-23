---
id: BUG-0018
title: load() descarta dos campos que el fichero lleva y el motor conserva: la frecuencia  (serie sin fechar) y cbands
status: fixed
severity: medium
component: inp-reader
found_in: fue-python
fixed_in: 0.1.16
reported: 2026-09-16
reporter: David (estudio del contrato de ficheros, bateria de conformidad)
tags:
  - inp
  - contrato
references: []
---

## Summary

Dos campos del `.inp` se leen y se tiran:

**La frecuencia `number`.** El formato admite la palabra `number` en el campo
de frecuencia: una serie **sin fechar**, numerada 1, 2, 3… en vez de datada.
`inp.py:243-247` la reconoce:

```python
if freq_str == "number":
    freq, numbering = 1, True
```

`numbering` es una **variable local que no se guarda en ningún sitio**. Fuera
de ese método no existe. El modelo que sale es indistinguible de una serie
ANUAL, así que todo escritor de Python la reescribe como ` 1` y la serie
adquiere un año de inicio que nunca tuvo.

**`cbands`.** El primer token de la línea de bandas ni siquiera se lee
(`inp.py:487`: `refactor = float(rf_toks[1])`). Todo escritor de Python emite
un `0` fijo en su lugar.

El motor en C conserva los dos: reescribe `" number "` (`fue.c:1748-1749`) y
`" %.2f %.2f"` con `cbands` delante (`fue.c:2329-2330`). Es la única de las
seis implementaciones que lo hace — el GUI también los pierde.

## Impact

Medio, y con un agravante de método.

Sobre los datos: una serie sin fechar que pase por Python vuelve como anual
empezando en el año que hubiera en ese campo. Las fechas de las intervenciones
se interpretan entonces contra un origen inventado.

El agravante: **la batería de conformidad usa `fue.load()` como juez**, así
que estos dos campos son un punto ciego del banco entero. Se añadieron al
corpus dos ficheros para ejercitarlos —`hueco_number.inp` y
`hueco_cbands.inp`— y los dos salen «iguales» pase lo que pase con el campo,
porque el juez no lo ve. Un lector que descarta un campo hace invisible la
pérdida de ese campo para todo lo que se apoye en él.

## Reproduction

```python
import fue

ts, m = fue.load("hueco_number.inp")   # su línea de frecuencia dice " number"
print(ts.freq)          # 1     — indistinguible de una serie anual
print(hasattr(ts, "numbering"), hasattr(m, "numbering"))   # False False
```

```python
# la línea de bandas del fichero dice " 2.50 50.00"
ts, m = fue.load("hueco_cbands.inp")
print(m.refactor)       # 50.0  — bien
# cbands (2.50) no está en ningún atributo del modelo ni de la serie
```

Y el ida y vuelta, con cualquier escritor de Python:

```
  entra:  ** Frequency of time series: either 1(A), 4(Q) or 12(M):
           number
  sale:   ** Frequency of time series: either 1(A), 4(Q) or 12(M):
           1
```

## Root cause

`src/fue/inp.py:243-247` (la frecuencia) y `src/fue/inp.py:484-490` (las
bandas). En los dos casos el campo se reconoce y se descarta; no es que el
parser no lo entienda.

## Fix

Guardarlos, para que los escritores puedan devolverlos:

```python
# inp.py, al construir la serie
ts.numbering = numbering          # o freq = "number" como valor sentinela

# inp.py, en la sección de bandas
cbands = float(rf_toks[0]) if rf_toks else 0.0
model.cbands = cbands
```

Y que `report.py::write_pre` y `art/pipeline.py::_write_inp` los emitan, con
el valor por defecto de hoy cuando el atributo no esté — así los ficheros
existentes no cambian.

Conviene decidir además si `numbering` merece ser un valor de `freq` (la
cadena `"number"`) en vez de un booleano aparte: hoy `freq == 1` significa dos
cosas distintas y sólo el fichero sabe cuál.

## Validation

`atws/conformidad/corpus/hueco_number.inp` y `hueco_cbands.inp` están puestos
para esto. Hoy pasan la batería sin ejercitar nada; con el arreglo, tienen que
empezar a acusar a los escritores que pierden los campos.

## Resolución (0.1.16)

Confirmado. Decisión sobre `number`: **no** como valor de `freq`. Ya existía
`TimeSeries.numbering` (propiedad, `freq == 0`); ahora también es asignable, y
el lector la pone a `True` en una serie `number` manteniendo `freq=1`, que es
como la estima el motor (`sper=1`) y lo que todo el código de arriba (art) sabe
manejar. De paso el motor recibe `numbering=1` como en C.

`cbands` se guarda en `model.cbands` (0.0 por defecto en `Model`).
`write_pre`/`write_fuf` y `art._write_inp` emiten ` number` y `cbands`. Los
clones (`add_intervention`, los de art) no los arrastran: vuelven al valor por
defecto, que es el mismo fichero de siempre.

Validación: `tests/test_bug_0017_0022_contrato_ficheros.py`.
