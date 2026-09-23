---
id: BUG-0022
title: En una serie ANUAL con nombre numerico, load() toma el nombre por el ano: la serie se lee empezando 254 anios mas tarde y pierde el nombre
status: fixed
severity: high
component: inp-reader
found_in: fue-python
fixed_in: 0.1.16
reported: 2026-09-16
reporter: David (conformidad, acuerdo entre lectores)
tags:
  - inp
  - anual
  - contrato
references: []
---

## Summary

La línea de fecha de una serie anual tiene tres formas históricas, y el lector
de Python las desambigua con una heurística: **el año es el último token
numérico antes del nombre** (`src/fue/inp.py:257-278`, con el comentario de
BUG-0018 que las enumera).

Si el nombre de la serie **es** un número, la heurística se lo come:

```
 258  1  1766 2020
```

El C lee `nobs=258, año=1766, nombre="2020"` — posicionalmente, que es como lo
lee el motor. Python lee `año=2020` y se queda **sin nombre**:

    C:      nobs=258 freq=1 year=1766 period=1 name=2020
    Python: nobs=258 freq=1 year=2020 period=1 name=series

254 años de desplazamiento, y el nombre sustituido por el genérico `series`.

## Impact

Alto, y de una clase que ningún otro control ve.

El año de inicio es el origen contra el que se interpretan **todas** las
fechas de las intervenciones. Con la serie corrida 254 años, cada `step`,
`impulse` o `ramp` fechada cae fuera de muestra — y una intervención fuera de
muestra produce un regresor idénticamente nulo con coeficiente **libre**, sin
que nada lo diga (regla 22 de `atws/REGLAS-NO-ESCRITAS.md`).

Y el nombre no es decoración: alimenta la inferencia de dominio de art.

Lo peor es que **es invisible para la batería de conformidad**. El fichero lo
aceptan los dos lectores, y el escritor en C lo reescribe byte a byte igual,
así que el juez —que es este mismo lector— compara su propia lectura
equivocada consigo misma y dice «iguales». Hizo falta una herramienta
distinta, `atws/conformidad/acuerdo.sh`, que compara **lo que cada lector
entiende** en vez de lo que cada escritor escribe.

Un nombre numérico no es rebuscado: las series de un panel se llaman por
código, y un año como nombre de columna sale de cualquier hoja de cálculo
pivotada.

## Reproduction

`hueco_nombre_num.inp` del corpus de conformidad es `en4_ar18.inp` con el
nombre cambiado a `2020`:

```
** Number of observations and starting date of time series:
 258  1  1766 2020
```

```python
import fue
ts, m = fue.load("hueco_nombre_num.inp")
print(ts.name, ts.start)      # series (2020, 1)

ts, m = fue.load("en4_ar18.inp")
print(ts.name, ts.start)      # EN.1   (1766, 1)
```

Y el lector en C del GUI sobre el mismo fichero:

```sh
$ gtkfue_roundtrip --cabecera hueco_nombre_num.inp
nobs=258 freq=1 year=1766 period=1 name=2020
```

## Root cause

`src/fue/inp.py:257-278`. Para `freq == 1` el lector no puede confiar en la
posición porque el segundo campo —el «outyear»— aparece en unos ficheros y no
en otros, así que busca el año como el último número de la línea. Un nombre
numérico es indistinguible de un año para esa regla.

El C no tiene el problema porque **no desambigua**: lee posicionalmente
`nobs`, el segundo campo, el año y el nombre, y si el fichero trae otra cosa
se rompe de otra manera. La tolerancia del puerto es real y valiosa —lee
ficheros de DRVUS y de drvec que el C no— pero aquí le cuesta un dato.

## Fix

La heurística sólo hace falta cuando hay **ambigüedad de cuenta**. Con cuatro
tokens no la hay: son `nobs`, el segundo campo, el año y el nombre, en ese
orden, igual que en el caso estacional. Reservar la búsqueda de `nums[-1]`
para las líneas de **tres** tokens, que son las que de verdad no dicen si el
segundo campo está:

```python
if freq == 1:
    if len(toks) >= 4:
        nobs, _outyear, begyear, name = toks[0], toks[1], toks[2], toks[3]
    else:
        ...      # la heuristica de hoy, para las formas cortas
```

Y avisar cuando la heurística se aplique y el nombre resultante quede vacío:
que un `.inp` se lea sin nombre de serie es siempre síntoma, no normalidad.

## Validation

`atws/conformidad/acuerdo.sh`, que compara la cabecera leída por los dos
lectores. Hoy acusa `hueco_nombre_num.inp` con las dos lecturas al lado; con
el arreglo tiene que salir en el grupo «de acuerdo». La batería de
conformidad, en cambio, no sirve para validar esto — y ésa es la mitad
interesante del informe.

## Resolución (0.1.16)

Confirmado. Con cuatro o más tokens y los dos primeros campos numéricos, la
cabecera anual se lee por posición (nobs, segundo campo, año, nombre). La
heurística del último número queda para las formas cortas y, si deja la serie
sin nombre, avisa con `RuntimeWarning`.

Validación: `tests/test_bug_0017_0022_contrato_ficheros.py` (las cuatro formas
de cabecera). Lo validaría también `atws/conformidad/acuerdo.sh`.
