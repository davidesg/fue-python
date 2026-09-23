---
id: BUG-0019
title: El puerto rechaza un .pre que escribe el propio motor: phi2 = -0.0 falla el test 'coef must be negative', que en C es '> 0.0' y pasa
status: fixed
severity: high
component: inp-reader
found_in: fue-python
fixed_in: 0.1.16
reported: 2026-09-16
reporter: David (bateria de conformidad, actor motor)
tags:
  - inp
  - contrato
  - frontera
references: []
---

## Summary

Un operador de frecuencia fija guarda un solo coeficiente, φ₂, y tiene que ser
negativo: φ₁ se deriva como `2·cos(2πk/s)·√(−φ₂)`. Los dos lados comprueban la
condición, y **no comprueban lo mismo**:

| | prueba | qué hace con −0.0 |
|---|---|---|
| motor (C) | `if (Ar1f[k][2] > 0.0) ifault = 1` (`fue-1.14/src/fue.c:2924`) | **pasa** (−0.0 > 0.0 es falso) |
| puerto (Python) | `if coef >= 0: raise ValueError("coef must be negative (phi2 < 0)")` (`fue/model.py:30`) | **rechaza** (−0.0 < 0 es False en Python) |

Y ese −0.0 no es hipotético: lo produce el propio motor. El `.pre` escribe los
coeficientes con `%.4f`, así que cualquier φ₂ menor que 5·10⁻⁵ en valor
absoluto sale como **`-0.0000`**.

Resultado: **`fue.load()` no puede leer un `.pre` que `fue` acaba de escribir.**

## Impact

Alto. Rompe el escalón de la escalera: un modelo estimado con el motor en C no
se puede continuar en Python, que es justo lo que el convenio `.pre` existe
para permitir. Y el mensaje no orienta — dice «coef must be negative» sobre un
número que en el fichero se lee `-0.0000` y que efectivamente lleva el signo
menos delante.

El fichero, además, ya venía dañado de antes por otra razón (ver más abajo), y
esta excepción es lo único que lo delata. Al arreglar la frontera conviene no
perder esa señal.

## Reproduction

`syn_ARF.inp` del corpus de conformidad declara un AR(2) de frecuencia fija con
φ₂ = −0.5 de semilla:

```
** Number and frequencies of regular AR(2) operators with fixed frequency:
1 2
**
-0.500000 1
```

```sh
$ cp syn_ARF.inp x.inp && fue x          # sale con 0, converge en 13 iteraciones
$ grep -A3 'AR(2) operators with fixed' x.pre
** Number and frequencies of regular AR(2) operators with fixed frequency:
1 2
**
-0.0000 1
```

El `.out` dice el valor de verdad:

```
  phi[ 2]   =   -0.0000056161
```

Y entonces:

```python
>>> import fue
>>> fue.load("x.pre")
ValueError: coef must be negative (phi2 < 0)
```

mientras que el motor relee ese mismo `.pre` sin protestar (salida 0).

## Root cause

Dos cosas, y conviene separarlas:

1. **La frontera** (este informe). `model.py:30` usa `>= 0` para rechazar, que
   excluye −0.0; el C usa `> 0.0`, que lo admite. Es la diferencia entre
   «negativo estricto» y «no positivo», y en el único punto donde el signo del
   cero importa.

2. **La cuantización del `.pre`** (del lado del C, no de este registro):
   `fue.c` escribe los coeficientes con `%.4f`. Un φ₂ de −5.6·10⁻⁶ vuelve como
   cero, y con él φ₁ = 2·cos(2πk/s)·√(−φ₂) = 0: **el operador entero
   desaparece**. O sea que ese `.pre` no es el óptimo que dice ser ni siquiera
   para el motor — reejecutarlo arranca de otro modelo. Es la grieta de
   precisión del invariante, documentada en `atws/CONTRATO.md` §1.2, con aquí
   su consecuencia estructural.

## Fix

Del lado del puerto, aceptar lo mismo que el motor:

```python
# fue/model.py, en el constructor del factor de frecuencia fija
if coef > 0.0:
    raise ValueError("coef must not be positive (phi2 <= 0)")
```

Y, ya que el caso tiene nombre, avisar en vez de callar cuando llegue un cero:

```python
if coef == 0.0:
    warnings.warn(
        "phi2 = 0 en un operador de frecuencia fija: phi1 = 2·cos(2πk/s)·√(−phi2) "
        "sale 0 y el factor desaparece. Si viene de un .pre, es la cuantización "
        "a %.4f del escritor: el valor real está en el .out.",
        RuntimeWarning)
```

Así el fichero sigue siendo legible —que es lo que pide el contrato— y la
pérdida deja de ser silenciosa, que es lo que importa.

## Validation

`atws/conformidad/corpus/syn_ARF.inp`, con el actor `motor` de la batería
(`sh bateria.sh motor`), que estima cada fichero con `fue` y relee su `.pre`.
Hoy ese fichero sale como fallo con el ValueError; con el arreglo tiene que
pasar, y el aviso tiene que aparecer.

## Resolución (0.1.16)

Confirmado. `FixedFreqFactor` rechaza `coef > 0.0` (la frontera del motor) y
avisa con `RuntimeWarning` cuando `coef == 0.0` (factor que desaparece). El
`.pre` del motor sobre `syn_ARF.inp` se lee y avisa. La batería `motor` ya no
acusa `syn_ARF.inp`.

La cuantización a `%.4f` del lado del C (segundo punto de la causa) sigue
pendiente en fue-1.14.
