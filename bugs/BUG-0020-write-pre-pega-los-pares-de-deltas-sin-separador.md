---
id: BUG-0020
title: write_pre pega los pares de deltas sin separador cuando hay dos o mas: el .pre que escribe no lo relee ni fue.load ni el motor
status: fixed
severity: high
component: inp-writer
found_in: fue-python
fixed_in: 0.1.16
reported: 2026-09-16
reporter: David (bateria de conformidad, actor pypre)
tags:
  - pre
  - contrato
  - deterministas
references: []
---

## Summary

En la sección de δ de los deterministas, `write_pre` emite cada par
`valor bandera` sin separador, y el salto de línea sólo al terminar el bucle
(`src/fue/report.py:1164-1170`):

```python
for val, free in zip(f['delta_vals'][i], itv.delta_free):
    flag = "1" if free else "0"
    out += f"{val:.4f}  {flag}"
out += "\n"
```

Con un solo δ no se nota. Con dos o más sale:

    1.8400  1-0.8631  1

El fichero deja de ser legible. `fue.load()` muere con
`ValueError: invalid literal for int() with base 10: '1-0.8631'`, y el motor
en C lo rechaza con salida 2: *«the flag (0 fixed, 1 estimated) of a delta of
a deterministic variable must be 0 or 1»*.

## Impact

Alto: **el `.pre` deja de ser reejecutable**, que es el invariante entero del
convenio de ficheros. Un modelo con un determinista de δ(B) de orden ≥ 2 —una
respuesta de segundo orden, que es un caso corriente en intervención— no se
puede continuar.

Es copia fiel de un defecto del motor: `fue-1.14/src/fue.c:1975-1983` tenía
exactamente la misma forma, con el `\n` fuera del bucle. **Ya está arreglado
ahí** (commit `0d09abb`), moviendo el salto dentro del bucle como hace la
sección de ω doce líneas más arriba. Este informe es el mismo arreglo en el
puerto.

Ha sobrevivido tanto porque el caso no estaba ejercitado en ningún banco: de
los 115 ficheros de `fue-1.14/tests/` y los 108 del corpus de conformidad,
**ninguno** tenía un determinista con dos δ.

## Reproduction

```python
import fue
from fue.report import write_pre

ts, m = fue.load("RIPC.3.delta2.inp")   # corpus de conformidad
m.fit()
write_pre(m, "x.pre")
```

```sh
$ grep -A2 '^0 0 0 0 0 0 0 0 0 0 0 2' x.pre
0 0 0 0 0 0 0 0 0 0 0 2
**
1.8400  1-0.8631  1
```

```python
>>> fue.load("x.pre")
ValueError: invalid literal for int() with base 10: '1-0.8631'
```

Y el motor sobre ese mismo fichero:

```
Error in the input file x.inp, line 54: the flag (0 fixed, 1 estimated) of a
delta of a deterministic variable must be 0 or 1
```

## Root cause

`src/fue/report.py:1169`: el `out += f"{val:.4f}  {flag}"` no termina en
espacio ni en salto, y el `\n` de la línea 1170 está fuera del bucle. La
sección de ω del mismo módulo (`report.py:1155-1159`) lo hace bien.

## Fix

Un par por línea, como en ω y como ya hace el motor:

```python
for val, free in zip(f['delta_vals'][i], itv.delta_free):
    flag = "1" if free else "0"
    out += f"{val:.4f}  {flag}\n"
```

(y quitar el `out += "\n"` de después del bucle).

## Validation

`atws/conformidad/corpus/RIPC.3.delta2.inp` está puesto en el corpus para
esto. Con el actor `pypre` de la batería (`sh bateria.sh pypre`), que estima
cada fichero y escribe su `.pre` con `write_pre`, ese fichero sale hoy como
fallo con el ValueError; con el arreglo tiene que pasar.

## Resolución (0.1.16)

Confirmado y arreglado como en el motor: un par `valor bandera` por línea.
**`write_fuf` era copia de `write_pre` y tenía el mismo defecto**: arreglado en
los dos.

Validación: `tests/test_bug_0017_0022_contrato_ficheros.py`; batería `pypre`:
`RIPC.3.delta2.inp` pasa.
