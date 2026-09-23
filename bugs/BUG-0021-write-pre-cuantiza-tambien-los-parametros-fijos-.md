---
id: BUG-0021
title: write_pre cuantiza tambien los parametros FIJOS: un AR fijado en 0.941176 vuelve fijado en 0.9412, que es otro modelo
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
  - precision
references: []
---

## Summary

`write_pre` escribe todos los coeficientes con un número fijo de decimales:
ARMA y δ con `{:.4f}`, ω y μ con `{:.6f}`, λ con `{:.2f}`, el armónico de un
cos/sin con `{:.0f}`.

Para un parámetro **libre** eso es discutible pero defendible: el valor es una
semilla, y volver a arrancar desde el valor redondeado a cuatro decimales casi
siempre lleva al mismo óptimo. Es la grieta de precisión conocida del
invariante, y está medida.

Para un parámetro **FIJO no lo es en absoluto.** Un parámetro fijo no es una
estimación: es parte de la **especificación**. El analista lo pinchó en un
valor concreto y el modelo es ese. Redondearlo no re-siembra nada — cambia el
modelo.

Medido sobre el corpus de conformidad:

    R.1_4.inp     AR(1) FIJO en 0.941176   ->   0.9412
    R.2_2.inp     AR(1) FIJO en 0.941176   ->   0.9412
    RIPC.3.1.inp  AR(1) FIJO en 0.944444   ->   0.9444

`0.941176` es 16/17 y `0.944444` es 17/18: valores puestos a mano, con
intención. El `.pre` los devuelve truncados, y reejecutar desde ahí ajusta
otro modelo.

Lo mismo con los demás campos que son especificación y no semilla:

    hueco_lambda.inp     lambda 1/3 = 0.333333  ->  0.33
    hueco_mu_fijo.inp    mu FIJO en -88.717979  ->  0.0  (BUG aparte: se tira)
    hueco_armonico.inp   cos 1.5                ->  cos 2.0   (redondea)

## Impact

Alto y silencioso. El fichero sale bien formado y se relee sin error; lo único
que ha cambiado es el modelo.

Afecta a todo uso legítimo de los parámetros fijos, que en esta escuela no son
raros: un AR fijado en un valor teórico, un operador de frecuencia fija
sembrado a mano, un testigo MEG con su coeficiente pinchado, un δ pinchado en
una respuesta conocida.

El motor en C tiene el mismo defecto —`write_pre` es copia fiel suya— así que
esto no es una divergencia del puerto, es un defecto compartido. Pero el puerto
es el que puede arreglarlo sin romper 30 años de ficheros, porque nadie compara
sus bytes con nada.

## Reproduction

```python
import fue
from fue.report import write_pre

ts, m = fue.load("R.1_4.inp")      # corpus de conformidad
print(m.ar, m.ar_free)             # [[0.941176]] [[False]]  <- FIJO

m.fit()
write_pre(m, "x.pre")

ts2, m2 = fue.load("x.pre")
print(m2.ar, m2.ar_free)           # [[0.9412]] [[False]]
```

El fichero de partida, para que se vea que la bandera dice fijo (ausente ⇒ 0):

```
**Number and orders of regular AR operators:
1 1
**
0.941176
```

## Root cause

`src/fue/report.py`, en todos los emisores de coeficientes: `_arma_body`
(línea 1077, `{:.4f}`), la sección de δ (1169, `{:.4f}`), la de ω (1157,
`{:.6f}`), μ (1195, `{:.6f}`), λ (1200, `{:.2f}`) y el armónico en
`_itv_name_line` (1545, `{:.0f}`).

Ninguno mira la bandera. La distinción semilla/especificación no existe en el
escritor.

## Fix

Dos reglas, y la segunda es la que importa:

1. Los valores **libres** pueden seguir con su formato de hoy, o subir a
   `.10f`. Es una decisión de gusto.
2. Los valores **FIJOS**, y todo campo que sea especificación y no estimación
   —λ, el armónico de un cos/sin, la frecuencia de un operador de frecuencia
   fija— se escriben con la representación más corta que **relea idéntica**:

```python
def _exacto(v):
    """El minimo de decimales que vuelve a leer el mismo float."""
    for d in range(6, 18):
        s = f"{v:.{d}f}"
        if float(s) == v:
            return s
    return repr(v)
```

Es lo que hace el GUI en C desde `a0d3509` (`gtk_fue/src/utils.c::inp_format`,
copiado a su vez de `fug/src/inpfile.c`), y por eso el GUI pasa hoy la batería
con 0 fallos.

## Validation

`sh bateria.sh pypre` en `atws/conformidad`. La batería compara los `.pre` en
modo estructura: se salta los valores que el fichero declaraba **libres** y
compara los **fijos**, que es exactamente esta distinción. Hoy acusa
`R.1_4.inp`, `R.2_2.inp`, `RIPC.3.1.inp`, `hueco_lambda.inp` y
`hueco_armonico.inp`; con el arreglo tienen que pasar los cinco.

## Resolución (0.1.16)

Confirmado. Nuevo `report._exacto(v, decimales)`: el formato de siempre si ya
relee idéntico (mismo byte), si no los mínimos decimales que lo consigan. Se
aplica a los valores FIJOS (ARMA, ω, δ, φ₂ de frecuencia fija), y siempre a λ,
al armónico cos/sin y a la frecuencia de un operador de frecuencia fija. Las
semillas libres conservan su formato.

Además, en `write_pre` **y** `write_fuf`:
- la μ FIJA no nula se tiraba (`0`): ahora sale `<μ exacta> 0`;
- los DATOS de la serie salen exactos (antes `.10f`): un `.inp` exacto y su
  `.pre` estimaban sobre datos distintos (art/BUG-0188).

Validación: `tests/test_bug_0017_0022_contrato_ficheros.py`; batería `pypre`:
pasan `R.1_4`, `R.2_2`, `RIPC.3.1`, `hueco_lambda`, `hueco_armonico` y
`hueco_mu_fijo`. Lo que queda en `pypre` (12) es la sección de horizonte de los
ficheros fuf, que `write_pre` no escribe por diseño.
