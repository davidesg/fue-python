---
id: BUG-0011
title: A DRVUS-era .inp does not load — the bands/refactor section did not exist then, and the format carries no version
status: open
severity: medium
component: inp
found_in: 0.1.9
fixed_in:
reported: 2026-08-12
reporter: al montar el banco de Box-Jenkins de DRVUS 1.2.01
tags:
  - inp
  - backward-compatibility
  - format-version
references:
  - src/fue/inp.py:442 (`[3.7] cbands + refactor`, leído incondicionalmente)
  - drvus-source/1.2.01/drvus/src/Box_y_Jenkings/ (los ficheros que no cargan)
  - docs/FILE_CONTRACT.md §2.7
  - BUG-0010 (la otra barrera con los mismos ficheros: la codificación)
---

## Summary

Un `.inp` de la época DRVUS **no carga**:

```
ValueError: Unexpected end of .inp file
```

La causa es que la sección

```
** ACF/PACF bands (0 Automatic) and reescaling factor:
 0 100.00
```

**no existía en DRVUS**. Sus ficheros van del bloque de factores individuales de
la diferencia anual directamente a la serie. El parser de hoy la lee
incondicionalmente (`inp.py:442`), así que consume la primera observación como si
fuera la línea de bandas, y se queda sin datos al final.

Insertando esas dos líneas el fichero carga y estima sin más cambios, así que el
resto del formato **sí** es compatible: lo que falta es una sección, no una
reinterpretación.

## Impact

Medio, y de la misma clase que BUG-0010 (codificación): **rompe la
compatibilidad con el antepasado del propio formato.** Los dos juntos dejan
inaccesible el banco de Box-Jenkins que viaja con DRVUS 1.2.01 —Series A a E,
con sus `.inp`, sus `.out` de referencia y los modelos que Mauricio ajustó—, que
es exactamente el material de validación externa que la documentación necesita.

Se encontró montando el caso `tests/test_published_benchmark_series_a.py`. Con
las dos líneas insertadas a mano, fue reproduce el `.out` de DRVUS a 2.7e-07.

Y el mensaje de error no ayuda: «Unexpected end of .inp file» no dice qué
sección falta.

## Root cause

`inp.py:442`:

```python
# [3.7] cbands + refactor
self._skip_sep()
rf_toks = self._next_data()
```

Sin comprobar que la sección esté. **Y el formato no lleva número de versión**,
así que el parser no puede saber con qué generación de fichero está tratando —
que es la recomendación D5 de la revisión externa, aquí con un caso concreto.

## Fix

Dos piezas:

1. **Tolerar su ausencia.** La sección se detecta por clave, como ya se hace con
   la de `fuf` y con las cuatro secciones opcionales comentadas de
   `fue-1.13.1`. Si la siguiente es la de la serie, `refactor = 1.0` — que es
   además lo que el parser ya hace cuando el campo vale cero.
2. **Un mensaje útil.** «Unexpected end of .inp file» debería decir qué sección
   se esperaba y en qué punto.

Y como mejora aparte, la que la revisión externa pedía: **una línea de versión de
formato** en la cabecera, para que la generación del fichero sea explícita en
vez de inferida.

## Validation

Cargar los `.inp` de `Box_y_Jenkings/Series{A..E}` sin tocarlos y comprobar que
`nobs`, la especificación y las estimaciones coinciden con los `.out` que los
acompañan. Series A ya está medida: θ̂ = 0.699384 contra 0.699384 de DRVUS.
