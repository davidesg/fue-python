---
id: BUG-0016
title: fue.__version__ estaba escrito a mano y se quedó TRES versiones atrás — la 0.1.14 publicada declara 0.1.11
status: fixed
severity: low
component: packaging
found_in: 0.1.11
fixed_in: 0.1.16
reported: 2026-09-07
reporter: David — al verificar la 0.1.14 recién publicada en PyPI
tags:
  - packaging
  - version
  - duplicacion
references:
  - src/fue/__init__.py (_version)
  - tools/gen_api_reference.py (_version — el mismo problema, resuelto sólo ahí)
  - bugs/BUG-0016-repro/repro.py
---

## Summary

`src/fue/__init__.py` tenía `__version__ = "0.1.11"` **escrito a mano**. Al
verificar la 0.1.14 recién publicada, instalándola de PyPI en un entorno limpio:

    pip show fue                  Version: 0.1.14   ✓
    importlib.metadata            0.1.14            ✓
    fue.__version__               0.1.11            ✗

Tres versiones de retraso, publicadas.

## Por qué nadie lo vio

Había una prueba —`test_smoke`— y sólo exigía que `__version__` fuese **una
cadena no vacía**. `"0.1.11"` lo es.

Y el mismo defecto ya había mordido por otro lado: *«la página de API publicada
de la 0.1.10 decía "fue 0.1.9"»*. Aquello se arregló **en el generador**
(`tools/gen_api_reference.py` lee el `pyproject`) sin tocar la raíz — de modo que
el dato siguió estando escrito dos veces, y la copia rezagada siguió ahí
esperando a que alguien la leyera.

Es la duplicación de concepto de siempre: **un dato escrito dos veces y la
segunda copia quedándose atrás**.

## Impact

Bajo, y conviene decirlo con precisión en vez de inflarlo: **`art` NO sella la
versión de `fue`** en sus guiones —comprobado—, así que ningún registro
científico lleva la cifra falsa. Lo que queda mal es lo que un usuario ve al
preguntar por la versión desde Python, y cualquier herramienta de terceros que
se fíe de `__version__` en vez de la metadata.

## Fix

`__version__` se deriva, con el mismo orden que usa `art.version_instrumento` y
por la misma razón:

  1. el `pyproject.toml` del árbol, **si existe** — en una instalación editable
     la metadata se escribió al instalar y no se regenera al subir la versión,
     así que sería ella la que mentiría;
  2. la metadata del paquete instalado, correcta para una rueda;
  3. `"0+desconocida"` — honesto: dice que no consta, no un número.

Nunca levanta: corre en cada `import fue`.

## Validation

Dos pruebas, y hacen falta las dos:

  - `test_dunder_version_matches_pyproject` — que coincidan;
  - `test_dunder_version_is_not_a_literal` — que se **derive**. Sin ésta, la
    primera pasaría también con el número correcto escrito a mano, y volvería a
    quedarse atrás en la siguiente subida.

## Nota — 2026-09-25

Arreglado en el árbol de la 0.1.15, que no llegó a publicarse; la primera versión
publicada que lo lleva es la **0.1.16**, y es la que figura en `fixed_in`.
