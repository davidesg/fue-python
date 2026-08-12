---
id: BUG-0010
title: The .inp parser assumes UTF-8, so files written by the original C program on a Latin-1 system cannot be read at all
status: open
severity: medium
component: inp
found_in: 0.1.9
fixed_in:
reported: 2026-08-12
reporter: David / archivo de la tesis (IPC Chile y Colombia)
tags:
  - inp
  - encoding
  - interop
  - c-parity
references:
  - fue/inp.py:96 (`for line in f` sobre un fichero abierto sin encoding explícito)
  - ~/Documents/Documentos/Tesis/Analisis/{Chile,Colombia}/ipc/mensuales/... (los ficheros)
---

## Summary

`fue.inp.load` abre el `.inp` con la codificación por defecto de la plataforma
—utf-8 en Linux— y revienta con `UnicodeDecodeError` sobre cualquier fichero
escrito en Latin-1:

```
UnicodeDecodeError: 'utf-8' codec can't decode byte 0x82 in position 122
```

Los `.inp` del archivo de la tesis (IPC de Chile y Colombia, 1986-2001) son
justo eso: los escribió el programa C original en un sistema Latin-1, y llevan
acentos en los nombres de las variables deterministas y en los comentarios.

## Impact

Medio, y de una clase concreta: **rompe la interoperabilidad con el propio
programa que definió el formato.** El `.inp` es el contrato central de la suite,
y un fichero legítimo, producido por el C, no lo lee el port en Python. No es un
caso hipotético — bloqueó el uso de las dos series I(2) de la tesis como banco de
pruebas del DCD regular y del MEG, que es donde se encontró.

Afecta a todo el archivo histórico de análisis, que es material valioso: modelos
ya identificados por un analista, con sus `.out` de referencia.

## Reproduction

```python
import fue
fue.inp.load(".../Chile/ipc/mensuales/analisis/muestra_1.86_12.01/PC.inp")
# UnicodeDecodeError
```

Workaround usado en `tests/test_thesis_i2_chile_colombia.py` de art: convertir a
un temporal antes de leer.

```python
open(dst, "w", encoding="utf-8").write(open(src, encoding="latin-1").read())
```

## Root cause

`inp.py:96` itera sobre un fichero abierto sin `encoding=`, así que Python usa
`locale.getpreferredencoding()`. El formato `.inp` es de los años noventa y su
codificación de facto es la de la plataforma que lo escribió.

## Fix

Abrir con `encoding="utf-8", errors="replace"`, o mejor: intentar utf-8 y caer a
latin-1 al fallar. El contenido que importa —números y palabras clave— es ASCII;
lo que trae bytes altos son nombres y comentarios, así que una lectura tolerante
no arriesga nada numérico.

Conviene hacerlo también en el ESCRITOR, declarando utf-8 explícitamente, para
que los ficheros nuevos tengan codificación conocida en vez de la del sistema.

## Validation

Cargar los `.inp` del archivo de la tesis sin conversión previa y comprobar que
`nobs`, `start` y los parámetros salen idénticos a los de la conversión manual.
Los ficheros de Chile y Colombia sirven de fixture y ya están en uso desde art.
