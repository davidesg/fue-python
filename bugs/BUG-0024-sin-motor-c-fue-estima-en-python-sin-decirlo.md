---
id: BUG-0024
title: Si la extensión C no carga, fue estima con el motor en Python sin decirlo — mismo modelo, otro optimizador, otra velocidad, y nada en la salida lo delata
status: open
severity: medium
component: engine
found_in: 0.1.16
fixed_in:
reported: 2026-09-26
reporter: David — prueba de instalación en frío de atsw 1.5.0
tags:
  - engine
  - fallback
  - silencioso
references:
  - src/fue/_engine.py:67-71
  - BUG-0013 (el respaldo DELIBERADO para modelos sin ARMA)
  - BUG-0005 (los óptimos que dependen del camino del optimizador)
---

## Summary

`_engine.estimate` intenta `from fue._fue_engine import ffi, lib` y, si falla con
`ImportError`, estima con `cast_us.estimate_py` —el puerto en Python— **sin
ningún aviso**:

```python
try:
    from fue._fue_engine import ffi, lib
except ImportError:
    from .cast_us import estimate_py
    return estimate_py(model)
```

## Impact

Medio. No es un número mal calculado: el motor en Python es un puerto
homologado. Pero no es el mismo camino: otro optimizador, otra velocidad (el
motor C es el que se usa en todas las pruebas de rendimiento), y en modelos con
verosimilitud poco nítida puede llegar a otro punto (BUG-0005). Quien tenga la
extensión rota —una rueda mal empaquetada, una GSL ausente, una plataforma sin
rueda que instale la rueda sin C— trabaja en otro modo sin saberlo, y un
resultado que no reproduce en otra máquina no tendría explicación visible.

Visto en la prueba en frío de atsw 1.5.0: el motor C cargó en todas las
versiones y plataformas probadas; el defecto es que, si no cargara, nada lo
diría.

## Root cause

El respaldo tiene dos usos y uno solo se anuncia: el de los modelos sin ARMA
(BUG-0013) es deliberado y se explica en el docstring; el del motor ausente se
trató igual, y en silencio.

## Fix (propuesto)

Un `RuntimeWarning` —una vez por proceso— cuando el motor C no carga: qué error
dio la importación y que se estima en Python. Y exponer el modo
(`fue.engine_backend()` → `"c"` o `"python"`) para que art pueda sellarlo en el
guion junto a la versión del instrumento. La rueda sin C (`FUE_SKIP_C=1`) es un
caso legítimo: el aviso lo dice, no lo impide.

## Validation

Simular la importación fallida (monkeypatch de `fue._fue_engine`) y exigir el
aviso con la causa; y que `engine_backend()` diga `"python"`.
