---
id: BUG-0014
title: La previsión de Python daba efecto NULO a compimp, easter y trend — un segundo generador del mismo regresor que se quedó atrás
status: fixed
severity: high
component: forecast
found_in: 0.1.12
fixed_in: 0.1.13
reported: 2026-09-06
reporter: David — «revisa si el bug de easter toca al motor de previsión, fuf»
tags:
  - silent
  - forecast
  - deterministic
  - easter
  - compimp
  - trend
  - duplicacion
references:
  - src/fue/forecast.py (_build_xi)
  - src/fue/cast_us.py (_build_indicator, el generador bueno)
  - bugs/BUG-0014-repro/repro.py
  - BUG-0007 (misma familia: un determinista que se pierde en silencio)
---

## Summary

`forecast._build_xi` construía el indicador de cada determinista **por segunda
vez**, indexado por `type_code` en vez de por nombre, con ramas para 8 de los 11
tipos y **sin `else` final**:

    TYPES     impulse 0 · step 1 · ramp 2 · seasonal 3 · cos 4 · sin 5
              alter 6 · custom 7 · compimp 8 · easter 9 · trend 10
    _build_xi 0 1 2 3 4 5 6 7 ................. ✗ ..... ✗ ...... ✗

`compimp`, `easter` y `trend` caían fuera, `D` quedaba a ceros y su efecto valía
**exactamente cero**. Sin aviso, sin error: el ω se estima, aparece en el `.out`
con su error típico, y la previsión lo ignora.

## El daño es doble

`xi` se usa en dos sitios de `forecast()`:

    nt[i] - xi[i]      limpia la HISTORIA para obtener el ruido
    f1[l] += xi[nobs+l]  añade el efecto al FUTURO

Con `xi ≡ 0` fallan los dos: el ruido que alimenta la recursión queda
contaminado por un determinista que nadie quitó —y eso sesga toda la trayectoria,
no sólo los meses afectados— y la previsión además no lleva el efecto.

## Medición

Serie sintética mensual con un efecto de Semana Santa de **+4%** inyectado; el
motor lo estima como ω = 399,4 (centésimas de log = 3,99%). Previsión a 18
meses, Pascua de 2020 en abril:

    ruta            03/2020    04/2020    05/2020
    fuf (C)         100.183    104.178    100.185     ← correcto
    Python (antes)  100.184    100.184    100.186     ← el efecto no está
    Python (ahora)  100.184    104.178    100.186

Y el síntoma que lo delataba a simple vista: la variación interanual de 04/2020
salía **−411,60%**, comparando una previsión sin Semana Santa contra un abril
observado que sí la tenía. Ahora sale −12,20%.

Comprobado también que `fuf` (el C) **no** tiene este defecto: escribe el
`easter` en el fichero de previsión y extiende el calendario al futuro.

## Root cause

**Dos generadores del mismo regresor.** `cast_us._build_indicator` los construye
por nombre y cubre los once; `forecast._build_xi` los reimplementaba por código y
cubría ocho. El segundo se quedó atrás cuando el primero creció — que es
exactamente lo que pasa con toda duplicación de concepto, y por eso una
duplicación no es nunca inocente.

## Fix

Se borra el duplicado: `_build_xi` llama a `_build_indicator(itv, T, freq, …)`
con `T = nobs + horizonte`. El generador acepta la longitud, así que cada tipo se
extiende **por su propia regla** —el easter por el calendario, el step por su
definición— sin repetir ninguna.

Dos cosas más, para que la clase quede cerrada:

  - `_build_indicator` gana un `else` que **levanta** `ValueError` ante un tipo
    desconocido. Devolver ceros es prever un modelo distinto del pedido sin
    decirlo, y ése era el fallo.
  - `custom` deja de asignar `itv.data[:nobs]` a una ventana que puede ser más
    larga que los datos: rellena lo que hay y el resto queda a cero. Sin esto,
    pedir el indicador hasta `nobs+horizonte` reventaba.

## Validación

`bugs/BUG-0014-repro/repro.py`: ningún tipo da regresor nulo, y un tipo
inventado se rechaza. Y de extremo a extremo, las dos rutas coinciden:
**104,1782 (Python) frente a 104,1779 (C)**, diferencia atribuible a los
decimales con que el fichero de previsión guarda los parámetros.
