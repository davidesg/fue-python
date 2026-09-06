---
id: BUG-0013
title: The C engine segfaults on a model with deterministic inputs and NO ARMA factors — no exception, no message, the interpreter dies
status: fixed
severity: critical
component: engine
found_in: 0.1.11
fixed_in: 0.1.12
reported: 2026-08-14
reporter: al construir el repro de art/BUG-0019
tags:
  - crash
  - binding
  - cast
references:
  - bugs/BUG-0013-repro/repro.py
  - src/fue/_engine.py (la rama que llama a lib.fue_estimate)
  - csrc/fue_api.c (cast_us / count_npar_build_par)
  - BUG-0008 (misma familia: el C muere en vez de informar)
  - art/bugs/BUG-0019 (el camino que construye exactamente este modelo)
---

## Resumen

Un modelo con **regresores deterministas y sin ningún factor ARMA** —
`ar=[]`, `ma=[]`, es decir una regresión con errores ARMA(0,0) — **mata el
intérprete** con SIGSEGV al llamar a `fit()`. No hay excepción, ni mensaje, ni
traza: el proceso desaparece.

Es una especificación legítima y no rebuscada: **es el primer peldaño de
cualquier escalera estacional** — armónicos deterministas, ruido blanco — y es
lo que construye el nodo guiado de `art` cuando la estacionalidad parece
determinista.

## Reproducción

```bash
python bugs/BUG-0013-repro/repro.py --safe   # los dos caminos que sí funcionan
python bugs/BUG-0013-repro/repro.py          # el tercero: SIGSEGV
```

Tres formas del **mismo modelo** sobre la misma serie:

| especificación | motor | resultado |
|---|---|---|
| `ar=[[0.0]]`, `ar_free=[[False]]` — un factor clavado en cero | C | OK, logL = −311.9352 |
| `ar=[]` | Python puro | OK, logL = **−311.9352** |
| `ar=[]` | **C** | **SIGSEGV** |

Los dos que funcionan coinciden hasta el último dígito, así que no hay duda de
que es el mismo modelo: lo que decide no son las matemáticas, es **si la
especificación lleva algún factor ARMA declarado**.

Medido también con `d=0` y `d=1`, con μ estimada y sin ella, y con uno o con
once armónicos: **las cuatro combinaciones revientan**. Lo único que lo evita es
que exista un factor.

## Impacto

**Crítico**, por tres razones y no por una:

1. **Mata el proceso.** Un crash no es un resultado equivocado que se pueda
   revisar: en un servidor MCP se lleva por delante la sesión, y en un guion por
   lotes no deja ni registro de dónde iba. Es la carencia de `BUG-0008` otra
   vez: *el C muere en vez de informar*.
2. **Está en el camino principal.** No hay que buscarlo: es el modelo de
   estacionalidad puramente determinista, el que la metodología propone como
   punto de partida.
3. **No lo tapa la batería.** Los casos reales del banco vienen de ficheros
   `.inp`, y **un `.inp` siempre declara la sección AR** — `art` y el propio
   escritor de `fue` emiten `1 1` con el coeficiente a cero y bandera de fijo.
   Por eso el `m00` del repro de `art/BUG-0019` estima sin problema: lleva
   `ar=[[0.0]]` con `ar_free=[[False]]`. **El crash sólo se alcanza construyendo
   el modelo por la API de Python**, que es lo que hace cualquiera que use el
   paquete como biblioteca.

## Dónde mirar

El motor Python puro trata el caso sin quejarse, así que el defecto está en el
puente o en el `cast`:

* `src/fue/_engine.py` — la rama que llama a `lib.fue_estimate`. Ya trata un
  caso límite parecido (`npar == 0` va por `eval_at_params` «porque el backend
  de C revienta en ese caso»), lo que sugiere que **este no es el primer sitio
  donde el C no admite una estructura vacía**.
* `csrc/fue_api.c` — `cast_us()` y el conteo de parámetros: con `p = q = 0` hay
  vectores que se dimensionan a cero y se escriben igual.

## Arreglado — 14-ago-2026 (desviación)

`_engine.estimate` desvía al motor Python cuando la especificación no declara
**ningún** factor ARMA (`_sin_estructura_arma`, que mira `ar`, `ma`, `ar_s`,
`ma_s`, `ar_f`, `ma_f`). Es la misma forma que ya tenía el caso `npar == 0`, y
por la misma razón: hay una respuesta correcta disponible y desviar es mejor
que lanzar — quien llama obtiene su ajuste y no se entera del agujero.

Un factor **clavado en cero sigue siendo un factor** y no se desvía: eso es lo
que escriben los ficheros, y va por el C como siempre.

Verificado en `tests/test_bug_0013_no_arma_factors.py` (6 pruebas): las cuatro
combinaciones que reventaban —d ∈ {0,1} × μ ∈ {sí,no}— estiman; y las dos
escrituras del mismo modelo coinciden a **1.9e-07**, que es el acuerdo entre
motores que documenta `docs/PERFORMANCE.md`, no epsilon de máquina.

⚠ **El defecto del C sigue ahí.** Esto es un rodeo, no una reparación: la
escritura fuera de rango con `p = q = 0` no se ha buscado. Lo que se ha quitado
es la muerte del proceso.

## Arreglo, el de fondo

Dos niveles, y conviene el primero aunque se haga el segundo:

1. **Que no muera.** El puente debe rechazar o desviar la especificación que el
   C no soporta, igual que hace con `npar == 0`. Desviar es preferible a
   rechazar: el motor Python puro **ya da la respuesta correcta**, así que
   `_engine.estimate` puede caer a `estimate_py` cuando no hay factores ARMA, y
   el usuario obtiene su ajuste sin enterarse.
2. **Que el C lo soporte.** Encontrar la escritura fuera de rango con `p=q=0`.
   Es lo correcto, pero es trabajo sobre el C y no debe bloquear lo anterior:
   entre un crash y un ajuste correcto por otra ruta, la elección no es difícil.

⚠ Mientras tanto, el rodeo que ya usan los ficheros: **declarar un factor AR(1)
fijado en cero**. Es el mismo modelo y no cambia ni un dígito.

---

## TODO — lo que queda abierto: el binario C y las wheels (2026-09-06)

El rodeo del puente resolvió **la API de Python**. No resolvió las otras dos
puertas, y una premisa escrita aquí resultó ser falsa.

`_engine.estimate` dice, justificando por qué bastaba con desviar:

> *«a file always declares the AR section, with the coefficient at zero and the
> flag fixed. **Only the Python API reaches here.**»*

**Medido sobre el ecosistema: 178 de 4.505 `.inp` declaran CERO factores ARMA de
todos los tipos.** No es un caso de laboratorio, y no son ficheros ajenos:

    ~/Dropbox/TFM_UCM/Tesis_Michael/replica/ITCER.inp        segfault
    ~/Dropbox/TFM_UCM/Tesis_Michael/replica/PGAS.inp
    ~/Dropbox/TFM_UCM/Tesis_Michael/replica/RATIO.inp
    ~/…/autonomo2/ITCER/ITCER_m10.inp                        segfault
    ~/Dropbox/Article_Multivariate Convergence/…/CINY.inp
    ~/Dropbox/Taller de ST Master/…/D.inp

Las tres series del TFM en curso y los modelos de una de ellas. Comprobado
contra `/usr/local/bin/fue` (FUE 1.13, construido de 1.13.1): **SIGSEGV**.

Y los escribe `art`: cuando el modelo no lleva ARMA, `pipeline._write_inp` emite
la sección con un `0` y no el `1 1 / 0.0 0` que el rodeo presupone. La creencia
de que sólo la API de Python llegaba aquí es lo que dejó la puerta del ejecutable
sin vigilar.

### Por qué sigue ahí

Porque el desvío funciona **demasiado bien**. Quien usa `fue` desde Python
obtiene su ajuste y no se entera del agujero; el `.inp` que queda escrito
—sin sección ARMA— sólo mata al proceso cuando alguien lo pasa por el binario o
por una wheel sin el desvío. Las dos rutas que fallan son justamente las que no
tienen tests en este paquete.

### El alcance real, en tres puertas

| puerta | estado |
|---|---|
| API de Python (`_engine.estimate`) | **desviada**, correcta |
| binario C `/usr/local/bin/fue` | **muere**, 178 ficheros expuestos |
| wheels sin el desvío | igual que el binario |

### Lo que hay que decidir

El arreglo de fondo —la escritura fuera de rango con `p=q=0` en el C— sigue sin
buscarse y **no es trivial**. Pero hay una mitigación barata que no depende de
él y cierra el caso para todo fichero futuro: **que `art` escriba siempre la
sección AR con `1 1 / 0.0 0`** —un AR(1) fijado en cero, que es el mismo modelo
y no mueve un dígito— en lugar de `0`. Es el rodeo que los ficheros antiguos ya
usaban sin saber por qué.

Eso deja tres cosas separadas: el C (de fondo, no trivial), el escritor de `art`
(barato, cierra el futuro) y los 178 ficheros ya escritos (recuperables uno a
uno con el mismo cambio de una línea).

### Estado (2026-09-06): la mitigación de `art`, aplicada

`art` **ya no escribe la sección AR vacía**: cuando el modelo no lleva ningún
factor ARMA, `pipeline._write_inp` emite `1 1 / 0.0 0` —un AR(1) fijado en
cero—, que es la convención que el resto del ecosistema ya seguía. Comprobado:

    antes   fue W eml → Segmentation fault (core dumped)
    después fue W eml → exit 0
    ℓ = −869.3408898 en los dos casos, npar = 2

Eso cierra el caso **para todo fichero que se escriba desde ahora**, y sólo eso.
Siguen abiertos, y no dependen de este paquete:

  1. **El C**, que es el arreglo de fondo: la escritura fuera de rango con
     `p=q=0` sigue sin buscarse. No es trivial.
  2. **Los 178 `.inp` ya escritos**, que hay que reescribir uno a uno —o
     regenerar— para que dejen de matar al binario.
  3. **Las wheels** sin el desvío del puente, que se comportan como el binario.
