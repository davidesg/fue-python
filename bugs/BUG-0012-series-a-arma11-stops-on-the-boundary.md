---
id: BUG-0012
title: On Box-Jenkins Series A the ARMA(1,1) stops on the AR boundary — and so does Mauricio's own C when compiled today: the 2001 reference was run at 80-bit x87 precision
status: fixed
severity: medium
component: estimation
found_in: 0.1.9
fixed_in: 0.1.10
reported: 2026-08-12
reporter: al montar el banco de Box-Jenkins de DRVUS 1.2.01
tags:
  - optimizer
  - boundary
  - benchmark
references:
  - drvus-source/1.2.01/drvus/src/Box_y_Jenkings/SeriesA/a1.inp (el caso)
  - drvus-source/1.2.01/drvus/src/Box_y_Jenkings/SeriesA/a1.out (el recorrido de DRVUS)
  - tests/test_optimizer_termcode.py
  - tests/test_box_jenkins_series.py
  - BUG-0005 (óptimo espurio según el build; misma familia)
---

## Resuelto — 12-ago-2026: no era fue

**El propio DRVUS, compilado hoy, falla exactamente igual.** Misma fuente de
Mauricio, sin una línea tocada, `gcc -O2` en esta máquina:

| | criterio | iteraciones | ‖g‖ | logelf |
|---|---|---|---|---|
| `a1.out` de ~2001 | GRADIENTE | 64 | 0.0000000 | **−50.7450915148** |
| **DRVUS compilado hoy** | **PASO** | **23** | **0.0100895** | **−57.6038617504** |
| fue (el puerto) | PASO | 23 | 0.0100924 | −57.6038617 |

Ninguna función cambió. El puerto reproduce el C original **incluido su fallo**,
que es lo que se le pide a un puerto.

### El testigo directo: los binarios de época corren hoy

No hizo falta creerse la reconstrucción, porque los binarios se conservan:

| binario | qué es | resultado hoy |
|---|---|---|
| `1.2.01/drvus/src/drvus` | ELF 32-bit, **octubre de 2006** | GRADIENTE, 64 iter, −50.7450915148 — **idéntico al `.out` archivado**, línea por línea en las 322 primeras |
| `1.0/Drvus/Drvus.exe` | PE32 Borland, **mayo de 2001**, bajo `wine` | PASO, 25 iter, ‖g‖=0.0002, −57.6037540996 |

O sea que **DRVUS 1.0 tampoco alcanza el óptimo en este caso**, y su `.exe` de
2001 coincide **bit a bit** con la recompilación de su fuente a 80 bits que se
hizo aquí — lo que valida el método antes de usarlo.

### El bisect entre versiones, con el mismo compilador y las mismas opciones

| versión | `-m32 -O0` | logelf |
|---|---|---|
| 1.0 | PASO, 25 iter | −57.6037540996 |
| **1.01** | **GRADIENTE, 64 iter** | **−50.7450915148** |

Y el `diff` de código entre las dos, con los comentarios fuera:

| módulo | líneas de código que difieren |
|---|---|
| `qnewtopt.c` → `qnewtop.c` (el optimizador) | **0** |
| `usmelard.c` (AS 197) | **0** |
| `elfvarma.c` (AS 311) | **0** |
| `drvmlest.c` (el driver) | **0** |
| `nlatools.c` | 2 — un `round()` renombrado |
| `drvus.c` | 32 — 30 son `strcmpi`→`strcmp` y el número de versión |

Lo que queda, y es **todo** lo funcional que cambió en el camino de estimación:

```c
   macheps = cmacheps();
-  gradtol = 1.0e-7;
-  steptol = 1.0e-5;
+  gradtol = pow( macheps, 1.1 / 3.0 );      /* ~ 1.98e-06 */
+  steptol = pow( macheps, 2.0 / 3.0 );      /* ~ 3.67e-11 */
```

`steptol` se hace **seis órdenes de magnitud más estricto**. Prueba: poniendo
esas dos líneas —y nada más— en la fuente de 1.0, la 1.0 pasa a dar
**GRADIENTE, 64 iteraciones, −50.7450915148**, exacto. fue lleva los valores de
1.01+ (`fue_api.c:70-71`), idénticos.

### La causa, localizada por barrido de compilación

Misma fuente, mismos ficheros, sólo cambian las opciones del compilador:

| opciones | criterio | iter | logelf |
|---|---|---|---|
| `-O0` | PASO | 23 | −57.6038617504 |
| `-O1` | PASO | 23 | −57.6038617504 |
| `-O2` | PASO | 23 | −57.6038617504 |
| `-O2 -ffloat-store` | PASO | 23 | −57.6038617504 |
| `-O2 -mfpmath=387` | PASO | 23 | −57.6038617465 |
| **`-O0 -mfpmath=387`** | **GRADIENTE** | **64** | **−50.7450915148** ← el de 2001, exacto |

Es la **precisión de los intermedios**: la x87 de 32 bits guardaba los cálculos
en registros de **80 bits**; hoy x86-64 usa SSE2 y trabaja en 64. Con `-O0
-mfpmath=387` los intermedios vuelven a ser de 80 bits y el recorrido de hoy
sigue al de 2001 **iteración por iteración** —mismos puntos en la 10, 20, 30,
40, 50 y 60, misma parada en la 64, φ=0.908683, θ=0.575839, μ=17.065276 frente a
17.065277—. Con `-O2 -mfpmath=387` el optimizador vuelca a memoria y redondea a
64, y el resultado vuelve al de SSE2: por eso hace falta `-O0`.

Once bits más de mantisa en los intermedios deciden, en un valle así de plano, si
la búsqueda puede seguir bajando o no.

Y **no es la regla de parada**, que era lo primero que había que descartar: con
`steptol = 0` —el criterio de paso desactivado— el build de 64 bits agota las 500
iteraciones y termina en el mismo sitio, −57.6038617504. En 64 bits la búsqueda
está genuinamente atascada ahí; el criterio de paso sólo lo anuncia.

Los dos ingredientes hacen falta a la vez: las tolerancias de 1.01 **y** los 80
bits. Con las tolerancias viejas ni a 80 bits se llega (1.0 para en la 25); con
las nuevas pero a 64 bits, tampoco (fue, y DRVUS 1.2.01 de hoy, paran en la 23).

### Y el óptimo de 2001 es el bueno

`statsmodels` —otro algoritmo, otros autores, ningún ancestro común— sobre la
misma ARMA(1,1) con media:

    statsmodels     phi=0.908685  theta=0.575841  media=17.0653  logL=-50.745092
    DRVUS 2001      phi=0.908683  theta=0.575839  mu   =17.065277 logL=-50.745092
    fue, mu0=17     phi=0.908685  theta=0.575841  mu   =17.065277 logL=-50.745092

Así que la verosimilitud de fue es la correcta y su optimizador llega: lo que
falla desde μ₀=2.5 es el camino, y falla igual en el C de Mauricio.

### Qué queda de esta ficha

Lo que se arregló —porque hacía falta y porque sin ello nada de lo anterior se
podía ver—:

1. **El motor devuelve su veredicto**: `termcode`, `niter`, `gnorm`.
2. **`converged` ya no significa «no reventó»**, y un alto que no es máximo sale
   por `RuntimeWarning`. Sin esto, fue devolvía −57.60 como bueno.
3. **La semilla de μ, documentada** (`docs/FILE_CONTRACT.md` §2.4): sembrar μ de
   la media de la variable diferenciada, que es lo que hace `art`.

Y una consecuencia para todo el archivo histórico, en §Lo que apareció.

---

## Summary

Sobre la **Serie A de Box-Jenkins**, ARMA(1,1) en el nivel (`a1.inp`, n=197,
semillas φ=0.87 θ=0.48), fue se detiene en

    φ = 0.999978   θ = 0.699565   logL = -57.603862

y el C original de Mauricio, **desde las mismas semillas y sobre los mismos
datos**, alcanza

    φ = 0.908683   θ = 0.575839   logL = -50.745092

**6.86 de verosimilitud mejor**.

⚠ Esta ficha decía aquí «y además es el valor publicado: Box y Jenkins dan
φ≈0.92, θ≈0.58». **Eso no está verificado** y no debe leerse como si lo
estuviera: el archivo que acompaña a las series (`Box_y_Jenkings/index.html`)
sólo trae la descripción de los datos, no la tabla de estimaciones del libro. Lo
único fijado contra el libro en toda la suite es θ=0.70 de la IMA(0,1,1) de la
Serie A (`tests/test_published_benchmark_series_a.py`). Y hay una razón de fondo
para no esperar coincidencia fina: **los valores del libro son de mínimos
cuadrados no lineales con retropredicción, no de ML exacta** —la misma distancia
que separa a TASTE de este motor, documentada en `docs/PROVENANCE.md` §3.6—, así
que lo razonable es acuerdo en dos decimales, no en siete.

Las otras ocho comparaciones del banco coinciden a **1e-11 o mejor**, así que no
es un problema general del motor: es este caso.

## Impact

Medio. El punto donde fue se detiene no es un óptimo local cualquiera: es la
**frontera del AR**, φ→1, donde la ARMA(1,1) degenera en la IMA(1,1) —que es
justamente `a2.inp`, y ahí fue clava el resultado de DRVUS a 1.5e-11—. Así que
el modelo que devuelve es interpretable, pero es el equivocado, y lo devuelve sin
avisar.

Importa más de lo que su severidad sugiere porque **la frontera es donde esta
suite trabaja**: el MEG, el DCD y los contrastes de no invertibilidad viven ahí.

## Reproduction

```python
import fue
ts, m = fue.load(".../Box_y_Jenkings/SeriesA/a1.inp")   # requiere BUG-0011 arreglado
m.fit()
m.ar[0][0], m.ma[0][0], m.loglik      # 0.999978, 0.699565, -57.603862
```

El `.out` de DRVUS conserva el recorrido y es lo más informativo del caso:

```
iter  0   φ=0.870000  θ=0.480000  |g|=14.684
iter 10   φ=0.999550  θ=0.526354  |g|= 0.666
iter 30   φ=0.999979  θ=0.700636  |g|= 0.047     <- donde fue se queda
iter 50   φ=0.995054  θ=0.715333  |g|= 0.006
iter 64   φ=0.908683  θ=0.575839  |g|= 0.000     <- DRVUS escapa
```

DRVUS **pasa por el punto de fue** hacia la iteración 30 y sale de él.

## Lo medido — 12-ago-2026

El motor C no devolvía el veredicto del optimizador: `raxopt` lo anuncia por
`outputv`, que el binding manda a `/dev/null`, así que la diagnosis estaba
inerte. Propagado (`qn_last_termcode`, `qn_last_nit`, `qn_last_gnorm`, que sólo
**registran** lo que raxopt ya calculó), el caso se lee de inmediato:

| ajuste | termcode | iteraciones | ‖g‖ | logL |
|---|---|---|---|---|
| `a1` tal cual | **2 — criterio de PASO** | 23 de 500 | 0.0101 | −57.603862 |
| `a1` desde μ₀=17 | 1 — gradiente | 7 | 0.0000004 | **−50.745092** |
| `a2` (IMA(0,1,1)) | 1 — gradiente | 7 | 0.0000000 | −53.508690 |
| DRVUS, el C original | gradiente | 64 | 0.0000000 | −50.745092 |

**fue no agota nada: se para.** El criterio que dispara no es el del gradiente
—que vale 0.0101, lejos de cero— sino el de paso: los iterados dejaron de
moverse. Y aun así `fit()` devolvió `converged = True`, porque `converged` era
`ifault == 0`, es decir «el motor no reventó», no «esto es un máximo».

### El disparador es la semilla de μ, y sólo ella

`a1.inp` trae μ₀ = 2.5 mientras la media de la serie es **17.0624**. Barriendo
μ₀ con todo lo demás intacto:

| μ₀ | 0 | 0.5 | 1 | 1.7 | 2.5 | 4 | 6 | 8 | 12 | 17 | 20 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| | ✔ | ✘ | ✘ | ✘ | ✘ | ✘ | ✔ | ✔ | ✔ | ✔ | ✔ | ✘ |

✔ = φ=0.908685, θ=0.575841, μ=17.0653, logL=−50.745092, termcode 1 — **el óptimo
de DRVUS a seis decimales**. ✘ = frontera, termcode 2.

Es decir: **el optimizador de fue alcanza el óptimo publicado**, desde μ₀=0 (86
iteraciones), desde μ₀=17 (7) y desde casi todo el intervalo intermedio. Lo que
hay no es un motor que no sabe llegar, sino una franja de semillas de μ desde la
que el camino pasa por una región plana pegada a φ=1 donde el criterio de paso
dispara antes de salir. El recorrido de DRVUS enseña esa salida: hacia la
iteración 40 μ salta de 2.59 a **16.84** —la media— y sólo entonces el AR vuelve
de 0.99995 a 0.91.

## Lo descartado, ahora con medida

- **No es el presupuesto de iteraciones** (era la hipótesis 1). fue se detiene en
  la 23 de 500. Subirlo no cambia nada: el criterio que dispara es el de paso.
- **No es `nlatools`** (hipótesis 2), ni la reescritura sobre GSL: el mismo
  binario alcanza el óptimo de DRVUS a seis decimales desde otras semillas de μ.
- **No son las tolerancias.** `cmacheps()` es idéntica en ambos árboles y
  devuelve `DBL_EPSILON`; DRVUS usa `gradtol = macheps^(1.1/3)`,
  `steptol = macheps^(2/3)`, `maxits=500`, `nrits=10` (`drvus.c:832-839`) y
  `fue_api.c:69-71` pasa exactamente los mismos valores.
- **No es la búsqueda restringida.** Con `chkma=False` el resultado no cambia en
  ningún dígito, así que el chequeo de invertibilidad —donde vive la única
  sustitución funcional, `gsl_eigenqr`— no interviene aquí.
- **No es el optimizador.** `drvus/src/qnewtop.c` contra
  `csrc/internal/qnewtopt.c`: 47 líneas, todas cabecera GPL, `#include`,
  `printf` → `if (outputv) fprintf(…)` y el registro nuevo. Ninguna numérica.
- **No es la verosimilitud.** Evaluada en el punto de DRVUS da exactamente
  −50.745092.
- **No es el porte a Python.** Los dos motores de fue se paran juntos
  (−57.603862 el C, −57.649946 `raxopt` en Python puro).

### Lo que queda: los dos binarios recorren caminos distintos, y ahora se ve

Con el veredicto propagado, el banco de Box-Jenkins compara **recorrido** y no
sólo destino, contra el `.out` que dejó DRVUS en ~2001:

| caso | logL fue | logL DRVUS | dif | termcode | iteraciones |
|---|---|---|---|---|---|
| a1 | −57.603862 | −50.745092 | **−6.86** | **2 / 1** | **23 / 64** |
| a2 | −53.508690 | −53.508690 | 1.5e-11 | 1 / 1 | 7 / 7 |
| b | −1249.974933 | −1249.974933 | −4.0e-11 | 1 / 1 | 3 / 3 |
| c | 131.668147 | 131.668147 | −4.7e-12 | 1 / 1 | 2 / 3 |
| c2 | 123.399306 | 123.399306 | 3.9e-11 | 1 / 1 | 6 / 6 |
| d | −67.751685 | −67.751685 | −2.6e-10 | 1 / 1 | **50 / 47** |
| d1 | −76.691867 | −76.691867 | −2.8e-11 | 1 / 1 | 3 / 3 |
| e1 | −414.617409 | −414.617409 | −3.4e-12 | 1 / 1 | 20 / 21 |
| e2 | −412.494817 | −412.494817 | 8.2e-12 | 1 / 1 | 24 / 25 |

Léase junto al otro contraste: contra los 28 `.out` de **`fue-1.13.1`** —su
propia línea de C— el puerto reproduce el número de iteraciones **exactamente,
28 de 28**. Contra el binario de **2001** el destino sigue siendo el mismo a
1e-10, pero el recorrido **deriva** donde el problema cuesta: 50 iteraciones
frente a 47 en `d`, 2 frente a 3 en `c`.

Eso ya no es conjetura: los dos binarios siguen caminos ligeramente distintos, y
está medido. En ocho casos la deriva no cambia nada porque el óptimo es nítido.
En `a1` —valle plano pegado a φ=1, μ sembrado en 2.5 contra una media de 17—
decide entre salir y no salir. **Lo que sigue sin saberse es por qué el binario
de 2001 cae del lado bueno**, y eso no se puede medir: ese binario no se
reproduce. Perseguirlo no rinde; lo que rinde es lo de abajo.

## Fix

Tres piezas, y ninguna toca el algoritmo publicado:

1. **Hecho: el motor devuelve su veredicto.** `FitResult.termcode`, `.niter` y
   `.gnorm` vienen del C (`qn_last_*` en `csrc/internal/qnewtopt.c`, que sólo
   registran). El informe `.out` ya tenía el hueco —estaba siempre vacío— y
   ahora escribe lo mismo que escribía el C.
2. **Hecho: `converged` ya no significa «no reventó».** Es
   `ifault == 0 and termcode in (0, 1)`, y un alto que no es máximo sale por
   `RuntimeWarning` con el motivo, las iteraciones y ‖g‖. El fallo del motor
   sigue siendo excepción; un ajuste que existe pero no es máximo, no: es algo
   sobre lo que hay que poder decidir. Era la carencia de BUG-0005.
3. **La semilla de μ es responsabilidad de quien escribe el `.inp`.** μ₀=2.5
   contra una media de 17 es una semilla obsoleta, no una especificación; `art`
   ya siembra μ de la media de la variable diferenciada, y por eso este caso no
   aparece por su camino. Documentar en `docs/FILE_CONTRACT.md` §2.4.

El caso se queda **abierto** por (3) y por la divergencia en sí, y el banco lo
sostiene por los dos lados: `test_box_jenkins_series.py` falla si `a1` empieza a
coincidir con DRVUS por su cuenta, y `test_optimizer_termcode.py` falla si deja
de avisar.

## Lo que apareció de camino

Propagar el veredicto convirtió los `.out` preservados en **oráculo**, porque
registran termcode, iteraciones y ‖g‖ de ejecuciones del C que ya nadie puede
repetir. Sobre los 28 del banco de casos reales, el puerto reproduce **los tres,
exactamente** — mismo código de parada, mismo número de iteraciones y la misma
norma del gradiente a 1e-4 relativo.

Queda una discrepancia entre paquetes, anotada y no resuelta aquí: `drtran`
cuenta el termcode 2 **como convergencia** (`estimate.py:203`,
`converged = termcode in (1, 2)`) y lo matiza aparte en `convergence_note`,
mientras que `fue` a partir de ahora no. La razón de `drtran` es buena para su
caso —siembra del `.pre`, donde el paso es diminuto por construcción— pero
`Coint/R.4`, con ‖g‖ = 1e5, muestra que el 2 también tapa lo contrario. Es
decisión de arquitectura de la suite, no de esta ficha.

Y de paso deshizo un malentendido documentado en `tests/test_real_cases.py`:
`Coint/R.4` no es «distinto BLAS, distinto óptimo». Es un ajuste que **nunca
convergió** —termcode 2, ‖g‖ = 116330.0394, y así lo dice su propio `.out` del
C—. Que las plataformas discrepen sobre él es la consecuencia, no la causa.

## Validation

`tests/test_box_jenkins_series.py` lleva el caso marcado como divergencia
conocida y las otras ocho como igualdad a 1e-9. Si `a1` empieza a coincidir, el
test lo dirá.

## Consecuencia para el archivo histórico

Los `.out` de los años 2000 son ejecuciones **a 80 bits**. Eso no los invalida
—en 8 de los 9 casos del banco de Box-Jenkins el destino coincide a 1e-10, y en
28 de 28 casos reales coinciden hasta el número de iteraciones— pero sí obliga a
decir qué son: la referencia de `a1` no es «lo que da el programa», es «lo que da
el programa con intermedios de 80 bits».

La buena noticia es que **son reproducibles**:

```bash
cp -r ~/Dropbox/SRC/drvus-source/1.2.01/drvus/src /tmp/drvus && cd /tmp/drvus
sed -i 's/\bround\s*(/drvus_round(/g' nlatools.c diagnose.c drvus.h   # choca con round() de C99
gcc -O0 -mfpmath=387 -o drvus drvus.c drvmlest.c elfvarma.c usmelard.c \
    qnewtopt.c nlatools.c diagnose.c -lm
./drvus a1 eml chk         # reproduce a1.out de 2001, iteración por iteración
```

El único cambio en la fuente es el nombre de `round()`, que en 2001 no chocaba
con `math.h` y hoy sí. Nada más.
