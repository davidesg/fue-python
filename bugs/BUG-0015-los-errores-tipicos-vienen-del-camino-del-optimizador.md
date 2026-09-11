---
id: BUG-0015
title: Los errores típicos vienen de la matriz que BFGS acumula por el CAMINO, no del hessiano en el óptimo — dos ejecuciones del mismo modelo dan SE distintos
status: open
severity: high
component: estimation
found_in: 0.1.7
fixed_in:
reported: 2026-07-12
reporter: David — al homologar drtran (puente fue → drvarma) contra fue
tags:
  - inferencia
  - errores-tipicos
  - covarianza
  - bfgs
  - limitacion-conocida
references:
  - csrc/internal/drvmlest.c:112 (la línea comentada: fdhess + choldcp)
  - fue-1.13.1/ERRORES_ESTANDAR.md (el estudio completo, con las tablas)
  - art-python/bugs/BUG-0090, BUG-0091 (el mismo hecho, visto desde art)
  - BUG-0012 (parar sin anular el gradiente: el caso extremo de esto)
---

## Summary

**Este informe existe porque la limitación no estaba en el índice.** El estudio
lleva desde julio en `ERRORES_ESTANDAR.md` y en dos commits del repo del C, pero
quien consultara `bugs/` no encontraba nada — y es la limitación más importante
que tiene el paquete hoy.

`fue` produce **errores típicos distintos en ejecuciones distintas del mismo
modelo**, con las **mismas** estimaciones puntuales.

| `ES_CPI_m10` | SE(μ) | SE(φ) |
|---|---|---|
| run A | **0.073304** | 0.062333 |
| run B | **0.028316** | 0.061577 |

μ = 0.154472 y φ = 0.402839 en las dos. Contra el **GLS exacto** sobre los mismos
datos, la run A se equivoca por un factor de **4 a 5** en los deterministas 3 a 6.
Y no es sesgo —la run A los infla, la run B infla unos y desinfla otros—: es
**ruido**.

**No afecta** a las estimaciones puntuales ni a la log-verosimilitud, que son
correctas. Afecta a la inferencia: SE, t, p-valores.

## La causa

`drvmlest.c` invierte el hessiano que **BFGS acumula a lo largo de la trayectoria
del optimizador**, no el hessiano evaluado en el óptimo. Esa matriz es un
subproducto del CAMINO: depende de por dónde se pasó y de cuántas iteraciones se
dieron, no sólo de dónde se llegó.

De ahí se sigue algo que no es evidente y que conviene tener presente al usar el
paquete: **arrancar en el óptimo produce los PEORES errores típicos**, porque el
optimizador no se mueve y la matriz se queda en la semilla. Es el motivo por el
que un `.pre` reestimado no puede publicar inferencia (art/BUG-0090, BUG-0091), y
el caso extremo es BUG-0012, parar sin anular el gradiente.

## El arreglo está identificado y NO se aplica

Descomentar `fdhess` en `drvmlest.c:112` —el hessiano por diferencias finitas— y
`drtran` demuestra que funciona: sus 17 SE del caso canónico clavan al binario y
no se mueven al perturbar el arranque.

**No se aplica porque `fue` es de uso general.** En un punto que no es el óptimo
el hessiano por diferencias finitas puede no ser definido positivo —en m6, 2 de
55 autovalores en las semillas— y `choldcp`, que es la Cholesky **modificada**,
parchearía los pivotes y publicaría números de aspecto impecable. La matriz del
BFGS no puede fallar así.

Es decir: el arreglo cambia un error **ruidoso y detectable** por uno **limpio e
indetectable**, salvo que se decida antes qué hacer cuando el hessiano no sea
definido positivo. El coste y el truncamiento de `xitol` se midieron y se
descartaron como explicación.

## Qué hacer mientras tanto

  - Los SE de un modelo estimado **desde su `.inp`**, con iteraciones de verdad,
    son utilizables con la reserva de arriba.
  - Los de un modelo reestimado **desde su `.pre`** no lo son: el optimizador no
    se mueve. `art` avisa de esto y separa `estimar` de `mirar` por esta razón.
  - El `.out` guarda la covarianza completa: es la constancia de lo que hubo,
    no una garantía de que sea el hessiano.

## Addendum medido — 11-sep-2026: el puerto Python YA tiene la vía, y el paso está mal

Medido al preguntar el analista si esto se puede arreglar de raíz. Tres hechos
nuevos, los tres reproducibles sobre `ES_CPI_m10`.

**1 · La vía del hessiano ya está cableada en Python.** `cast_us.py:551-566`
elige según el optimizador:

```python
if B_hess is not None:       # raxopt  → factor de Cholesky del BFGS (el CAMINO)
    cov[:, i] = 2.0 * obj_opt * _cholsol(B_hess, e) / n_eff
else:                        # lbfgsb  → _fdhess, el hessiano EN EL ÓPTIMO
    H = _fdhess(objective, x_opt.copy(), obj_opt, _SQRT_EPS)
    cov = 2.0 * obj_opt * np.linalg.inv(H) / n_eff
```

No hay que escribir `fdhess`: está, y se usa. Sólo que por esa rama sale
**cero**.

**2 · El paso es el equivocado, y eso explica parte de la no-definición
positiva.** `_SQRT_EPS = √ε ≈ 1,5e-08` es el paso de una derivada PRIMERA. Un
hessiano por diferencias necesita **ε^(1/4) ≈ 1,2e-04**: con √ε el error de
redondeo se amplifica por 1/η² ≈ 4,5e15 y lo que se mide es ruido.

Medido **en el óptimo** de `ES_CPI_m10` (13 parámetros, cond(H)=138):

| paso | diagonal negativa | SE de los 3 primeros |
|---|---:|---|
| ε^1/2 (el actual) | **2 de 13** | [0, 0, 0.03258] |
| ε^1/3 | 0 de 13 | [0.06833, 0.06829, 0.02769] |
| ε^1/4 | 0 de 13 | [0.06833, 0.06829, 0.02769] |
| ε^1/5 | 0 de 13 | [0.06833, 0.06829, 0.02769] |

**Estable en tres órdenes de magnitud del paso** — que es la firma de una
derivada bien calculada— y con la diagonal entera positiva.

Esto NO contradice el obstáculo del informe: aquél se midió **en las semillas**
(m6, 2 de 55 autovalores) y esto es **en el óptimo**. Pero sí acota la pregunta:
al menos parte de la no-definición positiva observada es un artefacto del paso y
no una propiedad del problema. La decisión pendiente hace falta menos veces de lo
que parecía.

**3 · Y el respaldo Python enmascara el fallo peor que el C.** Donde `drvmlest.c`
usaría `choldcp` —Cholesky modificada, que parchea pivotes—, aquí hay:

```python
std_errors = np.sqrt(np.maximum(diag, 0.0))
```

Una varianza negativa se convierte en **un error típico de 0.0**, publicado como
un número cualquiera. No es «impecable pero indetectable»: es cero, y un cero en
el denominador de una razón t no es un valor sospechoso, es un sinsentido.

### Y lo que contesta la pregunta de raíz

Con el paso corregido, la SE **deja de depender de dónde arrancó la optimización**:

| vía | origen | niter | SE de los 3 primeros |
|---|---|---:|---|
| BFGS (el camino) | `.inp` | 21 | [0.05667, 0.05619, 0.02705] |
| BFGS (el camino) | `.pre` | 9 | **[0.09554, 0.09621, 0.07887]** |
| hessiano ε^1/4 | `.inp` | 42 | [0.06833, 0.06829, 0.02769] |
| hessiano ε^1/4 | `.pre` | 12 | **[0.06833, 0.06829, 0.02769]** |

    diferencia máxima hessiano `.inp` vs `.pre` = 6,5e-09

Eso es exactamente lo que el paquete no tiene hoy. Y con ello caerían, por
innecesarias, media docena de defensas construidas alrededor del síntoma: el
detector de covarianza-semilla, el de casi-semilla, la cláusula de las SE en el
convenio de ficheros de `art`, y el aviso que publica factores de 4,23×.

**Lo que NO resuelve, y sigue siendo la decisión pendiente:** qué hacer cuando el
hessiano no sea definido positivo con el paso correcto. Y una consecuencia que
hay que mirar de frente antes de adoptarlo: **las SE cambian ~20 % respecto a las
que el paquete publica hoy** —[0.0683, 0.0683] frente a [0.0567, 0.0562] desde el
`.inp`— así que no es un arreglo transparente: reescribe la inferencia de todo lo
ya calculado.

## Lo que falta para cerrarlo

Una sesión propia: barrido empírico sobre la batería y **decisión de qué hacer
cuando el hessiano no sea definido positivo** —rechazar, avisar, o caer a la
matriz del BFGS diciéndolo—. Sin esa decisión, aplicar `fdhess` empeora el modo
de fallo aunque mejore el número.
