---
id: BUG-0005
title: ML optimizer converges to a spurious optimum on multimodal (seasonal-AR) likelihoods and reports converged=True with no diagnostic; the basin is platform-dependent (Windows vs Linux)
status: in-progress
severity: medium
component: estimation
found_in: 0.1.7
fixed_in:
reported: 2026-07-22
reporter: mtgp2 (DVR / GP-Note replication)
tags:
  - optimizer
  - convergence
  - multimodal-likelihood
  - seasonal-ar
  - reproducibility
  - platform-dependent
references:
  - csrc/internal/qnewtopt.c, drvmlest.c, elfvarma.c (the ML optimizer + likelihood)
  - src/fue/_build_cffi.py (Windows /O2 + vcpkg GSL vs Linux -O2 + system GSL)
  - ART BUG-0006 (art-python/bugs/BUG-0006-*): the seed contamination that first triggered it
  - art-python/bugs/BUG-0006-repro/ (US_CPI.pre + repro.py — reproduces on Windows)
  - Garcia-Hiernaux, Gonzalez-Perez & Guerrero (2026), Econ. Modelling 157, Table 2 (US CPI)
---

## Summary

On a **multimodal** likelihood — an AR(2)×AR(2)_12 model whose seasonal AR has
complex conjugate roots (US CPI, monthly, 2002–2026) — fue's C optimizer
(`qnewtopt`/BFGS) converges to a **different basin depending on the build**: from the
*same* source, *same* starting values and *same* data, the **Windows** wheel (MSVC
`/O2` + vcpkg GSL) settles on a **spurious optimum** (`μ̂=−0.144`, `σ_a=0.305`,
`AIC=−2511`, ~52 lower log-likelihood) while the **Linux** wheel (manylinux, GCC
`-O2` + system GSL) reaches the correct one (`μ̂=0.0021`, `σ_a=0.261`, `AIC=−2613`,
paper Table 2). In **both** cases fue reports **`converged=True`, `ifault=0` with no
diagnostic** — the absurd mean and the worse `σ_a`/AIC are not flagged.

Two things, then:
1. **Non-reproducibility** of the optimizer across builds on multimodal surfaces
   (the basin depends on last-ULP floating-point differences).
2. **No guard** against a clearly-degenerate optimum being reported as success.

## Impact

Silent, wrong "estimated" model on the affected platform — the estimate is not
reproducible across OSes for strongly-multimodal specs, and there is no signal that
anything went wrong. Surfaced in the DVR replication of Garcia-Hiernaux et al. (2026):
the US CPI row of Table 2 came out degenerate on Windows, corrupting everything
downstream of the US residuals (DVR, GARCH/GJR, Beveridge–Nelson). The consumer had
to hard-code per-country starting values, and drtran's exact-VARMA engine was used to
cross-check.

## Reproduction

Self-contained in `art-python/bugs/BUG-0006-repro/` (`US_CPI.pre`, n=293; `repro.py`):

| build | default seed → | verdict |
|---|---|---|
| Windows — fue 0.1.7 wheel (MSVC/vcpkg) | μ=−0.144, σ=0.305, AIC=−2511 | **SPURIOUS** (converged=True) |
| Linux — fue 0.1.7 wheel (manylinux)    | μ=+0.0021, σ=0.261, AIC=−2613 | correct |
| Linux — fue 0.1.7 editable (source)    | μ=+0.0021, σ=0.261, AIC=−2613 | correct |

Verified on Linux in an isolated venv with the exact PyPI stack (`fue==0.1.7`,
`art-tseries==0.1.2`): the bug does NOT reproduce — same source, same seed, correct
result. It reproduces on the reporter's Windows wheel. That a repro labelled
"fue≥0.1.7" fails to reproduce with fue 0.1.7 on Linux **is itself the finding**: it
is a build/platform-dependent numerical issue of the optimizer, not a logic bug.

## Root cause

The ML optimizer (`qnewtopt`/BFGS in `csrc/internal/`) is a gradient-based local
search. On a **unimodal** surface it lands on the same optimum regardless of last-bit
differences; on the US-CPI **multimodal** AR(2)×AR(2) surface the basin is decided by
tiny floating-point differences in the log-likelihood evaluation. Those differences
come from the **build**: MSVC vs GCC FP semantics (operation reordering / FMA
contraction under `/O2` vs `-O2`), the C runtime `libm` (MSVC CRT vs glibc:
`exp`/`log`/`cos`/`sin`/`pow`, used heavily by the log-likelihood, Box-Cox and the
harmonics), and the GSL build (vcpkg static vs system). No `long double` is involved.
Separately, fue never checks whether the "converged" optimum is sane, so a degenerate
basin is reported as success.

A wrong-sign seasonal-AR *seed* (ART BUG-0006, now fixed on the ART side) is what
first pushed the search toward the wrong basin; a correct seed makes this case robust.
But the underlying fue fragility remains: a sufficiently multimodal spec could still
diverge on some build, and would still be reported as `converged=True`.

## Fix

Not yet applied. Robustness, not bit-identical FP (which is infeasible across
compilers):
- **Guard against absurd optima.** After a "converged" fit, check that `μ̂` lies
  within a few sample-std of the differenced-series mean, and/or that the objective
  is not materially worse than at the starting values; if it fails, downgrade
  `converged`, raise `ifault`, or emit a warning instead of reporting silent success.
- **Multi-start for multimodal blocks.** When a seasonal AR/MA block is present,
  optimise from a small set of starting points (e.g. ±seed, HR/YW seed) and keep the
  best optimum — platform-independent (the best basin wins on any build).
- *(Optional, palliative)* homogenise FP across builds: try `/fp:precise` (MSVC) and
  `-ffp-contract=off` (GCC) so the wheels agree; does not remove the multimodality.

## Validation

When fixed, `repro.py` must reach the correct optimum (`σ_a≈0.261`, `μ̂≈0.0021`,
`AIC≈−2613`) on **every** platform from the default seed, or fail loudly rather than
report `converged=True` on the spurious basin. The other seven DVR economies (unimodal
or non-seasonal-AR) must be unchanged.


---

## Segundo intento (2026-09-07): la mitad 2 cerrada, la mitad 1 no reproduce

El informe declaraba **dos** cosas. Se ha cerrado una y se ha vuelto a probar la
otra.

### La mitad 1 — la no-reproducibilidad — sigue sin reproducir en Linux

`repro.py` con el `fue` de hoy, desde las **dos** semillas:

    default seed     → μ̂ = +0.002149   σ_a = 0.2608   converged=True
    identified seed  → μ̂ = +0.002149   σ_a = 0.2608   converged=True

Las dos llegan al óptimo **correcto** (Tabla 2 del artículo), y son idénticas
entre sí. La divergencia entre semillas que el informe describía **ya no
existe** aquí: la contaminación de la semilla estacional (art/BUG-0006) está
arreglada.

Lo que **no** se puede verificar desde esta máquina es la divergencia entre
PLATAFORMAS, que es la que da nombre al informe. Sigue abierta por falta de
Windows, no por falta de intento — y es aritmética de coma flotante de dos
compiladores, no algo que se arregle en el código.

### La mitad 2 — el disparate reportado como éxito — está cerrada

Esa sí se podía cerrar aquí, y es la que de verdad muerde: **un óptimo absurdo
se reportaba como `converged=True, ifault=0` sin un solo aviso.**

`model._avisa_si_la_media_es_absurda` compara μ̂ con la media de **`w`** —la
serie ya filtrada de deterministas y diferenciada, que es de la que μ es la
media— y avisa si se aleja más de `UMBRAL_MEDIA_ABSURDA = 5` desviaciones
típicas.

**El umbral está medido, no conjeturado.** Sobre los 1.570 modelos del
ecosistema con μ estimada:

    percentil 50                    0.0024
    percentil 99                    0.3313
    percentil 99.9                  0.7281
    MÁXIMO observado                1.3792
    por encima de 2·sd                   0

    el óptimo espurio de US CPI        47.2   ← lo que se quiere cazar

Entre lo peor legítimo (1,4) y el disparate (47,2) hay un factor de **34**.
Cualquier umbral entre 3 y 10 sirve; el 5 está cómodo en medio.

**Y sobre el corpus entero la guarda dispara 0 veces** — 1.570 óptimos juzgados,
ningún falso positivo — mientras que sobre el óptimo espurio dispara.

### Dos cosas que la medición enseñó y que no estaban en el informe

**μ es la media de `w`, no de la serie observada.** Con intervenciones o
armónicos la diferencia es enorme: los deterministas se llevan parte del nivel.
Comparando contra la serie cruda, 23 modelos salían «absurdos» siendo correctos.
`w` no se reconstruye: la da `cast_us_py`, la misma que usa el motor —
reimplementarla sería repetir BUG-0014.

**μ puede no estar identificada, y entonces no hay nada que juzgar.** El término
de deriva es μ·φ(1); con un AR de raíz unitaria escrito con `d=0`, φ(1) = 0 y μ
no entra en la verosimilitud: puede valer cualquier cosa. Son **47 de los 1.570**
modelos del ecosistema —los VIX con λ extrema y φ = 1 exacto— y avisar de ellos
serían 47 falsos positivos. La guarda los excluye explícitamente.

### Lo que queda

  - la no-reproducibilidad entre plataformas, **sin verificar** por falta de
    Windows;
  - el **multi-arranque** para bloques multimodales, que el informe propone y
    sigue sin hacerse. Con la guarda puesta, el caso al menos ya no pasa
    callando: avisa y dice qué hacer.

Por eso el informe queda `in-progress` y no `fixed`: la mitad que se podía cerrar desde aquí está cerrada, y la que no —la divergencia entre plataformas— sigue esperando una máquina Windows.


---

## La causa raíz de este informe está mal, y hoy se ha medido (2026-09-07)

El informe atribuye el caso a una **superficie multimodal** cuya cuenca decide la
aritmética del compilador. **La evidencia disponible hoy no sostiene esa
explicación**, y sí sostiene otra que el propio informe menciona de pasada:

> *«A wrong-sign seasonal-AR seed (ART BUG-0006, now fixed on the ART side) is
> what first pushed the search toward the wrong basin.»*

Era la línea importante y quedó como nota al margen.

### Prueba 1 — nueve arranques, un solo óptimo

Multi-arranque sobre el AR estacional de US_CPI, cubriendo la región y los dos
signos:

    (−0.109,−0.093)  ℓ=1322.513      (+0.800,−0.500)  ℓ=1322.513
    (+0.109,+0.093)  ℓ=1322.513      (+0.900, 0.000)  ℓ=1322.513
    (+0.500,+0.300)  ℓ=1322.513      ( 0.000,+0.900)  ℓ=1322.513
    (−0.500,−0.300)  ℓ=1322.513      ( 0.000, 0.000)  ℓ=1322.513
                                     (+0.950,−0.900)  ℓ=1322.513

**Nueve de nueve al mismo sitio**, incluido el signo cambiado. Dos semillas
fuera de la región estacionaria las rechaza el motor, que es lo correcto.

Si la superficie tuviera cuencas separadas, un barrido así las encontraría.

### Prueba 2 — el punto de Windows no es un punto estacionario

Objetivo a lo largo de μ, con el resto en el óptimo:

    μ = +0.002149    −268.677073    ← el óptimo
    μ = −0.050000    −268.859130         −0.18
    μ = −0.144000    −270.106980         −1.43
    μ = −0.300000    −274.788726         −6.11

    derivada del objetivo en μ = −0.144:  +19.57

**Decrece monótonamente y la derivada no se anula.** En la dirección de μ no hay
segunda cima: μ = −0,144 **no es un óptimo**, así que llamarlo «otra cuenca»
describe mal lo que ocurrió.

### Lo que se puede afirmar y lo que no

**Sostenido:** en esta plataforma y con este `fue`, el caso llega al mismo óptimo
desde nueve arranques distintos, y a lo largo de μ la verosimilitud es unimodal.
Arreglada la semilla, el caso es **robusto**.

**No sostenido, y hay que decirlo:** sólo se conocen μ̂ y σ̂ₐ del resultado de
Windows, no el vector completo de 16 parámetros. No se puede descartar que aquel
punto fuera un óptimo local en el espacio conjunto. Y sigue sin haber máquina
Windows para comprobarlo.

### La lectura que encaja con todo

La convención de signo de `fue` es la de Box y Jenkins para **todo** operador
—ω(B) = ω₀ − ω₁B − …— y eso confundió al consumidor, que sembró el AR estacional
con el signo cambiado. **Una semilla con el signo al revés no es «un punto de
partida algo peor»: es un punto de partida en la región equivocada**, y desde
ahí cualquier diferencia de último bit puede decidir adónde se va. Con la
semilla bien puesta, el caso no se mueve.

Eso reordena las conclusiones del informe:

  - el **multi-arranque** que proponía como arreglo **no hace falta para este
    caso** — nueve arranques dan lo mismo. Seguiría siendo prudente para
    superficies genuinamente multimodales, pero este no lo es;
  - lo que quedaba de defecto real era **el silencio**, y ése está cerrado;
  - la **no-reproducibilidad entre plataformas** no se ha refutado, pero tampoco
    tiene hoy ninguna evidencia a favor más allá de aquella observación única.

Queda `in-progress` porque no se puede cerrar lo que no se puede comprobar. Pero
el motivo por el que sigue abierto ya no es «el optimizador es frágil»: es
**«falta una máquina Windows para verificar una observación de julio que el
resto de la evidencia no acompaña»**, que es una razón muy distinta y mucho más
pequeña.
