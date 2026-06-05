# FUE — Arquitectura y rendimiento del estimador

---

## Arquitectura actual (2026-06-05)

```
model.fit()  /  fue model.inp  (CLI)
    └── _engine.estimate(model)
            │
            ├── _fue_engine.abi3.so presente?  ──YES──► lib.fue_estimate()
            │                                            C puro (cffi):
            │                                            cast_us + raxopt + flikam + elf
            │                                            (fue_api.c + GSL)
            │
            └── NO (ImportError) ──────────────────────► estimate_py(model)
                                                         Python puro:
                                                         cast_us_py + L-BFGS-B
                                                         + flikam_scalar + elf_scalar
```

**Módulos siempre en Python puro** (sin C, sin GSL):

| Módulo | Función |
|--------|---------|
| `series.py` | `TimeSeries` — datos y metadatos |
| `model.py` | `Model` — especificación, `fit()`, `forecast()` |
| `intervention.py` | `Intervention`, `FixedFreqFactor` |
| `inp.py` | Parser `.inp` — `fue.load()` |
| `report.py` | Informe ASCII `.out` — `write_out()`, `write_pre()` |
| `forecast.py` | Pronóstico ARMAX — `Model.forecast()` |
| `diagnostics.py` | ACF, PACF, Jarque-Bera, Ljung-Box |
| `plots.py` | Gráficos (matplotlib) |
| `cli.py` | Entrada CLI — `fue model.inp` |
| `elfvarma.py` | `elf_scalar` (Mauricio 1995), `flikam_scalar` (Melard 1984) |
| `cast_us.py` | `cast_us_py`, `estimate_py` — estimador Python puro |

**Módulo condicional (C cuando disponible, Python como fallback):**

| Módulo | Con C | Sin C |
|--------|-------|-------|
| `_engine.py` | llama `lib.fue_estimate()` vía `_fue_engine.abi3.so` | llama `estimate_py()` |

**Cómo forzar el camino Python puro:**

```bash
FUE_SKIP_C=1 pip install fue     # no compila .so; Python puro permanente
fue modelo                        # usa C si .so presente, Python si no
```

---

## Registro de rendimiento del estimador a lo largo de la migración

Este documento registra el coste de cada etapa de la migración C→Python
y establece la base para la decisión sobre el modelo híbrido definitivo.

Referencia hardware: Linux x86-64, CPython 3.12, gcc -O2.  
Script reproducible: `benchmarks/bench_estimator.py`.

---

## Casos de prueba canónicos

| Caso | Descripción | nobs | npar |
|------|-------------|-----:|-----:|
| AR(1) | AR(1) anual, sin intervenciones | 30 | 1 |
| SFNY.2 | step+delta + AR(1)×AR(2) + mu, boxlam=0 | 62 | 6 |
| RIPC.1 | 5 cos/sin + alter + step+delta + AR(1) fijo + mu, mensual | 72 | 14 |

---

## Etapas de la migración

### Etapa 0 — Híbrido Python-API + motor C  *(Fases 1–3)*

**Arquitectura**: Python define el modelo; `_engine.estimate()` delega la
estimación completa a la extensión cffi (`fue_estimate()`), que ejecuta en C:
`cast_us` + bucle del optimizador (`raxopt`) + `flikam` + `elf`.

Python no interviene en ningún bucle numérico.  
**Rendimiento = C puro.**

| Caso | Tiempo (ms) |
|------|------------:|
| AR(1) | 0.1 |
| SFNY.2 | 2.7 |
| RIPC.1 | 37 |

Esta es la línea base de rendimiento máximo. Requiere GSL y compilador C.

---

### Etapa 1 — Python puro  *(Fase 4, estado actual)*

**Arquitectura**: `estimate_py()` ejecuta todo en Python:
- `cast_us_py` desempaqueta parámetros y diferencia la serie
- `flikam_scalar` (Melard 1984) evalúa la función objetivo
- L-BFGS-B de scipy gestiona la optimización (gradiente numérico)
- `elf_scalar` (Mauricio 1995) evalúa el loglik final exacto

No requiere C ni GSL. Activo como fallback automático cuando la extensión C
no está disponible.

#### Coste por evaluación de la función objetivo (ms, mejor de 20)

| Caso | `cast_us_py` | `flikam_scalar` | Total |
|------|-------------:|----------------:|------:|
| AR(1) | 0.06 | 0.08 | 0.14 |
| SFNY.2 | 0.74 | 0.23 | 0.96 |
| RIPC.1 | 1.46 | 0.12 | 1.58 |

#### Evaluaciones del optimizador (L-BFGS-B, gradiente por diferencias finitas)

| Caso | npar | nfev | pasos de gradiente |
|------|-----:|-----:|-------------------:|
| AR(1) | 1 | 32 | 16 |
| SFNY.2 | 6 | 504 | 72 |
| RIPC.1 | 14 | 2985 | 199 |

`nfev ≈ (npar + 1) × grad_steps` porque cada gradiente requiere `npar+1`
evaluaciones (diferencias finitas hacia delante).

#### Tiempos totales y factor de ralentización

| Caso | C (ms) | Python puro (ms) | Factor |
|------|-------:|-----------------:|-------:|
| AR(1) | 0.1 | 7.8 | ×90 |
| SFNY.2 | 2.7 | 653 | ×240 |
| RIPC.1 | 37 | 6577 | ×174 |

**Causa**: los bucles Python puros en `cast_us_py` y `flikam_scalar` son
~100× más lentos que el equivalente C compilado, y el optimizador los llama
cientos de veces. No hay ninguna diferencia algorítmica — es únicamente
velocidad del intérprete.

Desglose del tiempo para SFNY.2 (por perfil `cProfile`):

| Función | Tiempo propio | % del total |
|---------|-------------:|------------:|
| `cast_us_py` | 472 ms | 63% |
| `flikam_scalar` | 142 ms | 19% |
| `calcnu_py` | 24 ms | 3% |
| `_unscramble` | 22 ms | 3% |
| `_twacf_scalar` | 20 ms | 3% |
| scipy overhead | ~17 ms | 2% |

---

## Opciones de arquitectura híbrida

Una vez completada la migración, hay tres arquitecturas posibles:

### Opción A — Híbrido: optimizador Python + bucle interno C

El optimizador L-BFGS-B corre en Python; cada evaluación de la función
objetivo llama a una función C expuesta vía cffi que ejecuta
`cast_us(x) + flikam(result)` de forma compilada.

**Requiere**: exponer una nueva función en `fue_api.c` y `fue_api.h`
(algo como `fue_objcfunc(x, spec, pi10x, pi20x)` → objetivo normalizado).

**Proyección** (asumiendo ~5 µs por evaluación en C):

| Caso | nfev | Proyección (ms) | vs. C puro |
|------|-----:|----------------:|-----------:|
| AR(1) | 32 | 0.2 | 1.2× |
| SFNY.2 | 504 | 2.5 | 0.9× |
| RIPC.1 | 2985 | 14.9 | 0.4× |

Rendimiento prácticamente igual al C puro. Requiere GSL.

**Trade-off**: igual dependencia que la Etapa 0, pero con el optimizador en
Python (más fácil de extender: criterios alternativos, restricciones, etc.).

### Opción B — Híbrido: optimizador Python + `cast_us_py` Python + C `flikam`

El bucle de `flikam` corre en C (cffi directo), pero `cast_us_py` sigue en
Python. Mejora parcial: `flikam` es solo el 19% del tiempo en SFNY.2 pero
`cast_us_py` (63%) permanece en Python.

**Proyección para SFNY.2**: 0.74ms × 504 + ~0 = ~373ms. Solo ×2.5 mejor
que Python puro, todavía ×138 peor que C.

No es una opción interesante: coste de implementación sin beneficio real.

### Opción C — Python puro optimizado (sin C)

Vectorizar los bucles de `cast_us_py` y `flikam_scalar` con NumPy/scipy.
Posibles herramientas: `np.convolve`, `scipy.signal.lfilter`.

No cambia los algoritmos publicados (convolucion es convolucion), pero
requiere reescribir los bucles. Posible ganancia estimada: ×5–20 en las
secciones vectorizadas → Python total quizá ×20–50 más lento que C.

---

## Decisión pendiente

| Opción | Rendimiento | Dependencias | Complejidad |
|--------|-------------|--------------|-------------|
| **Etapa 0** (actual C puro) | C puro | GSL + compilador | ya implementado |
| **Opción A** (Py-opt + C inner) | ≈ C puro | GSL + compilador | nueva cffi API |
| **Etapa 1** (Python puro) | ×90–240 vs C | ninguna | ya implementado |
| **Opción C** (Py vectorizado) | ×5–20 vs C puro | ninguna | refactor código |

La decisión híbrida se tomará una vez completada la migración, con la
Opción A como candidata principal si se necesita rendimiento C sin sacrificar
la flexibilidad del optimizador Python.

---

## Historial de mediciones

Reproducir: `pytest tests/test_performance.py::TestCvsPython::test_summary_table -v -s`

| Fecha | Commit | Caso | C (ms) | Python (ms) | Factor | Notas |
|-------|--------|------|-------:|------------:|-------:|-------|
| 2026-06-04 | 53e1699 | AR(1) | 0.1 | 7.8 | ×90 | Etapa 1 completada |
| 2026-06-04 | 53e1699 | SFNY.2 | 2.7 | 653 | ×240 | |
| 2026-06-04 | 53e1699 | RIPC.1 | 37 | 6577 | ×174 | |
| 2026-06-05 | f8d22bb | AR(1) | 0.1 | 6.2 | ×71 | test_performance.py añadido |
| 2026-06-05 | f8d22bb | SFNY.2 | 2.7 | 773 | ×289 | |
| 2026-06-05 | f8d22bb | RIPC.1 | 42 | 7864 | ×186 | |
