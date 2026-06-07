# FUE Python — Estado de la migración C → Python

Referencia: `fue-1.13.1` es el código C fuente de verdad.  
Última actualización: 2026-06-07 (rev 8)

---

## Resumen ejecutivo

| Módulo | Estado | Motor |
|--------|--------|-------|
| Lectura de `.inp` | ✅ completo | Python puro |
| Estimación ML exacta | ✅ completo | Python puro (raxopt) |
| Evaluación de verosimilitud | ✅ completo | Python puro (elf_scalar + flikam_scalar) |
| Pronóstico (fuf) | ✅ completo | Python puro |
| Generación de `.out` | ✅ completo | Python puro |
| Generación de `.pre` | ✅ completo | Python puro |
| Generación de `.fuf` | ✅ completo | Python puro |
| Gráficos diagnósticos | ✅ completo | Matplotlib |
| Gráfico de pronóstico | ✅ completo | Matplotlib |
| Informe HTML de pronóstico | ✅ completo | Jinja2 + SVG |
| Extensión C cffi | ✅ opcional | GSL + cffi |
| CLI `fue` | ✅ completo | Python puro |
| CLI `fuf` | ✅ completo | Python puro |
| Tests | ✅ 629 passing | — |

---

## FASE 1 — Núcleo C: extracción y API ✅

### 1.1 fue_api.c
- [x] `populate_globals()` — FueModelSpec → Tm + Ts + DataMat
- [x] `count_npar_build_par()` — contar parámetros libres, extraer par[] iniciales
- [x] `cast_us()`, `unscramble()`, `CalcNonsOp()`, `calcnu()`
- [x] `fue_estimate()` — entry point completo

### 1.2 Compilación y test
- [x] Extensión cffi: `python src/fue/_build_cffi.py`
- [x] Tests de smoke: test_api.py
- [x] Bugs corregidos:
  - `_build_cffi.py`: path incorrecto; cdef con macros no expandidas
  - `_engine.py`: dangling pointer (numpy array GC'd antes de `fue_estimate`)
  - `fue_api.c`: outputv=NULL segfault → fopen(/dev/null); `printf` siempre a stdout
  - AR/MA factor arrays no inicializados a -1.0 (requerido por `unscramble`)

### 1.3 Tests de equivalencia numérica
- [x] 4 casos: AR(1), IMA(1,1), SFNY.2, RIPC.1 mensual (14 params)
- [x] Tipos de intervención: cos/sin/alter; parámetros ar_free/ma_free; refactor
- [x] Bug: xitol sign en `fue_api.c` (ML exacto vs aproximado)

---

## FASE 2 — API Python iterativa ✅

- [x] `TimeSeries.plot()`, `plot_acf()`, `plot_pacf()`
- [x] `Model.compare(*others)` — tabla AIC/BIC/loglik
- [x] `Model.forecast(horizon)` — mirrors usfo.c/fuf.c
- [x] `TimeSeries.from_pandas(series)` — inferencia de freq/start

---

## FASE 3 — Distribución ✅

- [x] conda-recipe/meta.yaml con compiler('c'), gsl como host dep
- [x] pyproject.toml: build backend, cffi fuera de runtime deps, extra `c-engine`
- [x] cibuildwheel Linux (manylinux_2_28) + macOS (x86_64 + arm64) + Windows (vcpkg)
- [x] Job `build_pure_wheel` (FUE_SKIP_C=1) → `fue-*-py3-none-any.whl`
- [x] MANIFEST.in incluye csrc/**; cadena sdist→wheel verificada
- [x] README.md + COPYING (GPL-2.0) añadidos al repositorio
- [x] pyproject.toml: license, authors, keywords, classifiers; setuptools<69 + license-files=[] → Metadata 2.1 (sin rechazo PyPI)
- [x] datasets.py: `sfny()` y `ripc()` embebidos; quickstart.py workflow completo
- [x] **Publicado en PyPI: https://pypi.org/project/fue/0.1.0/** (2026-06-07)
  - `pip install fue` verificado ✅
  - `pip install "fue[report]"` verificado ✅

---

## FASE 4 — Migración a Python puro ✅

**Objetivo**: eliminar dependencia de GSL y del compilador C.

### 4.1 Evaluación de verosimilitud

**`elfvarma.py`** — dos algoritmos:

- `flikam_scalar` — Mélard (1984) AS 197. Kalman filter con quick recursions.
  Usado en el bucle interno del BFGS (rápido para n grande).
- `elf_scalar` — Mauricio (1997) AS 311 / Mauricio (2002) JTSA.
  Forma de innovaciones de Ansley (1979), adaptación multivariante a m=1.
  Cómputo sin inversión de matrices (Cholesky + sustitución).
  Usado para evaluación final exacta y cálculo de residuos.
- `_cgamma_scalar` — Subroutina CGAMMA de Mauricio (1997) AS 311.
  Versión corregida (Mauricio 1995b) del algoritmo de Kohn & Ansley (1982).

Bugs corregidos:
- Paso (f): `solve_triangular(M, h)` → `M.T @ h` (multiplicación, no resolución)
- Paso (e): no se restaba `mu` de `w` → añadido parámetro `mu=0.0`
- `elf_scalar` g=0 crash (p=q=0, modelo solo con media)
- `elf_scalar` Cholesky failure para AR órdenes grandes (p=26, UK.3): eigenvalor
  mínimo de V₁ΩV₁ᵀ ≈ −3.7e-20 por redondeo; reintento con shift diagonal mínimo

### 4.2 Estimador ML

**`cast_us.py`** — implementa la cadena completa de Mauricio (1995) [JASA95 §3]:

1. `build_est_spec(model)` — pre-computa DataMat, rnsop (equivalente a `populate_globals`)
2. `cast_us_py(x, spec)` — mapea vector de parámetros → (p, q, phi, theta, mu, w)
3. `_estimate_core(model)` — procedimiento de estimación:
   - Objetivo escalado F(x) = Π(x)/Π₀ ∈ (0,1) [JASA95 eq.3.5]
   - Inner loop: `flikam_scalar` (AS197, rápido)
   - Optimizador: `raxopt` (Dennis-Schnabel BFGS)
   - Evaluación final: `elf_scalar` (AS311, exacta)
   - Errores estándar: de la factor Cholesky B del BFGS en convergencia

### 4.3 Optimizador BFGS

**`qnewtopt.py`** — port Python de `qnewtopt.c`:

- `raxopt` — Algorithm A9.4.1 de Dennis & Schnabel (1983) ch.9
- `_lnsrch` — Algorithm A6.3.1 (búsqueda lineal cubica con backtracking)
- `_bfgsfac` — Algorithm A9.4.2 (actualización QR del factor Cholesky)
- `cdgrad` — gradiente por diferencias centrales
- Retorna 6-tupla: `(x, f, B, termcode, niter, gnorm)`
- La matriz de covarianza se calcula como `C[:,i] = 2·f·cholsol(B, e_i)/n`

### 4.4 Precisión numérica

Tests de fiabilidad (`tests/test_reliability*.py`): 282 tests en 4 baterías:

| Batería | Tests | Cobertura |
|---------|-------|-----------|
| reliability.py | 52 | Modelos básicos ARMA, equivalencia C/Python |
| reliability2.py | 37 | Modelos con intervenciones, boxcox |
| reliability3.py | 35 | Modelos estacionales, fixed-freq |
| reliability4.py | 39 | Casos reales completos, parser .inp |

Tolerancias: loglik < 1e-4, params < 1e-4 (C vs Python).

---

## FASE 5 — Pronóstico (fuf) ✅

**`forecast.py`** — mirrors `usfo.c` / `fuf.c`:

Bugs corregidos:
- `_unscramble()`: `new_p[0] = -1.0` en los dos bucles internos (factores AR regulares
  y estacionales) doblaba todos los coeficientes AR cuando el modelo combina AR×SAR;
  el término líder del polinomio se contaba dos veces al convolucionarse



- `ForecastResult` — nivel, diff1, diff estacional + desviaciones estándar
- `forecast(model, result, horizon)` — L pasos adelante en niveles originales
- `eval_at_params(model)` — evaluación sin reoptimizar (workflow fuf)
- `Model.forecast_fuf(horizon, sigma2)` — usa parámetros fijos del `.fuf`

**`model.py`** — `Model.write_fuf()` y `Model.forecast_fuf()`

**`report.py`** — `write_fuf_out()`: informe de pronóstico formato fuf

**`cli.py`** — CLI `fue`:
- `fue model [eml|aml] [chk|nochk] [-f [horizon]]`
- Con `-f`: genera `forecast_model.inp` (fuf format) con sección "Forecast horizon/sigma2"
- Cabecera del `.out` incluye `Input file` y `Output file`

**`fuf_cli.py`** — CLI `fuf`:
- `fuf forecast_model` — lee `forecast_model.inp`, escribe `forecast_model.out`
- No reestima; usa parámetros y sigma2 del fichero fuf

---

## FASE 6 — Informes `.out` / `.pre` ✅

**`report.py`** (~1600 líneas) genera informes ASCII que coinciden con `fue.c` al byte:

### Secciones del `.out`
1. Cabecera: ficheros inp/out, método, observaciones, convergencia, iteraciones, norma gradiente
2. Parámetros omega/delta por intervención con errores estándar
3. Operadores AR/MA (polinomios desarrollados `phi[k]` / `theta[k]`)
4. BoxCox: coeficiente y jacobiano
5. Sigma: σ, σ², log-verosimilitud, AIC, BIC
6. Matriz de correlación de parámetros (si npar > 1)
7. Estadísticos de residuos: media, mín, máx, error estándar, percentiles
8. Tabla de outliers (residuos > 2σ) con fecha y z-score
9. ACF de residuos con calibración (port exacto de `PlotCalibACF` en `diagnose.c`)
10. Gráfico ASCII de residuos (port de `PlotAsciiSer`)

### Secciones del `.pre`
Parámetros estimados como valores iniciales para reiniciar la estimación.

### Tests de regresión
- `test_write_out_ripc1`: compara `.out` generado con referencia C byte a byte
- `test_write_pre_ripc1`: idem para `.pre`

---

## FASE 7 — Gráficos diagnósticos ✅

**`plots.py`** — `plot_model_diagnostics(model)`:

Dos figuras con proporciones exactas del gnuplot de `fue.c`:

**Figura 1** (3 paneles: residuos + ACF/PACF apilados):
- Título `A.<stem>` centrado con `fig.suptitle()` (stem del fichero `.inp`)
- Residuos estandarizados en el tiempo
- ACF de residuos con banda 95% (título `acf`, igual que C)
- PACF de residuos con banda 95% (título `pacf`, igual que C)

**Figura 2** (histograma):
- Histograma de residuos con curva normal superpuesta
- Mismo título `A.<stem>`

`model._inp_stem` se guarda al cargar con `fue.load()` y `fue.load_fuf()`.

Función auxiliar `_snap_series_max` / `_snap_cmax` replica la elección de escala de gnuplot.

`plot_forecast(model, fr)` — figura con historia + pronóstico + bandas:

**Figura 1** (dos paneles):
- Panel superior: variación estacional histórica (linespoints) + pronóstico + bandas dashed; separador vertical
- Panel inferior: ERR residuos como impulsos con ±2σ; rango y snap igual que C (`prevcmax`)
- X-ticks: etiquetas de año a intervalos `freq` (12→anual, 4→bienal)
- Generado automáticamente por `fuf_cli.py` como `<base>_forecast.png`

---

## FASE 9 — Informe HTML de pronóstico SPS ✅

**`report_forecast.py`** — `write_forecast_report(model, fr, path, title=None, source=None, sps_name=None, narrative=None, pdf=False)`:

Genera un fichero `.html` auto-contenido con CSS y SVG embebidos.

### Diseño SPS (Sistema de Previsión y Seguimiento)

**Layout**: dos columnas — tabla izquierda, gráficos derecha.

**Columna izquierda — tabla**:
- `freq+1` filas históricas (ciclo completo + período actual)
- `freq` filas de previsión (un año adelante), fondo azul claro
- Fila en blanco (sin sombrear) + fila `H=horizonte` separadas
- Columnas: DATE | LEVEL (Value, Std%) | Monthly % (Std%) | Annual % (Std%) | ERR
- Separador (borde grueso) entre histórico y previsiones
- Detalles del modelo (σ², AIC, BIC, npar, muestra) en sección colapsable

**Columna derecha — gráfico único SVG (dos paneles, GridSpec)**:
- Panel superior: variación anual — linespoints (línea entrecortada + círculos) para toda la serie (histórico+previsión), círculos mayores para histórico; bandas ±1σ en línea discontinua; separador vertical en origen de previsión; grid en años
- Panel inferior: ERR — impulsos históricos; bandas ±2σ y línea cero terminan en el origen de previsión (spine truncado, `hlines` con xmax); grid en años históricos
- Ambos paneles comparten misma figura → ejes x alineados pixel a pixel
- Nota al pie: "Forecast bands ±1σ · ERR bands ±2σ"

**Parámetros de presentación**:
- `title`: título descriptivo de la serie (no del modelo), e.g. "Spain CPI Inflation"
- `source`: fuente de datos, e.g. "INE"
- `sps_name`: etiqueta SPS, e.g. "Spain Inflation"

**CLI**: `fuf --title TEXT --source TEXT --sps TEXT forecast_model` genera `<base>.html` automáticamente.  
**PDF**: `pdf=True` + `pip install "fue[pdf]"` vía weasyprint.  
**Narrativa**: parámetro `narrative` (HTML) para texto generado por LLM.

### Caso de referencia: España IPC
- Fichero: `/Inflation Volatility/Analisis/Spain/forecast_b2025/forecast_S.2.inp`
- Informe: `forecast_S.2.html`
- Documentación interna SPS: `Spain_S2.md`

---

## FASE 8 — Documentación y licencia ✅

### Literatura de referencia (`../literature/`)
Cuatro papers que cubren toda la implementación:

| Fichero | Referencia | Implementa |
|---------|-----------|-----------|
| `as197.pdf` | Mélard (1984) AS 197 | `flikam_scalar` |
| `518-2013-11-11-JAM197.pdf` | Mauricio (1997) AS 311 | `elf_scalar`, `_cgamma_scalar` |
| `9316.pdf` | Mauricio (1995) JASA/WP | `_estimate_core`, `raxopt` |
| `518-2013-11-11-JAM102.pdf` | Mauricio (2002) JTSA | `elf_scalar` (fórmulas compactas) |

### Documentación del código
- `elfvarma.py`: pasos (a)-(k) anotados con ecuaciones de AS311/JTSA02
- `cast_us.py`: procedimiento de Mauricio (1995) §3 documentado
- `qnewtopt.py`: algoritmos D&S83 A9.4.1, A6.3.1, A9.4.2 referenciados

### Licencia
- `pyproject.toml`: `license = { text = "GPL-2.0-or-later" }`
- Los tres módulos matemáticos llevan cabecera GPL con copyright de Mauricio, Treadway y Guerrero
- **Acuerdo**: Mauricio, Treadway y Guerrero liberan el código bajo GPL-2.0-or-later

---

## Pendiente

### Alta prioridad
- [x] **Francia**: ifault=6 resuelto (era consecuencia del bug `_unscramble`); F.3.inp actualizado a origin=12/2025 (INSEE serie 001759970, base 2015); añadida al SPS como séptimo país

### Media prioridad
- [ ] **conda recipe**: actualizar para builds sin extensión C (`FUE_SKIP_C=1`)
- [ ] **`pyproject.toml`**: marcar cffi/GSL como opcionales en las dependencias de build
- [ ] **`Model.write_out` sin ajuste**: actualmente requiere `.fit()` previo;
      considerar modo "evaluar en parámetros iniciales"

### Baja prioridad
- [ ] **Test `test_write_out_ripc1`**: actualmente `@requires_c`;
      generar una referencia Python equivalente para test sin C
- [ ] **Documentación de API**: docstrings Sphinx / mkdocs
- [ ] **Notebook actualizado**: demo con `estimate_py` y gráficos

---

## Suite de tests (2026-06-06)

**Total: 629 passing**

| Fichero | Tests | Qué cubre |
|---------|-------|-----------|
| test_api.py | 36 | Smoke tests: series, modelo, carga .inp, estimación |
| test_cast_us.py | 13 | cast_us_py, build_est_spec, calcnu_py |
| test_elfvarma.py | 11 | flikam_scalar, elf_scalar, _cgamma_scalar |
| test_estimation.py | 18 | Equivalencia C/Python: params, loglik, residuos |
| test_forecast.py | 13 | ForecastResult, fuf workflow, write_fuf_out |
| test_performance.py | 257 | Benchmarks C vs Python |
| test_qnewtopt.py | 29 | raxopt, cdgrad, _lnsrch, _bfgsfac |
| test_real_cases.py | 58 | Regresión .out/.pre byte-a-byte |
| test_reliability.py | 67 | Fiabilidad numérica: modelos básicos |
| test_reliability2.py | 37 | Intervenciones, BoxCox |
| test_reliability3.py | 44 | Estacionales, fixed-freq |
| test_reliability4.py | 46 | Casos reales, parser .inp |
