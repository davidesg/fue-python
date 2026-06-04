# FUE Python — Tareas pendientes

Estado: en progreso  
Referencia: fue-1.13.1 es el código C fuente de verdad.

---

## FASE 1 — Núcleo C: extracción y API (en progreso)

### 1.1 fue_api.c — implementación completa
- [x] `populate_globals()` — FueModelSpec → Tm + Ts + DataMat
  - [x] Generación de DataMat[0] (serie transformada, BoxCox)
  - [x] Generación de DataMat[i] para cada intervención (impulse, step, ramp, seasonal)
  - [x] Asignación de omega/delta/flags en Tusmodel (Tm)
  - [x] Asignación de factores AR/MA regulares y estacionales
- [x] `count_npar_build_par()` — contar parámetros libres y extraer par[] iniciales
- [x] `cast_us()` — extraído de fue.c:3645
- [x] `unscramble()` — extraído de fue.c:3939
- [x] `CalcNonsOp()` — extraído de fue.c:4303
- [x] `calcnu()` — extraído de fue.c:4489
- [x] `fue_estimate()` — implementación completa:
  1. populate_globals(spec)
  2. count_npar() → npar
  3. build_par_vector() → par[]
  4. Allocar dev[], cov[][], a[][]
  5. Llamar est(cast_us, npar, par, dev, cov, ...)
  6. Empaquetar FueResult y return

### 1.2 Compilación y test
- [x] Verificar que compila con gcc -c fue_api.c (sin gnuplot, sin pdflatex)
- [x] Compilar extensión cffi: python src/fue/_build_cffi.py
- [x] Activar cffi_modules en pyproject.toml
- [x] Test de smoke: 10/10 tests passing (test_api.py)
- [x] Bugs encontrados y corregidos:
  - _build_cffi.py: path incorrecto (3 niveles en lugar de 2)
  - _build_cffi.py: cdef con macros no expandidas → cdef manual con literales
  - _engine.py: dangling pointer (numpy array GC'd antes de fue_estimate) → ffi.from_buffer
  - fue_api.c: outputv=NULL causaba segfault en qnewtopt.c:report() → fopen(/dev/null)
  - qnewtopt.c: printf() siempre a stdout → if(outputv) fprintf(outputv,...)

### 1.3 Tests de equivalencia numérica
- [x] Crear tests/test_estimation.py con casos de prueba de fue-1.13.1
- [x] Cada test: mismo .inp → fue C vs fue Python → comparar params, loglik, residuos
- [x] Tolerancia: params dentro de 1e-6, loglik dentro de 1e-4
- [x] Bugs encontrados y corregidos:
  - test_ar1.inp tenía secciones extra (seasonal AR/MA, anual f-fixed) que el parser no lee → reescrito con estructura exacta de SFNY.2.inp
  - fue_api.c: xitol sign bug → fue.c usa -xitol para exact ML (fuerza recursiones Melard exactas en cxi()), fue_api.c usaba +xitol (ML aproximado) → corregido con eml flag
  - fue_api.c: tolerancias del optimizer no coincidían con fue.c → grtol/sptol = pow(DBL_EPSILON, ...), maxits=500
- [x] Casos de prueba:
  - Case 1: AR(1) anual, 30 obs SFNY30 → phi=0.9747519833, sigma2=0.2482607765, logelf=-23.1683049163
  - Case 2: IMA(1,1) anual, 30 obs SFNY30 → theta=-0.4228241648, sigma2=0.2060862798, logelf=-18.3455244692
  - Case 3: SFNY.2 completo (step + AR(1)×AR(2) + mu, boxlam=0) → 6 params, logelf=13.9573576937, sigma2=0.0370593261
- [x] Bug adicional corregido: fue_api.c: Ar1/Ar2/Ma1/Ma2[i][0] no se inicializaba a -1.0 como requiere unscramble() → con modelos de varios factores, los valores individuales por factor eran incorrectos aunque el polinomio combinado (y el loglik) eran correctos
- [x] Case 4: RIPC.1 mensual (cos/sin/alter + step + AR(1) fijo + mu, boxlam=0, refactor=100) → 14 params, logelf=-100.9274828448, sigma2=0.9662469111
  - [x] Nuevos tipos de intervención: FUE_ITV_COS, FUE_ITV_SIN, FUE_ITV_ALTER con campo `harmonic` en FueIntervention
  - [x] Soporte para coeficientes AR/MA fijos (coef_free=0) en Model: parámetros ar_free/ma_free/ar_s_free/ma_s_free
  - [x] Parámetro refactor en Model (default=1.0)

---

## FASE 2 — API Python iterativa

- [x] `TimeSeries.plot()` — gráfico con eje x en años decimales (plot_series mejorado)
- [x] `TimeSeries.plot_acf()` y `plot_pacf()` como métodos de instancia
- [x] `Model.compare(*others)` — tabla AIC/BIC/loglik comparativa
- [x] `Model.forecast(horizon)` — implementación Python pura en forecast.py (mirrors usfo.c/fuf.c)
- [x] `TimeSeries.from_pandas(series)` — inferencia automática de freq y start desde DatetimeIndex/PeriodIndex

---

## FASE 3 — Distribución

- [x] conda-forge recipe (entorno Linux prioritario)
      - `conda-recipe/meta.yaml`: recipe completa con compiler('c'), gsl como host dep, pip install
      - `_build_cffi.py`: detecta PREFIX (conda Linux/macOS) y LIBRARY_INC/LIBRARY_LIB (conda Windows)
- [x] pyproject.toml: build backend corregido (setuptools.build_meta + setup.py con cffi_modules)
      - `setuptools.backends.legacy:build` no disponible en setuptools del sistema → cambiado a `build_meta`
      - `cffi_modules` en `[tool.setuptools]` rechazado → movido a `setup.py`
      - `_build_cffi.py`: rutas con `../../` confundían distutils → `os.path.relpath` desde raíz del proyecto
      - `python setup.py bdist_wheel` genera `fue-0.1.0-cp312-cp312-linux_x86_64.whl` con `_fue_engine.abi3.so`
- [x] cibuildwheel para Linux/macOS — `.github/workflows/wheels.yml` + `[tool.cibuildwheel]` en pyproject.toml
      - Linux (manylinux_2_28): `dnf install -y gsl-devel` + auditwheel bundle libgsl.so
      - macOS: `brew install gsl` + delocate bundle; x86_64 + arm64
      - Python 3.10–3.13; sdist + wheels en el mismo job; publish a PyPI en tags v*
      - MANIFEST.in incluye csrc/**; sdist→wheel chain verificada localmente
- [x] Soporte Windows (AMD64, GSL estática vía vcpkg x64-windows-static-md)
      - fue_api.c: /dev/null → NUL en Windows (#ifdef _WIN32)
      - _build_cffi.py: flags MSVC (/O2 /W2), sin -lm, rutas GSL desde GSL_ROOT o vcpkg default
      - pyproject.toml: [tool.cibuildwheel.windows] con before-all=vcpkg install + environment GSL_ROOT
      - wheels.yml: windows-latest añadido a la matriz
- [x] Notebook de ejemplo: `notebooks/fue_example.ipynb`
      - AR(1) simple, SFNY.2 completo (step + AR(1)×AR(2) + mu + boxlam), RIPC mensual (cos/sin/alter)

---

## FASE 4 — Migración a Python puro (completada)

**Objetivo**: reemplazar la extensión cffi (`_fue_engine`) con una implementación
Python pura para eliminar la dependencia de GSL y el compilador C.

**Estado actual**:
- `elfvarma.py`: Python puro — `elf_scalar` y `flikam_scalar` ✓
- `cast_us.py`: Python puro — `calcnu_py`, `cast_us_py`, `estimate_py` ✓
- `_engine.estimate()`: usa C si disponible, cae a `estimate_py` si no ✓
- 140 tests pasan ✓

**Bugs corregidos en esta fase**:
- `elf_scalar` paso [6.2]: `solve_triangular(M, h)` (inversa) → `M.T @ h` (producto)
  — la C hacía multiplicación por Mᵀ, no resolución del sistema
- `elf_scalar`: no restaba `mu` de `w` en paso [5.1] → añadido parámetro `mu=0.0`
- `calcnu_py` test: expectativa incorrecta (convenio de signo ω(B) = ω₀−ω₁B−…)
- `test_sfny2`: modelo incorrecto (d=1, ar=[[0.9],[0.6]]) → modelo SFNY.2 real
- `estimate_py` Hessiano: `np.where(diag>0, sqrt(diag), 0)` → `sqrt(maximum(diag,0))`

### 4.1 cast_us_py — Python puro ✓
- [x] `calcnu_py()` — port de calcnu() en fue_api.c:4489
- [x] `cast_us_py(x, spec)` — ensambla `(phi, theta, mu, w)` desde vector de parámetros
- [x] Tests de `cast_us_py`: calcnu, build_est_spec, w correcto para AR(1) y d=1

### 4.2 optimizer_py — Python puro ✓
- [x] `estimate_py(model)` en `cast_us.py` (scipy L-BFGS-B + flikam_scalar + elf_scalar)
- [x] Tests de equivalencia: loglik 1e-3, params 1e-4/5e-4 para 3 casos de prueba

### 4.3 Limpieza final ✓
- [x] `_engine.py`: fallback a `estimate_py` cuando C no disponible
- [x] `model.fit()`: mensaje de error sin requerir C extension
- [ ] Actualizar `pyproject.toml`: marcar cffi/GSL como opcionales
- [ ] Actualizar conda recipe y cibuildwheel para builds sin C

---

## Notas de sesión

### 2026-05-26
- Creada estructura del proyecto en ../fue/
- fue_api.h completo (interfaz pública cffi)
- fue_api.c: stub, cast_us() pendiente de extraer
- Python package skeleton completo (series, intervention, model, diagnostics, plots)
- tests/test_api.py: smoke tests sin extensión C (todos deben pasar)
- Leído fue.c: cast_us en línea 3645, unscramble en 3939, CalcNonsOp en 4303,
  calcnu en 4489, BoxCox en 4511
- Globales en fue.c: Tm (l.38), Ts (l.39), DataMat (l.40) — replicar como static en fue_api.c
- DataMat construction: fue.c:280-436 (impulse/step/ramp/seasonal/easter/trend/cos/sin/alter)
- npar counting: fue.c:921-1048 (segunda pasada, post-lectura)
- SIGUIENTE: implementar fue_api.c completo (Fase 1.1)
