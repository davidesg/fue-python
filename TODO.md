# FUE Python — TODO

## PARA 0.3 — la puerta de entrada de `fuf` (2026-09-12)

- [ ] **`fuf` debería aceptar un `.pre` + horizonte, sin fichero intermedio.**
      Decisión del analista, 2026-09-12, después de ver el coste en el uso real
      —propio y de colegas—. Anotado también en `art-python/TODO.md`: el cambio
      toca a los dos.

      **El estado de hoy.** La previsión tiene su propio trío, paralelo al de
      `fue`:

        fue <modelo> -f <H>    →  forecast_<modelo>.inp    el fuf
        fuf <forecast_modelo>  →  forecast_<modelo>.out + _forecast.png + .html

      Y el fuf **es un `.inp` más una sección** —`** Forecast horizon and
      estimated innovation variance`—, así que lleva extensión `.inp`
      (`src/fue/inp.py:64-66`). Lo único que lo distingue por el nombre es el
      prefijo `forecast_`, que `load_fuf` quita para recuperar `_inp_stem`
      (`inp.py:72-76`): es convenio, no formato.

      **Por qué molesta.** En una carpeta de trabajo conviven dos clases de
      `.inp` que sólo se distinguen abriéndolas. Un humano lo aprende una vez;
      un asistente lo paga cada vez, y el analista lo ha medido como «un coste
      altísimo de tokens». Un convenio que hay que deducir leyendo el contenido
      es un acertijo.

      **La forma propuesta: entrar por el `.pre`.** Es como entra `drtran`
      (`load_pre`) y encaja con el resto de la escalera. El `.pre` ya es el
      óptimo en forma reejecutable, que es justo lo que `forecast_fuf`
      necesita: parámetros fijos, sin reestimar. Lo que el `.pre` no lleva
      —horizonte y σ²— son cosas de la LLAMADA y no del modelo.

      **Lo que hay que conservar, y es lo fácil de perder.** El fuf guarda σ²
      DENTRO. Eso es lo que hace comparables dos previsiones del mismo modelo
      hechas en momentos distintos. Si σ² pasa a recalcularse en cada llamada,
      la propiedad se va en silencio. Si se entra por `.pre`, σ² tiene que
      venir de algún sitio declarado —el `.out`, o un argumento explícito— y no
      por defecto.

      **Compatibilidad.** `load_fuf` debe seguir leyendo los fuf existentes:
      hay ejercicios publicados que dependen de ellos
      (`SF_MEG/empirical/sps/forecast_compare.py` exporta con `fue -f 77`).

- [ ] **`chkma` no está en el formato.** Un modelo estimado sin restricción de
      invertibilidad (`fue -e`) se relee como restringido: ni el `.inp` ni el
      `.pre` guardan ese bit, y `Model.__init__` lo pone a `True`. Hoy eso
      SALVA a todo el mundo —es lo que evita prever con la raíz no invertible
      del testigo MA_f, el error nº 1 documentado en
      `SF_MEG/empirical/FORECAST_COMPARISON.md`— pero es una propiedad que se
      sostiene porque el valor por defecto coincide con lo que hace falta, no
      porque el fichero la declare.

## LO QUE FALTA PARA UNA VERSIÓN ESTABLE (2026-09-06)

El motor está maduro y los números lo respaldan: **5.215 estimaciones** en el
ecosistema, **1.579** desde julio, **656 iteraciones** registradas en guiones, y
la suite de `art` en verde con 89 de sus 126 ficheros de prueba ejercitando el
motor.

Y sobre todo el ritmo de aparición de defectos:

    18 jul   BUG-0001 0002 0003
    22 jul   BUG-0004 0005
    30 jul   BUG-0006 0007 0008 0009
    12 ago   BUG-0010 0011 0012
    14 ago   BUG-0013
     6 sep   BUG-0014      ← y éste NO salió del uso: salió de ir a buscarlo

**13 defectos en 27 días, y luego 23 días sin ninguno**, cubriendo el período de
uso más intenso (la réplica del TFM y los casos UEM). Eso es lo que sostiene la
impresión de que el paquete está estable.

Faltan **tres puertas**, y ninguna es de código nuevo. Tras la revisión del
7 de septiembre, **la que parecía bloqueante ya no lo es** — ver abajo.

### 1 · BUG-0005 — YA NO ES EL QUE DECIDE (revisado 2026-09-07)

Se volvió a probar, y el resultado encoge la puerta hasta casi cerrarla.

**El silencio, que era la mitad grave, está cerrado** (0.1.14): un μ̂ absurdo ya
no se reporta como `converged=True`. Umbral medido sobre 1.570 modelos —máximo
legítimo 1,4 sd, el espurio 47,2— con 0 falsos positivos en el corpus entero.

**Y la causa raíz del informe resultó estar mal.** Nueve arranques cubriendo la
región y los dos signos llegan al MISMO óptimo, y el punto de Windows no es
siquiera estacionario (derivada +19,57). No hay evidencia de multimodalidad en
ese caso: lo que había era una semilla con el signo cambiado —la convención de
Box y Jenkins— y eso está arreglado en el consumidor.

**Consecuencia para la estable:** el multi-arranque que el informe proponía no
hace falta, y no hay fragilidad del optimizador que documentar. Lo único que
queda es **verificar en Windows** una observación de julio que el resto de la
evidencia no acompaña. Eso es una limitación conocida y documentable, no un
defecto oculto — y una estable puede salir con ello dicho en una línea:

> *en superficies multimodales el óptimo puede depender de la plataforma; `fue`
> avisa cuando el resultado es implausible, y la recomendación es reestimar
> desde varias semillas.*

### 2 · BUG-0015 — los errores típicos

No es «arreglarlo»: es **decidir**. El arreglo está identificado —descomentar
`fdhess` en `drvmlest.c:112`— y `drtran` demuestra que funciona. No se aplica
porque cambia un error **ruidoso y detectable** por uno **limpio e
indetectable**: fuera del óptimo el hessiano por diferencias finitas puede no ser
definido positivo, y `choldcp` parchearía los pivotes publicando números de
aspecto impecable.

Hace falta una sesión propia: barrido empírico sobre la batería y decisión de
qué hacer cuando el hessiano no sea definido positivo —rechazar, avisar, o caer
a la matriz del BFGS diciéndolo—.

Una versión estable puede salir con esto **abierto y documentado**; lo que no
puede es salir con ello sin decirlo. Por eso se levantó el informe.

### 3 · BUG-0013 — el AR(1) fijado en cero

El puente de Python desvía y `art` ya escribe siempre `1 1 / 0.0 0`, así que los
ficheros nuevos están cubiertos. Quedan:

  - **el C**, que es el arreglo de fondo: la escritura fuera de rango con
    `p=q=0` sigue sin buscarse. No es trivial;
  - **178 `.inp` ya escritos** que matan al binario, entre ellos ITCER, PGAS y
    RATIO del TFM en curso;
  - **las wheels** sin el desvío del puente, que se comportan como el binario.

### Y una señal de cobertura que conviene no ignorar

BUG-0014 mostró que **tres de once tipos deterministas daban regresor nulo en la
previsión y ninguna prueba lo cazó**. No es alarma —el defecto era de una ruta
menos usada que la de estimación— pero dice dónde mirar antes de sellar: la
previsión está menos ejercitada que la estimación.

---

# FUE Python — Estado de la migración C → Python

Referencia: `fue-1.13.1` es el código C fuente de verdad.  
Última actualización: 2026-06-07 (rev 9)

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
- [x] **Publicado en PyPI: https://pypi.org/project/fue/0.1.1/** (2026-06-07)
  - 0.1.0: publicación inicial
  - 0.1.1: David E. Guerrero como autor de contacto
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

### Los valores publicados de Box & Jenkins — falta la tabla (2026-08-13)

**Lo único fijado contra el libro en toda la suite es θ = 0.70** de la Serie A,
IMA(0,1,1) (`tests/test_published_benchmark_series_a.py:43`). El resto de las
nueve especificaciones del banco no está contrastado contra el libro **porque no
tenemos la tabla**: el archivo que acompaña a las series
(`Box_y_Jenkings/index.html`) sólo trae la descripción de los datos.

Lo que sí está medido, y es lo que hace que valga la pena cerrarlo: TASTE es el
algoritmo de Jenkins —mínimos cuadrados no lineales con retroprevisión— y por
tanto **da el estimador del libro**. Sus valores sobre el banco, ya obtenidos
(ver `Taste/oracle/bugs/BUG-0001`):

| caso | modelo | TASTE (el estimador del libro) | fue (ML exacta) |
|---|---|---|---|
| a1 | ARMA(1,1)+media | 0.91368  0.58034  17.07284 | 0.90868  0.57584  17.06528 |
| a2 | IMA(0,1,1) | 0.70503 | 0.69938 |
| b | IMA(0,1,1) | −0.08659 | −0.08636 |
| c | ARI(1,1,0) | 0.82382 | 0.82016 |
| c2 | IMA(0,2,2) | 0.12606  0.12104 | 0.12501  0.11939 |
| d | AR(1)+media | 0.87043  9.12052 | 0.86858  9.10876 |
| d1 | IMA(0,1,1) | 0.05913 | 0.05891 |
| e1 | AR(2)+media | 1.42175  −0.72622  47.91280 | 1.40757  −0.71281  48.19126 |
| e2 | AR(3)+media | 1.56876  −1.02205  0.21254  48.11996 | 1.55312  −1.00175  0.20634  48.44344 |

**Qué hace falta**: la tabla 7.13 de Box & Jenkins (1976). Con ella se fijan los
ocho valores restantes como prueba, al lado del 0.70 que ya está.

⚠ Y hay que decir con qué tolerancia: el libro estima por MCNL y `fue` por ML
exacta, así que lo razonable son **dos decimales**, no siete. Un test que exija
más estaría comparando estimadores distintos como si fueran el mismo — el mismo
error que ya se cometió una vez con TASTE y el modelo airline
(`tests/test_taste_nls_criterion.py`).

### Bugs pendientes

#### ~~`csrc/internal/nlatools.c` — `tensor()` crash con nrl < 0~~ ✅ CORREGIDO (2026-06-15)

**Síntoma resuelto**: segfault / "double free or corruption (out)" al estimar cualquier
modelo con AR (regular o AR_f) **más** MA_f.

**Causa**: `elfvarma.c` llamaba `tensor(-q+1, 0, 1, m, 1, m)` para alojar `gamwa`.
Con q ≥ 2 (MA_f aporta orden 2) resultaba `nrl = -1`.  La asignación `calloc(nrh+1, ...)`
solo reservaba 1 slot; el bucle de inicialización escribía en `t[-1]` → corrupción del heap.

**Nota histórica**: el binario standalone fue-1.13.1 no sufría este crash porque las
secciones Ar2f/Ma2f de `cast_us()` en `fue.c` están comentadas — los factores de
frecuencia fija se leen pero no se incluyen en phi/theta, por lo que q nunca supera
al orden del MA regular.  El híbrido Python (`fue_api.c`) sí expande Ma1f en theta
(`q1 += 2`), lo que activaba el bug.

**Fix aplicado** en `csrc/internal/nlatools.c`:
```c
/* tensor(): asignar nrh-nrl+1 slots y desplazar el puntero */
t = (double ***)calloc( (size_t)(nrh - nrl + 1), sizeof(double **) );
t -= nrl;   /* shift so t[nrl..nrh] are valid when nrl < 0 */

/* free_tensor(): deshacer el desplazamiento antes de free() */
if ( t ) { free( t[nrl][ncl] ); free( t[nrl] ); free( t + nrl ); }
```

**Impacto del fix**: MEG (AR + MA_f testigo) ahora usa el backend C → 0.13s/frecuencia
(vs ~215s con el estimador Python puro).  `art.formal_tests.dcd_f()` y `meg()`
actualizados para usar `model.fit()` en lugar del workaround `_fit_py()`.

---

#### Backend C — crash con AR_s + MA_s simultáneos (modelos SARIMA híbridos)

**Síntoma**: segfault / "Aborted (core dumped)" al estimar cualquier modelo que
combine operadores estacionales AR (ar_s, P≥1) **y** MA (ma_s, Q≥1) al mismo
tiempo mediante el backend C.

**Ejemplo mínimo que reproduce el crash**:
```python
import fue
ts, _ = fue.inp.load("R.1.inp")   # serie trimestral, s=4
m = fue.Model(ts, d=1, D=1, boxlam=0.0,
              ma=[[-0.3]], ma_free=[[True]],
              ar_s=[[0.0]], ar_s_free=[[True]],
              ma_s=[[-0.3]], ma_s_free=[[True]])
m.fit()   # → Aborted (core dumped)
```

**Causa probable**: misma raíz que el bug AR+MA_f anterior — el tensor `gamwa`
en `elfvarma.c` se asigna con tamaño incorrecto cuando se combina el orden
del factor AR estacional (s·P) con el MA estacional (s·Q), generando `nrl < 0`
en la llamada a `tensor()`.

**Impacto**: modelos SARIMA(p,d,q)(P,D,Q) con P≥1 y Q≥1 no se pueden estimar
con el backend C.  El estimador Python puro (`estimate_py`) no se ve afectado.

**Workaround activo**: `art.mcp_server._make_model` y `_build_inp` admiten
parámetros `P, Q`, pero los tests de art solo cubren P=0,Q=1 (o P=1,Q=0).
Hasta que se corrija, usar solo uno de los dos (AR_s OR MA_s, no ambos).

**Fix pendiente**: mismo que el bug tensor() anterior — corregir la asignación
en `nlatools.c:tensor()`.  Verificar que la corrección también elimina este crash.

---

#### Backend C — crash con p=0, q=0 (sin parámetros ARMA)

**Síntoma**: `Segmentation fault` al llamar `.fit()` en un modelo sin parámetros ARMA
(solo armónicos/intervenciones o modelo puro de ruido blanco diferenciado).

**Ejemplo mínimo**:
```python
ts, m = fue.inp.load("serie.inp")  # serie cualquiera
m.ar = []; m.ma = []               # sin ARMA
m.fit()                            # → Segfault
```

**Causa probable**: `elf_scalar` no inicializa correctamente el vector de parámetros
cuando g=0 (sin ARMA libre). El array de parámetros queda vacío o con puntero nulo,
y la función de verosimilitud intenta indexar fuera del rango.

**Impacto**: impide estimar modelos de solo armónicos sin ARMA, lo que es válido
metodológicamente (modelo con solo estacionalidad determinista, paso previo a la
identificación ARMA).

**Workaround aplicado** (2026-06-13):
- **Caso npar=0** (sin armónicos, p=q=0 fijo): `fue/_engine.py` detecta `len(_build_initial_x(model)) == 0` y enruta a `eval_at_params` (camino Python puro) sin llamar al backend C.
- **Caso con armónicos** (p=q=0 pero n_harmonics>0): `art/mcp_server.py` `_build_inp` y `_make_model` añaden `ar=[[0.0]], ar_free=[[False]]` (AR(1) con φ=0 fijo) cuando `p=0 AND q=0`. Esto aporta nar=1 al backend C y evita el crash, sin alterar la verosimilitud (φ fijo en 0 no se estima).

**Fix pendiente en C**: en `elf_scalar`, verificar g>0 antes de acceder al vector de
parámetros; si g=0, devolver solo el valor de verosimilitud del ruido blanco puro.

---

### Alta prioridad
- [x] **Francia**: ifault=6 resuelto (era consecuencia del bug `_unscramble`); F.3.inp actualizado a origin=12/2025 (INSEE serie 001759970, base 2015); añadida al SPS como séptimo país

### Media prioridad — Usar pyfug para gráficos diagnósticos

**Contexto** (jun-2026): fue Python es usado directamente por Claude como herramienta
de estimación en el flujo Box-Jenkins. Claude llama a `fue.Model` + `m.fit()` y luego
usa pyfug para los gráficos diagnósticos de residuos. No hace falta un ART complejo.

- [ ] **Reemplazar `plots.py` con pyfug** — **DECIDIDO 2026-09-25, para la 0.3 de
  la suite** (el paso 2 de BUG-0023; plan en art-python `TODO.md`, «PARA 0.3 — un
  solo dibujo para la figura de un modelo»). Forma decidida: **serie → pyfug;
  modelo → fue**. `plot_model_diagnostics` sigue siendo la única entrada para la
  figura de un modelo, delega en `pyfug.graphics.plot_combined` si pyfug está
  instalado (extra opcional `fue[graficos]`, importado de forma perezosa) y
  conserva `plots.py` como respaldo con los MISMOS números. fue no pasa a
  depender de pyfug. Condiciones previas, en pyfug: que deje de reescribir los
  `rcParams` globales al importarse y que sustituya statsmodels por ACF/PACF/Q
  propias (las de `fue.diagnostics` dan lo mismo a 1e-15), que es lo que le
  impone el tope de retardos n/2 − 1 y 1,4 s de importación.

  Original: sustituir los gráficos matplotlib actuales
  de `plot_model_diagnostics` y `plot_forecast` por los gráficos Jenkins-Treadway
  de `pyfug.graphics`, para coherencia visual con FUG y ART.
  - `plot_combined` → residuos + ACF/PACF en layout combinado
  - `plot_acf_pacf` → ACF/PACF standalone de residuos
  - `plot_histogram` → histograma de residuos con JB y p-valor
  - Referencia pyfug: `/home/david/Dropbox/SRC/atws/fug/pyfug`

  **Patrón de uso actual** (Claude directo, sin plots.py):
  ```python
  m = fue.Model(ts, boxlam=0.0, d=1, ar=[[...]], ar_free=[[True]])
  for at in [81, 83, 218, 219]:
      m = m.add_intervention('step', at=at)
  m.fit()
  resids = np.array(m.residuals.data)  # TimeSeries → ndarray
  # sigma = np.std(resids, ddof=0)
  # df de la Q: NO len(m.params) — eso cuenta armónicos e intervenciones.
  # Usa fue.diagnostics.free_arma_count(m) (BUG-0023).
  ```

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

## Suite de tests (2026-06-07, fue 0.1.1)

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

## Renderer determinista diagnóstico+ecuación (migración desde fue-C) — 2026-07-05

**Contexto (del MCP `art` / estudio SF_MEG).** El MCP presenta el modelo estimado como el
bloque canónico de dos ecuaciones (`art.describe.model_equation`) y le pide al LLM que lo
muestre "TAL CUAL". Eso funciona pero **depende de que el LLM obedezca**: a veces reconstruye
una tabla, cambia precisión o transcribe mal. La solución de fondo es un **renderer
determinista**, y el fue-C ya produce esa salida combinada (instrumentos diagnósticos +
ecuación). **Pendiente: migrar esa salida del C a fue-Python.**

- [ ] Portar del fue-C el **reporte diagnóstico+ecuación** (una salida unificada, no solo el `.out`).
- [ ] Consumir un **`model_equation_data`** estructurado (operadores, coefs, SE, μ, σ², ℓ, AIC/BIC)
      construido con el unpacker canónico `fue.forecast._reconstruct_params` — sin re-parsear texto.
- [ ] En el lado ART: poblar `Description.data` con ese estructurado y transportarlo en `_result`
      (hoy se descarta); entonces el renderer/host presenta sin pasar por el LLM.

Ver `~/Dropbox/SF_MEG/empirical/ART_MCP_REVIEW.md §2` (decisión: §2 diferido a esta migración).

## Covarianza: factor BFGS acumulado vs Hessiano en el óptimo — EVALUAR (no cambiar de inicio) — 2026-07-05

**Síntoma.** Re-estimar un modelo desde un `.pre` **sin modificar** (arranque ~en el óptimo) da
SE **inflados** frente al ajuste original. Ej. FR_CPI (estudio SF_MEG): μ SE = 0.0717 (`.pre`-reload)
vs 0.0156 (build fresco); el Hessiano de **diferencias finitas** en el MISMO óptimo da 0.0145
(correcto, independiente del arranque). Mismo punto estimado, misma verosimilitud.

**Causa.** `cast_us._estimate_core` (y `qnewtopt.c`) derivan la covarianza del **factor Cholesky del
Hessiano BFGS acumulado** (`B_hess`; B parte de I y se refina con cada iteración), no de un Hessiano
fresco en el óptimo. Con **pocas iteraciones** (arranque en el óptimo) B≈I → covarianza cruda. Es
**path-dependiente**: el SE reportado depende de la trayectoria del optimizador, no solo del óptimo.

**Por qué NO es un problema en la práctica.** El `.pre` se recarga para **modificar** el modelo, lo
que aleja el nuevo óptimo del arranque y da iteraciones suficientes para acumular B bien. Solo muerde
al **re-estimar idéntico** (tools de display: `model_equation_display`, `get_out_report`).

**Por qué NO cambiar a la ligera (trade-off real).** El factor BFGS es **siempre definido-positivo**
por construcción (Cholesky mantenido); un Hessiano de diferencias finitas puede salir **indefinido /
mal condicionado cerca de fronteras** (varianzas negativas). Probablemente por eso Mauricio lo
construyó así. Cambiar a `fdhess` afecta a TODOS los usuarios y podría alterar SE publicados.

**Acción — EVALUAR (no fix de inicio):**
- [ ] Revisar el **C original** (`qnewtopt.c`: cómo obtiene la covarianza) y los **papers/documentación
      del motor** (Mauricio) para entender la decisión de diseño.
- [ ] Evaluar con cuidado: covarianza vía `fdhess` en el óptimo (correcta, path-independiente) vs el
      factor BFGS (siempre PD, barato). Posible híbrido: `fdhess` con salvaguarda de PD (fallback a BFGS).
- [ ] Mitigación a nivel de **workflow** (sin tocar el motor): que los tools de display NO re-estimen un
      `.pre` sin cambios; que carguen los SE del ajuste original (`.out`) o que el `.pre` guarde los SE.

---

## Los errores típicos por el hessiano — para la 0.3 (11-sep-2026)

Sale de la corrida de ES_CPI del run 3 y de la pregunta del analista: *«sé que
fue tiene una función que los calcula desde el hessiano; me gustaría saber si se
puede solucionar de raíz, aunque sería tocar fue»*. Sí se puede. Está medido en
`bugs/BUG-0015` (addendum del 11-sep). Lo que falta es esto.

### 1 · El paso del `_fdhess` del puerto está MAL TRADUCIDO

No es una elección discutible: es un error de traducción del C, y explica por qué
esa rama devuelve ceros.

```c
/* qnewtopt.c — como lo escribió Mauricio */
cubreta = pow( eta, 1.0/3.0 );            /* ← la raíz cúbica va DENTRO */
step[i] = cubreta * rmax( x[i], 1.0 );

/* drvmlest.c:112 — y se le llama con macheps */
/* fdhess( objcfunc, npar, par, pi1, macheps, mtmp ); */
```

```python
# cast_us.py — el puerto
dx = eta * np.maximum(np.abs(x), 1.0)     # ← SIN raíz cúbica
...
H = _fdhess(objective, x_opt.copy(), obj_opt, _SQRT_EPS)   # ← y con √ε
```

    C:       ∛(2.2e-16)  =  6.04e-06   = ε^(1/3)     ✓
    Python:     √(2.2e-16) = 1.48e-08   = ε^(1/2)     ✗   407× más pequeño

El puerto **quitó la raíz cúbica y compensó con el exponente equivocado**. La
compensación correcta era pasar `eta = ε^(1/3)` directamente, no `√ε`.

Medido en el óptimo de `ES_CPI_m10` (13 parámetros, cond(H)=138): con ε^(1/2) la
diagonal del hessiano sale con **2 valores negativos** y
`np.sqrt(np.maximum(diag, 0.0))` los convierte en **errores típicos de 0.0**;
con ε^(1/3), ε^(1/4) o ε^(1/5) la diagonal es entera positiva y las SE son
**idénticas en los tres** — y coinciden con el **GLS exacto** de
`fue-1.13.1/ERRORES_ESTANDAR.md` (0.068328, 0.027692).

**Lo que el puerto SÍ mejoró** y conviene conservar: la diagonal la calcula con
diferencia CENTRADA (f₊ − 2f₀ + f₋) donde Mauricio usa adelantada, y la cruzada
con la fórmula de 4 puntos. Son más precisas. El defecto es sólo el paso.

### 2 · Verificado: está IMPLEMENTADA, es de Mauricio, y lleva comentada desde el origen

Se comprobó porque el analista lo preguntó bien: *«fdhess en Mauricio puede estar
comentada y no implementada. Verificar si es así.»* No lo está.

* **`fdhess` tiene cuerpo completo** en `qnewtopt.c` (~54 líneas: reserva
  vectores, toma ∛eta, escala el paso por |x|, llena diagonal y cruzadas, libera).
  No es una declaración vacía ni un esqueleto.
* **`choldcp` también existe** (`nlatools.c:124`) y no como adorno: `elfvarma.c`
  la usa en tres sitios. Las dos se compilan y se enlazan.
* **Está comentada desde el `drv` ORIGINAL de Mauricio** —`drvmlest.c:95`— y
  sigue comentada, igual, en las **siete** versiones de fue revisadas (1.01,
  1.05, 1.09, 1.11, 1.12.03, 1.13, 1.13.1). Nadie la añadió después ni la
  desactivó: nació así.

Y la dejó **lista para usar**, no a medias: la llamada que escribió pasa
`macheps`, que con la ∛ interna da ε^(1/3) = 6,0e-06 — el paso correcto. Lo que
está mal es la traducción al Python (§1), no el original.

> **NOTA DE MÉTODO, porque casi se escribe aquí lo contrario.** El primer censo
> dijo «no aparece en 1.01-1.13; la línea se añadió en 1.13.1». Era **falso**:
> `grep` trataba esos ficheros como BINARIOS —por el byte Latin-1 de «José» en
> la cabecera de copyright— y **suprimía la salida sin avisar**, exit code 1
> como si no hubiera coincidencias. Con `grep -a` aparece en todas. Cualquier
> censo sobre este árbol necesita `-a`.

### El porqué sigue sin saberse

El comentario que la precede es **neutro**: «This is an alternative way of
computing the second derivative matrix». No dice que esté rota ni que sea peor.

La hipótesis que encaja con todo lo que hay —y es hipótesis, no hallazgo:

* `drvmlest.c` es el driver de **VARMA multivariante** (`varmax.m` = nº de
  series). `fue` lo hereda siendo univariante.
* La línea de al lado es `choldcp`, la Cholesky **modificada**: parchea pivotes
  no definidos positivos y publica números de aspecto impecable.
* Y `ERRORES_ESTANDAR.md` documenta que en `drtran` —multivariante— activar
  `fdhess` destapó un hessiano **singular**, porque su cast metía dos varianzas
  libres en `x[]` y concentraba `sigma2`.

O sea: en el caso general multivariante `fdhess` + `choldcp` cambia un error
**ruidoso y detectable** (la matriz del BFGS) por uno **limpio e indetectable**.
En univariante ese problema no existe: `fue` no mete ninguna varianza en `x[]`.

**Qué hay que hacer antes de tocarlo:** buscar si Mauricio dejó la razón escrita
en alguna parte —artículo, notas, correspondencia— en vez de inferirla. El
recorrido por versiones ya está hecho y no dice nada: la línea es idéntica en las
siete, así que no hay un «antes y después» del que deducir el motivo.

Una función que un autor implementa entera, deja lista con el argumento correcto,
y aun así comenta —y mantiene comentada durante treinta años de versiones— tuvo
una razón. Descomentarla sin conocerla es repetir el experimento sin saber qué
falló la primera vez.

### 3 · La batería: hay dos, y apuntan a lados opuestos

* **`tests/test_reliability*.py`** fija `std_errors` contra la referencia
  **fue-1.13.1** (el binario en C) y contra el acuerdo C↔Python. Guarda la
  CONFORMIDAD del puerto, que es una propiedad declarada — y el arreglo de raíz
  la rompe a propósito, porque el C también está mal.
* **`fue-1.13.1/ERRORES_ESTANDAR.md`** tiene el **GLS exacto**, que es la verdad.
  Ahí el hessiano coincide al sexto decimal y deja de depender del arranque.

Hay que **decidir cuál se guarda**. Hoy la suite guarda la conformidad sin decir
que son cosas distintas.

Y falta **extender el GLS a más casos**: hoy hay uno (`ES_CPI_m10`, 13
parámetros). Para una batería hacen falta varios —con y sin estacionalidad, con
MA, con intervenciones—, y el GLS exacto sólo es fácil de construir cuando la
parte no determinista es tratable.

### 4 · Lo que hay que decidir además

* Qué hacer cuando el hessiano NO sea definido positivo con el paso correcto.
  Medido con ε^(1/4) la diagonal sale entera positiva **en el óptimo de este
  caso**; es un caso.
* Y que **las SE cambian ~20 %** respecto a lo que el paquete publica hoy. No es
  un arreglo transparente: reescribe la inferencia de todo lo ya calculado.

Si se adopta, caen por innecesarias media docena de defensas que `art` construyó
alrededor del síntoma: los dos detectores de covarianza-semilla (BUG-0027,
BUG-0041, BUG-0124), la cláusula de las SE del convenio de ficheros (BUG-0090,
BUG-0159) y el aviso de BUG-0168.
