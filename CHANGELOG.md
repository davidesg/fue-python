# Changelog — fue

Exact maximum-likelihood estimation of univariate time series (ARMAX with
transfer functions). Semantic-ish versioning; see `bugs/` for the full reports.

## 0.1.16 — 2026-09-25

La 0.1.15 no llegó a publicarse: su contenido (BUG-0016) sale en ésta. Tres
bloques, más una lección de empaquetado.

### El contrato de ficheros, sin mentir — BUG-0017 a BUG-0022

Seis defectos del estudio del contrato `.inp`/`.pre`/`.out` (el banco de
conformidad de atws), los seis confirmados contra el código:

* **BUG-0017** — `load()` rellenaba con ceros la columna de un determinista no
  estándar que faltaba: un fichero mutilado se leía como un modelo con un
  regresor nulo. Ahora falla, como el motor.
* **BUG-0018** — `number` (serie sin fechar) y `cbands` se perdían al leer; se
  conservan y los escritores los devuelven.
* **BUG-0019** — el puerto rechazaba un `.pre` escrito por el propio motor
  (`phi2 = -0.0` en un factor de frecuencia fija).
* **BUG-0020** — `write_pre` pegaba los pares de δ sin separador: con dos o más,
  el `.pre` no lo releía nadie.
* **BUG-0021** — `write_pre` cuantizaba también los parámetros FIJOS (un AR
  fijado en 0.941176 volvía en 0.9412, que es otro modelo), y tiraba la μ fija
  no nula. Lo que es especificación relee idéntico.
* **BUG-0022** — en una serie anual con nombre numérico, `load()` tomaba el
  nombre por el año.

### La figura de residuos de un modelo sabe del modelo — BUG-0023

`plot_model_diagnostics` rotulaba `Q(lags − npar)` con `npar` = TODOS los
parámetros, armónicos e intervenciones incluidos: **Q(28) donde el test tiene 39
grados de libertad** (IPC de España, AR(1) fijo en 0). Fechaba el primer residuo
en el primer dato, sin el desfase de la diferenciación. Y sus retardos no eran
la regla de fug C.

Lo que sólo sabe el modelo vive ahora en `fue.diagnostics`, una vez, y es API
pública: **`default_lags`** (la regla de fug C, `diagnose.c`),
**`free_arma_count`** (los ARMA estimados, con los factores de frecuencia fija),
**`differencing_offset`** (d + D·s + las raíces de `ifadf`) y
**`residuals_start`**. `plots.py` y el `.out` las usan. art las usa desde su
0.2.2.

### `__version__` sale de la metadata — BUG-0016

Era la entrada de la 0.1.15; ver abajo.

### Los tests, con numpy 2 y pandas 3

Probadas las ruedas en un entorno limpio (numpy 2.5, pandas 3.0), 40 tests
fallaban por construir sus ficheros con `repr()` de escalares de numpy
(`np.float64(10.0)`), por el alias de frecuencia `"A"` que pandas 3 retiró, o por
necesitar `drvarma` sin saltarse cuando falta. **La librería no:** escribe,
relee y reestima `.pre` con numpy 2 sin diferencias.

### Lo que sigue abierto

* **BUG-0015** — los errores típicos salen de la matriz que el BFGS acumula por
  el camino, no del hessiano en el óptimo. Esta versión trae el diagnóstico
  (la vía del hessiano existe y el paso está mal traducido), no el arreglo.
* **BUG-0005** — el caso R.4 llega a otro óptimo con numpy 2 que con numpy 1
  (logL 212.06 frente a 211.21 de C): la misma sensibilidad al entorno que el
  informe describe entre plataformas.

## 0.1.15 — 2026-09-07  ·  no publicada: su contenido sale en la 0.1.16

**`fue.__version__` estaba escrito a mano y se quedó tres versiones atrás**
(BUG-0016). La 0.1.14 recién publicada en PyPI declaraba `0.1.11` mientras su
metadata decía `0.1.14`.

Nadie lo vio porque la única prueba que lo tocaba —`test_smoke`— exigía que
fuese *una cadena no vacía*, y `"0.1.11"` lo es. Y el mismo defecto ya había
mordido por otro lado: la página de API de la 0.1.10 decía «fue 0.1.9». Aquello
se arregló **en el generador**, sin tocar la raíz, de modo que el dato siguió
escrito dos veces.

Ahora se deriva: el `pyproject.toml` del árbol si existe —en editable la
metadata se escribió al instalar y no se regenera—, si no la metadata del
paquete, y si no `"0+desconocida"`, que dice que no consta en vez de inventar un
número.

Dos pruebas, y hacen falta las dos: que coincida con el `pyproject`, y que **se
derive**. Sin la segunda, la primera pasaría con el número correcto escrito a
mano y volvería a quedarse atrás en la siguiente subida.

*Alcance, sin inflarlo:* `art` no sella la versión de `fue` en sus guiones, así
que ningún registro llevaba la cifra falsa. Lo que estaba mal es lo que ve quien
pregunta la versión desde Python.

## 0.1.14 — 2026-09-07

**Un óptimo absurdo ya no se reporta como éxito** (BUG-0005, la mitad que se
podía cerrar desde Linux).

- El optimizador es una búsqueda LOCAL sobre una superficie que puede ser
  multimodal, y no comprobaba si el óptimo al que llegaba tenía sentido. En el
  caso del informe la rueda de Windows se quedaba en una cuenca espuria con
  μ̂ = −0,144 —una inflación mensual del −14,4% para el IPC de EE.UU.— y lo daba
  por bueno con `converged=True, ifault=0` y ni un aviso.

  `_avisa_si_la_media_es_absurda` compara μ̂ con la media de `w` —la serie ya
  filtrada de deterministas y diferenciada, que es de la que μ es la media— y
  avisa si se aleja más de 5 desviaciones típicas.

- **El umbral está medido.** Sobre los 1.570 modelos con μ estimada del
  ecosistema:

      percentil 99.9                  0.73
      MÁXIMO observado                1.38
      por encima de 2·sd                 0
      el óptimo espurio de US CPI      47.2

  Factor 34 entre lo peor legítimo y el disparate. La guarda dispara **0 veces**
  sobre el corpus entero y sí sobre el espurio.

- Dos cosas que la medición enseñó y que el informe no decía, y que sin ellas la
  guarda sería inservible:

  **μ es la media de `w`, no de la serie observada.** Con intervenciones o
  armónicos los deterministas se llevan parte del nivel: comparando contra la
  serie cruda, 23 modelos correctos salían «absurdos». `w` la da `cast_us_py`,
  no se reconstruye — reimplementarla sería repetir BUG-0014.

  **μ puede no estar identificada.** La deriva es μ·φ(1); con un AR de raíz
  unitaria escrito con `d=0`, φ(1)=0 y μ no entra en la verosimilitud. Son 47 de
  los 1.570 modelos —los VIX con λ extrema y φ=1 exacto— y avisar de ellos serían
  47 falsos positivos.

- **Y la causa raíz del informe está mal.** BUG-0005 lo atribuía a una
  superficie multimodal cuya cuenca decide la aritmética del compilador. Medido
  hoy, la evidencia no lo sostiene:

  **Nueve arranques, un solo óptimo.** Multi-arranque sobre el AR estacional
  cubriendo la región y los dos signos: las nueve semillas válidas llegan a
  ℓ = 1322,513. Si hubiera cuencas separadas, un barrido así las encontraría.

  **Y el punto de Windows no es un punto estacionario.** A lo largo de μ el
  objetivo decrece monótonamente desde el óptimo (−268,677 → −270,107 en
  μ = −0,144) y la derivada allí vale **+19,57**. Un punto con gradiente no nulo
  no es una cuenca, así que llamarlo «otra cuenca» describe mal lo ocurrido.

  Lo que sí encaja es la línea que el propio informe dejó al margen: **una
  semilla con el signo cambiado**. La convención de `fue` es la de Box y Jenkins
  para todo operador —ω(B) = ω₀ − ω₁B − …— y sembrar al revés no es empezar
  «algo peor»: es empezar en la región equivocada, y desde ahí cualquier
  diferencia de último bit decide adónde se va. Arreglado eso en el consumidor
  (art/BUG-0006), el caso no se mueve desde ningún arranque.

  Consecuencia: el **multi-arranque** que el informe proponía como arreglo **no
  hace falta para este caso**, y lo que quedaba de defecto real era el silencio.

- Lo que **no** se puede afirmar, y por eso BUG-0005 queda `in-progress`: sólo se
  conocen μ̂ y σ̂ₐ del resultado de Windows, no el vector completo de 16
  parámetros, así que no se descarta que aquel punto fuera un óptimo local en el
  espacio conjunto. Y sigue sin haber máquina Windows. El motivo de que siga
  abierto ya no es «el optimizador es frágil» sino «falta verificar una
  observación de julio que el resto de la evidencia no acompaña».

## 0.1.13 — 2026-09-06

**La previsión daba efecto NULO a tres deterministas** (BUG-0014), y era una
duplicación de concepto la que lo causaba.

- `forecast._build_xi` construía el indicador **por segunda vez** —indexado por
  `type_code` en vez de por nombre— con ramas para 8 de los 11 tipos y sin
  `else` final. `compimp`, `easter` y `trend` caían fuera, su indicador salía
  idénticamente nulo y su efecto valía **exactamente cero**: sin aviso, sin
  error. El ω se estima, sale en el `.out` con su error típico, y la previsión
  lo ignora.

  El daño era doble, porque `xi` se usa dos veces: `nt − xi` limpia la HISTORIA
  para obtener el ruido y `f1 += xi` añade el efecto al FUTURO. Con `xi ≡ 0` el
  ruido que alimenta la recursión queda contaminado por un determinista que
  nadie quitó —lo que sesga toda la trayectoria, no sólo los meses afectados— y
  además la previsión no lleva el efecto.

  Medido con un efecto de Semana Santa de +4% (ω = 399,4 en centésimas de log):

      ruta            03/2020    04/2020    05/2020
      fuf (C)         100.183    104.178    100.185   ← correcto
      Python (antes)  100.184    100.184    100.186   ← el efecto no está
      Python (ahora)  100.184    104.178    100.186

  El síntoma estaba a la vista: la variación interanual de 04/2020 salía
  **−411,60%**, comparando una previsión sin Semana Santa contra un abril
  observado que sí la tenía. Ahora sale −12,20%.

  **`fuf` (el C) no tiene este defecto**: escribe el determinista en el fichero
  de previsión y extiende el calendario al futuro.

- **Se borra el duplicado**, que es el arreglo de fondo: `_build_xi` llama al
  generador único `cast_us._build_indicator` pidiéndole `nobs + horizonte`. Cada
  tipo se extiende **por su propia regla** —el easter por el calendario, el step
  por su definición— sin repetir ninguna. Dos generadores del mismo regresor no
  son una duplicación inocente: el segundo se queda atrás cuando el primero
  crece, y eso fue exactamente lo que pasó.

- **Un tipo desconocido ya no vale cero**: `_build_indicator` levanta
  `ValueError`. Devolver ceros es estimar —o prever— un modelo distinto del
  pedido sin decirlo.

- **`custom` deja de reventar** al pedirle una ventana más larga que sus datos:
  rellena lo que hay y el resto queda a cero. Sin esto no se podía pedir el
  indicador hasta `nobs + horizonte`.

- **BUG-0015 levantado**, que no es código sino un hueco del registro: la
  limitación de los errores típicos —vienen de la matriz que BFGS acumula por el
  CAMINO, no del hessiano en el óptimo— llevaba desde julio documentada en el
  repo del C y **no estaba en el índice de defectos**. Es la limitación más
  importante que tiene el paquete y quien consultara `bugs/` no la encontraba.

- La versión de `pyproject.toml` iba por detrás de sus propios informes: decía
  0.1.11 mientras BUG-0013 se cerraba «fixed_in 0.1.12».

## 0.1.11 — 2026-08-13

Lo que 0.1.10 dejó a medias, encontrado verificando la publicación.

- **La página de PyPI no enlazaba a ninguna parte**: `fue` no tenía
  `[project.urls]`. Ahora tiene Homepage, **Documentation** —el campo que PyPI
  muestra arriba del todo, apuntando al sitio— Repository, Issues, Changelog y
  el registro de defectos, que es público a propósito. El README enlaza también,
  porque es lo que PyPI muestra como cuerpo.

- **`fue.datasets.ripc` describía mal sus propios datos.** Decía «the series is
  the log of the Spanish CPI rescaled by 100», y no lo es: son los valores tal
  como fue los lee (~0.41-0.44), los mismos de `RIPC.1.inp`, y el modelo
  canónico aplica la transformación **él** (`boxlam=0`, `refactor=100`, `d=0`).
  Creerse el docstring significaba aplicar 100·log dos veces — un error que
  estima limpiamente y se lee plausible, que es la clase que sobrevive.

- **La referencia de API documentaba el paquete instalado, no el repositorio.**
  Primero fue la versión (0.1.10 publicó una página que decía «fue 0.1.9»), y
  después los docstrings: un arreglo en `src/` no llegaba al documento, y
  `--check` estaba de acuerdo porque las dos mitades leían la misma fuente
  equivocada. El generador antepone ahora `src/` y toma la versión de
  `pyproject.toml`.

- Y los **datos que viajan con el paquete entran en la referencia**: los cinco
  ejemplos empiezan con uno, y su docstring es donde se dice qué son los números.

## 0.1.10 — 2026-08-13

Documentation release, and one engine change that the documentation made
unavoidable.

### The engine says why it stopped

`raxopt` announces its verdict through `outputv`, which the binding sends to
`/dev/null`: it was computed and thrown away. `FitResult.converged` meant
`ifault == 0` — "nothing crashed" — so a fit that stopped because the iterates
froze, with a gradient of 0.01, came back as good.

- `FitResult.termcode`, `.niter`, `.gnorm` and `.termination` now come from the
  C. The three globals added to `qnewtopt.c` only **record** what raxopt already
  computed: no criterion, no announcement, no numerical behaviour changed.
- `converged` is `ifault == 0 and termcode in (0, 1)`, and anything else raises
  a `RuntimeWarning` naming the reason. Engine faults are still exceptions.
- The `.out` writes the convergence block the C wrote, wording included.

**BUG-0012 closed, and it was not the port**: Mauricio's own C, rebuilt today,
fails identically on Box-Jenkins Series A. The archived reference is a run at
80-bit x87 precision, reproducible with `-m32 -O0`
(`tools/reproduce_drvus_reference.sh`).

**BUG-0010 and BUG-0011 closed** — they had been fixed in 0.1.9 and the reports
were never updated.

### Verified against the publications, not only against ourselves

Everything that checked the likelihood descended from one implementation.

- **AS 197 executed from the article.** Melard's FORTRAN is printed in full;
  it is transcribed in `tests/fortran/as197.f`, compiled, and run: nine
  Box-Jenkins specifications, agreement to **5e-08**.
- **AS 311 by its published identities.** Equations (2)-(4) on the engine's own
  outputs, and — the genuinely external part — its quadratic form against
  Melard's, to **1e-14**.
- `qnewtopt.c` enters the verbatim invariant, with the two stopping criteria
  compared character for character.

### Documentation

From an empty `docs/` to 3.100 lines: what the model is, the file contract, the
formal tests with the critical values the 2011 manual left blank, convergence,
provenance, the port, migration from the C, a generated API reference, and
**why the wheel** — the two engines measured, in speed (median ×90) and in
answers (largest difference over 23 real models: 0.0002). Five graded examples,
checked in the battery.

### Wheels

Linux (x86_64 and aarch64, glibc and musl), macOS on Apple Silicon, and Windows
AMD64 — 26 files, as in 0.1.9. **Intel macOS remains the one target not built**,
because GitHub's Intel-mac runners are chronically starved and every current Mac
is arm64; those users get the sdist or the pure-Python wheel, which needs no
compiler and gives the same answers (largest difference over 23 real models:
0.0002 in log-likelihood).

## 0.1.9 — 2026-07-30

Deterministic-variables release. fue C builds **nine** deterministic regressors;
this package built six of them right. Found while porting drtran's transfer
network, which could not reproduce the m6 targets.

- **BUG-0006** (inp, **silent**): `compimp` — the COMPENSATED impulse, +1 at the
  date and **−1 the next period** — was read as a plain impulse, dropping the −1.
  Nothing failed: the file loaded, the model converged, the report looked healthy.
  It simply was not the model the `.pre` asked for. On `M6_EI.pre`: −292.495
  instead of −290.613, **1.88 of log-likelihood**. `easter` and `trend` were
  missing outright (those did fail, loudly). The three are now built in **both**
  backends — `cast_us.py` and `csrc/fue_api.c` — with `easter_date` and
  `obs_to_date` ported verbatim from fue C rather than taken from a calendar
  library: the indicator has to be the one fue C builds.
- **Shared vocabulary with fue C.** `impulse` is now the canonical type name —
  the school's word, and the format's — with `pulse` kept as a deprecated alias
  that is normalised away. This matters because **fue C does not reject a keyword
  it does not know**: it takes it for a non-standard variable and estimates
  something else, quietly. So writing a `.pre` now refuses to emit a type with no
  representation in the format (today only `seasonal`, which has no fue C
  regressor because deterministic seasonality goes in harmonics, not dummies).
- **BUG-0007** (interop, **silent**, *in fue C*): its `.pre` writer omits `easter`,
  tests `"time"` for `trend` — writing to the LaTeX file — and has no branch for
  **non-standard** variables either, so the type's line comes out empty and **fue C
  cannot re-read its own `.pre`**. A sweep of the ecosystem (5636 `.pre`/`.inp`)
  finds zero files with `easter`/`trend` and **98 corrupt**, all of them the
  non-standard case, where re-reading does not give a wrong number — it
  **segfaults**. The data columns were always written; only the name line was
  missing, so one word restores the file, and 97 of the 98 have an intact sibling
  `.inp`. Present since 1.01. Fixed in the fue C repo; guarded from here, since fue
  C has no battery of its own.

- **BUG-0008** (**crash**, *in fue C*): found cleaning up after BUG-0007. Two
  independent defects of the same shape — a numerical edge case **kills the
  process instead of being reported**. (a) The reporting plots segfault on a
  degenerate series: with zero-variance residuals `AbsMax` is `0/0` = NaN, which
  *neither* guard catches because every comparison against NaN is false, and the
  band writes run off the buffer; `PlotCor` repeats it with NaN correlations.
  (b) `gsl_eigenqr` calls GSL with no error handler installed, so when the QR
  iteration fails to converge **GSL's default handler aborts the process** — and
  not converging just means stationarity cannot be certified, which this interface
  already knows how to express (roots outside the unit circle → the caller sets
  `ifault` → the estimator moves away). Fixed in the fue C repo; normal output is
  byte-identical.
- **BUG-0009** (binding): defect (b) above is **in this package too**, because
  `csrc/internal/` embeds a copy of those C sources. In the standalone program an
  `abort()` kills a run; inside an extension module it kills the **interpreter** —
  notebook and all, with no traceback and no exception to catch. Latent rather
  than observed (the model that reliably aborts fue C fits fine here), fixed
  anyway: the code was byte-identical to the code that does crash, and a library
  may fail but may not take the interpreter with it.

- **nlatools (robustness):** `vector`/`ivector` now return the **offset** pointer
  (`v - nl`), so `v[nl..nh]` is addressable for any `nl` — the contract every
  caller assumes, and what the Numerical-Recipes cleanup dropped. It does not
  bite here today (fue allocates with `nl < 0` only in `elf`'s
  `gamwa = tensor(-q+1, 0, ...)`, already fixed), but it is a latent waiting for
  the first `vector(-k, k)`: that is exactly what bit drtran, whose
  identification allocates `vector(-nlags, nlags)`. `matrix`/`imatrix` are left
  alone on purpose — the same change breaks fue, so their layout is not
  interchangeable with the copy shared by drtran/drvarma.

Also: the package's console scripts are now `fue-py`/`fuf-py`. Declaring them as
`fue`/`fuf` shadowed the C programs, because `~/.local/bin` comes before
`/usr/local/bin` in the PATH — so `fue` ran the port while the user believed they
were running the original.

Regression baseline taken in a separate worktree at `a56677c` with the extension
rebuilt there: **651 passed** before, **651 passed** after, plus 25 new tests.

## 0.1.8 — 2026-07-23

Rescaling-consistency release. Traced with ART (`docs/RESCALING_ARCHITECTURE.md`):
the `refactor` (×100 conditioning) is a single per-model value, and every
attribute-consumer must see the *fit*, not the pre-fit seed.

- **BUG-0004** (forecast): `forecast_fuf` forecast from the **stale pre-fit seed
  attributes** — `eval_at_params`/`_build_initial_x` rebuilt `x0` from `ar/ar_s/mu0`,
  which `fit()` never overwrote. With ART's ×100 `μ0` seed the level exploded (euro
  HICP 103→136 in six months). Fix: `eval_at_params` reads `_result.params` when
  present. The fit itself (and the written `.pre`) were always correct.
- **fit sync (rescaling P4):** `Model.fit()` now calls `sync_params_to_attrs()` —
  the invertible-normalised `_result.params` are written back into
  `ar/ar_s/ma/ma_s/mu0` and the interventions (single scale, a plain copy). *The
  model IS the fit after fitting*, so forecast/`.pre`/reports all agree.
- **BUG (plots):** `plot_residuals_ts`'s percent header used `×refactor`; with the
  now-consistent `refactor=100` (residuals in the ×100 space) it double-counted and
  showed `σ̂_w = 25.30%` instead of `0.25%`. Fixed to `×100/refactor`, matching the
  rest of `plots.py` / `report_forecast.py`.
- **BUG-0005** (filed, open): optimizer can land in a spurious optimum on multimodal
  surfaces from a bad seed (guard + multi-start pending). Orthogonal to the rescale.

## 0.1.7 — 2026-07-19

**First binary wheels on PyPI.** cibuildwheel builds cp310–cp313 wheels for
Windows (amd64), macOS (arm64), and Linux — manylinux **and** musllinux, both
x86_64 and aarch64 — with GSL bundled inside the extension, plus the pure-Python
wheel and the sdist. `pip install fue` no longer needs a C compiler or GSL on
those platforms.

- **BUG-0003** (plots): `plot_residuals_ts` drew no year ticks/dividers for annual
  series (`freq==1`), so the decimal-year x-axis was unreadable. Added a `freq==1`
  branch replicating fue-C `gnuplot_File_PlotSer_CorrSer` (labels every 20 years
  anchored at the begin year, `tsby + 20·i`).
- CI (`wheels.yml`): fixed Windows GSL discovery (`$VCPKG_INSTALLATION_ROOT` bash
  expansion + forward slashes), macOS GSL discovery (`_discover_gsl_dirs` via
  `gsl-config`/Homebrew) + `MACOSX_DEPLOYMENT_TARGET` pinned to the runner, and
  the Linux `before-all` made portable (dnf on manylinux / apk on musllinux).
  Per-wheel test narrowed to a fast `test_smoke.py` (the golden battery is
  platform/BLAS-sensitive — e.g. the multimodal cointegration case R.4 — and stays
  a dev-only test). Intel macOS (macos-13) dropped from the matrix (runners
  chronically starved; Intel Macs are legacy — sdist/pure cover them).

## 0.1.6 — 2026-07-18

- **BUG-0002** (binding): the cffi `FueModelSpec` capped AR/MA blocks at 8 factors
  (`FueFactor[8]`) and each factor at order 16 (`coefs[16]`), so unfactored
  order ≥17 and ≥9-factor models crashed with `IndexError` in the Python binding
  where fue-C runs. The engine (Tusmodel) allocates factors dynamically — these
  were transport-buffer caps only. Raised to `FUE_MAX_FACTORS=32`,
  `FUE_MAX_POLYORD=64` (header + cdef in sync) with a clear `ValueError` guard.
  Validated vs fue-C on England: AR(18) and 9×AR(2) now match to 10–11 digits.

## 0.1.5 — 2026-07-18

- **BUG-0001** (forecast): the level forecast over-shot by `μ·φ/(1−φ)` (AR(1)) —
  the mean drift was double-counted (accumulated `l·μ` on top of the initial
  conditions). Catastrophic for `d=0` (the level exploded). Fixed to the mean
  form: seed the intercept `c = μ·(1−Σφ)` inside the level recursion. The same fix
  was applied to the C reference (fuf 1.08.2). `drtran`/`drvarma` were already
  correct.
