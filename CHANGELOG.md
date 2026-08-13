# Changelog — fue

Exact maximum-likelihood estimation of univariate time series (ARMAX with
transfer functions). Semantic-ish versioning; see `bugs/` for the full reports.

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
