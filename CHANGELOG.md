# Changelog — fue

Exact maximum-likelihood estimation of univariate time series (ARMAX with
transfer functions). Semantic-ish versioning; see `bugs/` for the full reports.

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
