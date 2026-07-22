---
id: BUG-0004
title: fit() does not write fitted params back to the model's ar/ar_s/mu0 attributes; forecast_fuf then forecasts from the stale pre-fit seeds and the level explodes
status: open
severity: high
component: forecast
found_in: 0.1.7
fixed_in:
reported: 2026-07-22
reporter: mtgp2 (DVR / GP-Note replication)
tags:
  - forecast
  - fit
  - state-sync
  - fuf
  - explosion
  - regression
references:
  - src/fue/*.py (Model.fit: does not sync fitted params -> self.ar/ar_s/ma/ma_s/mu0)
  - src/fue/cast_us.py (eval_at_params -> _build_initial_x(model): reads model ATTRIBUTES)
  - src/fue/forecast.py (forecast()), Model.forecast_fuf
  - bugs/BUG-0001-forecast-mean-drift.md (a DIFFERENT, smaller forecast bug, fixed in 0.1.5)
  - bugs/BUG-0004-repro/ (EA_CPI.pre + repro.py: self-contained reproduction)
---

## Summary

`Model.fit()` estimates correctly (the fitted values land in `_result.params`) but
**does not write them back into the model's own attributes** — `self.ar`, `self.ar_s`,
`self.ma`, `self.ma_s`, `self.mu0` keep their **pre-fit seed values**.
`forecast_fuf` forecasts via `eval_at_params(self)`, which rebuilds the parameter
vector with `_build_initial_x(model)` — i.e. from those **attributes**, not from
`_result.params`. So the forecast is computed with the *initial guess*, not the
fitted model. When the seeds differ materially from the fit — which is the normal
case, and is severe for models built by `art._make_model`, whose `mu0` seed is in
the **×100-rescaled** space (≈100× the raw-log drift) — the level forecast **runs
away**: euro-area HICP goes `107.6 → 112.7 → … → 136.3` over six months from a last
observation of `103.15`, instead of the correct mean-reverting `≈103.4`. The fit is
fine; only the forecast is wrong, and it warns nothing.

This is a **regression**: the identical `_make_model(...).fit().forecast_fuf()` path
produced correct forecasts under **fue 0.1.3**. It is *not* the drift double-count of
BUG-0001 (a small constant `μφ/(1−φ)`, fixed in 0.1.5).

## Impact

**Silent, catastrophic** for the standard workflow *build → fit → forecast* whenever
the seeds differ from the fit — i.e. essentially every real model, and every model
built with `art._make_model` (all the DVR-replication country CPI models). Consumers
of `forecast_fuf` (`report_forecast`, the fuf pipeline, forecast-path exercises) get
a runaway level and absurd YoY forecasts (Spanish CPI index `91 → 260` in six
months). Discovered while replicating Garcia-Hiernaux et al. (2026): the euro-area /
Spain / France / Germany inflation "spaghetti" charts had to be produced with
`drtran -0` (the exact-VARMA engine, which forecasts correctly from `_result`).

## Reproduction

Self-contained, in `bugs/BUG-0004-repro/` (`EA_CPI.pre` holds the euro-area HICP
data; `repro.py` rebuilds and fits the model):

```text
$ python repro.py
fitted (correct, _result.params):  phi=0.3133  Phi=0.2853  mu=0.00176
model attributes AFTER .fit (STALE): phi=0.1508  Phi=0.7372  mu=0.17945   <-- seeds, not fit
forecast_fuf level h1..6 (uses stale attrs -> EXPLODES): [107.6 112.7 118.3 124.1 130.3 136.3] | last obs 103.15
forecast_fuf level h1..6 (after sync -> CORRECT):        [103.4 103.5 103.7 103.8 104.  104.1]
fue.load('EA_CPI.pre').forecast_fuf (CORRECT):           [103.4 103.1 103.3 103.6 103.8 103.7]
```

Note `mu` attribute `0.17945 ≈ 100 × 0.00176` — the ×100 rescale seed — which is what
makes the `_make_model` case blow up rather than merely bias. Two things confirm the
mechanism: (i) writing the fitted model to `.pre` (from `_result`) and **reloading**
forecasts correctly — so the `.pre` alone reproduces nothing, the *fit-then-forecast*
step is essential; (ii) manually copying `_result.params` into `m.ar/m.ar_s/m.mu0`
before `forecast_fuf` also fixes it. (A model whose seeds already sit near the fit —
e.g. a plain AR(1) on log-Brent, `P=0` — forecasts fine, which is why the bug can be
missed.)

## Root cause

`Model.fit()` stores the estimate only in `self._result` and leaves the parameter
**attributes** (`ar`, `ar_s`, `ma`, `ma_s`, `mu0`, and intervention `omega`/`delta`)
at their pre-fit seed values. `forecast_fuf` → `cast_us.eval_at_params(self)` builds
`x0 = _build_initial_x(self)` **from those attributes** and evaluates the model at
`x0` (no optimisation, by design). Hence the forecast uses the seeds. The seeds are
whatever the constructor/`_make_model` planted; for `_make_model` the `mu0` seed is
carried in the rescaled (×100) space, so it is ~100× the drift the level recursion
expects, and the level diverges. (This is the state-sync analogue of BUG-0001, which
touched the same drift term in `forecast.py`.)

## Fix

Not yet applied in fue. Two clean options:
- **Sync on fit:** at the end of `Model.fit()`, write `_result.params` back into
  `self.ar/ar_s/ma/ma_s/mu0` (and the interventions), in the *same* units the
  attributes use, so a subsequent `forecast_fuf`/`eval_at_params` sees the fit; or
- **Read the fit:** have `forecast_fuf`/`eval_at_params` build `x0` from
  `self._result.params` when a fit is present, falling back to the attributes only
  when it is not.
Add a guard/warning when the forecast one-step difference does not approach `μ`.

Consumer workaround (DVR replication): forecast via the exact-VARMA engine —
`drtran <CPI>.pre <X>.pre -0 -f H` — which forecasts from `_result` and is correct;
or sync `_result.params` into the attributes before calling `forecast_fuf`.

## Regression vs 0.1.3

`_make_model(...).fit().forecast_fuf()` was **correct in fue 0.1.3**: the
DVR-replication code (`cases/DVR_2026/ecb_vs_arima.py`, `spaghetti.py`) called
`m.forecast_fuf(...).level` on this exact euro-area model and produced sensible,
mean-reverting inflation forecasts (annual paths of 2–8% converging to ~2%, saved in
`ecb_vs_arima.json`); the identical code under **0.1.7 explodes**. So either `fit()`
synced the attributes in 0.1.3, or `forecast_fuf` read `_result` there — and that was
lost by 0.1.7. Timeline: BUG-0001 (small drift bias) was 0.1.4→fixed 0.1.5; this
state-sync blow-up is a separate, larger regression present in 0.1.7 (0.1.4–0.1.6 not
bisected). Definitive check: `pip install fue==0.1.3`, run `repro.py` → correct;
`==0.1.7` → explodes.

## Validation

When fixed, `repro.py` must show `forecast_fuf` returning the mean-reverting level
(`≈103.4` at h1, matching `_result`/`.pre`/drtran) directly after `.fit()`, with no
manual attribute sync; and the euro-area / Spain / France / Germany forecast
spaghettis should be reproducible from `forecast_fuf` (they currently require
`drtran -0`).
