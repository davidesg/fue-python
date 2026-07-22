---
id: BUG-0004
title: forecast_fuf forecasts from the model's pre-fit SEED attributes (mu0/ar/...) instead of the fit; with a mis-scaled mu0 seed the level explodes
status: fixed
severity: high
component: forecast
found_in: 0.1.7
fixed_in: 0.1.8
reported: 2026-07-22
reporter: mtgp2 (DVR / GP-Note replication)
tags:
  - forecast
  - fit
  - state-sync
  - fuf
  - explosion
  - rescaling
references:
  - src/fue/cast_us.py (eval_at_params -> _build_initial_x(model): read model ATTRIBUTES)
  - src/fue/model.py (Model.fit: does not sync fitted params -> self.ar/ar_s/ma/ma_s/mu0)
  - tests/test_bug_0004_forecast_stale_attrs.py (regression test)
  - bugs/BUG-0004-repro/ (EA_CPI.pre + repro.py: self-contained reproduction)
  - ART BUG-0006 / _mu_seed: the ×100 mu0 seed that makes it explode (rescaling mismatch)
---

## Summary

`Model.fit()` estimates correctly (the fit lands in `_result.params`) but **does not
write the estimate back into the model's own attributes** — `self.ar`, `self.ar_s`,
`self.ma`, `self.ma_s`, `self.mu0` keep their **pre-fit SEED values**. `forecast_fuf`
forecasts via `eval_at_params(self)`, which rebuilds the parameter vector with
`_build_initial_x(model)` **from those attributes**, not from `_result.params`. So the
forecast is computed with the *initial guess*, not the fitted model.

When a seed differs materially from the fit — the normal case — the forecast is wrong.
It is **catastrophic** for models built by `art._make_model`, whose `mu0` seed is
`100 × drift` (ART's ×100 rescale seed) while the in-memory model carries
`refactor = 1`: the stale `mu0 ≈ 0.18` is ~100× the drift the level recursion expects,
so the level **runs away** — euro-area HICP goes `103.15 → 107.6 → … → 136.3` over six
months instead of the correct mean-reverting `≈103.4`. The fit is fine; only the
forecast is wrong, and it warns nothing.

## Impact

**Silent, catastrophic** for the standard *build → fit → forecast* workflow whenever a
seed differs from the fit — i.e. essentially every model built with `art._make_model`
(all the DVR-replication country CPI models). Consumers of `forecast_fuf`
(`report_forecast`, the fuf pipeline) got a runaway level and absurd YoY forecasts
(Spanish CPI index `91 → 260` in six months). Discovered replicating Garcia-Hiernaux
et al. (2026); the inflation "spaghetti" charts had to be produced with `drtran -0`
(the exact-VARMA engine, which forecasts from `_result`).

## NOT a fue version regression (corrected)

The original report claimed this worked in fue 0.1.3 and broke by 0.1.7. **Verified
false:** `repro.py` explodes on fue **0.1.3 too** (installed the PyPI wheel in an
isolated venv → EA level 103 → 303 in six months). `fit()` never synced the attributes,
and `_build_initial_x` always read `model.mu0`, in every fue version.

What actually changed is upstream, in **ART**: `_make_model` seeds `mu0 = _RESCALE_FACTOR
(100) × drift`. Before **ART BUG-0001**'s fix (`f3de4a3`, "seed mu at the rescaled
transformed mean"), `mu0` was seeded at **0**, so reading the stale seed was harmless
(mu=0 → no runaway drift). That ART fix — verified as the trigger: with `mu0 = 0` the
same repro forecasts correctly (`≈103.3`) — turned this latent fue read-from-seed bug
into a catastrophic explosion. So it is a *latent fue bug exposed by an ART change*, not
a fue 0.1.x→0.1.y regression.

## Reproduction

Self-contained in `bugs/BUG-0004-repro/` (`EA_CPI.pre` + `repro.py`):

```text
fitted (correct, _result.params):   phi=0.3133  Phi=0.2854  mu=0.00176
model attributes AFTER .fit (STALE): phi=0.3175  Phi=0.2831  mu=0.17945   <- seeds
forecast_fuf (pre-fix, uses stale attrs): [107.6 112.7 118.3 124.1 130.3 136.3]  EXPLODES
forecast_fuf (post-fix, uses _result):    [103.4 103.1 103.3 103.6 103.8 103.7]  correct
fue.load('EA_CPI.pre').forecast_fuf:      [103.4 103.1 103.3 103.6 103.8 103.7]  correct
```

`mu` attribute `0.17945 ≈ 100 × 0.00176` — the ×100 rescale seed. Two confirmations of
the mechanism: (i) writing the fitted model to `.pre` (from `_result`) and reloading
forecasts correctly — so the *fit-then-forecast* step is essential; (ii) manually
copying `_result.params` into the attributes before `forecast_fuf` also fixes it.

## Root cause

`Model.fit()` stores the estimate only in `self._result` and leaves the parameter
attributes at their seed values. `forecast_fuf` → `cast_us.eval_at_params(self)` built
`x0 = _build_initial_x(self)` **from those attributes** (`cast_us.py`, the mu term is
`x.append(float(model.mu0))`) and evaluated at `x0` (no optimisation, by design). Hence
the forecast used the seeds. The `mu0` seed is mis-scaled (×100 vs the in-memory model's
`refactor=1`), so the level diverges — see the rescaling-architecture note: ART's
hardcoded ×100 seed is decoupled from `model.refactor`, and `fit()` not syncing the
attributes lets that mis-scaled seed leak to every attribute consumer (forecast here,
also `_write_inp`).

## Fix

Applied in `src/fue/cast_us.py::eval_at_params`: when a fit is present, build `x0` from
`model._result.params` (the fit), falling back to `_build_initial_x` (the attributes)
only when there is no fit — e.g. a model just loaded from a `.pre`, whose attributes ARE
the fitted values. `_result.params` is the same flat free-parameter vector in the same
canonical order, already normalised to the invertible MA root, so it is a drop-in.

```python
x0   = _build_initial_x(model)
_res = getattr(model, "_result", None)
_rp  = getattr(_res, "params", None) if _res is not None else None
if _rp is not None and len(_rp) == len(x0):
    x0 = np.asarray(_rp, dtype=float)
```

This is a correct point fix (the forecast always reads the fit). The **systemic** fixes,
tracked separately, are: (a) fue — `Model.fit()` should sync `_result.params` into the
attributes, so *every* attribute consumer (forecast, `_write_inp`, …) sees the fit; and
(b) ART — make `model.refactor` the single source of truth for the ×100 rescale (seed
mu on the already-rescaled series, set the Model's refactor), removing the hardcoded
×100 mismatch.

## Validation

`repro.py` now returns the mean-reverting level (`≈103.4` at h1, matching
`_result`/`.pre`/drtran) directly after `.fit()`, with no manual sync. Regression test
`tests/test_bug_0004_forecast_stale_attrs.py`: a fitted model forecasts near the last
observation, and the forecast is INVARIANT to corrupting the (unsynced) parameter
attributes — proving it reads the fit, not the attributes. Full fue suite: 651 passed.
