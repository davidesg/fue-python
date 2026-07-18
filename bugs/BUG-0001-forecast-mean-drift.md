---
id: BUG-0001
title: Forecast level over-shoots by mu*phi(1)^-1 in the mean drift (drift double-counted)
status: fixed
severity: high
component: forecast
found_in: 0.1.4
fixed_in: 0.1.5
reported: 2026-07-18
reporter: D. E. Guerrero
tags:
  - forecast
  - mean
  - drift
references:
  - src/fue/forecast.py (forecast(), step [7])
  - fuf-1.08.1 src/usfo.c forecast() [1]-[2]
  - drtran docs/FUF_FORECAST_BUG.md
  - fuf-1.08.1-fix (patched C reference)
---

## Summary

For a model with a non-zero mean on the differenced series
`phi(B)[∇Y_t − μ] = a_t`, the point forecast of the **level** is systematically
too high by a constant that converges to `μ·φ/(1−φ)` (AR(1) case).  The
differenced/annual **rates** converge to the right value, so the error is easy to
miss; only the level is biased.  `fue.forecast.forecast()` reproduces the same
composition as the C `fuf` (`usfo.c`), so both carry the bug.  `drtran` and
`drvarma` (which forecast in first differences with the mean form and integrate)
are correct.

## Impact

Every level forecast of an inflation series with drift (μ≠0) — i.e. essentially
all the country CPI reports generated via `report_forecast` — has a small
positive bias in the level (≈0.1 % for the Spanish CPI, up to ~2.8 index points
for high-drift 1970s samples).  Growth-rate forecasts (month-on-month,
year-on-year) are correct asymptotically but biased by the same term for the
first `s` horizons (where the rate is forecast-minus-observed).  Models with μ=0
are unaffected.

## Reproduction

Forecast a differenced model with an estimated non-zero mean and compare the
level to the closed form `E[Y_{n+l}] = Y_n + l·μ + (∇Y_n − μ)·Σ_{j≤l} φ^j`
(or to drtran).  Hand example (φ=0.5, μ=1, ∇Y_n=3, Y_n=100):

| l | correct | fue/fuf |
|---|---------|---------|
| 1 | 102.00  | 102.50 (+0.5) |
| 2 | 103.50  | 104.25 (+0.75) |
| ∞ | —       | +μφ/(1−φ)=1.0 |

Real case (Spanish CPI, μ=0.154, φ=0.403): fue/fuf level is +0.104 (in the
100·log space) above drtran, exactly `μφ/(1−φ)`.

## Root cause

`forecast.py` step [6]–[7] (mirroring `usfo.c`) folds the differencing into the
AR (`phi0 = varphi(phi, rnsop)`), runs the **homogeneous** level recursion on the
raw de-seasonalised level `w = nt − xi`, and then adds the mean as an
**accumulated drift** `l·μ`:

```python
s2 += mu; f1[l] += s2            # step [7] — accumulated l·μ  ← the bug
```

The homogeneous recursion is seeded with the *actual* initial conditions, which
already contain the drift `μ·t`.  Adding `l·μ` on top counts the drift twice in
the transient, leaving a residual `+φ^l·μ` per step.  Formally the level obeys
`phi0(B) w_t = μ·φ(1) + a_t`, so the correct intercept is `c = μ·φ(1)` **inside**
the recursion, not `l·μ` afterwards.

## Fix

**Applied** in `src/fue/forecast.py` (step [6]–[7]): add the intercept
`c = μ·(1 − Σ φ_i)` at each step of the level recursion and drop the `l·μ`
accumulation.

```python
drift = mu * (1.0 - float(np.sum(phi_coefs))) if mu else 0.0   # μ·φ(1); d-independent
for l in range(1, L + 1):
    ...
    f1[l] = vtmp1 - vtmp2 + drift          # step [6]: + intercept
# step [7]: keep only the deterministic xi[nobs+l]; the l·μ block is gone
```

The same patch is applied and validated in the C reference (`fuf-1.08.1-fix`,
`usfo.c`: `drift[k] = mu[k]*(1 − Σ phi[i][k][k])`, added in [1]; `l·μ` removed
in [2]).

**d=0 is catastrophic, not just biased.**  For a *stationary* model (no
differencing) the accumulated `l·μ` adds a linear trend in the transformed space,
so the level *explodes* (e.g. the SFNY.2 golden values were 1.30 → … → 100.32
over 5 steps).  The fix converges to the mean.

## Validation

- Fixed `forecast()` equals the mean-form closed form to machine precision
  (≤1e-13) for d=1 ARIMA(1,1,0)+drift, AR(2), and d=0 stationary AR(1).
- Closed form and the C fix agree with drtran to display precision.
- Battery of 71 country models (8 countries, μ from 0 to 0.55): the patched C
  `fuf` matches drtran on every case where fuf runs (max |fix − drtran| ≤ 0.005,
  i.e. rounding), with over-shoots up to 2.8 index points removed; μ=0 models are
  a no-op (patched == original).  See drtran `docs/FUF_FORECAST_BUG.md`.
- Regression tests: `tests/test_bug_0001_forecast_drift.py` (closed-form match,
  d=0 convergence, μ=0 no-op).  The fuf-homologation golden values in
  `tests/test_forecast.py` (SFNY.2, Spain S.2) were regenerated with the fixed
  fuf; the Spain reference `.out` lives in `tests/real_cases/`.  Full suite:
  642 passed.
