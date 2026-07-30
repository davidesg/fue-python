---
id: BUG-0006
title: compimp is read as a plain impulse, and easter/trend are not supported — a .pre estimates a different model in silence
status: fixed
severity: high
component: inp
found_in: 0.1.8
fixed_in: 0.1.9
reported: 2026-07-30
reporter: David E. Guerrero
tags: [silent, interop, deterministic, compimp, easter, trend]
references: [BUG-0007]
---

## Summary

fue C builds **nine** deterministic regressors (`fue.c:295-435`): `impulse`,
`compimp`, `step`, `ramp`, `easter`, `trend`, `cos`, `sin`, `alter`. This package
built six of them correctly and got three wrong:

* **`compimp`** — the COMPENSATED impulse (+1 at the date, **−1 at the next
  period**) was mapped to a plain `pulse` by the reader (`inp.py:276`), dropping
  the −1. **This is a silent failure**: the file loads, the model converges, the
  report looks healthy — it is simply not the model the `.pre` asked for.
* **`easter`** and **`trend`** — not in `Intervention.TYPES` at all. These fail
  loudly, but with an unhelpful message, and one of the two paths fails somewhere
  else entirely (see *Impact*).

A fourth, related defect is about shared vocabulary rather than arithmetic. fue C
**does not reject a keyword it does not know**: it takes it for a *non-standard*
variable whose data comes as extra columns of the series block. So every word the
two programs do not share is a silent misreading waiting to happen — and
`.pre`/`.inp` files are meant to be interchangeable between the two interpreters.
This package used `pulse` where the school and the file format say `impulse`, and
it has a `seasonal` type that fue C has no regressor for.

## Impact

Any model with a compensated impulse was **estimated wrong without any
diagnostic**. Measured on a synthetic monthly AR(1): −299.844021 (fue C) against
−299.712394 (this package). On `M6_EI.pre` (Relloso 1997, Table 4, the m6 system):
−290.613205 against −292.495149, a difference of **1.88 in log-likelihood**, which
propagates to any multivariate work built on those univariate models — it is what
blocked reproducing the m6 targets in the drtran port.

`easter`/`trend` raise, so nothing silent there, but the two entry points fail
differently: loading an `.inp` raises `ValueError: type must be one of [...]`,
while loading a `.pre` written by fue C dies earlier and confusingly at
`inp.py:300` with `invalid literal for int()`. That second failure is not this
package's fault — the `.pre` is already corrupt when it arrives (BUG-0007).

The vocabulary mismatch bites in the direction Python → C: a file that says
`pulse` (accepted here, unknown there) or `seasonal` (a type with no keyword) is
read by fue C as a non-standard variable, quietly. Measured: −263.317088 for any
unknown keyword, against −299.445725 for the real thing.

## Reproduction

```python
import fue
from fue.cast_us import _build_indicator
from fue.intervention import Intervention

ind = _build_indicator(Intervention("compimp", at=4), 12, 12, 1, 2002)
# before: [.., 1, 0, ..]  — a plain pulse
# after:  [.., 1, -1, ..] — sums to zero
```

End to end, with fue C as the arbiter: build a `.pre` with a single deterministic,
let fue C fit it and rewrite the `.pre` with its estimates, then evaluate this
package **at that same point**. If the regressor matches, the two
log-likelihoods match. Before the fix `compimp` did not, and `easter`/`trend`
raised. See `tests/test_bug_0006_deterministic_types.py`.

## Root cause

* `inp.py:276` folded `impulse`, `pulse` and `compimp` into a single `pulse`
  type. The compensated impulse is a *different regressor*, not a spelling.
* `Intervention.TYPES` had eight entries; `easter` and `trend` were missing, so
  `inp.py` handed the reader a type the constructor rejects.
* `_build_indicator` (`cast_us.py`) and `populate_globals` (`csrc/fue_api.c`) —
  the two places that build the indicator, one per backend — had no branch for
  any of the three.
* `report.py:_itv_name_line` ended in `else: return t`, emitting the internal
  type name for anything unknown; for `seasonal` that writes a word fue C
  misreads.

Why it survived: the canonical test cases use `impulse`, `step`, `cos`, `sin` and
`alter`. Of the six m6 series only `M6_EI.pre` carries a `compimp`, and the
difference does not look like a bug — it looks like a slightly worse fit.

## Fix

* `intervention.py` — `impulse` is now the canonical name (`pulse` kept as a
  deprecated alias, normalised on construction) and `compimp`, `easter`, `trend`
  are added. Codes are only ever appended: they are the contract with the C
  engine's `FUE_ITV_*`.
* `cast_us.py` — indicators for the three, plus `easter_date` and `obs_to_date`,
  ported verbatim from `nlatools.c:693` and `diagnose.c:40`. Ported rather than
  taken from a calendar library on purpose: the indicator has to be the one fue C
  builds, and a different-but-more-correct rule would silently change every
  estimate that uses it.
* `csrc/fue_api.{h,c}` — the same three regressors in the compiled backend, so
  both paths agree.
* `inp.py` — `compimp` is no longer folded into `impulse`.
* `report.py` — writes the format's own keywords, and **refuses** to write a type
  with no `.pre` representation (today only `seasonal`) instead of emitting a
  word fue C would take for a non-standard variable. Deterministic seasonality is
  written with harmonics (`cos`/`sin` plus the Nyquist `alter`), not with dummies,
  which is why fue C has no such regressor.

## Validation

`tests/test_bug_0006_deterministic_types.py` (25 tests). The end-to-end
parametrisation runs the fue C binary for **all nine** types and compares at fue
C's own optimum; `easter_date` is checked against the compiled `Easter()` for
2002-2020.

Regression baseline established in a separate worktree at `a56677c` with the
extension rebuilt there: **651 passed** before, **651 passed** after, then 25 new
tests on top.
