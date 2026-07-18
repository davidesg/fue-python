---
id: BUG-0002
title: Python binding caps AR/MA at 8 factors and factor order at 16 (fixed cdata arrays) — long-order models crash with IndexError
status: open
severity: high
component: binding
found_in: 0.1.5
fixed_in: 
reported: 2026-07-18
reporter: D. E. Guerrero
tags:
  - binding
  - cffi
  - factors
references:
  - src/fue/_build_cffi.py (FueFactor.coefs[16], FueModelSpec.ar1/ar2/ma1/ma2[8])
  - src/fue/_engine.py (_fill_factors)
  - /home/david/Dropbox/Cycles/bugs_art_fue.md (BUG #6, England EN.1-EN.4)
---

## Summary

The cffi `FueModelSpec` struct declares every ARMA block as a **fixed-size array**
of factors, and every factor as a **fixed-size array** of coefficients:

- `FueFactor.coefs[16]` / `coef_free[16]` — a single unfactored operator of order
  **≥17** overflows.
- `FueModelSpec.ar1[8]`, `ar2[8]`, `ma1[8]`, `ma2[8]` — a factorised operator with
  **≥9 factors** overflows.

`_fill_factors` (`_engine.py`) writes `spec_arr[i].coefs[j] = v` straight into
these arrays, so a model past either cap raises `IndexError` in cffi and the
estimation aborts. The underlying C engine (`fue` C / `usmelard`) has **no such
fixed limit** — it estimates these models fine — so fue-Python silently fails
where fue-C succeeds, breaking the documented invariant that the two backends
agree.

## Impact

Any long-order or many-factor model cannot be estimated through the Python
binding:

- **AR/MA order ≥ 17 unfactored** → `IndexError: index too large for cdata
  'double[16]'`. (AR(16) is the largest that works.)
- **≥ 9 AR(2)/MA(2) factors** → `IndexError: index too large for cdata
  'FueFactor[8]'`. (8 factors is the largest that works.)

Concretely this blocks the England annual models in the *Joseph's Cycles*
project — EN.1 = AR(16), EN.2 = AR(15), EN.3 = AR(14), **EN.4 = AR(18)** and its
factorised form (9 × AR(2)) — all of which fue-C estimates but fue-Python cannot.
The workaround has been to fall back to fue-C (`/usr/local/bin/fue`).

Where both backends *can* run, they agree exactly: AR(16) unfactored on England
gives `logelf = -1600.0924204792` in **both** backends (identical to 10 digits),
confirming the only defect is the artificial cap, not the numerics.

## Reproduction

No fit needed — the caps are visible at the struct level:

```python
from fue._fue_engine import ffi
spec = ffi.new("FueModelSpec *")
spec.ar1[0].coefs[16] = 1.0   # IndexError: index too large for cdata 'double[16]'
spec.ar1[8].order     = 2     # IndexError: index too large for cdata 'FueFactor[8]'
spec.ar1[0].coefs[15] = 1.0   # OK  (order 16 is the last that fits)
spec.ar1[7].order     = 2     # OK  (8 factors is the last that fits)
```

End-to-end: fit `Model(..., ar=[[0.0]*18])` (one AR(18) factor) or a factorised
model with 9 `AR(2)` factors — both abort in `_fill_factors`.

| Model (England, n=259)   | fue C                    | fue Python              |
|--------------------------|--------------------------|-------------------------|
| AR(16) unfactored        | logelf −1600.0924204792  | −1600.0924204792 ✓      |
| AR(17), AR(18) unfactored| OK (EN.4 = AR18)         | IndexError `double[16]` |
| 8 × AR(2) factors        | OK                       | OK                      |
| 9 × AR(2) factors (EN.4) | logelf −1597.8503        | IndexError `FueFactor[8]` |

## Root cause

`src/fue/_build_cffi.py`, the cffi `cdef`:

```c
typedef struct {
    int    order;
    double coefs[16];        /* <- caps a single factor at order 16 */
    int    coef_free[16];
} FueFactor;

typedef struct {
    ...
    int       nar1;  FueFactor ar1[8];   /* <- caps each block at 8 factors */
    int       nar2;  FueFactor ar2[8];
    int       nma1;  FueFactor ma1[8];
    int       nma2;  FueFactor ma2[8];
    ...
} FueModelSpec;
```

These sizes are a binding-side convenience, not a C-engine constraint. The engine
consumes the coefficients through `usmelard`/`cast_us` in dynamically sized form;
only the cffi marshalling struct imposes `[16]` and `[8]`.

(Related fixed sizes in the same struct — `interventions[64]`, `ar1f_*[8]`,
`ifadf[8]` — are not exercised by these cases but share the pattern and should be
reviewed together.)

## Fix

Two options, in order of preference:

1. **Size the marshalling to the model.** Compute the needed capacity from the
   model (`max_order = max(len(f) for f in all_factors)`, `max_factors =
   max(nar1, nar2, nma1, nma2)`) and either (a) generate the `cdef` with those
   sizes, or (b) pass coefficients through a flat, length-prefixed `double[]`
   buffer + counts (as `data` is already passed via `ffi.from_buffer`), removing
   the fixed factor/coef arrays from the struct entirely. Option (b) is the clean
   fix and mirrors how the C engine already consumes the data.

2. **Raise the caps** to a safe bound (e.g. `coefs[64]`, `ar1[32]`) as a stopgap.
   Cheap, but re-introduces the same failure mode for even larger models and
   wastes struct space. Only acceptable as an interim measure.

Whichever is chosen, `_fill_factors` should also **validate** `len(factor)` and
the factor count against the capacity and raise a clear `ValueError` ("AR order N
exceeds binding capacity M; rebuild with larger cdef") instead of a raw cffi
`IndexError`.

## Validation

- Add a regression test that fits (or at least marshals) an AR(18) unfactored
  model and a 9×AR(2) factorised model through the Python binding without
  `IndexError`, and asserts `logelf` matches the fue-C reference
  (EN.4 unfactored and factorised: −1600.09… / −1597.8503).
- Guard test: a model one past the (new) cap raises `ValueError`, not `IndexError`.
- Cross-check: for every model where both backends run, `logelf` agrees to ≥10
  digits (already true for AR(16); extend the England battery once the cap is
  lifted).
