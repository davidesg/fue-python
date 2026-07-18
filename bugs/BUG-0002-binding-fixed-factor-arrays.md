---
id: BUG-0002
title: Python binding caps AR/MA at 8 factors and factor order at 16 (fixed cdata arrays) — long-order models crash with IndexError
status: fixed
severity: high
component: binding
found_in: 0.1.5
fixed_in: 0.1.6
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

**Applied** in 0.1.6. Confirmed that `fue_api.c` copies each factor into
**dynamically allocated** engine structures (`Tm.Ar1[i] = vector(0, f->order)`,
`Tm.Ar1 = calloc(NumAr1+1, …)`) and that the macros are referenced **nowhere** in
the engine internals (`usmelard.c`, `elfvarma.c`, …) — the fixed arrays are a
pure transport buffer, so raising the caps is safe and needs no engine change.

Changes:

1. **Raise the transport caps** well beyond any realistic univariate model, in
   lockstep across the C header and the cffi cdef:
   - `csrc/fue_api.h`: `FUE_MAX_FACTORS 8 → 32`, `FUE_MAX_POLYORD 16 → 64`.
   - `src/fue/_build_cffi.py`: the mirrored cdef literals `[8] → [32]`,
     `[16] → [64]` (`coefs`, `coef_free`, `omega`, `delta`, `ar1..ma2`,
     `ar1f_*`, `ma1f_*`). `interventions[64]` and `ifadf[8]` are unrelated and
     left as-is.
2. **Python guard** in `src/fue/_engine.py` `_fill_factors`: validate the factor
   count against `_MAX_FACTORS` and each factor order against `_MAX_POLYORD`,
   raising a clear `ValueError` naming the limit instead of a raw cffi
   `IndexError`. Module constants `_MAX_FACTORS=32`, `_MAX_POLYORD=64` are kept
   in sync with the header.

A fully-dynamic marshalling (flat length-prefixed buffer, no fixed arrays)
remains a possible future refinement, but the raised caps already exceed any
realistic model and the guard makes the boundary explicit.

## Validation

Cross-backend on the real England annual series (n=258), the cases that used to
crash — fue-C 1.13 (`/usr/local/bin/fue`) vs fue-Python 0.1.6:

| Model                         | fue C            | fue Python       | Δ        |
|-------------------------------|------------------|------------------|----------|
| EN.4 = AR(18) unfactored      | −1590.8128831276 | −1590.8128831276 | 2.9e-11  |
| EN.4 = 9 × AR(2) factorised   | −1597.8503371512 | −1597.8503371512 | <1e-10   |

Both previously raised `IndexError` (`double[16]` / `FueFactor[8]`); both now
estimate with `ifault=0` and match fue-C to 10–11 digits (optimizer tolerance).

Regression tests `tests/test_bug_0002_binding_factor_caps.py` (fixtures
`tests/real_cases/en4_ar18.inp`, `en4_fac9x2.inp`):
- AR(18) unfactored and 9×AR(2) fit and match the fue-C references to 1e-6;
- marshalling at the new cap (`coefs[63]`, `ar1[31]`) succeeds, one past
  (`coefs[64]`, `ar1[32]`) still raises `IndexError` at the cffi layer;
- a model one past the cap (`_MAX_FACTORS+1` factors, order `_MAX_POLYORD+1`)
  raises a clear `ValueError`, not `IndexError`.

Full suite green after rebuilding the extension.
