---
id: BUG-0009
title: The embedded C engine calls GSL with no error handler, so a failed eigensolve aborts the Python interpreter
status: fixed
severity: high
component: binding
found_in: 0.1.8
fixed_in: 0.1.9
reported: 2026-07-30
reporter: David E. Guerrero
tags: [crash, abort, robustness, gsl, embedded-c]
references: [BUG-0008]
---

## Summary

This package does not merely *call* fue C — it **embeds a copy of it**
(`csrc/internal/`: `nlatools.c`, `elfvarma.c`, `drvmlest.c`, `qnewtopt.c`,
`usmelard.c`), compiled into `_fue_engine`. So a defect in those files is a defect
here too, and BUG-0008(b) is one of them:
`csrc/internal/nlatools.c:346` calls `gsl_eigen_nonsymm` with **no GSL error
handler installed**, and GSL's default handler calls `abort()`.

In the standalone program that kills a run. **Inside an extension module it kills
the interpreter** — the notebook, the batch script, the whole session — with no
Python traceback, no exception to catch, nothing but a signal. Same defect, worse
blast radius.

## Impact

Latent rather than observed: the model that reliably aborts fue C
(`drvus-source/.../PV4.11.inp`, where the QR iteration reports *maximum iterations
reached without finding all eigenvalues*) **does not** abort here — fue Python
fits it to `logL = -51.1208` and returns normally, because its optimizer takes a
different path and never hands GSL that matrix.

So this is a hazard, not a live failure. It is filed and fixed anyway on two
grounds: the code is identical to the code that does crash (verified by diffing
the solver region against `fue-1.13.1/src/nlatools.c` — byte-for-byte the same
before the fix), and the consequence of it ever triggering is disproportionate for
a library. A library may fail; it may not take the interpreter with it.

Note what does **not** apply: BUG-0008(a) — the plotting segfaults on a
zero-variance series — cannot happen here. `diagnose.c` is not among the embedded
sources; this package draws with matplotlib (`plots.py`), which at worst raises.

## Reproduction

Not reproducible from Python on the known input (see *Impact*). The defect is
established by inspection plus the C-side reproduction in BUG-0008:

```sh
grep -n 'gsl_eigen_nonsymm\|gsl_set_error_handler' csrc/internal/nlatools.c
# 346:   gsl_eigen_nonsymm( A, eval, w );      ← no handler anywhere
```

## Root cause

`gsl_set_error_handler` appears nowhere in the embedded sources, so GSL's default
— `abort()` — governs every call into the library.

## Fix

The same one applied to fue C in BUG-0008: turn the handler off around the call,
check the returned status, and on failure return roots **outside** the unit
circle. That is how "cannot certify stationarity" is expressed in this interface —
the caller (`elfvarma.c`) sets `ifault` when a modulus reaches 1.00005, and the
estimator moves away from the point. A rejected point is an answer; `abort()` is
not.

Keeping the two copies in step matters more than the fix itself: the embedded tree
is a *copy*, and a copy that drifts from its original is how a bug fixed in one
place survives in the other.

## Validation

Extension rebuilt and the suite re-run: **676 passed**, unchanged. The model that
aborts fue C still fits here to `logL = -51.1208`, so the fix costs nothing on the
normal path — which is the point: it only changes what happens where the program
used to die.
