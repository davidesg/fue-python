---
id: BUG-0025
title: The line search never returns when the objective is NaN — the C spins for ever, the Python port crashes on the NaN gradient
status: fixed
severity: high
component: estimation
found_in: 0.1.16
fixed_in: 0.1.17 (unreleased)
reported: 2026-09-27
reporter: David — found by comparing the wheel's C with the atsw-gui monorepo after BUG-0015
tags:
  - optimizer
  - lnsrch
  - hang
  - nan
references:
  - csrc/internal/qnewtopt.c (lnsrch)
  - src/fue/qnewtopt.py (_lnsrch, _cholsol)
  - atsw-gui lib/optim/lnsrch.c (the same fix, shared by drvarma, drtran, fue and fuf CLIs)
  - drvarma-python BUG-0006 (the same defect, fixed there on 2026-09-26)
  - BUG-0015 (whose objcfunc guard already closes the path through Model.fit)
---

## Summary

`lnsrch`, the backtracking line search of raxopt (Dennis & Schnabel A6.3.1),
decides with comparisons. With a NaN every comparison is false. In the C the
step was neither accepted nor abandoned, lambda became NaN, and the loop spun
for ever: no error, no output, until the process was killed. The Python port
survived only by accident, because one comparison is written the other way
round. It then took the NaN into the gradient and died with a scipy
`ValueError`.

The rest of the family had already fixed this: drvarma-python BUG-0006 and
the shared `lib/optim/lnsrch.c` of the atsw-gui monorepo (drtran met it on
real data and ran for an hour and a half). The wheel's copy of the C had not
been fixed.

## Impact

- **Through `Model.fit()`:** closed since the BUG-0015 commit (c8860e5). Its
  `objcfunc` turns a non-finite objective into the inadmissible value 1.0, so
  no NaN reaches the line search. Before that commit, a likelihood that
  overflowed in the direction of a step hung the fit.
- **Through the optimiser directly:** open until this fix.
  `fue.qnewtopt.raxopt` is public, and the C `raxopt` is reachable by any
  caller of the engine.

## Reproduction

`bugs/BUG-0025-repro/`: a quadratic that is NaN beyond x = 1, where the first
full step lands.

    timeout 20 python3 bugs/BUG-0025-repro/repro_python.py
      before: ValueError: array must not contain infs or NaNs
      after:  returned: [0.99999468 -0.33333156] 0.44444680... termcode 2 niter 25

    gcc -Icsrc/internal bugs/BUG-0025-repro/repro_c.c csrc/internal/qnewtopt.c \
        csrc/internal/nlatools.c $(pkg-config --cflags --libs gsl) -lm -o /tmp/r
    timeout 10 /tmp/r; echo exit=$?
      before: exit=124 (killed by the timeout)
      after:  returned: lambda 0.1 retcode 0 x 0.6 f 5.76

After the fix, the full `raxopt` in C and in Python stops at the same point,
with the same termcode and the same iteration count.

## Root cause

In the C, `if ( tlambda <= 0.1 * lambda ) lambda = 0.1 * lambda; else lambda =
tlambda;` with `tlambda` NaN takes the `else`, so lambda becomes NaN. From
then on `lambda < minlam` is false and the loop never ends. The first-backtrack
test `lambda == 1.0` also assumed that a smaller lambda always had a previous
point to interpolate from.

The Python port writes the same test as `tlambda > 0.1 * lam`. With NaN it
takes the other branch and shrinks, so it survived by accident. It then
interpolated with the NaN and, a few iterations later, met the NaN in the
central-difference gradient. scipy's `solve_triangular` refuses NaN, where
the C's `cholsol` passes it through.

## Fix

The same fix as `lib/optim/lnsrch.c`:

- **A non-finite trial value is an inadmissible point.** lambda shrinks by
  0.1 without interpolating. Once lambda is below `minlam`, the search gives
  up like any other failed search.
- **The first-backtrack test is a flag,** `haveprev`, not `lambda == 1.0`.

On a finite trajectory the arithmetic is the original's, operation for
operation. The optimiser's test files and every regression pin are unchanged.

- **C:** every added line is marked `BUG-0025`. The provenance test
  (`test_the_optimizer_is_still_mauricios`) undoes that exact block before
  comparing with fue-1.13.1. The block must be there exactly once, so it
  cannot hide anything else.
- **Python:**
  - `_lnsrch` has the same explicit branch.
  - `_cholsol` no longer asks scipy to check finiteness. Like the C, a NaN
    gradient gives a NaN direction, the line search gives up, and raxopt
    stops.
  - `minlam` is infinite for that NaN direction, as it is in the C.

## Validation

`tests/test_bug_0025_lnsrch_nan.py` runs each case in a subprocess with a
timeout, so a regression fails instead of hanging the suite. It covers the
Python `raxopt`, the Python `_lnsrch`, and the C `lnsrch` built from the
harness. The harness hangs on the previous C (exit 124).
