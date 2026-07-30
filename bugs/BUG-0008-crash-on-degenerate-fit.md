---
id: BUG-0008
title: fue C dies instead of reporting — the reporting plots segfault on a degenerate (zero-variance) series, and GSL's default handler aborts the process when the eigensolver fails
status: fixed
severity: high
component: interop
found_in: 0.1.8
fixed_in: 0.1.9
reported: 2026-07-30
reporter: David E. Guerrero
tags: [fue-c, crash, segfault, robustness, nan]
references: [BUG-0007, BUG-0009]
---

## Summary

Not in this package: in fue C, found from here while cleaning up after BUG-0007,
where 1 of the 95 regenerated files still crashed and two more aborted. Two
independent defects, same shape — **a numerical edge case kills the process
instead of being reported**:

* **A. The reporting plots segfault on a degenerate series.** When a fit collapses
  and leaves residuals that are all zero, `AbsMax` in `File_PlotSer` is computed as
  `0/0` = **NaN**, and *neither* of that routine's two guards catches it, because
  every comparison against NaN is false. From there `HorInc = 25/NaN`,
  `iround(NaN)` is garbage, and `Tmpstr[27 ± garbage]` runs off the buffer. The
  ACF/PACF plot (`PlotCor`) has the same defect from the other end: a NaN
  correlation gives a garbage run length, which the compiler turns into a `memset`
  past the end.
* **B. GSL aborts the process.** `gsl_eigenqr` (`nlatools.c`) calls
  `gsl_eigen_nonsymm` without installing an error handler. When the QR iteration
  does not converge (`francis.c:209: maximum iterations reached without finding
  all eigenvalues`), **GSL's default handler calls `abort()`**. Not converging is
  not a catastrophe: it is a point where stationarity cannot be certified.

## Impact

A crash, not a wrong number — but a crash in the *reporting* stage, after the
estimation, so what is lost is the run. On a batch (or an optimizer sweep, or a
teaching session) it takes the whole process down with no diagnostic beyond the
signal.

Reproduced on three real files:

| file | before |
|---|---|
| `dolarization/.../Mod/Coint/R.4.pre` | SIGSEGV in `File_PlotSer` |
| `SRC/drvus-source/{1.02,1.2.01}/.../PV4.11.inp` | SIGABRT from GSL |

## Reproduction

```sh
fue R4 eml     # R4.inp = the R.4.pre above
# Program received signal SIGSEGV
# 0x… File_PlotSer (ser=…) at src/diagnose.c:890
# 890    if ( Tmpstr[27 + iround( BandPos1 )] == ' ' )
# #1     main () at src/fue.c:1451        <- plotting the RESIDUALS
```

With gdb at the crash: `ser->var = 0`, `ser->mean = 0`, `nobs = 68` — the residual
series is identically zero (that model reports `logelf 0.0000000000`).

```sh
fue PV eml     # PV.inp = PV4.11.inp
# gsl: francis.c:209: ERROR: maximum iterations reached without finding all eigenvalues
# #5  gsl_error ()          <- the DEFAULT handler calls abort()
# #7  gsl_eigen_francis ()
# #8  gsl_eigen_nonsymm ()  <- nlatools.c:346
```

## Root cause

**A.** `diagnose.c:823-838`. `rtmp4 = sqrt(ser->var)`; with `var == 0`,
`AbsMax = |(x−mean)/0| = NaN`. The guards `if (AbsMax <= 2.0)` and
`if (AbsMax > 8.0)` are both false for NaN, so control falls through to
`HorInc = 25.0/AbsMax`. Telling detail: the `*` marker at `diagnose.c:887` already
carries a bounds check — the same defect was found and patched *there*, and the
four band writes right below it (890-897) were left unguarded. `PlotCor`
(`diagnose.c:1348-1360`) repeats it: `posi = abs(iround(corr[i] * 25.0))` with a
NaN correlation, then a `for` loop writing `posi` characters.

**B.** `nlatools.c:346`. No `gsl_set_error_handler` anywhere in the program, so
GSL's default — `abort()` — is in force for every call.

## Fix

**A.** In `File_PlotSer`, reject the degenerate series before computing anything
from it (`!(rtmp4 > 0.0) || !(AbsMax == AbsMax)` → warn and skip the plot, using
the routine's existing escape). In `PlotCor`, treat a non-finite correlation as
zero and clamp `posi` to 25 — a correlation lives in [−1, 1], so 25 is the only
legitimate maximum — and bounds-check the two band writes. Written as
`!(x > 0)` and `x == x` rather than `<=` and `isnan` so the NaN is caught without
depending on compiler flags.

**B.** Turn GSL's handler off around the call, check the status, and on failure
return roots **outside** the unit circle. That is how "no vale" is spelled in this
interface: the caller (`elfvarma.c:532`) sets `ifault` when a modulus reaches
1.00005, and the estimator moves away from the point. A rejected point is an
answer; `abort()` is not.

Applied in the fue C repo (`fue-1.13.1/src/diagnose.c`, `src/nlatools.c`).

## Validation

Both files now finish with `rc = 0`:

| file | after |
|---|---|
| `R.4.pre` | `rc=0`, `logelf 0.0000000000`, `Warning: series with zero variance; plot skipped` |
| `PV4.11.inp` | `rc=0`, `logelf -nan` |

Normal output is untouched — the same synthetic monthly case gives
`logelf -300.0622719126` before and after, with its ACF bands and markers
identical. The fue Python suite (which drives this binary for all nine
deterministic types) stays at **676 passed**.

Defect **B** also lives in this package, because `csrc/internal/` embeds a copy of
the same sources — filed and fixed as BUG-0009. Defect **A** does not:
`diagnose.c` is not among the embedded files.

## What this does NOT fix

The crashes are gone; the *diagnosis* is still poor. `R.4` reports a
log-likelihood of exactly 0 from a fit whose residuals are identically zero, and
`PV4.11` reports `logelf -nan`. Both are degenerate fits announced as if they were
results. A program that cannot fit a model should say so, not print a number —
worth its own report, and a different kind of work from this one.
