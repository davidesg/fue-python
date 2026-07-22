---
id: BUG-0005
title: ML optimizer converges to a spurious optimum on multimodal (seasonal-AR) likelihoods and reports converged=True with no diagnostic; the basin is platform-dependent (Windows vs Linux)
status: open
severity: medium
component: estimation
found_in: 0.1.7
fixed_in:
reported: 2026-07-22
reporter: mtgp2 (DVR / GP-Note replication)
tags:
  - optimizer
  - convergence
  - multimodal-likelihood
  - seasonal-ar
  - reproducibility
  - platform-dependent
references:
  - csrc/internal/qnewtopt.c, drvmlest.c, elfvarma.c (the ML optimizer + likelihood)
  - src/fue/_build_cffi.py (Windows /O2 + vcpkg GSL vs Linux -O2 + system GSL)
  - ART BUG-0006 (art-python/bugs/BUG-0006-*): the seed contamination that first triggered it
  - art-python/bugs/BUG-0006-repro/ (US_CPI.pre + repro.py — reproduces on Windows)
  - Garcia-Hiernaux, Gonzalez-Perez & Guerrero (2026), Econ. Modelling 157, Table 2 (US CPI)
---

## Summary

On a **multimodal** likelihood — an AR(2)×AR(2)_12 model whose seasonal AR has
complex conjugate roots (US CPI, monthly, 2002–2026) — fue's C optimizer
(`qnewtopt`/BFGS) converges to a **different basin depending on the build**: from the
*same* source, *same* starting values and *same* data, the **Windows** wheel (MSVC
`/O2` + vcpkg GSL) settles on a **spurious optimum** (`μ̂=−0.144`, `σ_a=0.305`,
`AIC=−2511`, ~52 lower log-likelihood) while the **Linux** wheel (manylinux, GCC
`-O2` + system GSL) reaches the correct one (`μ̂=0.0021`, `σ_a=0.261`, `AIC=−2613`,
paper Table 2). In **both** cases fue reports **`converged=True`, `ifault=0` with no
diagnostic** — the absurd mean and the worse `σ_a`/AIC are not flagged.

Two things, then:
1. **Non-reproducibility** of the optimizer across builds on multimodal surfaces
   (the basin depends on last-ULP floating-point differences).
2. **No guard** against a clearly-degenerate optimum being reported as success.

## Impact

Silent, wrong "estimated" model on the affected platform — the estimate is not
reproducible across OSes for strongly-multimodal specs, and there is no signal that
anything went wrong. Surfaced in the DVR replication of Garcia-Hiernaux et al. (2026):
the US CPI row of Table 2 came out degenerate on Windows, corrupting everything
downstream of the US residuals (DVR, GARCH/GJR, Beveridge–Nelson). The consumer had
to hard-code per-country starting values, and drtran's exact-VARMA engine was used to
cross-check.

## Reproduction

Self-contained in `art-python/bugs/BUG-0006-repro/` (`US_CPI.pre`, n=293; `repro.py`):

| build | default seed → | verdict |
|---|---|---|
| Windows — fue 0.1.7 wheel (MSVC/vcpkg) | μ=−0.144, σ=0.305, AIC=−2511 | **SPURIOUS** (converged=True) |
| Linux — fue 0.1.7 wheel (manylinux)    | μ=+0.0021, σ=0.261, AIC=−2613 | correct |
| Linux — fue 0.1.7 editable (source)    | μ=+0.0021, σ=0.261, AIC=−2613 | correct |

Verified on Linux in an isolated venv with the exact PyPI stack (`fue==0.1.7`,
`art-tseries==0.1.2`): the bug does NOT reproduce — same source, same seed, correct
result. It reproduces on the reporter's Windows wheel. That a repro labelled
"fue≥0.1.7" fails to reproduce with fue 0.1.7 on Linux **is itself the finding**: it
is a build/platform-dependent numerical issue of the optimizer, not a logic bug.

## Root cause

The ML optimizer (`qnewtopt`/BFGS in `csrc/internal/`) is a gradient-based local
search. On a **unimodal** surface it lands on the same optimum regardless of last-bit
differences; on the US-CPI **multimodal** AR(2)×AR(2) surface the basin is decided by
tiny floating-point differences in the log-likelihood evaluation. Those differences
come from the **build**: MSVC vs GCC FP semantics (operation reordering / FMA
contraction under `/O2` vs `-O2`), the C runtime `libm` (MSVC CRT vs glibc:
`exp`/`log`/`cos`/`sin`/`pow`, used heavily by the log-likelihood, Box-Cox and the
harmonics), and the GSL build (vcpkg static vs system). No `long double` is involved.
Separately, fue never checks whether the "converged" optimum is sane, so a degenerate
basin is reported as success.

A wrong-sign seasonal-AR *seed* (ART BUG-0006, now fixed on the ART side) is what
first pushed the search toward the wrong basin; a correct seed makes this case robust.
But the underlying fue fragility remains: a sufficiently multimodal spec could still
diverge on some build, and would still be reported as `converged=True`.

## Fix

Not yet applied. Robustness, not bit-identical FP (which is infeasible across
compilers):
- **Guard against absurd optima.** After a "converged" fit, check that `μ̂` lies
  within a few sample-std of the differenced-series mean, and/or that the objective
  is not materially worse than at the starting values; if it fails, downgrade
  `converged`, raise `ifault`, or emit a warning instead of reporting silent success.
- **Multi-start for multimodal blocks.** When a seasonal AR/MA block is present,
  optimise from a small set of starting points (e.g. ±seed, HR/YW seed) and keep the
  best optimum — platform-independent (the best basin wins on any build).
- *(Optional, palliative)* homogenise FP across builds: try `/fp:precise` (MSVC) and
  `-ffp-contract=off` (GCC) so the wheels agree; does not remove the multimodality.

## Validation

When fixed, `repro.py` must reach the correct optimum (`σ_a≈0.261`, `μ̂≈0.0021`,
`AIC≈−2613`) on **every** platform from the default seed, or fail loudly rather than
report `converged=True` on the spurious basin. The other seven DVR economies (unimodal
or non-seasonal-AR) must be unchanged.
