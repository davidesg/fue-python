---
id: BUG-0003
title: plot_residuals_ts draws no year ticks or vertical dividers for annual series (freq==1) — X-axis unreadable
status: fixed
severity: medium
component: plots
found_in: 0.1.6
fixed_in: 0.1.7
reported: 2026-07-18
reporter: D. E. Guerrero
tags:
  - plots
  - matplotlib
  - residuals
  - annual
references:
  - src/fue/plots.py:65 (plot_residuals_ts)
  - src/fue/plots.py:105-118 (tick/divider block gated on freq > 1)
  - /home/david/Dropbox/Cycles/bugs_art_fue.md (ISSUE #3)
---

## Summary

`plot_residuals_ts` builds the year ticks and vertical dividers **only** inside an
`if freq > 1:` branch (monthly/quarterly).  For **annual series (freq == 1)** the
block is skipped entirely: `set_xticks` is never called and no vertical gridlines
are drawn, so the decimal-year x-axis falls back to matplotlib's default locator
and the year labels run together into an unreadable strip
(e.g. `17707707807908008108208308408508608708808909…`).  The ACF/PACF panels are
unaffected — only the residual time-series panel.

## Impact

The residual figure is unusable for annual series — precisely the long climate
series in the *Joseph's Cycles* project (Geneva precipitation days n≈248, GEP,
Zurich).  This function is what the ART MCP delegates to for its residual panel,
so the defect surfaces there too.  The target is the fue (C) PDF layout
(`Analysis/idem/GEP/GEP.2.pdf`, `Analysis/idem/GE/GE.3.pdf`): year ticks every
~20 years with vertical gridlines, legible.

## Reproduction

```python
import numpy as np
from fue.plots import plot_residuals_ts
r = np.random.default_rng(0).normal(0, 1, 248)   # 248 annual obs
plot_residuals_ts(r, title="annual residuals")   # x-axis labels overlap
```

Compare with a `freq=12` model, whose branch produces clean 2-year dividers.

## Root cause

`src/fue/plots.py`, `plot_residuals_ts`: the tick-position / vertical-divider loop
(lines ~105–118) is guarded by `if freq > 1:`.  `freq == 1` produces no
`tick_pos`, so `set_xticks` is skipped and the axis is left to the default
auto-locator with per-observation decimal-year values.

## Fix

**Applied** in 0.1.7 (`src/fue/plots.py`, `plot_residuals_ts`).  Added an annual
(`freq == 1`) branch to the tick block, replicating fue (C)
`gnuplot_File_PlotSer_CorrSer`: year labels every **20 years anchored at the begin
year** (fue C step `2*f*10 = 20`, label `tsby + 20*i` → 1768, 1788, …, 2008), with
a vertical divider at each labelled year except the first (fue C draws `set arrow`
from `i=1` while labels start at `i=0`).  This deliberately anchors to the series
start, **not** to round centuries, so the panel matches fue (C) exactly.  The
`freq > 1` path is unchanged (2-year seasonal dividers, fue C step `2*f`).

## Validation

- fue (C) parity: for a Geneva-like annual series (start 1768, 248 obs) the ticks
  are `1768, 1788, …, 2008` — identical values and spacing to
  `gnuplot_File_PlotSer_CorrSer` (`tsby + 20*i`).
- Regression test `tests/test_bug_0003_residual_plot_annual.py`: for `freq == 1`
  the x-ticks are non-empty, anchored at the begin year, spaced exactly 20 years,
  and within the sample span (13 ticks, not one-per-observation); the `freq > 1`
  path still produces year dividers.
