"""
Regression test for BUG-0003 — annual residual plot has an unreadable x-axis.

plot_residuals_ts gated its year-tick / vertical-divider block on ``freq > 1``,
so annual series (freq == 1) got no ticks and the decimal-year labels merged into
a strip.  The fix sub-samples year ticks to a round decade/score step for
freq == 1.  See bugs/BUG-0003-residual-plot-annual-xaxis.md.
"""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt

from fue.plots import plot_residuals_ts


class _TS:
    def __init__(self, freq, start):
        self.freq = freq
        self.start = start


class _M:
    def __init__(self, freq=1, start=(1768, 1), refactor=100.0):
        self.series = _TS(freq, start)
        self.refactor = refactor


def test_annual_residual_plot_subsamples_year_ticks():
    # fue (C) scheme: labels every 20 years anchored at the begin year
    # (tsby + 20*i), not round centuries.  Geneva-like series 1768..2015.
    r = np.random.default_rng(0).normal(0, 1, 248)
    fig, ax = plt.subplots()
    plot_residuals_ts(r, model=_M(freq=1, start=(1768, 1)), title="annual", ax=ax)
    ticks = [int(t) for t in ax.get_xticks()]
    plt.close(fig)

    assert ticks                                       # not empty (the bug)
    assert ticks[0] == 1768                             # anchored at the begin year
    assert all((t - 1768) % 20 == 0 for t in ticks)     # exactly every 20 years
    assert all(1768 <= t <= 2015 for t in ticks)       # within the sample span
    assert len(ticks) == 13                             # 1768,1788,...,2008


def test_seasonal_residual_plot_still_has_year_dividers():
    r = np.random.default_rng(1).normal(0, 1, 120)     # 10 years monthly
    fig, ax = plt.subplots()
    plot_residuals_ts(r, model=_M(freq=12, start=(2000, 1)), title="monthly", ax=ax)
    ticks = list(ax.get_xticks())
    plt.close(fig)
    assert ticks                                       # freq>1 path unchanged
