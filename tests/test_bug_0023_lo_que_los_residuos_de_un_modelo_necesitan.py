"""
Regression tests for BUG-0023 — the residual figure of a MODEL got three things
wrong that are not in the residuals themselves: the Ljung-Box df (it subtracted
EVERY parameter, harmonics included), the date of the first residual (it ignored
the differencing) and the default number of lags (not fug C's rule).

The four .pre files are real models of the SF_MEG run 3 (Spanish CPI):

    A_m00   10 harmonics, AR(1) FIXED at 0    → 0 estimated ARMA parameters
    A_m02   10 harmonics, MA(1) free          → 1
    A_m06   AR(1) + fixed-frequency MA(2)     → 2, and one ifadf root
    B_m11   d=1 and two ifadf roots           → offset 5

See bugs/BUG-0023-*.md.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

import fue
from fue.diagnostics import (default_lags, free_arma_count,
                             differencing_offset, residuals_start)
from fue.plots import plot_model_diagnostics, _obs_to_decimal_year

DATA = Path(__file__).parent / "data" / "bug_0023"


def _load(stem):
    out = fue.load(str(DATA / f"ES_CPI_{stem}.pre"))
    m = next(x for x in (out if isinstance(out, tuple) else (out,))
             if hasattr(x, "fit"))
    m._inp_stem = f"ES_CPI_{stem}"
    return m


# ── fug C's default_lags, case by case (diagnose.c) ─────────────────────────

@pytest.mark.parametrize("nobs,freq,lags", [
    (216, 12, 39),    # monthly: 3·(f+1)
    (60,   4, 15),    # quarterly: 3·(f+1)
    (100,  1,  9),    # annual
    (250,  1, 45),    # annual, long: the branch fue's plots had lost
    (10,   4,  8),    # short: nobs − f/2
    (5,    1,  3),    # short, capped at nobs − 2
    (3,   12,  1),    # never below 1
])
def test_default_lags_is_fug_c(nobs, freq, lags):
    assert default_lags(nobs, freq) == lags


# ── the ARMA parameters that were ESTIMATED ──────────────────────────────────

@pytest.mark.parametrize("stem,n_arma", [
    ("A_m00", 0),     # the AR(1) is fixed at 0: it estimates nothing
    ("A_m02", 1),
    ("A_m06", 2),     # the fixed-frequency MA(2) factor counts
])
def test_free_arma_count(stem, n_arma):
    assert free_arma_count(_load(stem)) == n_arma


# ── the observations the differencing consumes ───────────────────────────────

@pytest.mark.parametrize("stem,offset", [
    ("A_m00", 1),     # d=1
    ("A_m06", 3),     # d=1 + one interior ifadf root (2)
    ("B_m11", 5),     # d=1 + two interior ifadf roots (2 + 2)
])
def test_differencing_offset_matches_the_residuals(stem, offset):
    m = _load(stem)
    m.fit()
    assert differencing_offset(m) == offset
    assert len(m._result.residuals) == len(m.series.data) - offset


# ── the figure ───────────────────────────────────────────────────────────────

@pytest.mark.parametrize("stem,label", [
    ("A_m00", "Q(39)"),   # was Q(28): 39 − 11 parameters, harmonics included
    ("A_m02", "Q(38)"),   # was Q(27)
])
def test_the_q_label_carries_the_arma_df(stem, label):
    m = _load(stem)
    m.fit()
    fig, fig_h = plot_model_diagnostics(m)
    try:
        ax_acf = fig.axes[1]
        assert ax_acf.get_xlabel().startswith(label + " = ")
    finally:
        plt.close(fig)
        plt.close(fig_h)


def test_the_first_residual_is_dated_after_the_differencing():
    m = _load("B_m11")
    m.fit()
    fig, fig_h = plot_model_diagnostics(m)
    try:
        x0 = fig.axes[0].lines[0].get_xdata()[0]
        y, p = residuals_start(m)
        assert x0 == pytest.approx(_obs_to_decimal_year(1, y, p, m.series.freq)[0])
        # and it is NOT the series start (the bug)
        y0, p0 = m.series.start
        assert x0 != pytest.approx(_obs_to_decimal_year(1, y0, p0, m.series.freq)[0])
    finally:
        plt.close(fig)
        plt.close(fig_h)


def test_df_le_zero_is_not_a_test():
    from fue.plots import plot_acf_pacf
    import numpy as np
    r = np.random.default_rng(0).standard_normal(120)
    fig = plot_acf_pacf(r, npar=12, freq=1, lags=9)   # df = 9 − 12 < 1
    try:
        lbl = fig.axes[0].get_xlabel()
        assert lbl.startswith("Q = ") and "no test" in lbl
    finally:
        plt.close(fig)
