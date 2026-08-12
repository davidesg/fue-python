"""The published benchmark: Box & Jenkins Series A, IMA(0,1,1).

This closes the link `docs/PROVENANCE.md` §6 listed as missing. The rest of the
suite verifies **port against binary** — `test_estimation.py` says so itself:
*"hard-coded reference values obtained from the reference binary"*. That checks
the port and not the algorithm. This file checks the algorithm, against the
canonical published benchmark of the field and against an implementation that
shares no ancestry with any of it.

    Series A — chemical process concentration readings, every two hours,
    n = 197. Box & Jenkins (1976), Series A; the IMA(0,1,1) fit is the
    textbook example, with θ ≈ 0.70.

Measured 12 August 2026:

    Box & Jenkins (book)          θ = 0.70
    DRVUS  (Mauricio's C, ~2001)  θ = 0.699384   σ = 0.317382   logL = −53.508690
    fue    (Python, 2026)         θ = 0.699384   σ = 0.317382   logL = −53.508690
    statsmodels (exact MLE)       θ = 0.699384   σ = 0.317382   logL = −53.508686

`statsmodels` is the one that matters here: different authors, a different
algorithm (state-space/innovations), a different decade, and no shared code with
DRVUS, FUE or this port. Agreement to **4.5e-09** on θ is not a coincidence of
lineage.

The data ships in `tests/data/bj_series_a.txt`, taken from Mauricio's own DRVUS
1.2.01 fixture. The DRVUS reference values are from the `.out` that shipped with
it — a run of the original C, preserved since ~2001.
"""
import os

import numpy as np
import pytest

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "data", "bj_series_a.txt")

#: What Mauricio's own C produced on this series, from `SeriesA/a2.out`.
DRVUS = dict(theta=0.699384, se_theta=0.064511,
             sigma2=0.1007314876, sigma=0.3173822421, loglik=-53.5086902793)

#: The textbook value.
BOOK_THETA = 0.70


def _series():
    if not os.path.exists(_DATA):
        pytest.skip("Series A fixture missing")
    with open(_DATA) as fh:
        return np.array([float(l) for l in fh if l.strip()
                         and not l.startswith("#")])


def _fit_fue(y):
    import fue

    ts = fue.TimeSeries(list(y), freq=12, start=(1, 1), name="BJ_SeriesA")
    m = fue.Model(ts, d=1, D=0, boxlam=1.0,
                  ar=[], ma=[[0.53]], ma_free=[[True]],
                  ar_s=[], ma_s=[], interventions=[], ifadf=[0] * 7,
                  mu=0.0, estimate_mu=False, refactor=1.0)
    m.fit()
    return m


# ── the data itself ────────────────────────────────────────────────────────

def test_the_fixture_is_series_a():
    y = _series()
    assert len(y) == 197
    assert 16.0 < y.min() and y.max() < 18.5     # concentration readings


# ── against the original C ─────────────────────────────────────────────────

def test_fue_reproduces_mauricios_own_c():
    """The 2026 port against a run of the 2001 C preserved in its fixtures."""
    m = _fit_fue(_series())
    assert m.ma[0][0] == pytest.approx(DRVUS["theta"], abs=1e-5)
    assert m.loglik == pytest.approx(DRVUS["loglik"], abs=1e-4)


# ── against the published value ────────────────────────────────────────────

def test_fue_reproduces_the_textbook_estimate():
    """θ ≈ 0.70 — the number in Box & Jenkins for Series A."""
    m = _fit_fue(_series())
    assert m.ma[0][0] == pytest.approx(BOOK_THETA, abs=0.005)


# ── against an implementation that shares no ancestry ──────────────────────

def test_fue_agrees_with_statsmodels():
    """The check that is genuinely external.

    `statsmodels` computes the exact likelihood by a different route and owes
    nothing to Mauricio's code. If both land on the same optimum to 1e-8, the
    likelihood being maximised is the same function.

    Note the sign convention: `statsmodels` writes the MA as (1 + θB) and fue as
    (1 − θB), so the estimates are negatives of each other.
    """
    sm = pytest.importorskip("statsmodels.tsa.arima.model")

    y = _series()
    m = _fit_fue(y)
    r = sm.ARIMA(y, order=(0, 1, 1), trend="n").fit(method="innovations_mle")

    assert m.ma[0][0] == pytest.approx(-r.params[0], abs=1e-6)
    assert m.loglik == pytest.approx(r.llf, abs=1e-4)


def test_the_three_agree_with_each_other():
    """The whole chain in one assertion, and the tolerances are not decoration:
    fue↔DRVUS is 2.7e-07 and fue↔statsmodels 4.5e-09."""
    sm = pytest.importorskip("statsmodels.tsa.arima.model")

    y = _series()
    theta_fue = _fit_fue(y).ma[0][0]
    theta_sm = -sm.ARIMA(y, order=(0, 1, 1), trend="n").fit(
        method="innovations_mle").params[0]

    assert abs(theta_fue - DRVUS["theta"]) < 1e-5
    assert abs(theta_fue - theta_sm) < 1e-6
    assert abs(theta_fue - BOOK_THETA) < 0.005
