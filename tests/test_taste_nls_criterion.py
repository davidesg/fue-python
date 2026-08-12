"""fue and TASTE disagree on the airline model, and the reason is the criterion.

TASTE estimates by **nonlinear least squares with backforecasting** — the method
of Box & Jenkins (1976) §7.1.4 — and `fue` by **exact unconditional maximum
likelihood** (Mauricio 1995/1997). On a series whose seasonal MA sits near the
non-invertibility boundary the two criteria pull apart:

    (0,1,1)(0,1,1)_12 on the Spanish CPI, 2002-2019, n=216, λ=0, ×100

        fue   (exact ML)        θ = -0.421156   Θ = +0.814706
        TASTE (BJ NLS)          θ = -0.426710   Θ = +0.906060

**Neither is wrong: each sits at its own optimum.** On fue's exact likelihood
TASTE's point is 1.494 worse; on the backforecasting sum of squares fue's point
is 14.418 against TASTE's 14.186.

This is worth a test because the difference was first read as a defect. What the
oracle case `fue_es_cpi_airline` can validate for `fue` is therefore **not the
estimate** — it is that both programs are looking at the *same series through the
same operator*. That is what these tests check: the Box-Jenkins criterion,
evaluated on the series `fue` produces after ∇∇₁₂, lands on TASTE's published
estimates. If `fue`'s differencing ever drifted, this would move.

Two things ruled out along the way, recorded so nobody repeats them:

  * **Not convergence.** TASTE gives identical values at 60 and at 300 iterations.
  * **Not the plain conditional sum of squares** (pre-sample residuals set to
    zero), which optimises at Θ = +0.784781 — *below* fue, in the opposite
    direction from TASTE.
"""
import numpy as np
import pytest

_S = 12
_N_BACKCASTS = 13          # the MA operator (1−θB)(1−ΘB¹²) has degree 13

#: What TASTE reports, from the oracle battery.
TASTE = dict(theta=-0.426710, Theta=0.906060)
#: What fue's exact ML gives on the same file.
FUE = dict(theta=-0.421156, Theta=0.814706)


def _series():
    """∇∇₁₂ of 100·log(CPI) — the stationary series both programs work on."""
    import fue as _fue

    import os
    for cand in (
        os.path.expanduser("~/Dropbox/SRC/atws/Taste/oracle/data/sfmeg/SF_ES_CPI_airline.pre"),
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "..", "..", "ART", "art-python", "bugs", "BUG-0010-repro",
                     "IPC_ES_m10.pre"),
    ):
        if os.path.exists(cand):
            ts, _m = _fue.load(cand)
            z = np.log(np.asarray(ts.data, float)) * 100.0
            w = np.diff(z)
            return w[_S:] - w[:-_S]
    pytest.skip("the Spanish CPI fixture is not available")


def _resid(w, th, TH):
    """Forward residuals of (1−θB)(1−ΘB¹²), pre-sample set to zero."""
    a = np.zeros(len(w))
    for t in range(len(w)):
        v = w[t]
        if t >= 1:
            v += th * a[t - 1]
        if t >= _S:
            v += TH * a[t - _S]
        if t >= _S + 1:
            v -= th * TH * a[t - _S - 1]
        a[t] = v
    return a


def _backforecasts(w, th, TH, q=_N_BACKCASTS):
    """Box-Jenkins §7.1.4 backforecasting: the SERIES is forecast backwards.

    Reverse w, take its residuals, and forecast forward on the reversed series —
    which is backwards on the original. For a pure MA of degree 13 only thirteen
    backforecasts are non-zero.

    Backforecasting the *residuals* instead is the natural-looking mistake and
    gives Θ ≈ 0.71, further from both programs, not closer.
    """
    ar = _resid(w[::-1], th, TH)
    ae = np.concatenate([ar, np.zeros(q + _S + 2)])
    m = len(w)
    return np.array([-th * ae[m + l - 2] - TH * ae[m + l - 1 - _S]
                     + th * TH * ae[m + l - 2 - _S] for l in range(1, q + 1)])


def uss(w, th, TH):
    """Unconditional (backforecast) sum of squares — TASTE's criterion."""
    wb = _backforecasts(w, th, TH)
    return float(np.sum(_resid(np.concatenate([wb[::-1], w]), th, TH) ** 2))


def css(w, th, TH):
    """Conditional sum of squares — pre-sample residuals at zero."""
    return float(np.sum(_resid(w, th, TH) ** 2))


# ── the criterion explains the difference ──────────────────────────────────

def test_the_backforecasting_criterion_reproduces_taste():
    """The claim: fue's differenced series + Box-Jenkins NLS = TASTE."""
    from scipy.optimize import minimize

    w = _series()
    r = minimize(lambda p: uss(w, *p), [-0.42, 0.81], method="Nelder-Mead",
                 options=dict(xatol=1e-9, fatol=1e-11, maxiter=8000))
    theta, Theta = r.x
    assert theta == pytest.approx(TASTE["theta"], abs=1e-4), theta
    assert Theta == pytest.approx(TASTE["Theta"], abs=1e-3), Theta


def test_each_estimator_sits_at_its_own_optimum():
    """Neither program is wrong — they optimise different things.

    On the backforecasting sum of squares TASTE's point wins; on fue's exact
    likelihood fue's point wins by 1.494. That symmetry is the whole finding.
    """
    w = _series()
    assert uss(w, TASTE["theta"], TASTE["Theta"]) < uss(w, FUE["theta"], FUE["Theta"])


def test_the_conditional_criterion_is_not_the_explanation():
    """Ruled out: plain CSS goes the other way, to Θ ≈ 0.785."""
    from scipy.optimize import minimize

    w = _series()
    r = minimize(lambda p: css(w, *p), [-0.42, 0.81], method="Nelder-Mead",
                 options=dict(xatol=1e-9, fatol=1e-11, maxiter=8000))
    assert r.x[1] < FUE["Theta"], "CSS should fall below fue, not rise to TASTE"
    assert abs(r.x[1] - TASTE["Theta"]) > 0.1


# ── and what fue itself gives ──────────────────────────────────────────────

def test_fue_exact_ml_is_where_the_oracle_note_says():
    import fue as _fue
    import os

    p = os.path.expanduser(
        "~/Dropbox/SRC/atws/Taste/oracle/data/sfmeg/SF_ES_CPI_airline.pre")
    if not os.path.exists(p):
        pytest.skip("oracle fixture not available")
    ts, m = _fue.load(p)
    m.fit()
    assert m.ma[0][0] == pytest.approx(FUE["theta"], abs=5e-4)
    assert m.ma_s[0][0] == pytest.approx(FUE["Theta"], abs=5e-4)
