"""fue.scale: the rescaling factor's one rule for the suite (BUG-0028 and the
scale study, atsw-gui docs/ESTUDIO-escala.md)."""
import math

import numpy as np
import pytest

import fue
from fue import scale
from fue.forecast import _boxcox, _inv_boxcox


def _serie(nivel=100.0, sd=0.01, n=160, seed=3):
    rng = np.random.default_rng(seed)
    e = np.zeros(n)
    for t in range(1, n):
        e[t] = 0.6 * e[t - 1] + rng.normal(0, sd)
    return nivel * np.exp(e)


@pytest.mark.parametrize("lam", [-1.0, -0.5, 0.0, 0.333333, 0.5, 1.0, 2.0])
def test_boxcox_is_the_engines(lam):
    y = _serie()
    a = scale.boxcox(y, lam, 100.0)
    assert a == pytest.approx([_boxcox(float(v), lam, 100.0) for v in y], rel=1e-14)
    assert scale.boxcox(float(y[5]), lam, 7.0) == _boxcox(float(y[5]), lam, 7.0)
    assert scale.inv_boxcox(a, lam, 100.0) == pytest.approx(y, rel=1e-10)
    assert scale.inv_boxcox(float(a[5]), lam, 100.0) == pytest.approx(
        _inv_boxcox(float(a[5]), lam, 100.0), rel=1e-14)


def test_lambda_one_is_y_not_y_minus_one():
    assert scale.boxcox(5.0, 1.0) == 5.0
    assert fue.boxcox is scale.boxcox and fue.inv_boxcox is scale.inv_boxcox


def test_display_scale():
    assert scale.display_scale(0.0, 100.0) == (1.0, "%", True)
    assert scale.display_scale(0.0, 1.0) == (100.0, "%", True)
    assert scale.display_scale(1.0, 100.0) == (0.01, "series", False)
    assert scale.display_scale(-1.0, 10.0) == (0.1, "transformed", False)
    assert scale.param_scale(100.0) == 0.01
    assert scale.cross_scale(100.0, 100.0) == 1.0
    assert scale.cross_scale(100.0, 1.0) == 0.01          # Y at 100, X at 1


def test_format_refactor_keeps_100_and_does_not_lose_small_ones():
    assert scale.format_refactor(100.0) == "100.00"
    assert scale.format_refactor(1.0) == "1.00"
    assert scale.format_refactor(0.01) == "0.01"
    assert float(scale.format_refactor(0.001)) == 0.001   # two decimals wrote 0.00
    assert float(scale.format_refactor(1e6)) == 1e6


def test_a_small_factor_survives_the_pre(tmp_path):
    y = 1e6 + 1e4 * (_serie(1.0) - 1.0) * 100
    ts = fue.TimeSeries(list(map(float, y)), freq=12, start=(2000, 1), name="Y")
    m = fue.Model(ts, d=1, D=0, boxlam=1.0, ar=[[0.1]], ar_free=[[True]], refactor=0.001)
    m.fit()
    p = tmp_path / "y.pre"
    m.write_pre(str(p))
    _, m2 = fue.load(str(p))
    assert m2.refactor == 0.001


@pytest.mark.parametrize("nivel,lam,d,ok,rec", [
    (100.0, 0.0, 0, True, 100.0),         # logs at 100: the norm
    (100.0, 0.0, 1, True, 100.0),
    (100.0, 1.0, 0, True, 100.0),         # levels about 100: 122, inside
    (1e4, 1.0, 0, False, 1.0),            # 1.2e4 -> 122 at factor 1
    (1e4, -1.0, 1, False, 1e6),           # 1e-4 -> ~1 at 1e6
])
def test_check_scale(nivel, lam, d, ok, rec):
    c = scale.check_scale(_serie(nivel), lam, d=d, refactor=100.0)
    assert c.ok is ok
    assert c.recommended == rec
    assert scale.SPREAD_BAND[0] <= c.spread / 100.0 * rec <= scale.SPREAD_BAND[1]


def test_the_base_files_factor_is_kept_while_it_works():
    y = _serie(100.0)
    assert scale.recommended_refactor(y, 0.0, d=1, current=10.0) == 10.0   # 0.1: inside
    assert scale.recommended_refactor(y, 0.0, d=1, current=1.0) == 100.0   # 0.01: outside


def test_the_ratio_is_a_problem_no_factor_mends():
    y = _serie(1e6)                        # lambda = -1, d = 0: ~1 +- 1e-8
    c = scale.check_scale(y, -1.0, d=0, refactor=1e6)
    assert not c.ok and "ratio" in [k for k, _ in c.problems]
    assert c.ratio > scale.RATIO_MAX
    assert scale.check_scale(y, -1.0, d=1, refactor=1e6).ratio == 0.0     # differenced


def test_no_admissible_factor_says_so():
    y = _serie(1e12)                       # levels of 1e12: 1e10 at factor 1
    c = scale.check_scale(y, 1.0, d=0, refactor=100.0)
    assert c.recommended is None
    assert "other units" in c.problems[0][1]


@pytest.mark.parametrize("lam", [0.0, 1.0, -0.5])
def test_loglik_y_does_not_depend_on_the_factor(lam):
    y = _serie(100.0, sd=0.05, n=120)
    ts = fue.TimeSeries(list(map(float, y)), freq=12, start=(2000, 1), name="Y")
    out = []
    for f in (1.0, 100.0, 1000.0):
        w = scale.boxcox(y, lam, f)
        m = fue.Model(ts, d=0, D=0, boxlam=lam, ar=[[0.1]], ar_free=[[True]],
                      mu=float(np.mean(w)), estimate_mu=True, refactor=f)
        m.fit()
        out.append(scale.loglik_y(m.loglik, y, lam, f, n=len(m._result.residuals)))
    assert out[1] == pytest.approx(out[0], abs=2e-3)
    assert out[2] == pytest.approx(out[0], abs=2e-3)


def test_loglik_y_ranks_lambdas():
    y = _serie(100.0, sd=0.25, n=240, seed=2)       # generated in logs
    ts = fue.TimeSeries(list(map(float, y)), freq=12, start=(2000, 1), name="Y")
    ly = {}
    for lam, f in ((0.0, 100.0), (1.0, 1.0)):
        w = scale.boxcox(y, lam, f)
        m = fue.Model(ts, d=0, D=0, boxlam=lam, ar=[[0.1]], ar_free=[[True]],
                      mu=float(np.mean(w)), estimate_mu=True, refactor=f)
        m.fit()
        ly[lam] = scale.loglik_y(m.loglik, y, lam, f, n=len(m._result.residuals))
    assert ly[0.0] > ly[1.0] + 5.0
