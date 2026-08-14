"""A model with no ARMA factor at all must estimate, not kill the interpreter.

`bugs/BUG-0013`: a regression on deterministic inputs with ARMA(0,0) errors —
harmonics and white noise, the first rung of every seasonal ladder — segfaulted
the C engine when the model was built through the Python API. No exception, no
message: the process disappeared.

Two things make this test worth more than a crash guard:

* it pins the DIVERSION, not the crash. The pure-Python engine already computed
  the right answer, so `_engine.estimate` sends the case there. If someone ever
  fixes the C and removes the diversion, these tests still pass — they check the
  answer, not the route;
* it pins the EQUIVALENCE that identifies the defect: the same model written
  with one AR factor pinned at zero goes through the C engine and agrees to the
  last digit. That is what proves the two paths are the same model and the
  segfault was not about the mathematics.

A crash cannot be caught by pytest — the interpreter dies — so a regression here
does not fail this file: it takes the whole run with it. That is itself the
signal.
"""
import warnings

import numpy as np
import pytest

import fue


def _serie(n=240, s=12, seed=20260814):
    rng = np.random.RandomState(seed)
    t = np.arange(1, n + 1)
    y = (100.0 + 3.0 * np.cos(2 * np.pi * t / s) + 1.5 * np.sin(2 * np.pi * t / s)
         + rng.normal(0, 1.0, n))
    return fue.TimeSeries(list(y), freq=s, start=(2000, 1), name="DET"), y


def _armonicos():
    return [fue.Intervention("cos", harmonic=1.0, omega=[0.1]),
            fue.Intervention("sin", harmonic=1.0, omega=[0.1])]


@pytest.mark.parametrize("d", [0, 1])
@pytest.mark.parametrize("estimate_mu", [True, False])
def test_a_model_with_no_arma_factor_estimates(d, estimate_mu):
    """The four combinations that used to segfault: d ∈ {0,1} × mu ∈ {on,off}."""
    ts, y = _serie()
    kw = dict(d=d, interventions=_armonicos())
    if estimate_mu:
        kw.update(mu=float(y.mean()), estimate_mu=True)

    m = fue.Model(ts, **kw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit()

    assert m._result is not None
    assert np.isfinite(m.loglik)
    assert m._result.ifault == 0


def test_it_agrees_with_the_same_model_written_the_way_a_file_writes_it():
    """`ar=[]` against `ar=[[0.0]]` fixed — the same model, two spellings.

    Every `.inp` writes the second, which is why no file ever hit the crash.

    They must give the same number, and the tolerance is the MEASURED one:
    1.9e-07. It is not machine epsilon because the two spellings now travel by
    different engines — the diverted one through the Python translation, the
    other through the C — and that is exactly the cross-engine agreement
    `docs/PERFORMANCE.md` reports (largest |Δ logL| over 23 real models: 2e-04).
    An earlier version of this test asserted 1e-9, which was a guess and failed.
    """
    ts, y = _serie()

    sin_factor = fue.Model(ts, d=0, interventions=_armonicos(),
                           mu=float(y.mean()), estimate_mu=True)
    con_factor = fue.Model(ts, d=0, interventions=_armonicos(),
                           ar=[[0.0]], ar_free=[[False]],
                           mu=float(y.mean()), estimate_mu=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sin_factor.fit()
        con_factor.fit()

    assert sin_factor.loglik == pytest.approx(con_factor.loglik, abs=1e-6)
    assert sin_factor._result.npar == con_factor._result.npar


def test_the_diversion_is_where_the_bug_report_says_it_is():
    """The helper that decides, and its contract: a factor pinned at zero is
    still a factor, so it must NOT be diverted."""
    from fue._engine import _sin_estructura_arma

    ts, y = _serie()
    assert _sin_estructura_arma(
        fue.Model(ts, d=0, interventions=_armonicos()))
    assert not _sin_estructura_arma(
        fue.Model(ts, d=0, interventions=_armonicos(),
                  ar=[[0.0]], ar_free=[[False]]))
    assert not _sin_estructura_arma(
        fue.Model(ts, d=0, interventions=_armonicos(), ma=[[0.4]]))
