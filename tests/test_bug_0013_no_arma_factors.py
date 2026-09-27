"""A model with no ARMA factor at all must estimate, not kill the interpreter.

`bugs/BUG-0013`: a regression on deterministic inputs with ARMA(0,0) errors —
harmonics and white noise, the first rung of every seasonal ladder — segfaulted
the C engine when the model was built through the Python API. No exception, no
message: the process disappeared.

Two things make this test worth more than a crash guard:

* it checks the answer, not only the absence of a crash. Until 0.1.16
  `_engine.estimate` diverted the case to the pure-Python engine; since 0.1.17
  the C itself handles it (nlatools.c accepts the 0 x 0 matrices elf() uses
  when max(p,q) = 0), and the diversion is gone;
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

    Both now travel through the C engine, so they must agree to the last digit
    (while the first was diverted to the Python engine the measured difference
    was 1.9e-07, the cross-engine agreement of `docs/PERFORMANCE.md`).
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

    assert sin_factor.loglik == pytest.approx(con_factor.loglik, abs=1e-9)
    assert sin_factor._result.npar == con_factor._result.npar
    np.testing.assert_allclose(sin_factor._result.params, con_factor._result.params,
                               rtol=0, atol=1e-9)


def test_the_c_engine_handles_it_without_the_python_engine(monkeypatch):
    """The root fix, not a route around it: with the Python engine made
    unavailable, `_engine.estimate` still fits the model with no ARMA factor."""
    pytest.importorskip("fue._fue_engine")
    import fue.cast_us

    def _no(*a, **k):
        raise AssertionError("diverted to the Python engine")
    monkeypatch.setattr(fue.cast_us, "estimate_py", _no)

    ts, y = _serie()
    m = fue.Model(ts, d=0, interventions=_armonicos(),
                  mu=float(y.mean()), estimate_mu=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit()
    assert m._result.ifault == 0 and np.isfinite(m.loglik)
    # with no ARMA, the SE of mu has the closed form sigma/sqrt(n)
    n = len(m._result.residuals)
    assert m._result.std_errors[-1] == pytest.approx(
        np.sqrt(m._result.sigma2 / n), rel=1e-3)
