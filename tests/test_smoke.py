"""
Fast smoke test for the built wheel — run by cibuildwheel's ``test-command`` on
every platform/arch and by the pure-Python wheel job.

It verifies that the installed package imports and estimates/forecasts a small
model with a finite, converged result.  It deliberately does NOT check exact
golden log-likelihoods: those are platform/BLAS-sensitive for ill-conditioned
models (e.g. the multimodal cointegration case R.4) and belong in the dev test
job on a fixed platform, not in a per-wheel gate on 16 build targets.
"""

import numpy as np

import fue


def test_import_and_version():
    assert isinstance(fue.__version__, str) and fue.__version__


def test_fit_and_forecast_small_model():
    # ARIMA(1,1,0) + drift on a synthetic level: must converge with a finite
    # log-likelihood and produce finite level forecasts (exercises the estimation
    # engine — C extension if the wheel is binary, pure-Python otherwise).
    rng = np.random.default_rng(0)
    y = np.cumsum(0.2 + rng.normal(0.0, 1.0, 150)) + 100.0
    ts = fue.TimeSeries(y, freq=1, start=(1, 1900))
    m = fue.Model(ts, ar=[[0.0]], d=1, estimate_mu=True, boxlam=1.0)
    m.fit()
    r = m._result
    assert r.ifault == 0
    assert np.isfinite(r.loglik)
    assert r.sigma2 > 0

    fr = m.forecast(6)
    assert len(fr.level) == 6
    assert np.all(np.isfinite(fr.level))
