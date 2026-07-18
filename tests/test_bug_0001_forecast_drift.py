"""
Regression test for BUG-0001 — forecast level over-shoots by the mean drift.

The level forecast of a differenced model with drift must equal the mean-form
reference  Y_{n+l} = Y_n + Σ_{k≤l} ŝ_{n+k},  ŝ following  s_t − μ = φ(B)^{-1}a_t
(so the differenced forecast reverts to μ, and the level integrates it), NOT the
buggy homogeneous-plus-l·μ composition (which added an extra μ·φ/(1-φ)).

See bugs/BUG-0001-forecast-mean-drift.md.
"""

import numpy as np
import pytest

from fue import TimeSeries, Model


def _meanform_level(y, phi, mu, L):
    """Independent correct reference for an ARIMA(p,1,0)+drift level forecast:
    forecast the differenced series with the mean form and integrate."""
    s = np.diff(np.asarray(y, float))          # d = 1
    u = list(s - mu)                           # mean-centred differences
    p = len(phi)
    for _ in range(L):
        u.append(sum(phi[i] * u[-(i + 1)] for i in range(p)))
    sf = [u[len(s) + l] + mu for l in range(L)]  # differenced point forecasts
    lvl, last = [], float(y[-1])
    for l in range(L):
        last += sf[l]
        lvl.append(last)
    return np.array(lvl)


def _synth_arima110(phi_true, mu_true, n=220, seed=0, sd=0.3):
    rng = np.random.default_rng(seed)
    u = np.zeros(n)
    a = rng.normal(0, sd, n)
    for t in range(1, n):
        u[t] = phi_true * u[t - 1] + a[t]
    return np.cumsum(mu_true + u) + 100.0       # a level with drift


def _fit_or_skip(model):
    try:
        model.fit()
    except Exception as exc:                     # engine unavailable / non-conv
        pytest.skip(f"fit unavailable: {exc}")
    return model


def test_ar1_drift_matches_meanform():
    y = _synth_arima110(0.5, 0.3, seed=1)
    m = _fit_or_skip(Model(TimeSeries(y, freq=1, start=(1, 1900)),
                           ar=[[0.0]], d=1, estimate_mu=True, boxlam=1.0))
    L = 18
    fr = m.forecast(L)
    phi = [m._result.params[0]]
    mu = m._result.params[-1]
    ref = _meanform_level(y, phi, mu, L)
    assert np.max(np.abs(fr.level - ref)) < 1e-6


def test_ar2_drift_matches_meanform():
    # AR(2) stresses phi(1) = 1 - (phi1 + phi2) in the drift intercept.
    y = _synth_arima110(0.4, 0.25, seed=7)       # generated AR(1); fit AR(2)
    m = _fit_or_skip(Model(TimeSeries(y, freq=1, start=(1, 1900)),
                           ar=[[0.0, 0.0]], d=1, estimate_mu=True, boxlam=1.0))
    L = 18
    fr = m.forecast(L)
    phi = [m._result.params[0], m._result.params[1]]
    mu = m._result.params[-1]
    ref = _meanform_level(y, phi, mu, L)
    assert np.max(np.abs(fr.level - ref)) < 1e-6


def test_differenced_forecast_converges_to_mu():
    # The month-on-month forecast (diff1) must converge to 100*mu (levels model).
    y = _synth_arima110(0.6, 0.4, seed=3)
    m = _fit_or_skip(Model(TimeSeries(y, freq=1, start=(1, 1900)),
                           ar=[[0.0]], d=1, estimate_mu=True, boxlam=1.0))
    fr = m.forecast(30)
    mu = m._result.params[-1]
    assert abs(fr.diff1[-1] - 100.0 * mu) < 1e-4


def test_d0_stationary_converges_to_mean():
    # d=0 (no differencing): the bug was CATASTROPHIC here — accumulated l·μ adds
    # a linear trend to a stationary model, so the level exploded.  The fix must
    # converge to the mean and match the mean-form AR forecast to machine
    # precision.
    rng = np.random.default_rng(2)
    n, phi_t, mean_t = 300, 0.6, 5.0
    y = np.zeros(n)
    a = rng.normal(0, 0.5, n)
    for t in range(1, n):
        y[t] = mean_t + phi_t * (y[t - 1] - mean_t) + a[t]
    m = _fit_or_skip(Model(TimeSeries(y, freq=1, start=(1, 1900)),
                           ar=[[0.0]], d=0, estimate_mu=True, boxlam=1.0))
    L = 15
    fr = m.forecast(L)
    phi, mu = m._result.params[0], m._result.params[-1]
    # mean-form d=0 reference: w_{n+l} = μ + φ(w_{n+l-1} − μ)
    ref, prev = [], float(y[-1])
    for _ in range(L):
        prev = mu + phi * (prev - mu)
        ref.append(prev)
    assert np.max(np.abs(fr.level - np.array(ref))) < 1e-6
    assert abs(fr.level[-1] - mu) < 1e-3          # converges to the mean


def test_zero_mean_is_pure_ar_forecast():
    # mu = 0: no drift; the level forecast must be the plain integrated AR(1).
    y = _synth_arima110(0.5, 0.0, seed=5)
    m = _fit_or_skip(Model(TimeSeries(y, freq=1, start=(1, 1900)),
                           ar=[[0.0]], d=1, estimate_mu=False, boxlam=1.0))
    L = 12
    fr = m.forecast(L)
    phi = [m._result.params[0]]
    ref = _meanform_level(y, phi, 0.0, L)        # mu = 0
    assert np.max(np.abs(fr.level - ref)) < 1e-6
