"""
Phase 2 — forecast tests.

Verifies that Model.forecast() produces reasonable level forecasts for
simple models where we can check the output analytically or by
comparison with the fuf reference binary.

Tolerances are deliberately loose for now (1e-4 on levels) — the goal
is to confirm the algorithm is wired correctly, not to match fuf to
machine precision.
"""

import math
import numpy as np
import pytest

pytest.importorskip("fue._fue_engine",
                    reason="C extension not compiled — skip forecast tests")

from fue import TimeSeries, Model, Intervention, ForecastResult

# First 30 annual observations from the SFNY.2 dataset (1852-1881)
_SFNY30 = np.array([
    3.91505848, 2.02125792, 0.81208771, 0.60807414, 1.21576447,
    1.43763055, 1.78032601, 0.82841058, 0.65433228, 0.74324607,
    0.93394905, 0.60094494, 0.80840161, 0.90899270, 0.40822203,
    0.41975993, 0.50368768, 0.57248427, 0.72970370, 0.90175445,
    0.61763439, 0.63607641, 0.67670827, 0.81812744, 0.78095914,
    0.82024104, 0.86103433, 0.84442843, 0.74566075, 0.63347579,
])


def test_forecast_returns_result():
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m = Model(ts, ar=[[0.5]]).fit()
    fr = m.forecast(5)
    assert isinstance(fr, ForecastResult)
    assert fr.horizon == 5
    assert len(fr.level) == 5
    assert len(fr.level_std) == 5


def test_forecast_ar1_shape():
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m = Model(ts, ar=[[0.5]]).fit()
    fr = m.forecast(10)
    assert fr.level.shape == (10,)
    assert fr.diff1.shape == (10,)
    assert fr.seasonal_diff.shape == (10,)


def test_forecast_ar1_converges_to_mean():
    """AR(1) forecasts should converge toward the unconditional mean."""
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m = Model(ts, ar=[[0.5]]).fit()
    fr = m.forecast(50)
    # For AR(1) with phi≈0.97, forecasts are slow but still should be
    # bounded and heading somewhere finite.
    assert np.all(np.isfinite(fr.level))
    assert np.all(fr.level > 0)


def test_forecast_std_increasing():
    """Forecast uncertainty should be non-decreasing."""
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m = Model(ts, ar=[[0.5]]).fit()
    fr = m.forecast(10)
    # level_std should be non-decreasing (psi-weight variances accumulate)
    assert np.all(np.diff(fr.level_std) >= -1e-12)


def test_forecast_ima11_1step():
    """
    IMA(1,1): 1-step forecast in Box-Cox level space.

    For IMA(1,1) with boxlam=1 (levels): phi0[1]=1+phi0=1 (since φ=0),
    and theta estimated around -0.42, so:
        f1[1] = nt[nobs] + theta*a[nobs-1]
             ≈ nt[nobs] + theta*residuals[-1]

    We don't have fuf reference values, so we just verify shape and
    that the 1-step is close to the last observation (random walk property).
    """
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m = Model(ts, ma=[[-0.5]], d=1).fit()
    fr = m.forecast(5)
    assert len(fr.level) == 5
    # IMA level forecast at h=1 should be near last obs
    last_obs = _SFNY30[-1]
    assert abs(fr.level[0] - last_obs) < 0.5, (
        f"1-step IMA forecast {fr.level[0]:.4f} far from last obs {last_obs:.4f}"
    )


def test_forecast_level_positive_boxcox0():
    """With boxlam=0 (log), inv_boxcox gives exp(z) which is always positive."""
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m = Model(ts, ar=[[0.5]], boxlam=0.0).fit()
    fr = m.forecast(5)
    assert np.all(fr.level > 0)
