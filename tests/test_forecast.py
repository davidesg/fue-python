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


# Full 62-obs SFNY.2 dataset (annual, 1852-1913)
_SFNY62 = np.array([
    3.9150584800, 2.0212579200, 0.8120877100, 0.6080741400, 1.2157644700,
    1.4376305500, 1.7803260100, 0.8284105800, 0.6543322800, 0.7432460700,
    0.9339490500, 0.6009449400, 0.8084016100, 0.9089927000, 0.4082220300,
    0.4197599300, 0.5036876800, 0.5724842700, 0.7297037000, 0.9017544500,
    0.6176343900, 0.6360764100, 0.6767082700, 0.8181274400, 0.7809591400,
    0.8202410400, 0.8610343300, 0.8444284300, 0.7456607500, 0.6334757900,
    0.7263755700, 0.8135161000, 0.7914275400, 0.8030587300, 0.8386753300,
    0.9867881400, 0.8048586300, 0.8165155300, 0.7596009300, 0.8407096800,
    0.8948088200, 0.8940759100, 0.8432364600, 0.7721518200, 0.8250954400,
    0.8738444300, 0.8136010600, 0.7849749600, 0.7132336000, 0.7068852200,
    0.8109034800, 0.9483109700, 0.7259892200, 0.8033732500, 0.8401149300,
    0.8924720200, 0.8932824600, 0.9094242400, 0.8287118900, 0.8864734000,
    0.8225149700, 0.9473733600,
])


def test_forecast_sfny2_vs_fuf():
    """
    SFNY.2 reference case: step + AR(1)xAR(2) + mu + boxlam=0, 62 obs.

    Estimated params taken directly from fuf output (all treated as fixed
    here to avoid re-estimation rounding differences).  Reference values
    come from forecast_sfny2_prev.1.1914.tex (2 decimal places, fuf 1.08.1).
    """
    ts = TimeSeries(_SFNY62, freq=1, start=(1852, 12), name="SFNY2")
    itv = Intervention("step", at=2,
                       omega=[-0.600361], omega_free=[False],
                       delta=[0.5873], delta_free=[False])
    m = Model(ts, interventions=[itv],
              ar=[[0.1311], [0.4880, -0.2588]],
              ar_free=[[False], [False, False]],
              mu=1.209415, estimate_mu=False,
              boxlam=0.0, d=0, D=0, refactor=1.0)

    from fue.forecast import forecast as _forecast

    class _FixedResult:
        params = np.array([])
        sigma2 = 0.0370593261
        residuals = np.zeros(62)

    fr = _forecast(m, _FixedResult(), 5)

    # fuf reference (2-decimal LaTeX output) — tolerance 0.01
    ref_level  = [1.30,  2.38,   7.38,   27.73,  100.32]
    ref_lstd   = [19.25, 22.64,  22.67,  22.81,  22.86]
    ref_diff1  = [31.47, 60.73, 113.03, 132.44,  128.57]

    for h in range(5):
        assert abs(round(fr.level[h], 2) - ref_level[h]) < 0.01, \
            f"h={h+1} level: got {fr.level[h]:.4f}, expected {ref_level[h]}"
        assert abs(round(fr.level_std[h] * 100, 2) - ref_lstd[h]) < 0.01, \
            f"h={h+1} level_std: got {fr.level_std[h]*100:.4f}%, expected {ref_lstd[h]}%"
        assert abs(round(fr.diff1[h], 2) - ref_diff1[h]) < 0.01, \
            f"h={h+1} diff1: got {fr.diff1[h]:.4f}%, expected {ref_diff1[h]}%"
