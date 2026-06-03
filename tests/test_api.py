"""
Smoke tests for the Python API layer (no C extension required).

Tests that verify the C engine results go in test_estimation.py
and are marked with pytest.mark.engine.
"""

import numpy as np
import pytest
from fue import TimeSeries, Intervention, Model


def _dummy_series(n=120, freq=12):
    rng = np.random.default_rng(42)
    return TimeSeries(rng.standard_normal(n), freq=freq, start=(1990, 1),
                      name="test")


# ── TimeSeries ────────────────────────────────────────────────────────────

def test_timeseries_basic():
    ts = _dummy_series()
    assert ts.nobs == 120
    assert ts.freq == 12
    assert ts.start == (1990, 1)


def test_timeseries_from_array():
    data = list(range(48))
    ts = TimeSeries.from_array(data, freq=4, start=(2000, 1))
    assert ts.nobs == 48
    assert ts.freq == 4


def test_describe_monthly():
    import math
    data = [1.2, 0.8, 1.5, 0.9, 1.1, 0.7, 1.3, 1.0, 0.6, 1.4, 1.1, 0.8]
    ts = TimeSeries.from_array(data, freq=12, start=(2020, 1), name="test")
    out = ts.describe()
    n = len(data)
    import numpy as np
    x = np.array(data)
    mu  = x.mean()
    std = x.std(ddof=0)
    skew = (((x - mu) / std) ** 3).mean()
    kurt = (((x - mu) / std) ** 4).mean() - 3.0
    jb   = n / 6.0 * (skew**2 + kurt**2 / 4.0)
    assert f"{mu:.6f}" in out
    assert f"{jb:.6f}" in out
    assert "9/2020" in out      # minimum at observation 9
    assert "3/2020" in out      # maximum at observation 3


def test_describe_annual():
    ts = TimeSeries.from_array([100, 105, 102], freq=1, start=(2000, 1))
    out = ts.describe()
    assert "2000" in out
    assert "2002" in out
    assert "from 2000 to 2002" in out


def test_describe_obs_to_date():
    ts = TimeSeries.from_array([0.0] * 13, freq=12, start=(2020, 1))
    assert ts._obs_to_date(1)  == (2020, 1)
    assert ts._obs_to_date(12) == (2020, 12)
    assert ts._obs_to_date(13) == (2021, 1)


# ── Intervention ──────────────────────────────────────────────────────────

def test_intervention_pulse():
    itv = Intervention("pulse", at=10)
    assert itv.type_code == 0
    assert itv.omega == [1.0]
    assert itv.delta == []


def test_intervention_step_with_delta():
    itv = Intervention("step", at=24, omega=[1.0], delta=[0.8])
    assert itv.type_code == 1
    assert len(itv.delta) == 1


def test_intervention_bad_type():
    with pytest.raises(ValueError):
        Intervention("invalid", at=5)


def test_intervention_mismatched_free():
    with pytest.raises(ValueError):
        Intervention("pulse", at=0, omega=[1.0, 0.5], omega_free=[True])


# ── Model (no C extension) ────────────────────────────────────────────────

def test_model_creation():
    ts  = _dummy_series()
    m   = Model(ts, ar=[[0.5]], ma=[[0.3]], d=1)
    assert m.d == 1
    assert len(m.ar) == 1
    assert len(m.interventions) == 0


def test_model_add_intervention():
    ts  = _dummy_series()
    m1  = Model(ts, ar=[[0.5]], ma=[[0.3]])
    m2  = m1.add_intervention("step", at=30, omega=[1.0], delta=[0.5])
    assert len(m1.interventions) == 0   # original unchanged
    assert len(m2.interventions) == 1


def test_model_fit_arma11():
    """ARMA(1,1) smoke fit: engine must converge and return sensible values."""
    ts = _dummy_series()
    m  = Model(ts, ar=[[0.5]], ma=[[0.3]])
    try:
        m.fit()
        assert m._result.converged
        assert m._result.npar == 2
        assert m._result.loglik < 0          # log-likelihood is always negative
        assert len(m._result.residuals) == ts.nobs
    except ImportError:
        pytest.skip("C extension not compiled")


def test_model_requires_fit_before_results():
    ts = _dummy_series()
    m  = Model(ts)
    with pytest.raises(RuntimeError, match="not been fitted"):
        _ = m.residuals


# ── TimeSeries.from_pandas ────────────────────────────────────────────────

def test_from_pandas_annual():
    pd = pytest.importorskip("pandas")
    s = pd.Series(range(20),
                  index=pd.period_range("1990", periods=20, freq="A"))
    ts = TimeSeries.from_pandas(s)
    assert ts.nobs == 20
    assert ts.freq == 1
    assert ts.start == (1990, 1)


def test_from_pandas_monthly():
    pd = pytest.importorskip("pandas")
    s = pd.Series(range(36),
                  index=pd.period_range("2000-01", periods=36, freq="M"))
    ts = TimeSeries.from_pandas(s, name="test_m")
    assert ts.nobs == 36
    assert ts.freq == 12
    assert ts.start == (2000, 1)
    assert ts.name == "test_m"


def test_from_pandas_quarterly():
    pd = pytest.importorskip("pandas")
    s = pd.Series(range(12),
                  index=pd.period_range("1995Q1", periods=12, freq="Q"))
    ts = TimeSeries.from_pandas(s)
    assert ts.freq == 4
    assert ts.start == (1995, 1)


def test_from_pandas_datetime_index():
    pd = pytest.importorskip("pandas")
    s = pd.Series(range(24),
                  index=pd.date_range("2010-01", periods=24, freq="MS"))
    ts = TimeSeries.from_pandas(s, freq=12)
    assert ts.freq == 12
    assert ts.start == (2010, 1)


# ── Model.compare ─────────────────────────────────────────────────────────

def test_model_compare_requires_fitted():
    ts = _dummy_series()
    m1 = Model(ts, ar=[[0.5]])
    m2 = Model(ts, ma=[[0.3]])
    with pytest.raises(RuntimeError):
        m1.compare(m2)


def test_model_compare_output():
    ts = _dummy_series()
    m1 = Model(ts, ar=[[0.5]], ma=[[0.3]])
    try:
        m1.fit()
    except ImportError:
        pytest.skip("C extension not compiled")
    table = m1.compare()
    assert "loglik" in table
    assert "AIC" in table
