"""
Smoke tests for the Python API layer (no C extension required).

Tests that verify the C engine results go in test_estimation.py
and are marked with pytest.mark.engine.
"""

import os
import numpy as np
import pytest
from fue import TimeSeries, Intervention, Model
import fue

_INP_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "fue-1.12.02_win", "src"
)
_SFNY2_INP = os.path.join(_INP_DIR, "SFNY.2.inp")
_V8_INP    = os.path.join(_INP_DIR, "V8.inp")
_PU1_INP   = os.path.join(
    os.path.dirname(__file__), "..", "..", "fue-1.13", "examples", "PU.1.inp"
)


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
    jb   = (n // 6) * (skew**2 + kurt**2 / 4.0)
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


# ── inp.load() ────────────────────────────────────────────────────────────────

def _skip_if_no_inp(path):
    if not os.path.exists(path):
        pytest.skip(f"inp file not found: {path}")


def test_load_sfny2_metadata():
    _skip_if_no_inp(_SFNY2_INP)
    ts, m = fue.load(_SFNY2_INP)
    assert ts.freq == 1
    assert ts.nobs == 62
    assert ts.start == (1852, 1)
    assert ts.name == "SFNY"
    assert m.d == 0
    assert m.D == 0
    assert m.boxlam == 0.0
    assert m.estimate_mu is True


def test_load_sfny2_ar_factors():
    _skip_if_no_inp(_SFNY2_INP)
    _, m = fue.load(_SFNY2_INP)
    assert len(m.ar) == 2
    assert m.ar[0] == pytest.approx([0.8])
    assert m.ar[1] == pytest.approx([-0.1, -0.1])


def test_load_sfny2_intervention():
    _skip_if_no_inp(_SFNY2_INP)
    _, m = fue.load(_SFNY2_INP)
    assert len(m.interventions) == 1
    itv = m.interventions[0]
    assert itv.type_code == 1                      # step
    assert itv.at == 1                             # 1853 → obs 2 → at=1 (0-based)
    assert itv.omega == pytest.approx([0.08])
    assert itv.delta == pytest.approx([0.6])


def test_load_pu1_metadata():
    _skip_if_no_inp(_PU1_INP)
    ts, m = fue.load(_PU1_INP)
    assert ts.freq == 12
    assert ts.nobs == 115
    assert ts.start == (2000, 1)
    assert m.d == 2
    assert m.boxlam == 0.0
    assert m.refactor == pytest.approx(100.0)


def test_load_pu1_interventions():
    _skip_if_no_inp(_PU1_INP)
    _, m = fue.load(_PU1_INP)
    assert len(m.interventions) == 11
    assert m.interventions[0].type_code == 4       # cos
    assert m.interventions[0].harmonic == pytest.approx(1.0)
    assert m.interventions[10].type_code == 6      # alter


def test_load_v8_metadata():
    _skip_if_no_inp(_V8_INP)
    ts, m = fue.load(_V8_INP)
    assert ts.freq == 1
    assert ts.nobs == 253
    assert m.d == 1
    assert m.boxlam == pytest.approx(-0.4)
    assert len(m.interventions) == 6


def test_load_v8_multi_omega():
    """Impulses 5 and 6 in V8.inp have 2-element omega (Nomega=1)."""
    _skip_if_no_inp(_V8_INP)
    _, m = fue.load(_V8_INP)
    itv4 = m.interventions[4]   # impulse 198
    assert len(itv4.omega) == 2
    assert itv4.omega == pytest.approx([-0.046, 0.046])


# ── report.write_out ──────────────────────────────────────────────────────────

def _make_fitted_model():
    """Return a Model with a mock FitResult matching SFNY.2 estimates."""
    from fue.model import FitResult
    from fue.report import write_out

    ts = TimeSeries.from_array(
        list(range(62)), freq=1, start=(1852, 1), name="SFNY"
    )
    m = Model(
        ts,
        ar=[[0.8], [-0.1, -0.1]],
        interventions=[Intervention("step", at=1, omega=[-0.08], delta=[0.6])],
        estimate_mu=True, mu=0.0, boxlam=0.0,
    )
    # Mock FitResult matching SFNY.2.out values
    params     = np.array([-0.600310, 0.587290, 0.131007, 0.488024, -0.258792, 1.209391])
    std_errors = np.array([ 0.268425, 0.142186, 0.640792, 0.624776,  0.233014, 0.249425])
    n = len(params)
    cov = np.diag(std_errors ** 2)
    residuals = np.zeros(62)

    m._result = FitResult({
        'ifault': 0, 'npar': n, 'nresiduals': 62,
        'sigma2': 0.0370592811, 'loglik': 13.9573881339,
        'aic': 0.0, 'bic': 0.0,
        'params': params, 'std_errors': std_errors,
        'cov_matrix': cov,
        'residuals': residuals,
    })
    return m


def test_write_out_requires_fit():
    ts = _dummy_series()
    m  = Model(ts, ar=[[0.5]])
    with pytest.raises(RuntimeError, match="not been fitted"):
        m.write_out()


def test_write_out_returns_string():
    m = _make_fitted_model()
    out = m.write_out()
    assert isinstance(out, str)
    assert len(out) > 100


def test_write_out_param_table():
    m   = _make_fitted_model()
    out = m.write_out()
    assert "Omegas for deterministic variable 1" in out
    assert "Deltas for deterministic variable 1" in out
    assert "Coefficients for regular AR factor 1" in out
    assert "Coefficients for regular AR factor 2" in out
    assert "Mean parameter (mu)" in out
    # Check fitted values appear
    assert "-0.600310" in out
    assert "0.587290" in out
    assert "1.209391" in out


def test_write_out_hq_schwarz():
    m   = _make_fitted_model()
    out = m.write_out()
    assert "Hannan-Quinn" in out
    assert "Schwarz" in out
    assert "-3.11" in out
    assert "-2.76" in out


def test_write_out_ar_operator():
    m   = _make_fitted_model()
    out = m.write_out()
    assert "Coefficients of the Autoregressive operator" in out
    assert "phi[ 1]" in out
    # SFNY.2: phi[1] ≈ 0.619031
    assert "0.6190" in out


def test_write_out_sigma():
    m   = _make_fitted_model()
    out = m.write_out()
    assert "sigma2:" in out
    assert "0.0370592811" in out
    assert "logelf:" in out
    assert "13.9573881339" in out


def test_write_out_residual_stats():
    m   = _make_fitted_model()
    out = m.write_out()
    assert "Unconditional residuals" in out
    assert "seasonal period: 1" in out


def test_write_out_acf_section():
    # Use a real random series so ACF/PACF have something to compute
    rng = np.random.default_rng(0)
    ts  = TimeSeries(rng.standard_normal(62), freq=1, start=(1852, 1), name="SFNY")
    from fue.model import FitResult
    m   = Model(ts, ar=[[0.5]], ma=[[0.3]], estimate_mu=True)
    n   = 3
    m._result = FitResult({
        'ifault': 0, 'npar': n, 'nresiduals': 62,
        'sigma2': 0.8, 'loglik': -80.0,
        'aic': 0.0, 'bic': 0.0,
        'params':     np.array([0.4, 0.2, 0.0]),
        'std_errors': np.array([0.1, 0.1, 0.1]),
        'cov_matrix': np.eye(n) * 0.01,
        'residuals':  rng.standard_normal(62),
    })
    out = m.write_out()
    assert "Autocorrelation function" in out
    assert "Partial autocorrelation function" in out


def test_write_out_to_file(tmp_path):
    m    = _make_fitted_model()
    path = str(tmp_path / "test.out")
    ret  = m.write_out(path=path)
    assert os.path.exists(path)
    with open(path) as f:
        content = f.read()
    assert "-0.600310" in content
    assert ret.rstrip("\n") == content.rstrip("\n")


def test_report_utilities():
    """Unit tests for internal report utilities."""
    from fue.report import _iround, _poly_mul, _count_nparma, _chi_test

    # _iround: round-half-away-from-zero
    assert _iround(0.5)  ==  1
    assert _iround(-0.5) == -1
    assert _iround(0.4)  ==  0
    assert _iround(-0.4) ==  0
    assert _iround(6.636) == 7

    # _poly_mul: (1-aB)(1-bB) = 1-(a+b)B+abB^2
    p = [1.0, -0.2]
    q = [1.0, -0.3]
    r = _poly_mul(p, q)
    assert r == pytest.approx([1.0, -0.5, 0.06])

    # _count_nparma for SFNY.2-style model
    ts = _dummy_series()
    m  = Model(ts, ar=[[0.5], [0.3, 0.1]], ma=[[0.2]])
    assert _count_nparma(m) == 4   # 1+2+1

    # _chi_test: Q = n*(n+2)*sum(r^2/(n-k))
    corr = np.array([0.1, 0.2])
    n    = 50
    q    = _chi_test(corr, 2, n)
    expected = n * (n + 2) * (0.01 / 49 + 0.04 / 48)
    assert q == pytest.approx(expected)
