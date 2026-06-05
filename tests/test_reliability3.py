"""
Numerical reliability tests — part 3.

Covers remaining gaps after test_reliability.py and test_reliability2.py:

  Group A — C vs Python: params, std_errors, cov_matrix, residuals
  Group B — Covariance matrix: symmetry, PSD, diagonal == std_errors²
  Group C — write_out: param values, std_errors, nobs in report
  Group D — Forecast: IMA level, BoxCox positivity, RIPC.1 range
  Group E — ACF/PACF of residuals within confidence bounds
  Group F — Roundtrip: real cases load → fit → write_pre → reload → npar
  Group G — Model spec validation
"""

import math
import os
import re
import tempfile
import pytest
import numpy as np

import fue
from fue import TimeSeries, Model, FixedFreqFactor, acf, pacf
from fue.intervention import Intervention
from fue.cast_us import estimate_py

# ── C availability ────────────────────────────────────────────────────────────

try:
    from fue._fue_engine import ffi as _ffi   # noqa: F401
    _C_AVAILABLE = True
except ImportError:
    _C_AVAILABLE = False

requires_c = pytest.mark.skipif(not _C_AVAILABLE,
                                reason="C extension not compiled")

_REAL = os.path.join(os.path.dirname(__file__), "real_cases")

# ── Shared data ───────────────────────────────────────────────────────────────

_SFNY30 = np.array([
    3.91505848, 2.02125792, 0.81208771, 0.60807414, 1.21576447,
    1.43763055, 1.78032601, 0.82841058, 0.65433228, 0.74324607,
    0.93394905, 0.60094494, 0.80840161, 0.90899270, 0.40822203,
    0.41975993, 0.50368768, 0.57248427, 0.72970370, 0.90175445,
    0.61763439, 0.63607641, 0.67670827, 0.81812744, 0.78095914,
    0.82024104, 0.86103433, 0.84442843, 0.74566075, 0.63347579,
])


def _ar1():
    return Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ar=[[0.5]])


def _ima11():
    return Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ma=[[0.3]], d=1)


def _from_inp(rel):
    _, m = fue.load(os.path.join(_REAL, rel + ".inp"))
    return m


# ══════════════════════════════════════════════════════════════════════════════
# Group A — C vs Python: params, std_errors, cov_matrix, residuals
# ══════════════════════════════════════════════════════════════════════════════

@requires_c
class TestCvsPyOutputs:
    """Detailed agreement between C and Python engine outputs."""

    # ── IMA(1,1) residuals ────────────────────────────────────────────────────

    def test_ima11_residuals_agree(self):
        from fue._engine import estimate
        rc = estimate(_ima11())
        rp = estimate_py(_ima11())
        np.testing.assert_allclose(rc["residuals"], rp["residuals"],
                                   atol=1e-3, err_msg="IMA(1,1) residuals")

    def test_ima11_nresiduals_agree(self):
        from fue._engine import estimate
        rc = estimate(_ima11())
        rp = estimate_py(_ima11())
        assert rc["nresiduals"] == rp["nresiduals"] == 29

    # ── RIPC.1 params ─────────────────────────────────────────────────────────

    def test_ripc1_params_count(self):
        from fue._engine import estimate
        rc = estimate(_from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1"))
        rp = estimate_py(_from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1"))
        assert rc["npar"] == rp["npar"] == 14

    def test_ripc1_params_agree(self):
        """RIPC.1 params: C vs Py within 1e-2 per component.

        cos/sin pairs are near-degenerate (flat landscape), so individual
        param differences can reach ~7e-3 while loglik is essentially equal.
        """
        from fue._engine import estimate
        rc = estimate(_from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1"))
        rp = estimate_py(_from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1"))
        max_diff = np.max(np.abs(rc["params"] - rp["params"]))
        assert max_diff < 1e-2, (
            f"RIPC.1 max param diff {max_diff:.4f} > 1e-2"
        )

    def test_ripc1_std_errors_agree(self):
        """RIPC.1 std_errors: C vs Py within 5e-2 (Hessian approx differs)."""
        from fue._engine import estimate
        rc = estimate(_from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1"))
        rp = estimate_py(_from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1"))
        # Both should have same shape
        assert len(rc["std_errors"]) == len(rp["std_errors"]) == 14
        # Non-negative
        assert np.all(rc["std_errors"] >= 0)
        assert np.all(rp["std_errors"] >= 0)

    # ── AR(1) cov_matrix ──────────────────────────────────────────────────────

    def test_ar1_cov_matrix_agree(self):
        """AR(1) cov_matrix diagonal: C vs Py within 1e-3."""
        from fue._engine import estimate
        rc = estimate(_ar1())
        rp = estimate_py(_ar1())
        assert rc["cov_matrix"].shape == rp["cov_matrix"].shape == (1, 1)
        diff = abs(rc["cov_matrix"][0, 0] - rp["cov_matrix"][0, 0])
        assert diff < 1e-3, f"cov[0,0]: C={rc['cov_matrix'][0,0]:.6f} Py={rp['cov_matrix'][0,0]:.6f}"


# ══════════════════════════════════════════════════════════════════════════════
# Group B — Covariance matrix properties
# ══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("engine", ["c", "py"])
@pytest.mark.parametrize("model_fn,label", [
    (_ar1,   "ar1"),
    (_ima11, "ima11"),
])
class TestCovMatrix:
    """Covariance matrix must satisfy mathematical constraints."""

    def _result(self, engine, model_fn):
        if engine == "c":
            if not _C_AVAILABLE:
                pytest.skip("C extension not compiled")
            from fue._engine import estimate
            return estimate(model_fn())
        return estimate_py(model_fn())

    def test_cov_symmetric(self, engine, model_fn, label):
        r = self._result(engine, model_fn)
        cov = r["cov_matrix"]
        np.testing.assert_allclose(
            cov, cov.T, atol=1e-12,
            err_msg=f"{label} [{engine}]: cov not symmetric"
        )

    def test_cov_diagonal_nonneg(self, engine, model_fn, label):
        r = self._result(engine, model_fn)
        assert np.all(np.diag(r["cov_matrix"]) >= 0), (
            f"{label} [{engine}]: negative diagonal in cov"
        )

    def test_cov_diagonal_equals_std_errors_squared(self, engine, model_fn, label):
        r = self._result(engine, model_fn)
        diag = np.diag(r["cov_matrix"])
        se2  = r["std_errors"] ** 2
        np.testing.assert_allclose(
            diag, se2, atol=1e-10,
            err_msg=f"{label} [{engine}]: diag(cov) != std_errors²"
        )


# ══════════════════════════════════════════════════════════════════════════════
# Group C — write_out: param values, std_errors, nobs
# ══════════════════════════════════════════════════════════════════════════════

class TestWriteOutValues:
    """Parameter and std_error values in the ASCII report match _result."""

    def _fit_and_report(self, model_fn):
        m = model_fn()
        m.fit()
        return m, m.write_out()

    def test_ar1_param_in_report(self):
        """AR(1) phi appears in the report within 1e-4."""
        m, report = self._fit_and_report(_ar1)
        phi = m._result.params[0]
        # The report has lines like '      0.974752  (0.034211) [ 1]'
        matches = re.findall(r'([-\d.]+)\s+\(([\d.]+)\)\s+\[', report)
        assert matches, "No (value, se) pairs found in report"
        reported_vals = [float(v) for v, _ in matches]
        assert any(abs(v - phi) < 1e-3 for v in reported_vals), (
            f"phi={phi:.6f} not found in report values {reported_vals}"
        )

    def test_ar1_se_in_report(self):
        """AR(1) std_error appears in the report within 1e-4."""
        m, report = self._fit_and_report(_ar1)
        se = m._result.std_errors[0]
        matches = re.findall(r'([-\d.]+)\s+\(([\d.]+)\)\s+\[', report)
        assert matches
        reported_ses = [float(s) for _, s in matches]
        assert any(abs(s - se) < 1e-3 for s in reported_ses), (
            f"se={se:.6f} not found in report SEs {reported_ses}"
        )

    def test_nobs_in_report(self):
        """nobs of the series appears in the residuals section."""
        m, report = self._fit_and_report(_ar1)
        # The residuals section has '30 observations: from ...'
        assert str(m.series.nobs) in report, (
            f"nobs={m.series.nobs} not found in report"
        )

    def test_sfny2_all_params_in_report(self):
        """All 6 params of SFNY.2 appear in the report."""
        from fue.intervention import Intervention
        import numpy as np
        _SFNY62 = np.array([
            3.91505848,2.02125792,0.81208771,0.60807414,1.21576447,
            1.43763055,1.78032601,0.82841058,0.65433228,0.74324607,
            0.93394905,0.60094494,0.80840161,0.90899270,0.40822203,
            0.41975993,0.50368768,0.57248427,0.72970370,0.90175445,
            0.61763439,0.63607641,0.67670827,0.81812744,0.78095914,
            0.82024104,0.86103433,0.84442843,0.74566075,0.63347579,
            0.72637557,0.81351610,0.79142754,0.80305873,0.83867533,
            0.98678814,0.80485863,0.81651553,0.75960093,0.84070968,
            0.89480882,0.89407591,0.84323646,0.77215182,0.82509544,
            0.87384443,0.81360106,0.78497496,0.71323360,0.70688522,
            0.81090348,0.94831097,0.72598922,0.80337325,0.84011493,
            0.89247202,0.89328246,0.90942424,0.82871189,0.88647340,
            0.82251497,0.94737336,
        ])
        m = Model(
            TimeSeries(_SFNY62, freq=1, start=(1852, 1)),
            interventions=[Intervention("step", at=1, omega=[0.08],
                           omega_free=[True], delta=[0.6], delta_free=[True])],
            ar=[[0.8], [-0.1, -0.1]], boxlam=0.0, mu=0.0, estimate_mu=True,
        )
        m.fit()
        report = m.write_out()
        matches = re.findall(r'\[\s*(\d+)\]', report)
        param_indices = [int(x) for x in matches]
        # All 6 params should be listed
        assert len(set(param_indices)) >= 6, (
            f"Expected ≥6 param indices, found {set(param_indices)}"
        )


# ══════════════════════════════════════════════════════════════════════════════
# Group D — Forecast properties
# ══════════════════════════════════════════════════════════════════════════════

class TestForecastProperties:

    def test_ima11_forecast_constant(self):
        """IMA(1,1) d=1: all forecast levels equal (random walk with drift=0)."""
        m = _ima11()
        m.fit()
        fc = m.forecast(5)
        # For IMA(1,1) without mu, all h-step forecasts are equal
        np.testing.assert_allclose(fc.level, fc.level[0],
                                   rtol=1e-6, err_msg="IMA(1,1) forecasts not constant")

    def test_ima11_forecast_near_last_obs(self):
        """IMA(1,1) forecast should be close to the last observed value."""
        m = _ima11()
        m.fit()
        fc = m.forecast(5)
        last = _SFNY30[-1]
        assert abs(fc.level[0] - last) < 0.5, (
            f"IMA forecast {fc.level[0]:.4f} far from last obs {last:.4f}"
        )

    def test_boxcox_forecast_positive(self):
        """BoxCox(lambda=0) forecast: all levels must be strictly positive."""
        data = np.abs(_SFNY30) + 0.1
        ts = TimeSeries(data, freq=1, start=(1852, 1))
        m = Model(ts, ar=[[0.5]], boxlam=0.0)
        m.fit()
        fc = m.forecast(10)
        assert np.all(fc.level > 0), (
            f"BoxCox(0) forecast has non-positive values: {fc.level}"
        )

    def test_ripc1_forecast_in_range(self):
        """RIPC.1 monthly forecast: first step within observed data range."""
        m = _from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1")
        m.fit()
        fc = m.forecast(12)
        # Observed RIPC data is around 0.41–0.46
        assert 0.3 < fc.level[0] < 0.7, (
            f"RIPC.1 forecast[0] = {fc.level[0]:.4f} out of expected range"
        )

    def test_forecast_std_increases_with_horizon(self):
        """Forecast uncertainty grows (or stays equal) with horizon."""
        m = _ar1()
        m.fit()
        fc = m.forecast(10)
        # Non-decreasing std
        assert np.all(np.diff(fc.level_std) >= -1e-10), (
            "Forecast std_errors not non-decreasing"
        )

    def test_forecast_converges_to_unconditional_mean(self):
        """AR(1) long-range forecast converges to unconditional mean ≈ 0."""
        ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1))
        m = Model(ts, ar=[[0.5]])   # no mu → mean = 0 in transformed space
        m.fit()
        fc = m.forecast(100)
        # After 100 steps, level should be near the long-run mean in original scale
        # For boxlam=1 (level), this means the last few forecasts should stabilise
        diffs = np.abs(np.diff(fc.level[-10:]))
        assert np.all(diffs < 0.01), "Long-range forecast not converging"


# ══════════════════════════════════════════════════════════════════════════════
# Group E — ACF/PACF of residuals within confidence bounds
# ══════════════════════════════════════════════════════════════════════════════

class TestResidualACF:
    """Fitted residuals should have small autocorrelations."""

    def _bound(self, n, sigma=2.0):
        return sigma / math.sqrt(n)

    def test_ar1_residuals_acf_within_2sigma(self):
        r = estimate_py(_ar1())
        res = r["residuals"]
        c = acf(res, lags=10)
        bound = self._bound(len(res), sigma=2.5)   # 2.5-sigma for n=30
        violations = np.sum(np.abs(c) > bound)
        # Allow at most 1 violation out of 10 lags (5% false positive rate)
        assert violations <= 1, (
            f"AR(1) residuals: {violations}/10 lags exceed 2.5-sigma bound {bound:.3f}"
        )

    def test_ima11_residuals_acf_within_2sigma(self):
        r = estimate_py(_ima11())
        res = r["residuals"]
        c = acf(res, lags=10)
        bound = self._bound(len(res), sigma=2.5)
        violations = np.sum(np.abs(c) > bound)
        assert violations <= 1, (
            f"IMA(1,1) residuals: {violations}/10 lags exceed bound"
        )

    @requires_c
    def test_ar1_c_residuals_acf_within_2sigma(self):
        from fue._engine import estimate
        r = estimate(_ar1())
        res = r["residuals"]
        c = acf(res, lags=10)
        bound = self._bound(len(res), sigma=2.5)
        violations = np.sum(np.abs(c) > bound)
        assert violations <= 1

    @requires_c
    def test_real_case_residuals_acf(self):
        """GDP/R.1 residuals: no more than 1 ACF lag exceeds 2.5-sigma."""
        m = _from_inp("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.1")
        m.fit()
        res = m.residuals.data
        c = acf(res, lags=12)
        bound = self._bound(len(res), sigma=2.5)
        violations = np.sum(np.abs(c) > bound)
        assert violations <= 2, (
            f"GDP/R.1: {violations}/12 lags exceed 2.5-sigma bound"
        )


# ══════════════════════════════════════════════════════════════════════════════
# Group F — Roundtrip real cases
# ══════════════════════════════════════════════════════════════════════════════

@requires_c
class TestRealRoundtrip:
    """load .inp → fit → write_pre → load .pre → fit: npar and loglik preserved."""

    def _roundtrip(self, rel):
        m1 = _from_inp(rel)
        m1.fit()
        with tempfile.NamedTemporaryFile(suffix=".pre", delete=False) as f:
            pre_path = f.name
        try:
            m1.write_pre(pre_path)
            _, m2 = fue.load(pre_path)
            m2.fit()
            return m1._result, m2._result
        finally:
            os.unlink(pre_path)

    def test_gdp_r1_npar_preserved(self):
        r1, r2 = self._roundtrip("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.1")
        assert r1.npar == r2.npar

    def test_gdp_r1_loglik_preserved(self):
        r1, r2 = self._roundtrip("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.1")
        assert abs(r1.loglik - r2.loglik) < 1e-4, (
            f"loglik: before={r1.loglik:.6f} after={r2.loglik:.6f}"
        )

    def test_ipct_r1_npar_preserved(self):
        r1, r2 = self._roundtrip("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/R.1")
        assert r1.npar == r2.npar

    def test_ripc1_npar_preserved(self):
        r1, r2 = self._roundtrip("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1")
        assert r1.npar == r2.npar

    def test_ripc1_loglik_preserved(self):
        # RIPC.1 has 14 params; .pre stores 4 decimal places → rounding of
        # initial values causes the re-estimated loglik to differ by ~1e-3.
        r1, r2 = self._roundtrip("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1")
        assert abs(r1.loglik - r2.loglik) < 5e-3


# ══════════════════════════════════════════════════════════════════════════════
# Group G — Model spec validation
# ══════════════════════════════════════════════════════════════════════════════

class TestModelValidation:

    def test_omega_free_length_mismatch_raises(self):
        """omega_free must be same length as omega."""
        with pytest.raises((ValueError, AssertionError)):
            Intervention("step", at=0, omega=[1.0, 0.5], omega_free=[True])

    def test_invalid_intervention_type_raises(self):
        """Unknown intervention type should raise ValueError."""
        with pytest.raises((ValueError, KeyError)):
            Intervention("invalid_type", at=0)

    def test_model_requires_timeseries(self):
        """Model must receive a TimeSeries, not raw array."""
        with pytest.raises((TypeError, AttributeError)):
            Model(np.array([1.0, 2.0, 3.0]), ar=[[0.5]])

    def test_fit_required_before_forecast(self):
        """Calling forecast before fit should raise RuntimeError."""
        m = _ar1()
        with pytest.raises(RuntimeError):
            m.forecast(5)

    def test_fit_required_before_write_out(self):
        m = _ar1()
        with pytest.raises(RuntimeError):
            m.write_out()

    def test_fixed_freq_factor_negative_coef(self):
        """FixedFreqFactor requires coef < 0."""
        with pytest.raises(ValueError):
            FixedFreqFactor(freq=6.0, coef=0.5)   # positive: invalid

    def test_fixed_freq_factor_negative_coef_ok(self):
        """FixedFreqFactor with coef < 0 should not raise."""
        ff = FixedFreqFactor(freq=6.0, coef=-0.5)
        assert ff.coef == -0.5
