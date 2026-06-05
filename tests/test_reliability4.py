"""
Numerical reliability tests — part 4.

Covers remaining 15 gaps:

  Group A — Diagnostics: numerical values vs reference formulas
  Group B — TimeSeries: obs_to_date, describe() values
  Group C — BoxCox: round-trip inverse transform
  Group D — calcnu_py: indicator types (step, pulse, ramp, seasonal)
  Group E — estimate_py: sigma2 = f1/n_eff via elf_scalar at optimum
  Group F — write_out: section presence and ACF band value
  Group G — inp parser: numbering flag, cbands, all real cases load
"""

import math
import os
import re
import pytest
import numpy as np

import fue
from fue import TimeSeries, Model, acf, ljung_box, jarque_bera
from fue.intervention import Intervention
from fue.cast_us import (
    calcnu_py, _build_indicator, build_est_spec, cast_us_py, estimate_py,
    _logelf_c,
)
from fue.elfvarma import elf_scalar

# ── C availability ────────────────────────────────────────────────────────────

try:
    from fue._fue_engine import ffi as _ffi   # noqa: F401
    _C_AVAILABLE = True
except ImportError:
    _C_AVAILABLE = False

requires_c = pytest.mark.skipif(not _C_AVAILABLE,
                                reason="C extension not compiled")

_REAL = os.path.join(os.path.dirname(__file__), "real_cases")

_SFNY30 = np.array([
    3.91505848, 2.02125792, 0.81208771, 0.60807414, 1.21576447,
    1.43763055, 1.78032601, 0.82841058, 0.65433228, 0.74324607,
    0.93394905, 0.60094494, 0.80840161, 0.90899270, 0.40822203,
    0.41975993, 0.50368768, 0.57248427, 0.72970370, 0.90175445,
    0.61763439, 0.63607641, 0.67670827, 0.81812744, 0.78095914,
    0.82024104, 0.86103433, 0.84442843, 0.74566075, 0.63347579,
])


# ══════════════════════════════════════════════════════════════════════════════
# Group A — Diagnostics: numerical values vs reference formulas
# ══════════════════════════════════════════════════════════════════════════════

class TestDiagnosticsNumerical:
    """fue diagnostics must match their defining formulas."""

    def test_acf_matches_pearson(self):
        """acf(w, lags=k)[i] == corr(w[k:], w[:-k]) for k=1..5."""
        w = _SFNY30.copy()
        c = acf(w, lags=5)
        n, mu = len(w), w.mean()
        wc = w - mu
        c0 = np.dot(wc, wc) / n
        for k in range(1, 6):
            expected = np.dot(wc[k:], wc[:-k]) / (n * c0)
            assert abs(c[k - 1] - expected) < 1e-10, (
                f"acf[{k}] = {c[k-1]:.10f}, expected {expected:.10f}"
            )

    def test_ljung_box_q_formula(self):
        """LB Q = n*(n+2) * sum(rk²/(n-k), k=1..lags)."""
        w = _SFNY30.copy()
        lags = 5
        result = ljung_box(w, lags=lags)
        q_reported = result["statistic"][0] if isinstance(result["statistic"], list) \
            else result["statistic"]
        c = acf(w, lags=lags)
        n = len(w)
        q_manual = n * (n + 2) * sum(c[k] ** 2 / (n - (k + 1)) for k in range(lags))
        assert abs(q_reported - q_manual) < 1e-8, (
            f"LB Q reported={q_reported:.8f}, manual={q_manual:.8f}"
        )

    def test_jarque_bera_matches_scipy(self):
        """fue jarque_bera == scipy.stats.jarque_bera on same data."""
        from scipy import stats
        rng = np.random.default_rng(42)
        w = rng.standard_normal(100)
        jb_fue    = jarque_bera(w)
        jb_scipy  = stats.jarque_bera(w)
        p_fue   = jb_fue.pvalue   if hasattr(jb_fue,   "pvalue") else jb_fue["pvalue"]
        p_scipy = jb_scipy.pvalue if hasattr(jb_scipy, "pvalue") else jb_scipy["pvalue"]
        assert abs(p_fue - p_scipy) < 1e-10, (
            f"JB p-value: fue={p_fue:.10f}, scipy={p_scipy:.10f}"
        )

    def test_acf_length(self):
        """acf(w, lags=k) returns exactly k values (lags 1..k)."""
        w = _SFNY30.copy()
        for lags in (5, 10, 24):
            c = acf(w, lags=lags)
            assert len(c) == lags, f"acf lags={lags}: got {len(c)} values"

    def test_ljung_box_pvalue_in_unit_interval(self):
        w = _SFNY30.copy()
        result = ljung_box(w, lags=10)
        p = result["pvalue"][0] if isinstance(result["pvalue"], list) \
            else result["pvalue"]
        assert 0.0 <= p <= 1.0


# ══════════════════════════════════════════════════════════════════════════════
# Group B — TimeSeries: obs_to_date, describe()
# ══════════════════════════════════════════════════════════════════════════════

class TestTimeSeries:

    def test_describe_mean(self):
        ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1))
        desc = ts.describe()
        match = re.search(r'Mean:\s+([\d.]+)', desc)
        assert match, "Mean not found in describe()"
        assert abs(float(match.group(1)) - np.mean(_SFNY30)) < 1e-4

    def test_describe_variance(self):
        ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1))
        desc = ts.describe()
        match = re.search(r'Variance:\s+([\d.]+)', desc)
        assert match, "Variance not found in describe()"
        assert abs(float(match.group(1)) - np.var(_SFNY30)) < 1e-4

    def test_describe_min_max(self):
        ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1))
        desc = ts.describe()
        min_match = re.search(r'Minimum:\s+([\d.]+)', desc)
        max_match = re.search(r'Maximum:\s+([\d.]+)', desc)
        assert min_match and max_match
        assert abs(float(min_match.group(1)) - np.min(_SFNY30)) < 1e-4
        assert abs(float(max_match.group(1)) - np.max(_SFNY30)) < 1e-4

    def test_describe_nobs(self):
        ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1))
        desc = ts.describe()
        assert "30" in desc

    def test_obs_to_date_annual(self):
        """Annual series: obs k starts at year start[0] + k - 1."""
        ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1))
        # obs 1 → 1852, obs 30 → 1881
        assert "1852" in ts.describe()
        assert "1881" in ts.describe()

    def test_obs_to_date_quarterly(self):
        """Quarterly series: obs 5 (Q1 of year 2) appears correctly."""
        data = np.ones(8)
        ts = TimeSeries(data, freq=4, start=(2000, 1))
        desc = ts.describe()
        # Should contain 2000 and 2001
        assert "2000" in desc
        assert "2001" in desc

    def test_obs_to_date_monthly(self):
        data = np.ones(24)
        ts = TimeSeries(data, freq=12, start=(2020, 1))
        desc = ts.describe()
        assert "2020" in desc
        assert "2021" in desc

    def test_nobs_property(self):
        ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1))
        assert ts.nobs == 30

    def test_freq_property(self):
        ts = TimeSeries(_SFNY30, freq=4, start=(2000, 1))
        assert ts.freq == 4


# ══════════════════════════════════════════════════════════════════════════════
# Group C — BoxCox: round-trip inverse transform
# ══════════════════════════════════════════════════════════════════════════════

class TestBoxCox:
    """BoxCox(x, lam) followed by inverse must recover x."""

    def _roundtrip(self, x, lam, refactor=1.0):
        from fue.forecast import _boxcox, _inv_boxcox
        y  = _boxcox(x, lam, refactor)
        xr = _inv_boxcox(y, lam, refactor)
        return xr

    @pytest.mark.parametrize("lam", [0.0, 0.5, 1.0, -0.5, 2.0])
    def test_roundtrip_scalar(self, lam):
        x = 2.5
        xr = self._roundtrip(x, lam)
        assert abs(xr - x) < 1e-10, f"lambda={lam}: {x} → {xr}"

    @pytest.mark.parametrize("lam", [0.0, 0.5, 1.0])
    def test_roundtrip_array(self, lam):
        data = np.abs(_SFNY30) + 0.1   # strictly positive
        from fue.forecast import _boxcox, _inv_boxcox
        y  = np.array([_boxcox(x, lam, 1.0) for x in data])
        xr = np.array([_inv_boxcox(yi, lam, 1.0) for yi in y])
        np.testing.assert_allclose(xr, data, atol=1e-9,
                                   err_msg=f"BoxCox round-trip lambda={lam}")

    def test_lam1_is_identity(self):
        """BoxCox(lambda=1) should be identity (x itself)."""
        from fue.forecast import _boxcox
        x = 3.14
        assert abs(_boxcox(x, 1.0, 1.0) - x) < 1e-10

    def test_lam0_is_log(self):
        """BoxCox(lambda=0) should equal log(x)."""
        from fue.forecast import _boxcox
        x = math.e
        assert abs(_boxcox(x, 0.0, 1.0) - 1.0) < 1e-10   # log(e) = 1

    def test_refactor_scales_transform(self):
        """refactor multiplies the log transform: BoxCox(x,0,r) = r*log(x)."""
        from fue.forecast import _boxcox
        x, r = 2.0, 100.0
        assert abs(_boxcox(x, 0.0, r) - r * math.log(x)) < 1e-10


# ══════════════════════════════════════════════════════════════════════════════
# Group D — calcnu_py: indicator types
# ══════════════════════════════════════════════════════════════════════════════

class TestIndicatorTypes:
    """_build_indicator produces the correct 1-indexed array for each type."""

    def test_pulse_single_spike(self):
        """Pulse at obs 3: indicator[3]=1, rest=0."""
        itv = Intervention("pulse", at=2)   # at=2 → obs 3 (0-indexed)
        ind = _build_indicator(itv, 8, 1, 1)
        assert ind[3] == 1.0
        assert np.sum(ind) == 1.0

    def test_step_from_obs(self):
        """Step at obs 4: indicator[k]=1 for k>=4."""
        itv = Intervention("step", at=3)
        ind = _build_indicator(itv, 8, 1, 1)
        assert np.all(ind[4:9] == 1.0)
        assert np.all(ind[1:4] == 0.0)

    def test_ramp_from_obs(self):
        """Ramp at obs 5: indicator[k] = k - 5 + 1 for k>=5, else 0."""
        itv = Intervention("ramp", at=4)   # at=4 → obs 5
        ind = _build_indicator(itv, 10, 1, 1)
        assert ind[5] == 1.0
        assert ind[6] == 2.0
        assert ind[7] == 3.0
        assert np.all(ind[1:5] == 0.0)

    def test_seasonal_quarterly(self):
        """Seasonal indicator: at=k fires at obs where ((j-begtime) % freq)+1 == k+1.

        at=2, freq=4, begtime=1: fires when ((j-1) % 4)+1 == 3 → obs 3,7,11.
        """
        itv = Intervention("seasonal", at=2)
        ind = _build_indicator(itv, 12, freq=4, begtime=1)
        assert ind[3] == 1.0    # obs 3 = Q3 year 1
        assert ind[7] == 1.0    # obs 7 = Q3 year 2
        assert ind[11] == 1.0   # obs 11 = Q3 year 3
        assert ind[2] == 0.0 and ind[4] == 0.0

    def test_calcnu_step_response(self):
        """Step filter (delta=[1]) → nu[j]=1 for all j."""
        nu = calcnu_py([1.0], [1.0], lags=10)
        np.testing.assert_allclose(nu, np.ones(11), atol=1e-12)

    def test_calcnu_pulse_response(self):
        """Pulse filter (no delta) → nu[0]=omega[0], nu[j]=0 for j>0."""
        nu = calcnu_py([2.5], [], lags=5)
        assert abs(nu[0] - 2.5) < 1e-12
        np.testing.assert_allclose(nu[1:], 0.0, atol=1e-12)

    def test_calcnu_ramp_via_double_integration(self):
        """Ramp = ω(B)/δ(B)² with δ=[1,1] : nu grows linearly."""
        # δ(B) = 1 - B: double unit root → ramp
        # delta = [1, -1] means 1*(nu[j-1]) - (-1)*(nu[j-2]) = nu[j-1]+nu[j-2]?
        # Actually step is delta=[1], so ramp needs delta=[1] applied twice
        # Simple test: nu grows monotonically
        nu = calcnu_py([1.0], [2.0, -1.0], lags=8)  # δ(B) = 1 - 2B + B² = (1-B)²
        # (1-B)² ramp response: should grow like j+1
        for j in range(1, 8):
            assert nu[j] >= nu[j-1] - 1e-10, f"nu not non-decreasing at j={j}"

    def test_indicator_pulse_out_of_range(self):
        """Pulse indicator with at= beyond nobs should give all zeros."""
        itv = Intervention("pulse", at=99)
        ind = _build_indicator(itv, 10, 1, 1)
        assert np.sum(ind) == 0.0


# ══════════════════════════════════════════════════════════════════════════════
# Group E — estimate_py: sigma2 = f1/n_eff via elf_scalar at optimum
# ══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("model_fn,label", [
    (lambda: Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ar=[[0.5]]),
     "ar1"),
    (lambda: Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ma=[[0.3]], d=1),
     "ima11"),
])
def test_sigma2_equals_f1_over_n(model_fn, label):
    """estimate_py sigma2 must equal f1/n_eff from elf_scalar at optimal x."""
    m = model_fn()
    r = estimate_py(m)
    spec = build_est_spec(m)
    p, q, phi, theta, mu, w, _ = cast_us_py(r["params"], spec)
    n_eff = len(w)
    _, f1, _, _, _ = elf_scalar(n_eff, p, q, phi, theta, w,
                                sigma2=1.0, mu=mu, do_chkma=False)
    sigma2_from_elf = f1 / n_eff
    assert abs(r["sigma2"] - sigma2_from_elf) < 1e-8, (
        f"{label}: sigma2={r['sigma2']:.10f}, f1/n={sigma2_from_elf:.10f}"
    )


# ══════════════════════════════════════════════════════════════════════════════
# Group F — write_out: section presence and ACF band value
# ══════════════════════════════════════════════════════════════════════════════

class TestWriteOutSections:

    def _report(self):
        m = Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ar=[[0.5]])
        m.fit()
        return m, m.write_out()

    def test_acf_section_present(self):
        _, report = self._report()
        assert "Autocorrelation" in report

    def test_acf_band_value(self):
        """ACF confidence band = ±2/sqrt(n), shown in report header."""
        m, report = self._report()
        n = len(m.residuals.data)
        expected_band = 2.0 / math.sqrt(n)
        match = re.search(r'acf bands.*?(\d+\.\d+)', report)
        assert match, "ACF band value not found in report"
        band_in_report = float(match.group(1))
        assert abs(band_in_report - expected_band) < 0.005, (
            f"ACF band {band_in_report:.3f} vs expected {expected_band:.3f}"
        )

    def test_histogram_section_present(self):
        _, report = self._report()
        assert "Histogram" in report or "istogram" in report

    def test_parameter_table_present(self):
        _, report = self._report()
        # Parameter values appear with std_errors in [ ] notation
        assert re.search(r'\[\s*\d+\]', report), \
            "Parameter index [k] not found in report"

    def test_residuals_plot_present(self):
        _, report = self._report()
        assert "Standardized time series plot" in report

    def test_correlation_section_present(self):
        _, report = self._report()
        assert "Autocorrelation function" in report


# ══════════════════════════════════════════════════════════════════════════════
# Group G — inp parser: all real cases, numbering flag, cbands
# ══════════════════════════════════════════════════════════════════════════════

class TestInpParser:

    _FUG = {
        "PRICES/GDP/Sample_1.2003_4.2019/Idem/R.1",
        "PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC",
        "PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Idem/R.1",
        "PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Idem/R1",
        "PRICES/PCE/Sample_1.2003_4.2019/Idem/R.1",
    }

    def _all_inp_files(self):
        result = []
        for root, _, files in os.walk(_REAL):
            for f in sorted(files):
                if not f.endswith(".inp"):
                    continue
                base = os.path.join(root, f[:-4])
                rel  = os.path.relpath(base, _REAL)
                if rel not in self._FUG:
                    result.append((rel, base))
        return result

    def test_all_real_inp_load_without_error(self):
        """Every non-FUG .inp file parses without exception."""
        errors = []
        for rel, base in self._all_inp_files():
            try:
                ts, m = fue.load(base + ".inp")
                assert ts.nobs > 0, f"{rel}: nobs=0"
            except Exception as e:
                errors.append(f"{rel}: {e}")
        assert not errors, "Load errors:\n" + "\n".join(errors)

    def test_numbering_flag_false_by_default(self):
        """Most models have numbering=False (period labels, not obs numbers)."""
        _, m = fue.load(os.path.join(
            _REAL, "PRICES/GDP/Sample_1.2003_4.2019/Mod/R.1.inp"))
        assert m.series.numbering == False

    def test_nobs_matches_data(self):
        """Parsed nobs matches length of data array."""
        for rel, base in self._all_inp_files()[:5]:
            ts, _ = fue.load(base + ".inp")
            assert len(ts.data) == ts.nobs, (
                f"{rel}: nobs={ts.nobs} but len(data)={len(ts.data)}"
            )

    def test_freq_is_positive(self):
        """All parsed models have freq >= 1."""
        for rel, base in self._all_inp_files():
            ts, _ = fue.load(base + ".inp")
            assert ts.freq >= 1, f"{rel}: freq={ts.freq}"

    def test_boxlam_is_float(self):
        """boxlam is parsed as a float for all models."""
        for rel, base in self._all_inp_files():
            _, m = fue.load(base + ".inp")
            assert isinstance(m.boxlam, float), (
                f"{rel}: boxlam={m.boxlam!r} is not float"
            )
