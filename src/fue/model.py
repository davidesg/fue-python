"""Model: ARMAX specification with interventions, fitted by exact ML via FUE."""

from .series import TimeSeries
from .intervention import Intervention


class FixedFreqFactor:
    """Second-order AR or MA factor with fixed spectral frequency.

    Polynomial: 1 − phi1·B − phi2·B²
    where phi1 = 2·cos(2π·freq/sper)·√(−phi2) is derived from the fixed
    frequency, and only phi2 (equivalently the spectral radius r = √(−phi2))
    is estimated.

    Parameters
    ----------
    freq : float
        Fixed frequency in cycles per seasonal period (pfre1 in fue.c).
        For monthly data (sper=12): freq=6 → biennial cycle.
    coef : float
        Initial value for phi2 (AR) or theta2 (MA).  Must be < 0.
    free : bool
        Estimate *coef* by ML (default True).
    """

    def __init__(self, freq, coef=-0.5, free=True):
        if float(coef) >= 0:
            raise ValueError("coef must be negative (phi2 < 0)")
        self.freq = float(freq)
        self.coef = float(coef)
        self.free = bool(free)

    def __repr__(self):
        return f"FixedFreqFactor(freq={self.freq}, coef={self.coef}, free={self.free})"


class FitResult:
    """Container for estimation results returned by the C engine."""

    def __init__(self, data):
        self.ifault     = data['ifault']
        self.converged  = data['ifault'] == 0
        self.npar       = data['npar']
        self.sigma2     = data['sigma2']
        self.loglik     = data['loglik']
        self.aic        = data['aic']
        self.bic        = data['bic']
        self.params     = data['params']
        self.std_errors = data['std_errors']
        self.cov_matrix = data['cov_matrix']
        self.residuals  = data['residuals']


class Model:
    """
    ARMAX model with linear transfer function interventions.

    Parameters
    ----------
    series : TimeSeries
        The dependent variable.
    ar : list of list of float, optional
        Regular AR factors.  Each inner list is [φ₁, …, φ_p] for one factor.
        Example: ``[[0.7]]`` for AR(1), ``[[0.5, 0.2]]`` for AR(2).
    ma : list of list of float, optional
        Regular MA factors.  Same structure as *ar*.
    ar_s : list of list of float, optional
        Seasonal AR factors Φ(Bˢ).
    ma_s : list of list of float, optional
        Seasonal MA factors Θ(Bˢ).
    d : int
        Regular differencing order (default 0).
    D : int
        Seasonal differencing order (default 0).
    interventions : list of Intervention, optional
        Deterministic components with linear transfer functions.
    mu : float, optional
        Initial value for the mean parameter (default 0.0).
    estimate_mu : bool
        Whether to include μ in estimation (default False).
    boxlam : float
        Box-Cox parameter: 0.0 = log, 1.0 = levels (default 1.0).
    eml : bool
        True = exact ML (default), False = approximate ML.
    chkma : bool
        Enforce MA invertibility (default True).
    """

    def __init__(self, series, ar=None, ma=None, ar_s=None, ma_s=None,
                 ar_free=None, ma_free=None, ar_s_free=None, ma_s_free=None,
                 ar_f=None, ma_f=None,
                 d=0, D=0, ifadf=None, interventions=None, mu=0.0,
                 estimate_mu=False, boxlam=1.0, refactor=1.0,
                 eml=True, chkma=True):
        if not isinstance(series, TimeSeries):
            raise TypeError("series must be a TimeSeries instance")
        self.series        = series
        self.ar            = ar   or []
        self.ma            = ma   or []
        self.ar_s          = ar_s or []
        self.ma_s          = ma_s or []
        # ar_free/ma_free: list of lists of bool, same shape as ar/ma/ar_s/ma_s.
        # None means all coefficients are free.
        self.ar_free       = ar_free
        self.ma_free       = ma_free
        self.ar_s_free     = ar_s_free
        self.ma_s_free     = ma_s_free
        # Fixed-frequency second-order factors (list of FixedFreqFactor)
        self.ar_f          = list(ar_f or [])
        self.ma_f          = list(ma_f or [])
        self.d             = int(d)
        self.D             = int(D)
        self.ifadf         = list(ifadf) if ifadf else []
        self.interventions = list(interventions or [])
        self.mu0           = float(mu)
        self.estimate_mu   = bool(estimate_mu)
        self.boxlam        = float(boxlam)
        self.refactor      = float(refactor)
        self.eml           = bool(eml)
        self.chkma         = bool(chkma)
        self._result       = None

    # ── Model building helpers ────────────────────────────────────────────

    def add_intervention(self, type, at, omega=None, delta=None,
                         omega_free=None, delta_free=None):
        """Return a new Model with one extra intervention appended."""
        itv = Intervention(type, at, omega=omega, delta=delta,
                           omega_free=omega_free, delta_free=delta_free)
        new = Model(
            self.series, ar=self.ar, ma=self.ma,
            ar_s=self.ar_s, ma_s=self.ma_s,
            ar_free=self.ar_free, ma_free=self.ma_free,
            ar_s_free=self.ar_s_free, ma_s_free=self.ma_s_free,
            ar_f=self.ar_f, ma_f=self.ma_f,
            d=self.d, D=self.D, ifadf=self.ifadf,
            interventions=self.interventions + [itv],
            mu=self.mu0, estimate_mu=self.estimate_mu,
            boxlam=self.boxlam, refactor=self.refactor,
            eml=self.eml, chkma=self.chkma,
        )
        return new

    # ── Estimation ────────────────────────────────────────────────────────

    def fit(self):
        """
        Estimate model parameters by exact maximum likelihood.

        Sets self._result and returns self (for chaining).
        Raises RuntimeError if estimation returns a non-zero ifault.
        """
        from ._engine import estimate
        raw = estimate(self)
        self._result = FitResult(raw)
        if not self._result.converged:
            try:
                from fue._fue_engine import ffi, lib
                msg = ffi.string(lib.fue_strerror(self._result.ifault)).decode()
            except ImportError:
                msg = f"ifault={self._result.ifault}"
            raise RuntimeError(f"FUE estimation failed: {msg}")
        return self

    def forecast_fuf(self, horizon=None, sigma2=None):
        """
        Compute forecasts using the current parameter values as fixed estimates.

        This mirrors the fuf workflow: parameters are read from the model as-is
        (no re-estimation).  Residuals are computed in a single forward pass.

        If the model was loaded from a fuf file (via fue.load_fuf()), the
        horizon and sigma2 from the file are used when not explicitly provided.

        Parameters
        ----------
        horizon : int, optional
            Forecast horizon (number of steps ahead). Required if the model
            was not loaded from a fuf file.
        sigma2 : float, optional
            Innovation variance.  If None and the model has a stored fuf sigma2
            (from load_fuf), that value is used; otherwise it is estimated from
            the data at the provided parameter values.

        Returns
        -------
        ForecastResult
        """
        if horizon is None:
            horizon = getattr(self, "_fuf_horizon", None)
        if horizon is None:
            raise ValueError("forecast_fuf: horizon must be provided")
        if sigma2 is None:
            sigma2 = getattr(self, "_fuf_sigma2", None)

        from .cast_us import eval_at_params
        raw = eval_at_params(self)
        if raw["ifault"] != 0:
            raise RuntimeError(f"forecast_fuf: eval_at_params failed (ifault={raw['ifault']})")
        if sigma2 is None:
            sigma2 = raw["sigma2"]

        # Build a synthetic FitResult using the provided sigma2
        raw_with_sigma = {**raw, "sigma2": sigma2}
        result = FitResult(raw_with_sigma)

        from .forecast import forecast as _forecast
        return _forecast(self, result, int(horizon))

    # ── Results ───────────────────────────────────────────────────────────

    @property
    def residuals(self):
        self._require_fit()
        from .series import TimeSeries
        return TimeSeries(self._result.residuals,
                          freq=self.series.freq,
                          name="residuals")

    @property
    def params(self):
        self._require_fit()
        return self._result.params

    @property
    def std_errors(self):
        self._require_fit()
        return self._result.std_errors

    @property
    def loglik(self):
        self._require_fit()
        return self._result.loglik

    @property
    def aic(self):
        self._require_fit()
        return self._result.aic

    @property
    def bic(self):
        self._require_fit()
        return self._result.bic

    def forecast(self, horizon):
        """
        Compute L-step-ahead ARMAX forecasts.

        Parameters
        ----------
        horizon : int
            Number of periods ahead to forecast.

        Returns
        -------
        ForecastResult
            Dataclass with level, diff1, seasonal_diff arrays and their
            standard deviations (all length *horizon*).
        """
        self._require_fit()
        from .forecast import forecast as _forecast
        return _forecast(self, self._result, horizon)

    def compare(self, *others):
        """
        Print a comparison table of fitted models.

        Parameters
        ----------
        *others : Model
            Additional fitted models to compare against *self*.

        Returns
        -------
        str
            Formatted table (also printed to stdout).
        """
        models = [self] + list(others)
        for m in models:
            if m._result is None:
                raise RuntimeError("All models must be fitted before comparing.")

        header = f"{'Model':<20} {'npar':>5} {'loglik':>12} {'sigma2':>12} {'AIC':>10} {'BIC':>10}"
        sep    = "-" * len(header)
        rows   = [header, sep]
        for m in models:
            r = m._result
            label = m.series.name[:19]
            rows.append(
                f"{label:<20} {r.npar:>5d} {r.loglik:>12.4f} {r.sigma2:>12.6f}"
                f" {r.aic:>10.4f} {r.bic:>10.4f}"
            )
        table = "\n".join(rows)
        print(table)
        return table

    def summary(self):
        self._require_fit()
        r = self._result
        lines = [
            f"FUE Model — {self.series.name}",
            f"  nobs      : {self.series.nobs}",
            f"  freq      : {self.series.freq}",
            f"  d / D     : {self.d} / {self.D}",
            f"  AR factors: {len(self.ar)} regular, {len(self.ar_s)} seasonal, {len(self.ar_f)} f-fixed",
            f"  MA factors: {len(self.ma)} regular, {len(self.ma_s)} seasonal, {len(self.ma_f)} f-fixed",
            f"  Interv.   : {len(self.interventions)}",
            f"  npar      : {r.npar}",
            f"  loglik    : {r.loglik:.6f}",
            f"  sigma²    : {r.sigma2:.6f}",
            f"  AIC       : {r.aic:.4f}",
            f"  BIC       : {r.bic:.4f}",
            "",
            "  Parameters:",
        ]
        for i, (p, se) in enumerate(zip(r.params, r.std_errors)):
            lines.append(f"    [{i:2d}]  {p:12.6f}  (se {se:.6f})")
        return "\n".join(lines)

    # ── Diagnostics / plots ───────────────────────────────────────────────

    def write_out(self, path=None):
        """
        Generate an estimation report in fue .out format.

        Parameters
        ----------
        path : str or None
            Write to this file path, or return as a string if None.

        Returns
        -------
        str
        """
        self._require_fit()
        from .report import write_out as _write_out
        return _write_out(self, path=path)

    def write_pre(self, path):
        """
        Write a .pre file with estimated parameters as new initial values.

        Parameters
        ----------
        path : str
            Output path, e.g. "RIPC.1.pre".
        """
        self._require_fit()
        from .report import write_pre as _write_pre
        _write_pre(self, path=path)

    def write_fuf(self, horizon, sigma2=None, path=None):
        """
        Write a fuf forecast input file.

        The file contains the model's current parameter values (fitted if
        available, initial otherwise) plus the "Forecast horizon / sigma2"
        section that fuf/forecast_fuf require.

        Parameters
        ----------
        horizon : int
            Steps ahead to forecast.
        sigma2 : float, optional
            Innovation variance. Defaults to the fitted sigma2 (if the model
            has been fitted) or the fuf sigma2 stored on the model.
        path : str or None
            Write to file; return as string if None.
        """
        if sigma2 is None:
            if self._result is not None:
                sigma2 = self._result.sigma2
            elif hasattr(self, "_fuf_sigma2"):
                sigma2 = self._fuf_sigma2
            else:
                raise ValueError("write_fuf: sigma2 required (model not fitted)")
        from .report import write_fuf as _write_fuf
        return _write_fuf(self, horizon=horizon, sigma2=sigma2, path=path)

    def plot_residuals(self, lags=24):
        from .plots import plot_residual_diagnostics
        self._require_fit()
        plot_residual_diagnostics(self._result.residuals,
                                  lags=lags, title=self.series.name)

    # ── Internal ──────────────────────────────────────────────────────────

    def _require_fit(self):
        if self._result is None:
            raise RuntimeError("Model has not been fitted yet. Call .fit() first.")

    def __repr__(self):
        status = "fitted" if self._result else "unfitted"
        return (f"Model({self.series.name!r}, d={self.d}, D={self.D}, "
                f"ar={self.ar}, ma={self.ma}, "
                f"interventions={len(self.interventions)}, {status})")
