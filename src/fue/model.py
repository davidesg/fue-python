"""Model: ARMAX specification with interventions, fitted by exact ML via FUE."""

from .series import TimeSeries
from .intervention import Intervention


class FitResult:
    """Container for estimation results returned by the C engine."""

    def __init__(self, raw):
        self._raw       = raw          # FueResult* (cffi cdata, kept alive)
        self.ifault     = raw.ifault
        self.converged  = raw.ifault == 0
        self.npar       = raw.npar
        self.sigma2     = raw.sigma2
        self.loglik     = raw.loglik
        self.aic        = raw.aic
        self.bic        = raw.bic

        import numpy as np
        n = raw.npar
        self.params     = np.frombuffer(raw.params[0:n],      dtype=float).copy()
        self.std_errors = np.frombuffer(raw.std_errors[0:n],  dtype=float).copy()
        self.cov_matrix = np.frombuffer(raw.cov_matrix[0:n*n],dtype=float).reshape(n,n).copy()
        self.residuals  = np.frombuffer(raw.residuals[0:raw.nresiduals], dtype=float).copy()


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
                 d=0, D=0, interventions=None, mu=0.0, estimate_mu=False,
                 boxlam=1.0, eml=True, chkma=True):
        if not isinstance(series, TimeSeries):
            raise TypeError("series must be a TimeSeries instance")
        self.series        = series
        self.ar            = ar   or []
        self.ma            = ma   or []
        self.ar_s          = ar_s or []
        self.ma_s          = ma_s or []
        self.d             = int(d)
        self.D             = int(D)
        self.interventions = list(interventions or [])
        self.mu0           = float(mu)
        self.estimate_mu   = bool(estimate_mu)
        self.boxlam        = float(boxlam)
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
            d=self.d, D=self.D,
            interventions=self.interventions + [itv],
            mu=self.mu0, estimate_mu=self.estimate_mu,
            boxlam=self.boxlam, eml=self.eml, chkma=self.chkma,
        )
        return new

    # ── Estimation ────────────────────────────────────────────────────────

    def fit(self):
        """
        Estimate model parameters by exact maximum likelihood.

        Sets self._result and returns self (for chaining).
        Raises RuntimeError if the C engine returns a non-zero ifault.
        """
        from ._engine import estimate   # cffi bridge (built in Phase 1)
        raw = estimate(self)
        self._result = FitResult(raw)
        if not self._result.converged:
            from .._engine import fue_strerror   # noqa: F401 — will exist
            msg = f"ifault={self._result.ifault}"
            raise RuntimeError(f"FUE estimation failed: {msg}")
        return self

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

    def summary(self):
        self._require_fit()
        r = self._result
        lines = [
            f"FUE Model — {self.series.name}",
            f"  nobs      : {self.series.nobs}",
            f"  freq      : {self.series.freq}",
            f"  d / D     : {self.d} / {self.D}",
            f"  AR factors: {len(self.ar)} regular, {len(self.ar_s)} seasonal",
            f"  MA factors: {len(self.ma)} regular, {len(self.ma_s)} seasonal",
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
