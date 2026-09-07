"""Model: ARMAX specification with interventions, fitted by exact ML via FUE."""

import warnings

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


#: Los códigos de parada de `raxopt` (qnewtopt.c: report()), en castellano.
_TERMINATION = {
    1: "criterio del gradiente satisfecho",
    2: "criterio de paso: los iterados dejaron de moverse",
    3: "el último paso no encontró un punto mejor",
    4: "límite de iteraciones alcanzado",
    5: "cinco pasos consecutivos de longitud máxima",
}


class FitResult:
    """Container for estimation results returned by the C engine."""

    def __init__(self, data):
        self.ifault     = data['ifault']
        self.npar       = data['npar']
        self.sigma2     = data['sigma2']
        self.loglik     = data['loglik']
        self.aic        = data['aic']
        self.bic        = data['bic']
        self.params     = data['params']
        self.std_errors = data['std_errors']
        self.cov_matrix = data['cov_matrix']
        self.residuals  = data['residuals']
        self.niter      = data.get('niter')
        self.gnorm      = data.get('gnorm')
        self.termcode   = data.get('termcode')
        self.w          = data.get('w')
        # BUG-0012: `converged` significaba `ifault == 0`, es decir «el motor no
        # reventó» — y con eso devolvía como bueno un alto por criterio de PASO
        # con ‖g‖=0.01, que no es un máximo. El veredicto del optimizador ya
        # viaja (termcode), así que la afirmación puede ser la honesta.
        #   termcode 1 = el gradiente se anuló: esto sí es un máximo
        #   termcode 0 = no registrado (motores anteriores a 0.1.10)
        self.converged  = (self.ifault == 0
                           and (self.termcode is None or self.termcode in (0, 1)))

    @property
    def termination(self):
        """Por qué paró el optimizador, en una línea. `None` si no consta."""
        if self.termcode is None:
            return None
        return _TERMINATION.get(self.termcode, f"termination code {self.termcode}")


#: Cuántas desviaciones típicas de `w` puede alejarse μ̂ de la media de `w`
#: antes de que el ajuste se considere sospechoso. Medido, no conjeturado: sobre
#: **1.523 modelos ajustados** del ecosistema con μ identificada, el máximo
#: observado es **1,38** y el percentil 99,9 es **0,73** — ninguno pasa de 2. El
#: óptimo espurio de BUG-0005 (US CPI en Windows) está en **47,2**.
#:
#: El umbral va en 5: muy por encima de todo lo legítimo que se ha visto y muy
#: por debajo de lo que se quiere cazar. No hay nada mágico en el 5; lo que
#: importa es que entre 1,4 y 47 hay un factor de 34, y ahí cabe cualquier
#: umbral razonable.
UMBRAL_MEDIA_ABSURDA = 5.0


def _avisa_si_la_media_es_absurda(model) -> None:
    """¿Es sano el óptimo al que se ha llegado? (BUG-0005)

    El optimizador es una búsqueda local sobre una superficie que puede ser
    multimodal, y **no comprueba si el óptimo al que llega tiene sentido**: en
    el caso del informe, la rueda de Windows se quedaba en una cuenca espuria
    con μ̂ = −0,144 sobre una serie cuya media diferenciada es +0,0022, y lo
    reportaba como `converged=True, ifault=0` sin un solo aviso.

    La comprobación es la que propone el propio informe, hecha sobre la serie
    correcta: μ es la media del proceso ARMA, es decir de **`w`** —la serie ya
    filtrada de deterministas y diferenciada— y NO de la serie observada. Con
    intervenciones o armónicos la diferencia es enorme, porque los deterministas
    se llevan parte del nivel.

    `w` no se reconstruye aquí: la da `cast_us_py`, que es la que usa el motor.
    Reimplementarla sería repetir el error de BUG-0014.

    **μ sólo tiene sentido si está identificada.** El término de deriva es
    μ·φ(1), así que cuando φ(1) ≈ 0 —un AR con raíz unitaria escrito con d=0—
    μ no entra en la verosimilitud y puede valer cualquier cosa. Medido: 47 de
    los 1.570 modelos del ecosistema están en ese caso, y avisar de ellos serían
    47 falsos positivos. No se comprueban.

    Es un AVISO, no una excepción: el ajuste existe y puede ser el bueno; lo que
    no puede es pasar en silencio.
    """
    import numpy as _np
    try:
        r = getattr(model, "_result", None)
        if r is None or not getattr(model, "estimate_mu", False):
            return
        from .cast_us import build_est_spec, cast_us_py
        params = _np.asarray(getattr(r, "params", None), float)
        if params.size == 0:
            return
        _p, _q, phi, _th, mu, w, ifault = cast_us_py(params, build_est_spec(model))
        if ifault or w.size < 8:
            return
        phi1 = 1.0 - float(_np.sum(phi))
        if abs(phi1) <= 0.01:          # μ no identificada: no hay nada que juzgar
            return
        sd = float(w.std(ddof=1))
        if not _np.isfinite(sd) or sd <= 0:
            return
        d = abs(float(mu) - float(w.mean())) / sd
        if d > UMBRAL_MEDIA_ABSURDA:
            warnings.warn(
                f"fue: la media estimada NO es plausible — μ̂={float(mu):.6g} "
                f"está a {d:.1f} desviaciones típicas de la media de la serie "
                f"filtrada y diferenciada ({float(w.mean()):.6g}, sd={sd:.6g}). "
                f"Sobre 1.523 ajustes reales el máximo observado es 1,4. El "
                f"optimizador es una búsqueda LOCAL y puede haber caído en una "
                f"cuenca espuria: reestima desde otras semillas y quédate con la "
                f"mayor verosimilitud (fue/bugs/BUG-0005).",
                RuntimeWarning, stacklevel=3)
    except Exception:
        # Una comprobación de sanidad no puede tumbar una estimación válida.
        pass


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
        self._inp_stem     = ""

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
        # El fallo del motor sigue siendo excepción; un alto que no es máximo,
        # no: es un ajuste que existe y sobre el que hay que poder decidir.
        if self._result.ifault != 0:
            try:
                from fue._fue_engine import ffi, lib
                msg = ffi.string(lib.fue_strerror(self._result.ifault)).decode()
            except ImportError:
                msg = f"ifault={self._result.ifault}"
            raise RuntimeError(f"FUE estimation failed: {msg}")
        if not self._result.converged:
            warnings.warn(
                f"fue: la estimación paró sin anular el gradiente "
                f"({self._result.termination}; {self._result.niter} iteraciones, "
                f"‖g‖={self._result.gnorm:.4g}). Los valores devueltos NO son un "
                f"máximo verificado — revisa las semillas antes de usarlos "
                f"(fue/bugs/BUG-0012).",
                RuntimeWarning, stacklevel=2)
        _avisa_si_la_media_es_absurda(self)
        from .cast_us import normalize_ma_invertibility, sync_params_to_attrs
        normalize_ma_invertibility(self)
        # BUG-0004 / rescaling-architecture P4: after fitting, the model IS the fitted
        # model — write the estimate back into the parameter attributes so every
        # consumer (forecast_fuf, _write_inp, reports) reads the fit, not the seeds.
        sync_params_to_attrs(self)
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
        self._result = result

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

    def write_out(self, path=None, inp_name="", out_name=""):
        """
        Generate an estimation report in fue .out format.

        Parameters
        ----------
        path : str or None
            Write to this file path, or return as a string if None.
        inp_name : str
            Label for the "Input file" header line.
        out_name : str
            Label for the "Output file" header line.
            If empty and *path* is given, the basename of *path* is used.

        Returns
        -------
        str
        """
        self._require_fit()
        from .report import write_out as _write_out
        return _write_out(self, path=path, inp_name=inp_name, out_name=out_name)

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

    def write_fuf_out(self, fr, path=None, inp_name="", out_name=""):
        """
        Generate a forecast report in fuf .out format.

        Parameters
        ----------
        fr : ForecastResult  (from model.forecast_fuf())
        path : str or None
            Write to file; return as string if None.
        inp_name, out_name : str
            Optional file-name labels shown in the header.
        """
        from .report import write_fuf_out as _write_fuf_out
        return _write_fuf_out(self, fr, path=path,
                              inp_name=inp_name, out_name=out_name)

    def plot_residuals(self, lags=None):
        from .plots import plot_model_diagnostics
        self._require_fit()
        plot_model_diagnostics(self, lags=lags)

    # ── Internal ──────────────────────────────────────────────────────────

    def _require_fit(self):
        if self._result is None:
            raise RuntimeError("Model has not been fitted yet. Call .fit() first.")

    def __repr__(self):
        status = "fitted" if self._result else "unfitted"
        return (f"Model({self.series.name!r}, d={self.d}, D={self.D}, "
                f"ar={self.ar}, ma={self.ma}, "
                f"interventions={len(self.interventions)}, {status})")
