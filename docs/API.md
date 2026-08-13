# API reference

*Generated from the docstrings by `tools/gen_api_reference.py`. Do not edit: fix the docstring and regenerate. `tests/test_docs_match_the_code.py` checks that this file is current.*

`fue` 0.1.11

---

## Building a model

### `TimeSeries(data, freq=12, start=(1900, 1), name='series')`

A univariate time series with frequency and start-date metadata.

Parameters
----------
data : array-like
    Observations in chronological order.
freq : int
    Observations per year: 1 (annual), 4 (quarterly), 12 (monthly).
    Pass 0 to use plain observation numbering.
start : tuple (year, period)
    First observation.  Period is 1-based (Jan = 1 for monthly).
name : str, optional
    Series label used in summaries and plots.

#### `TimeSeries.plot(title=None, ax=None)`

Time-series line plot with calendar x-axis.

#### `TimeSeries.plot_acf(lags=24, confidence=0.95, ax=None)`

Autocorrelation function stem plot.

#### `TimeSeries.plot_pacf(lags=24, confidence=0.95, ax=None)`

Partial autocorrelation function stem plot.

#### `TimeSeries.describe()`

Sample statistics matching fue's File_StatSer output.

Uses population moments (divisor n) to match the C implementation.
Returns the formatted string (also printed to stdout).

### `Model(series, ar=None, ma=None, ar_s=None, ma_s=None, ar_free=None, ma_free=None, ar_s_free=None, ma_s_free=None, ar_f=None, ma_f=None, d=0, D=0, ifadf=None, interventions=None, mu=0.0, estimate_mu=False, boxlam=1.0, refactor=1.0, eml=True, chkma=True)`

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

#### `Model.add_intervention(type, at, omega=None, delta=None, omega_free=None, delta_free=None)`

Return a new Model with one extra intervention appended.

#### `Model.fit()`

Estimate model parameters by exact maximum likelihood.

Sets self._result and returns self (for chaining).
Raises RuntimeError if estimation returns a non-zero ifault.

#### `Model.forecast_fuf(horizon=None, sigma2=None)`

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

#### `Model.forecast(horizon)`

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

#### `Model.compare(*others)`

Print a comparison table of fitted models.

Parameters
----------
*others : Model
    Additional fitted models to compare against *self*.

Returns
-------
str
    Formatted table (also printed to stdout).

#### `Model.summary()`

*(no docstring)*

#### `Model.write_out(path=None, inp_name='', out_name='')`

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

#### `Model.write_pre(path)`

Write a .pre file with estimated parameters as new initial values.

Parameters
----------
path : str
    Output path, e.g. "RIPC.1.pre".

#### `Model.write_fuf(horizon, sigma2=None, path=None)`

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

#### `Model.write_fuf_out(fr, path=None, inp_name='', out_name='')`

Generate a forecast report in fuf .out format.

Parameters
----------
fr : ForecastResult  (from model.forecast_fuf())
path : str or None
    Write to file; return as string if None.
inp_name, out_name : str
    Optional file-name labels shown in the header.

#### `Model.plot_residuals(lags=None)`

*(no docstring)*

### `Intervention(type, at=0, omega=None, delta=None, omega_free=None, delta_free=None, harmonic=1.0, data=None)`

Deterministic component with linear transfer function.

The effect on the series is  ω(B)/δ(B) · x_t  where x_t is a binary
indicator determined by *type* and *at*.

Parameters
----------
type : str
    ``'impulse'``  — isolated impulse at *at*.  ``'pulse'`` is accepted as a
                     deprecated alias and normalised to ``'impulse'``: the
                     school's vocabulary — and the one in the `.pre`/`.inp`
                     format and in fue C — is ``impulse``, and having two
                     names for one thing is how a file ends up with a
                     keyword the other interpreter does not know.
    ``'compimp'``  — COMPENSATED impulse: +1 at *at* and **−1 at *at*+1**
                     (``compimp`` in the .pre/.inp format).  A pulse that is
                     undone the next period: the level returns to where it
                     was, so it does not shift the mean of a differenced
                     series the way a plain pulse does.
    ``'step'``     — permanent level shift starting at *at*
    ``'ramp'``     — linear ramp starting at *at*
    ``'seasonal'`` — periodic seasonal dummy (*at* = 0-based period within
                     year).  **Python-only, and outside the methodology:**
                     deterministic seasonality is parameterised with
                     HARMONICS (``cos``/``sin`` plus the Nyquist
                     ``alter``), not with dummies, which is why fue C has no
                     such regressor and the format has no keyword for it.
                     Writing one raises rather than emit a word fue C would
                     silently take for a non-standard variable.
    ``'easter'``   — Easter-holiday variable; **monthly series only**
                     (freq == 12), *at* unused.  See `_build_indicator`.
    ``'trend'``    — deterministic linear trend, 1, 2, …, n; *at* unused
    ``'cos'``      — cosine component cos(2π·harmonic/freq·j); *at* unused
    ``'sin'``      — sine component   sin(2π·harmonic/freq·j); *at* unused
    ``'alter'``    — alternating sign (-1)^j; *at* unused
    ``'custom'``   — external indicator supplied as *data* array
at : int
    0-based observation index for pulse/compimp/step/ramp (0 = first
    observation); 0-based period within year for seasonal.  Unused for
    easter/trend/cos/sin/alter/custom.
omega : list of float
    Numerator polynomial coefficients [ω₀, ω₁, …].  Default ``[1.0]``.
delta : list of float
    Denominator polynomial coefficients [δ₁, δ₂, …].  Default ``[]``
    (no denominator → pure FIR).
omega_free : list of bool, optional
    Which omega coefficients to estimate.  Defaults to all True.
delta_free : list of bool, optional
    Which delta coefficients to estimate.  Defaults to all True.
data : array-like, optional
    Pre-computed indicator values, length nobs.  Required for type='custom'.

### `FixedFreqFactor(freq, coef=-0.5, free=True)`

Second-order AR or MA factor with fixed spectral frequency.

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


## Reading and writing files

### `load(path)`

Parse a fue .inp file and return (TimeSeries, Model).

The returned Model is unfitted; call .fit() to estimate parameters.
If the file is in fuf format (contains the forecast horizon/sigma2 section),
the extra fields are stored in model._fuf_horizon and model._fuf_sigma2.

Parameters
----------
path : str or path-like
    Path to the .inp file (with or without the .inp extension).

Returns
-------
ts : TimeSeries
model : Model  (unfitted)

### `load_fuf(path)`

Parse a fuf forecast specification file and return (TimeSeries, Model).

fuf files are like fue .inp files but contain an extra section after the
observations line: "** Forecast horizon and estimated innovation variance"
with two values: L (forecast horizon) and sigma2 (estimated variance).
All parameter values in the file are treated as pre-estimated (fixed).

Parameters
----------
path : str or path-like
    Path to the fuf .inp file.

Returns
-------
ts : TimeSeries
model : Model  (unfitted; call model.forecast_fuf() to get forecasts)

### `write_out(model, path=None, inp_name='', out_name='')`

Generate an estimation report in fue .out format.

Parameters
----------
model : Model  (must be fitted)
path : str or None
    If given, write to this file path.  If None, return the text.
inp_name : str
    Label shown in "Input file" header line.
out_name : str
    Label shown in "Output file" header line.  If empty and *path* is
    given, the basename of *path* is used.

Returns
-------
str

### `write_fuf(model, horizon, sigma2, path=None)`

Write a fuf forecast input file.

Same format as the .pre file but with the FUF header and an extra
"Forecast horizon / sigma2" section after the observations line.

Works with both fitted and unfitted models: fitted params are written
for a fitted model; initial param values for an unfitted one.

Parameters
----------
model : Model
horizon : int
    Steps ahead to forecast.
sigma2 : float
    Estimated innovation variance.
path : str or None
    Write to file; return as string if None.

### `write_fuf_out(model, fr, path=None, inp_name='', out_name='')`

Generate a forecast report in fuf .out format.

The output mirrors fuf-1.08.1's output: a forecast table (observed history
+ future forecasts) followed by residual diagnostics, ACF, PACF, and the
calibration of distortions table.

Parameters
----------
model : Model  (fitted or unfitted; parameters read from initial values)
fr : ForecastResult  (from model.forecast_fuf())
path : str or None
    Write to file; return as string if None.
inp_name : str
    Name shown in the "Input file" header line.
out_name : str
    Name shown in the "Output file" header line.

Returns
-------
str

### `write_forecast_report(model, fr, path, title=None, source=None, sps_name=None, narrative=None, pdf=False)`

Write a self-contained HTML forecast report.


## Diagnostics

### `acf(data, lags=24)`

Sample autocorrelation function.

Returns array of length *lags* with r[k] = Corr(x_t, x_{t-k}).

### `pacf(data, lags=24)`

Partial autocorrelation function via Durbin-Levinson recursion.

Returns array of length *lags*.

### `ljung_box(data, lags=None, df_correction=0)`

Ljung-Box portmanteau test.

Parameters
----------
data : array-like
    Residuals.
lags : int or list of int
    Lag(s) at which to compute the test.  Defaults to min(10, nobs//5).
df_correction : int
    Number of estimated ARMA parameters (subtracted from degrees of freedom).

Returns
-------
dict with keys 'statistic', 'pvalue', 'lags'.

### `jarque_bera(data)`

Jarque-Bera normality test.

Returns (statistic, p-value).


## Results

### `ForecastResult(horizon: int, level: numpy.ndarray, level_std: numpy.ndarray, diff1: numpy.ndarray, diff1_std: numpy.ndarray, seasonal_diff: numpy.ndarray, seasonal_diff_std: numpy.ndarray, sigma2: float) -> None`

Point forecasts and standard errors from Model.forecast().


## Datasets

Shipped with the package; `from fue.datasets import ripc`.

### `fue.datasets.ripc() -> fue.series.TimeSeries`

RIPC, monthly, January 2002 – December 2007 (72 observations).

Standard monthly test case for fue, used to verify seasonal ARMAX
estimation with Fourier harmonics and alternator interventions. It is the
series of `tests/real_cases/.../RIPC.1.inp`, byte for byte.

⚠ The values are the series **as fue reads it**, around 0.41–0.44 — NOT a
transformed variable. The canonical model applies the transformation
itself, with `boxlam=0` and `refactor=100` (that is, 100·log) and `d=0`;
passing those to `fue.Model` is what reproduces `RIPC.1`.

This docstring said the opposite until 2026-08-13 — "the series is the log
of the Spanish CPI rescaled by 100" — which would have meant applying the
transformation twice. It is wrong in a way that estimates cleanly and reads
plausibly, which is the kind that survives.

Returns
-------
TimeSeries
    Monthly series (freq=12), start=(2002, 1), name="RIPC".

### `fue.datasets.sfny() -> fue.series.TimeSeries`

SFNY annual precipitation index, 1852–1913 (62 observations).

The series is a sunspot-New York precipitation proxy used as a
standard test case for the FUE estimation engine (Mauricio 1995,
JASA §4 example SFNY.2).

The recommended model is an ARMAX with a level shift at 1853:

    log(y_t) = ω/(1 − δB) · S_t  +  AR(1) × AR(2)  +  μ  +  ε_t

where S_t is a step function starting at t=2 (1853).

Returns
-------
TimeSeries
    Annual series (freq=1), start=(1852, 1), name="SFNY".


## Modules

Re-exported at package level; their contents are reached through the objects above.

* **`fue.datasets`** — Built-in datasets for fue examples and tests.
* **`fue.diagnostics`** — Diagnostic statistics: ACF, PACF, Jarque-Bera, Ljung-Box.
* **`fue.forecast`** — ARMAX forecast engine — pure Python implementation.
* **`fue.inp`** — Parser for fue .inp model-data files.
* **`fue.intervention`** — Intervention: linear transfer function ω(B)/δ(B) applied to an indicator.
* **`fue.model`** — Model: ARMAX specification with interventions, fitted by exact ML via FUE.
* **`fue.report`** — Generate .out-style estimation reports matching fue's ASCII output format.
* **`fue.report_forecast`** — HTML forecast report generator — SPS (Sistema de Previsión y Seguimiento).
* **`fue.series`** — TimeSeries: lightweight wrapper around a numpy array with date metadata.
