# fue — Exact Maximum Likelihood for Univariate Time Series

**fue** is a Python implementation of the FUE/FUF estimation engine originally
written in C by Arthur B. Treadway and David E. Guerrero, based on the
algorithms designed and coded by José Alberto Mauricio.
It fits ARMAX models with linear transfer-function interventions by exact
maximum likelihood using the Ansley (1979) innovations form and the
Mauricio (1997, AS 311) / Mauricio (1995, JASA) algorithms.

## Features

- Exact ML estimation of ARIMA/SARIMA/ARMAX models (Box–Jenkins family)
- Seasonal and non-seasonal AR/MA operators; regular and seasonal differencing
- Linear transfer function interventions: impulse, step, ramp, Fourier
  harmonics, alternator, and arbitrary user-defined regressors
- Fixed-frequency AR/MA factors (AR(2) with constrained spectral peak)
- Box–Cox transformation with automatic back-transformation of forecasts
- Multi-step forecasts with asymptotic prediction intervals
- ASCII `.out` and `.pre` report output compatible with the C FUE binary
- HTML forecast reports (requires `jinja2`)
- CLI tools `fue` and `fuf` mirroring the C binaries

## Installation

```bash
pip install fue                  # pure Python (numpy + scipy + matplotlib)
pip install "fue[report]"        # + HTML forecast reports (jinja2)
pip install "fue[pdf]"           # + PDF export (weasyprint)
```

Python 3.10+ required.

## Quick start

```python
import numpy as np
import fue

# Build a time series
data = np.array([...])   # monthly observations
ts = fue.TimeSeries(data, freq=12, start=(2002, 1), name="CPI")

# Specify an ARIMA(1,1,0)(1,0,0)_12 model with a log transform
m = fue.Model(ts,
              ar=[[0.3]],          # AR(1), initial value 0.3
              sar=[[0.2]],         # seasonal AR(1) at lag 12
              d=1,                 # one regular difference
              boxlam=0.0,          # log transform (Box–Cox λ=0)
              refactor=100.0)

result = m.fit()
print(result.sigma2, result.aic)
m.plot_residuals()
```

### Load from a `.inp` file (fue format)

```python
ts, m = fue.load("model.inp")
result = m.fit()
m.write_out("model.out")
m.write_pre("model.pre")
```

### Forecast from a pre-estimated `.inp` file (fuf format)

```python
ts, m = fue.load_fuf("forecast_model.inp")
fr = m.forecast_fuf(horizon=24)
print(fr.level)          # point forecasts in original units
print(fr.level_std)      # forecast standard deviations
```

### Command-line interface

```bash
# Estimate a model and write .out / .pre
fue model_name [eml|aml] [chk|nochk] [-f horizon]

# Generate forecasts from a pre-estimated .inp file
fuf forecast_model
```

## Numerical methods

| Algorithm | Reference | Used for |
|-----------|-----------|----------|
| Ansley (1979) innovations form | Mauricio (1997) AS 311 | Exact log-likelihood (`elf_scalar`) |
| Kalman filter (quick recursions) | Mélard (1984) AS 197 | Inner BFGS loop (`flikam_scalar`) |
| BFGS with Cholesky factor update | Dennis & Schnabel (1983) ch. 9 | Optimization (`raxopt`) |
| Scaled objective Π(x)/Π₀ | Mauricio (1995) JASA §3 | Numerical conditioning |

## Bug tracking

Bugs are tracked in-repo under [`bugs/`](bugs/README.md) — one Markdown file per
report (`BUG-NNNN-slug.md`) with a small frontmatter schema. A fix references the
id in its commit, e.g. `fix(forecast): BUG-0001 …`.

```bash
fue-bug list                       # list reports (open marked with *)
fue-bug show BUG-0001              # print a report
fue-bug new "title" --component forecast   # file a new report
fue-bug check                      # validate all reports (runs in CI, tests/test_bugs.py)
fue-bug index                      # regenerate bugs/README.md
```

## Authors and licence

**fue** is developed by Arthur B. Treadway and David E. Guerrero, based on the
algorithms designed and coded by José Alberto Mauricio.

Released under the **GNU General Public Licence v2.0 or later** (GPL-2.0-or-later).
See [COPYING](COPYING) for the full licence text.
