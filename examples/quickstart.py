"""
fue quickstart — full Box-Jenkins workflow with SFNY data
=========================================================

Demonstrates:
  1. Loading a built-in dataset
  2. ACF/PACF exploration
  3. Fitting a reference AR(1) model
  4. Fitting the full SFNY.2 ARMAX model (level shift + AR(1)×AR(2) + log)
  5. Model comparison (AIC/BIC)
  6. Residual diagnostics
  7. Multi-step forecasts with prediction intervals
  8. Writing an HTML forecast report  (requires: pip install "fue[report]")

Run
---
    python examples/quickstart.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")          # remove if running interactively
import matplotlib.pyplot as plt

import fue
from fue import Intervention, Model
from fue.datasets import sfny

# ── 1. Data ───────────────────────────────────────────────────────────────────

ts = sfny()
print(ts)
# TimeSeries('SFNY', nobs=62, freq=1, start=(1852, 1))

# ── 2. ACF / PACF ─────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
ts.plot_acf(lags=20, ax=axes[0])
ts.plot_pacf(lags=20, ax=axes[1])
fig.suptitle("SFNY — ACF / PACF of raw series", fontsize=10)
plt.tight_layout()
plt.savefig("sfny_acf.png", dpi=120, bbox_inches="tight")
print("Saved sfny_acf.png")

# ── 3. Reference model: AR(1) ─────────────────────────────────────────────────

m_ar1 = Model(ts, ar=[[0.5]])
m_ar1.fit()
print("\n--- AR(1) ---")
print(m_ar1.summary())

# ── 4. SFNY.2: ARMAX with level shift ─────────────────────────────────────────
#
# Model:  log(y_t) = ω/(1-δB)·S_t  +  φ₁(B)·φ₂(B)·y_t  +  ε_t
#
#   S_t       : step function starting at t=2 (1853), estimated transfer function
#   AR(1)×AR(2): product of a first- and second-order AR operator
#   boxlam=0.0: log transformation applied before fitting

step = Intervention(
    "step",
    at=1,                    # step starts at observation index 1 (= 1853)
    omega=[0.08],            # initial ω₀
    omega_free=[True],
    delta=[0.6],             # initial δ (decay)
    delta_free=[True],
)

m_sfny2 = Model(
    ts,
    interventions=[step],
    ar=[[0.8], [-0.1, -0.1]],   # AR(1) × AR(2) — two separate factors
    boxlam=0.0,                  # log transform
)
m_sfny2.fit()
print("\n--- SFNY.2 ---")
print(m_sfny2.summary())

# ── 5. Model comparison ───────────────────────────────────────────────────────

print("\n--- Model comparison ---")
m_ar1.compare(m_sfny2)

# ── 6. Residual diagnostics ───────────────────────────────────────────────────

from fue.plots import plot_model_diagnostics

fig_diag, fig_hist = plot_model_diagnostics(m_sfny2)
fig_diag.savefig("sfny_diag.png", dpi=120, bbox_inches="tight")
fig_hist.savefig("sfny_hist.png", dpi=120, bbox_inches="tight")
print("\nSaved sfny_diag.png, sfny_hist.png")

# ── 7. Forecast ───────────────────────────────────────────────────────────────

HORIZON = 10
fr = m_sfny2.forecast(horizon=HORIZON)

print(f"\n--- {HORIZON}-step forecast (original scale) ---")
print(f"{'h':>4}  {'level':>10}  {'±1σ':>10}")
for h in range(HORIZON):
    print(f"{h+1:>4}  {fr.level[h]:>10.5f}  {fr.level_std[h]:>10.5f}")

# ── 8. HTML report ────────────────────────────────────────────────────────────

try:
    from fue.report_forecast import write_forecast_report
    write_forecast_report(
        m_sfny2, fr,
        path="sfny_forecast.html",
        title="SFNY Annual Precipitation Index",
        source="Mauricio (1995)",
        sps_name="fue example",
    )
    print("\nSaved sfny_forecast.html")
except ImportError:
    print("\nSkipping HTML report — pip install 'fue[report]' to enable")
