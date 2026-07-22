"""
BUG-0004 reproduction — forecast_fuf forecasts from the model's STALE pre-fit
attributes (ar / ar_s / mu0), not the fitted parameters in _result.params.

Run from this folder:  python repro.py
Needs: fue (>=0.1.7 to see the bug) and art (only to seed a realistic model).

The euro-area HICP data is embedded in EA_CPI.pre; we load the series from it and
rebuild the model with the normal art seeds. .fit() estimates correctly, but the
model's ar/ar_s/mu0 attributes are left at their pre-fit seeds; forecast_fuf reads
those, so the level runs away. Writing the fitted model to .pre (from _result) and
reloading gives the correct forecast, which is why the .pre reproduces nothing.
"""
import numpy as np
import fue
from art.pipeline import _make_model               # only to seed a realistic model

ld = fue.load("EA_CPI.pre")
ts = ld[0] if isinstance(ld, tuple) else ld.series  # euro-area HICP series

m = _make_model(ts, lam=0.0, d=1, D=0, p=1, q=0, n_harmonics=5, P=1, Q=0, estimate_mu=True)
m.fit()
p = m._result.params
print("fitted (correct, _result.params):  phi=%.4f  Phi=%.4f  mu=%.5f" % (p[-3], p[-2], p[-1]))
print("model attributes AFTER .fit (STALE): phi=%.4f  Phi=%.4f  mu=%.5f" % (m.ar[0][0], m.ar_s[0][0], m.mu0))

lvl = np.asarray(m.forecast_fuf(horizon=6).level, float)
print("forecast_fuf level h1..6 (uses stale attrs -> EXPLODES):", np.round(lvl, 1),
      " | last obs %.2f" % float(ts.data[-1]))

m.ar = [[p[-3]]]; m.ar_s = [[p[-2]]]; m.mu0 = float(p[-1])   # sync fitted -> attributes
lvl2 = np.asarray(m.forecast_fuf(horizon=6).level, float)
print("forecast_fuf level h1..6 (after sync -> CORRECT):     ", np.round(lvl2, 1))

ld2 = fue.load("EA_CPI.pre"); ml = ld2[1] if isinstance(ld2, tuple) else ld2
print("fue.load('EA_CPI.pre').forecast_fuf (CORRECT):        ",
      np.round(np.asarray(ml.forecast_fuf(horizon=6).level, float), 1))
