"""Example 1 — the minimal flow: load, fit, read.

ARIMA(1,1,0) with a mean on a monthly price index. Everything the rest of the
examples do is this, with more model.

    python examples/01_minimal_arima.py
"""
import numpy as np

import fue
from fue.datasets import ripc

ts = ripc()                      # RIPC, monthly, 2002-2007, n = 72
                                 # the series as fue reads it (~0.41-0.44):
                                 # the transformation below is fue's, not the
                                 # data's -- see fue.datasets.ripc
print(f"{ts.name}: n={ts.nobs}, freq={ts.freq}, from {ts.start[1]}/{ts.start[0]}")

# 100·log(y) is the usual working scale for a price index: differences of the
# log are relative changes, and ×100 puts them in percent.  fue does the
# transformation itself -- boxlam=0 is the log -- and `refactor` is the ×100.
#
# d=1 because a price index is I(1): the level wanders, the monthly change does
# not.  mu is the mean of the DIFFERENCED variable, i.e. average inflation per
# month, so it is estimated rather than fixed at zero.
w = np.diff(np.log(np.asarray(ts.data)) * 100.0)

m = fue.Model(ts,
              d=1, boxlam=0.0, refactor=100.0,
              ar=[[0.3]], ar_free=[[True]],
              mu=float(w.mean()), estimate_mu=True)
m.fit()
r = m._result

print(f"\n  phi_1  = {m.ar[0][0]:9.6f}  ({r.std_errors[0]:.6f})")
print(f"  mu     = {r.params[-1]:9.6f}  ({r.std_errors[-1]:.6f})   % per month")
print(f"  sigma  = {np.sqrt(r.sigma2):9.6f}")
print(f"  logL   = {m.loglik:9.6f}    AIC {r.aic:.3f}   BIC {r.bic:.3f}")

# Always read the verdict, not just the numbers: `converged` is now the honest
# claim (gradient criterion), and anything else says why it stopped.
print(f"\n  stopped by: {r.termination}  after {r.niter} iterations, |g|={r.gnorm:.2e}")
print(f"  converged : {r.converged}")

lb = fue.ljung_box(r.residuals, lags=12, df_correction=r.npar)
q, p_val = lb["statistic"][0], lb["pvalue"][0]
print(f"\n  Ljung-Box(12) = {q:.2f}   p = {p_val:.3f}"
      f"   {'residuals look white' if p_val > 0.05 else 'STRUCTURE LEFT'}")
