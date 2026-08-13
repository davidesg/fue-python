"""Example 2 — the airline model, (0,1,1)(0,1,1)_12.

The canonical Box-Jenkins seasonal model, and the baseline that examples 4 and 5
are measured against: it treats seasonality as ONE operator with ONE parameter.

    python examples/02_airline.py
"""
import numpy as np

import fue
from fue.datasets import ripc

ts = ripc()

# (0,1,1)(0,1,1)_12 on 100·log: one regular difference, one annual difference,
# one regular MA and one seasonal MA.  No mean: after ∇∇_12 the drift is gone.
m = fue.Model(ts,
              d=1, D=1, boxlam=0.0, refactor=100.0,
              ma=[[0.4]], ma_free=[[True]],
              ma_s=[[0.6]], ma_s_free=[[True]],
              estimate_mu=False)
m.fit()
r = m._result

print("airline (0,1,1)(0,1,1)_12 on 100·log(CPI)")
print(f"  theta   = {m.ma[0][0]:9.6f}  ({r.std_errors[0]:.6f})")
print(f"  Theta   = {m.ma_s[0][0]:9.6f}  ({r.std_errors[1]:.6f})")
print(f"  sigma   = {np.sqrt(r.sigma2):9.6f}")
print(f"  logL    = {m.loglik:9.6f}    BIC {r.bic:.3f}")
print(f"  stopped by {r.termination}, {r.niter} iterations")

# What the annual difference costs, and what it assumes.
#
# ∇_12 = (1−B)·(1+B+…+B^11): besides the annual moving sum it applies a REGULAR
# difference, which has nothing to do with seasonality.  And the single Theta
# forces every seasonal frequency to share one root modulus.  Both are what
# examples 4 and 5 take apart.
perdidas = 1 + 12
print(f"\n  observations: {ts.nobs} − {perdidas} (∇ and ∇_12) = {ts.nobs - perdidas} usable")
print("  one parameter, Theta, is doing the work of six frequencies")

if abs(m.ma_s[0][0]) > 0.9:
    print(f"\n  ⚠ Theta = {m.ma_s[0][0]:.4f} is near the non-invertibility boundary.")
    print("    That is where the DCD test belongs — docs/FORMAL_TESTS.md §3 —")
    print("    and where the exact-ML and Box-Jenkins criteria pull apart.")
