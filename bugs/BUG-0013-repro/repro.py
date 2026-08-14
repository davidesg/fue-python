#!/usr/bin/env python3
"""BUG-0013 — the C engine segfaults on a model with no ARMA factors.

Three lines are enough. A regression on deterministic inputs with white-noise
errors — ARMA(0,0), which is a legitimate specification and the first rung of
every seasonal ladder — kills the interpreter: no exception, no message, no
traceback. `python repro.py` prints the diagnosis and then dies on the last fit
unless you pass --safe.

    python repro.py            # crashes, on purpose
    python repro.py --safe     # shows the two specifications that DO work
"""
import sys
import warnings

import numpy as np

warnings.simplefilter("ignore")
import fue

rng = np.random.RandomState(20260814)
n, s = 240, 12
t = np.arange(1, n + 1)
y = (100.0 + 3.0 * np.cos(2 * np.pi * t / s) + 1.5 * np.sin(2 * np.pi * t / s)
     + rng.normal(0, 1.0, n))
ts = fue.TimeSeries(list(y), freq=s, start=(2000, 1), name="DET")


def armonicos():
    return [fue.Intervention("cos", harmonic=1.0, omega=[0.1]),
            fue.Intervention("sin", harmonic=1.0, omega=[0.1])]


print(f"fue {fue.__version__}\n")

# ── what works, and is the same model ──────────────────────────────────────
m_fijo = fue.Model(ts, d=0, interventions=armonicos(),
                   ar=[[0.0]], ar_free=[[False]],      # one factor, pinned at 0
                   mu=float(y.mean()), estimate_mu=True)
m_fijo.fit()
print(f"1. ar=[[0.0]] fixed  ->  OK   npar={m_fijo._result.npar}  "
      f"logL={m_fijo.loglik:.4f}")

from fue.cast_us import estimate_py
from fue.model import FitResult

m_py = fue.Model(ts, d=0, interventions=armonicos(),
                 mu=float(y.mean()), estimate_mu=True)
r = FitResult(estimate_py(m_py))
print(f"2. ar=[] on the PYTHON engine  ->  OK   npar={r.npar}  "
      f"logL={r.loglik:.4f}")

print("\n   Same model three ways. The first two agree; the third is below.")
print("   The difference is not the mathematics — it is whether the")
print("   specification carries an ARMA factor at all.\n")

if "--safe" in sys.argv:
    print("3. ar=[] on the C engine  ->  NOT RUN (--safe)")
    sys.exit(0)

print("3. ar=[] on the C engine  ->  fitting now; the interpreter dies here")
sys.stdout.flush()
m_c = fue.Model(ts, d=0, interventions=armonicos(),
                mu=float(y.mean()), estimate_mu=True)
m_c.fit()
print("   ...if you are reading this, the bug is fixed.")
