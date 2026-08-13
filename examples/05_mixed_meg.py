"""Example 5 — a mixed MEG: some frequencies deterministic, others stochastic.

This is the class of model no other program in this family can specify. The
airline model of example 2 forces every seasonal frequency to share one root;
the annual difference of example 4 makes every frequency non-stationary at
once. Neither can say what is usually true of a real series: **some seasonal
frequencies evolve and others do not.**

The series here has stochastic seasonality at f=1 only, and a fixed harmonic at
f=2. The example specifies exactly that.

    python examples/05_mixed_meg.py
"""
import numpy as np

import fue

n, s = 360, 12
SIGMA = 1.0
rng = np.random.RandomState(20260813)

# f = 1 STOCHASTIC: a cycle whose amplitude wanders because its factor sits on
# the unit circle — (1 − 2cos(2π/12)B + B²) applied to noise, i.e. integrated at
# that frequency.
c1 = 2 * np.cos(2 * np.pi * 1 / s)
u = rng.normal(0.0, 0.30, n)
x1 = np.zeros(n)
for t in range(2, n):
    x1[t] = c1 * x1[t - 1] - x1[t - 2] + u[t]

# f = 2 DETERMINISTIC: fixed amplitude, forever.  Index from 1: fue's convention.
t = np.arange(1, n + 1)
x2 = 2.5 * np.cos(2 * np.pi * 2 * t / s) + 1.0 * np.sin(2 * np.pi * 2 * t / s)

a = rng.normal(0.0, SIGMA, n)
y = 100.0 + x1 + x2 + a
ts = fue.TimeSeries(list(y), freq=s, start=(2000, 1), name="MIXED")

# NOTE mu.  Slot 0 of ifadf is OFF, so nothing removes the level: with
# estimate_mu=False the engine would fit a series of level 100 as if it had mean
# zero and report sigma = 100.  Integrating at f=1 says nothing about f=0.
#
# The specification: ifadf slot 1 ON (integrate at f=1), everything else OFF,
# and f=2 carried by a deterministic harmonic pair.  An MA_f at f=1 is the
# over-differencing witness that docs/FORMAL_TESTS.md §5 asks for.
ifadf = [0, 1, 0, 0, 0, 0, 0]

m = fue.Model(ts, d=0, D=0, ifadf=ifadf,
              ma_f=[fue.FixedFreqFactor(freq=1.0, coef=-0.81)],
              interventions=[
                  fue.Intervention("cos", harmonic=2.0, omega=[0.1]),
                  fue.Intervention("sin", harmonic=2.0, omega=[0.1]),
              ],
              mu=float(np.mean(y)), estimate_mu=True)
m.fit()
r = m._result

print("MEG: stochastic at f=1, deterministic at f=2\n")
print(f"  ifadf = {ifadf}   (slot 1 on: integrate at f=1 only)")
print(f"  observations lost: 2  (one interior factor), against {s} for ∇_12\n")
print(f"  cos(f=2) = {m.interventions[0].omega[0]:8.4f}   true 2.50")
print(f"  sin(f=2) = {m.interventions[1].omega[0]:8.4f}   true 1.00")
print(f"  sigma    = {np.sqrt(r.sigma2):8.4f}   true {SIGMA:.2f}")
print(f"  logL     = {m.loglik:.4f},  stopped by {r.termination}")

# The witness: an MA_f at the integrated frequency.  If it comes out
# non-invertible, the integration was not needed — the seasonality at that
# frequency is deterministic after all.  That is the MEG decision rule, and the
# DCD test is what makes it formal.
# The parameter estimated is `coef` = −r², with r the spectral radius of the
# factor; r = 1 is the non-invertibility boundary.
coef = m.ma_f[0].coef
r_mod = np.sqrt(-coef) if coef < 0 else float("nan")
print(f"\n  MA_f witness at f=1: coef = {coef:.4f}  →  radius r = {r_mod:.4f}")
veredicto = ("NEAR the boundary (r→1): the witness says over-differencing"
             if r_mod > 0.95 else
             "well inside (r<0.95): the integration at f=1 is doing work")
print(f"    {veredicto}")
print("\n  The formal decision is the DCD test of docs/FORMAL_TESTS.md §3, with the")
print("  complex-pair critical values (1.11/2.04/4.52 asymptotically) — NOT the")
print("  real-root ones, because an interior frequency is a second-order factor.")
