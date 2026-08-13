"""Example 4 — deterministic harmonics against the annual difference.

Two ways of handling seasonality that are NOT the same hypothesis:

  * deterministic — the seasonal pattern is fixed: sines and cosines, estimated
    as regressors, nothing differenced;
  * stochastic  — the pattern evolves: the annual difference ∇_12, which is a
    non-stationarity assumption at every seasonal frequency at once.

The example fits both to a series whose seasonality is DETERMINISTIC by
construction, and shows what each costs.

    python examples/04_harmonics_vs_annual_difference.py
"""
import numpy as np

import fue

n, s = 240, 12
PHI, SIGMA = 0.6, 1.0

rng = np.random.RandomState(20260813)
a = rng.normal(0.0, SIGMA, n)
N = np.zeros(n)
N[0] = rng.normal(0.0, SIGMA / np.sqrt(1 - PHI ** 2))
for t in range(1, n):
    N[t] = PHI * N[t - 1] + a[t]

# A FIXED seasonal pattern: two harmonics that never change amplitude.
#
# NOTE THE INDEX.  fue builds its harmonics as cos(2π·k·j/s) with j = 1…n
# (`fue_api.c:806`), not from j = 0.  Simulating from 0 shifts the phase by one
# period: the amplitude survives but the (cos, sin) pair rotates into each
# other, and the example then appears to recover 1.78/3.02 where the truth was
# 3.00/1.50.  Same model, different origin.
t = np.arange(1, n + 1)
seasonal = (3.0 * np.cos(2 * np.pi * 1 * t / s) + 1.5 * np.sin(2 * np.pi * 1 * t / s)
            + 2.0 * np.cos(2 * np.pi * 2 * t / s))
y = 50.0 + seasonal + N

ts = fue.TimeSeries(list(y), freq=s, start=(2000, 1), name="SEAS")

# ── (a) deterministic: harmonics as regressors, no differencing ────────────
det = []
for f in (1, 2, 3, 4, 5):
    det.append(fue.Intervention("cos", harmonic=float(f), omega=[0.1]))
    det.append(fue.Intervention("sin", harmonic=float(f), omega=[0.1]))
det.append(fue.Intervention("alter", omega=[0.1]))          # f = s/2, Nyquist

m_det = fue.Model(ts, d=0, ar=[[0.4]], ar_free=[[True]],
                  interventions=det, mu=float(y.mean()), estimate_mu=True)
m_det.fit()

# ── (b) stochastic: the annual difference ─────────────────────────────────
m_dif = fue.Model(ts, d=0, D=1,
                  ma_s=[[0.5]], ma_s_free=[[True]],
                  ar=[[0.4]], ar_free=[[True]], estimate_mu=False)
m_dif.fit()

print("the same series, two hypotheses about its seasonality\n")
print("⚠ the two log-likelihoods are NOT comparable: (b) differences the series,")
print("  so it is a likelihood for 228 observations of a different variable.")
print("  Comparing them — or their BIC — is the commonest error in this business.\n")
print(f"{'':32} {'logL':>11} {'BIC':>10} {'npar':>5} {'obs used':>9}")
print("-" * 72)
print(f"{'(a) deterministic harmonics':32} {m_det.loglik:11.4f} "
      f"{m_det._result.bic:10.3f} {m_det._result.npar:5d} {n:9d}")
print(f"{'(b) annual difference ∇_12':32} {m_dif.loglik:11.4f} "
      f"{m_dif._result.bic:10.3f} {m_dif._result.npar:5d} {n - s:9d}")

print("\nthe two harmonics that are really there, as fue recovers them:")
for iv, verdad in zip(m_det.interventions[:4], [3.0, 1.5, 2.0, 0.0]):
    print(f"  {iv.type}(f={iv.harmonic:.0f}) = {iv.omega[0]:8.4f}   true {verdad:5.2f}")

print(f"\n  ∇_12 costs {s} observations and assumes non-stationarity at SIX")
print(f"  frequencies at once.  Here the truth is deterministic, so it is the")
print(f"  wrong hypothesis — but no likelihood comparison will tell you that,")
print(f"  because the two models do not describe the same variable.")
print(f"  What decides it is the MEG sweep of docs/FORMAL_TESTS.md §5, frequency")
print(f"  by frequency, which is example 05.")

# The equivalence worth knowing, and the distinction it rests on:
#
#   ifadf = [1,1,1,1,1,1,1]  →  (1−B)·∏(seasonal factors)  =  ∇_12   IDENTICAL
#   ifadf = [0,1,1,1,1,1,1]  →  the annual moving sum alone, degree 11
#
# The f=0 slot IS the regular difference.  Turning it off is what separates
# "seasonality" from "the trend that ∇_12 drags along with it".
m_all = fue.Model(ts, d=0, D=0, ifadf=[1, 1, 1, 1, 1, 1, 1],
                  ma_s=[[0.5]], ma_s_free=[[True]],
                  ar=[[0.4]], ar_free=[[True]], estimate_mu=False)
m_all.fit()
m_sum = fue.Model(ts, d=0, D=0, ifadf=[0, 1, 1, 1, 1, 1, 1],
                  ma_s=[[0.5]], ma_s_free=[[True]],
                  ar=[[0.4]], ar_free=[[True]], estimate_mu=False)
m_sum.fit()
print(f"\n  ifadf all set   logL {m_all.loglik:.6f}   vs  D=1  {m_dif.loglik:.6f}"
      f"   difference {abs(m_all.loglik - m_dif.loglik):.1e}   ← the same operator")
print(f"  moving sum only logL {m_sum.loglik:.6f}"
      f"   ← a DIFFERENT operator: no (1−B), one observation more")
