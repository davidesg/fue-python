"""Example 3 — ARMAX: interventions with a date. THE HALF THAT IS THE THESIS.

Exact unconditional ML with intervention inputs is what `fue` has and the
alternatives do not. To make it checkable rather than merely illustrative, the
series is SIMULATED with known effects and the example reports what fue
recovers.

    y_t = 12·step(t ≥ 60) + (−8)·impulse(t = 100) + N_t ,   N_t ~ AR(1), φ=0.7

    python examples/03_armax_interventions.py
"""
import numpy as np

import fue

n, t_step, t_pulse = 600, 200, 400
OMEGA_STEP, OMEGA_PULSE, PHI, SIGMA = 12.0, -8.0, 0.7, 2.0

rng = np.random.RandomState(20260813)
a = rng.normal(0.0, SIGMA, n)
N = np.zeros(n)
N[0] = rng.normal(0.0, SIGMA / np.sqrt(1.0 - PHI ** 2))   # start stationary
for t in range(1, n):
    N[t] = PHI * N[t - 1] + a[t]

y = 100.0 + N.copy()
y[t_step:] += OMEGA_STEP                        # step: permanent level shift
y[t_pulse] += OMEGA_PULSE                       # impulse: one period only

ts = fue.TimeSeries(list(y), freq=12, start=(2000, 1), name="SIM")

# `at` is a 0-BASED observation index — the same index used to build the series
# above.  Getting this off by one is not a small error: it moves the regressor
# one period away from the event, and the estimate absorbs whatever it finds
# there.  Written the first time with a 1-based index, this example recovered
# omega = +1.86 for a true −8.
#
# The dates are part of the SPECIFICATION: fue does not search for them.  An
# intervention is extra-sample information about something that happened, not a
# feature to be mined.
m = fue.Model(ts,
              d=0,
              ar=[[0.5]], ar_free=[[True]],
              interventions=[
                  fue.Intervention("step",    at=t_step,  omega=[1.0]),
                  fue.Intervention("impulse", at=t_pulse, omega=[1.0]),
              ],
              mu=float(np.mean(y)), estimate_mu=True)
m.fit()
r = m._result

print("ARMAX: AR(1) noise + a step and an impulse, both at known dates\n")
print(f"{'parameter':>12} {'true':>9} {'estimated':>11} {'s.e.':>9}")
print("-" * 45)
omegas = [iv.omega[0] for iv in m.interventions]
print(f"{'omega step':>12} {OMEGA_STEP:9.3f} {omegas[0]:11.3f} {r.std_errors[0]:9.3f}")
print(f"{'omega pulse':>12} {OMEGA_PULSE:9.3f} {omegas[1]:11.3f} {r.std_errors[1]:9.3f}")
print(f"{'phi':>12} {PHI:9.3f} {m.ar[0][0]:11.3f} {r.std_errors[2]:9.3f}")
print(f"{'sigma':>12} {SIGMA:9.3f} {np.sqrt(r.sigma2):11.3f}")

print(f"\n  logL = {m.loglik:.4f}, stopped by {r.termination}")

# A note on reading this table, learned by writing it wrong.  With a SHORT
# sample the estimates can sit several standard errors from the parameters and
# the engine is not to blame: a persistent AR(1) wanders, so the realisation
# genuinely has a different level.  At n=180 this example recovered omega=9.66
# for a true 12 — and ordinary least squares on the same design gave 9.55, i.e.
# the data said 9.6, not 12.  The check that means something is fue against
# another estimator on the SAME realisation, not against the parameters that
# generated it.

# The long-run gain g is what the simplification test of docs/FORMAL_TESTS.md §6
# is about.  With a single omega and no delta, g = omega: the step moves the
# level permanently, the impulse does not move it at all.
print(f"\n  long-run gain of the step   g = {omegas[0]:.3f}  (permanent)")
print(f"  long-run gain of the impulse g = 0.000  (by construction: it is one period)")

# What happens if the intervention is left out — the reason it is in the model.
m0 = fue.Model(ts, d=0, ar=[[0.5]], ar_free=[[True]],
               mu=float(np.mean(y)), estimate_mu=True)
m0.fit()
print(f"\n  without the interventions: phi = {m0.ar[0][0]:.3f} (true {PHI}), "
      f"sigma = {np.sqrt(m0._result.sigma2):.3f} (true {SIGMA})")
print("  the level shift is absorbed by the AR, which drifts towards a unit root")
