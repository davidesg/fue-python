"""
benchmarks/bench_estimator.py — reproducible performance benchmark.

Measures wall-clock time and per-call costs for the three canonical test cases
at each architectural configuration available in the current build.

Usage:
    python benchmarks/bench_estimator.py

Output: tab-separated table to stdout; suitable for pasting into PERFORMANCE.md.
"""

import time
import numpy as np
from scipy.optimize import minimize

from fue import TimeSeries, Model
from fue.intervention import Intervention
from fue.cast_us import build_est_spec, cast_us_py, _build_initial_x
from fue.elfvarma import flikam_scalar, elf_scalar

# ── Test cases ────────────────────────────────────────────────────────────────

_SFNY30 = np.array([
    3.91505848, 2.02125792, 0.81208771, 0.60807414, 1.21576447,
    1.43763055, 1.78032601, 0.82841058, 0.65433228, 0.74324607,
    0.93394905, 0.60094494, 0.80840161, 0.90899270, 0.40822203,
    0.41975993, 0.50368768, 0.57248427, 0.72970370, 0.90175445,
    0.61763439, 0.63607641, 0.67670827, 0.81812744, 0.78095914,
    0.82024104, 0.86103433, 0.84442843, 0.74566075, 0.63347579,
])
_SFNY62 = np.array([
    3.91505848, 2.02125792, 0.81208771, 0.60807414, 1.21576447,
    1.43763055, 1.78032601, 0.82841058, 0.65433228, 0.74324607,
    0.93394905, 0.60094494, 0.80840161, 0.90899270, 0.40822203,
    0.41975993, 0.50368768, 0.57248427, 0.72970370, 0.90175445,
    0.61763439, 0.63607641, 0.67670827, 0.81812744, 0.78095914,
    0.82024104, 0.86103433, 0.84442843, 0.74566075, 0.63347579,
    0.72637557, 0.81351610, 0.79142754, 0.80305873, 0.83867533,
    0.98678814, 0.80485863, 0.81651553, 0.75960093, 0.84070968,
    0.89480882, 0.89407591, 0.84323646, 0.77215182, 0.82509544,
    0.87384443, 0.81360106, 0.78497496, 0.71323360, 0.70688522,
    0.81090348, 0.94831097, 0.72598922, 0.80337325, 0.84011493,
    0.89247202, 0.89328246, 0.90942424, 0.82871189, 0.88647340,
    0.82251497, 0.94737336,
])
_RIPC1 = np.array([
    0.413459, 0.416226, 0.418544, 0.422442, 0.424508, 0.425892,
    0.425137, 0.425577, 0.427322, 0.429367, 0.432350, 0.434795,
    0.443644, 0.443617, 0.443454, 0.448741, 0.450270, 0.448844,
    0.448505, 0.447079, 0.449154, 0.449650, 0.452374, 0.452679,
    0.452326, 0.452981, 0.453226, 0.454730, 0.449935, 0.447131,
    0.445082, 0.444966, 0.445047, 0.443958, 0.445585, 0.446936,
    0.447090, 0.445735, 0.443441, 0.444167, 0.445403, 0.445490,
    0.442748, 0.439849, 0.437653, 0.438303, 0.442595, 0.445719,
    0.444463, 0.446707, 0.447143, 0.443674, 0.440874, 0.438993,
    0.437829, 0.437908, 0.442587, 0.446552, 0.447956, 0.447148,
    0.447111, 0.445032, 0.441439, 0.438548, 0.436016, 0.436859,
    0.438797, 0.439924, 0.441830, 0.441481, 0.441057, 0.443876,
])


def _make_ar1():
    return Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ar=[[0.5]])


def _make_sfny2():
    return Model(
        TimeSeries(_SFNY62, freq=1, start=(1852, 1)),
        interventions=[Intervention(
            "step", at=1, omega=[0.08], omega_free=[True],
            delta=[0.6], delta_free=[True])],
        ar=[[0.8], [-0.1, -0.1]], boxlam=0.0, mu=0.0, estimate_mu=True,
    )


def _make_ripc1():
    return Model(
        TimeSeries(_RIPC1, freq=12, start=(2002, 1)),
        interventions=[
            Intervention("cos",   harmonic=1.0, omega=[0.0], omega_free=[True]),
            Intervention("sin",   harmonic=1.0, omega=[0.0], omega_free=[True]),
            Intervention("cos",   harmonic=2.0, omega=[0.0], omega_free=[True]),
            Intervention("sin",   harmonic=2.0, omega=[0.0], omega_free=[True]),
            Intervention("cos",   harmonic=3.0, omega=[0.0], omega_free=[True]),
            Intervention("sin",   harmonic=3.0, omega=[0.0], omega_free=[True]),
            Intervention("cos",   harmonic=4.0, omega=[0.0], omega_free=[True]),
            Intervention("sin",   harmonic=4.0, omega=[0.0], omega_free=[True]),
            Intervention("cos",   harmonic=5.0, omega=[0.0], omega_free=[True]),
            Intervention("sin",   harmonic=5.0, omega=[0.0], omega_free=[True]),
            Intervention("alter", omega=[0.0], omega_free=[True]),
            Intervention("step", at=1, omega=[0.014], omega_free=[True],
                         delta=[0.6], delta_free=[True]),
        ],
        ar=[[0.0]], ar_free=[[False]],
        boxlam=0.0, refactor=100.0, mu=0.0, estimate_mu=True,
    )


CASES = [
    ("AR(1)",   _make_ar1,   30,  1),
    ("SFNY.2",  _make_sfny2, 62,  6),
    ("RIPC.1",  _make_ripc1, 72, 14),
]


# ── Timing helpers ────────────────────────────────────────────────────────────

def _wall(fn, model_fn, reps=3):
    times = []
    for _ in range(reps):
        m = model_fn()
        t0 = time.perf_counter()
        fn(m)
        times.append(time.perf_counter() - t0)
    return min(times) * 1000


def _per_call(model_fn, reps=20):
    """Return (cast_us_py_ms, flikam_scalar_ms) per single evaluation."""
    m = model_fn()
    spec = build_est_spec(m)
    x0 = _build_initial_x(m)
    n_eff = m.series.nobs - spec.ornsop
    p, q, phi, theta, mu, w, _ = cast_us_py(x0, spec)

    for _ in range(3):
        cast_us_py(x0, spec)
    t0 = time.perf_counter()
    for _ in range(reps):
        cast_us_py(x0, spec)
    t_cast = (time.perf_counter() - t0) / reps * 1000

    for _ in range(3):
        flikam_scalar(n_eff, p, q, phi, theta, mu, w, xitol=1e-3)
    t0 = time.perf_counter()
    for _ in range(reps):
        flikam_scalar(n_eff, p, q, phi, theta, mu, w, xitol=1e-3)
    t_flik = (time.perf_counter() - t0) / reps * 1000

    return t_cast, t_flik


def _nfev(model_fn):
    """Count L-BFGS-B function evaluations for this model."""
    m = model_fn()
    spec = build_est_spec(m)
    x0 = _build_initial_x(m)
    n_eff = m.series.nobs - spec.ornsop
    p0, q0, phi0, theta0, mu0, w0, _ = cast_us_py(x0, spec)
    s0, f0, _, _, _ = flikam_scalar(n_eff, p0, q0, phi0, theta0, mu0, w0, xitol=1e-3)
    cnt = [0]

    def obj(x):
        cnt[0] += 1
        p, q, phi, theta, mu, w, fault = cast_us_py(x, spec)
        if fault or len(w) == 0:
            return 1.0
        s, f, _, _, iflt = flikam_scalar(n_eff, p, q, phi, theta, mu, w, xitol=1e-3)
        if iflt or s <= 0 or f <= 0:
            return 1.0
        return (s / s0) * (f / f0)

    minimize(obj, x0, method="L-BFGS-B",
             options={"maxiter": 500, "ftol": 1e-14, "gtol": 1e-7})
    return cnt[0]


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    try:
        from fue._engine import estimate as estimate_c
        has_c = True
    except Exception:
        has_c = False

    from fue.cast_us import estimate_py

    print("=" * 72)
    print("FUE estimator benchmark")
    print("=" * 72)

    print("\n--- Per-call cost (ms, best of 20) ---")
    print(f"{'Case':10s}  {'nobs':>5}  {'npar':>5}  "
          f"{'cast_us_py':>12}  {'flikam_py':>10}  {'total':>8}")
    for label, fn, nobs, npar in CASES:
        tc, tf = _per_call(fn)
        print(f"{label:10s}  {nobs:5d}  {npar:5d}  "
              f"{tc:12.3f}  {tf:10.3f}  {tc+tf:8.3f}")

    print("\n--- Optimizer evaluations (L-BFGS-B) ---")
    print(f"{'Case':10s}  {'npar':>5}  {'nfev':>6}  {'grad_steps':>10}")
    nfevs = {}
    for label, fn, nobs, npar in CASES:
        nf = _nfev(fn)
        nfevs[label] = nf
        print(f"{label:10s}  {npar:5d}  {nf:6d}  {nf//(npar+1):10d}")

    print("\n--- Wall-clock time (ms, best of 3) ---")
    reps = {"AR(1)": 3, "SFNY.2": 3, "RIPC.1": 1}
    print(f"{'Case':10s}  {'C (ms)':>10}  {'Py (ms)':>10}  {'factor':>8}")
    for label, fn, nobs, npar in CASES:
        t_py = _wall(estimate_py, fn, reps[label])
        if has_c:
            t_c = _wall(estimate_c, fn, reps[label])
            print(f"{label:10s}  {t_c:10.1f}  {t_py:10.1f}  {t_py/t_c:8.0f}x")
        else:
            print(f"{label:10s}  {'n/a':>10}  {t_py:10.1f}  {'n/a':>8}")

    # Projected hybrid (Python opt + C inner loop)
    if has_c:
        print("\n--- Projected: hybrid Py-optimizer + C inner loop ---")
        print("  (C inner loop = cast_us + flikam compiled; ~5 µs per eval)")
        C_EVAL_US = 5.0
        print(f"{'Case':10s}  {'nfev':>6}  {'proj_ms':>10}  {'vs_C':>8}")
        for label, fn, nobs, npar in CASES:
            nf = nfevs[label]
            t_c = _wall(estimate_c, fn, reps[label])
            proj = nf * C_EVAL_US / 1000
            print(f"{label:10s}  {nf:6d}  {proj:10.1f}  {proj/t_c:8.1f}x")


if __name__ == "__main__":
    main()
