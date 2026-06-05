"""
Performance tests: fue C engine vs fue Python engine.

Run all:          pytest tests/test_performance.py -v -s
C-only:           pytest tests/test_performance.py -v -s -k "c_engine"
Python-only:      pytest tests/test_performance.py -v -s -k "py_engine"

These tests verify:
  1. Numerical equivalence — Python results match C within tolerance.
  2. Timing ratio — Python is slower than C (structural sanity).
  3. Python absolute ceiling — no runaway regressions.

Timing assertions are intentionally loose (×1000 ceiling) because wall-clock
time is environment-dependent.  The printed table is the real deliverable.

Reference values (fue-1.13.1, Linux x86-64, gcc -O2, 2026-06-05):
  see PERFORMANCE.md — Historial de mediciones.
"""

import time
import math
import pytest
import numpy as np

from fue import TimeSeries, Model
from fue.intervention import Intervention
from fue.cast_us import estimate_py


# ── Fixtures — canonical test models ─────────────────────────────────────────

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


def _ar1():
    return Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ar=[[0.5]])


def _sfny2():
    return Model(
        TimeSeries(_SFNY62, freq=1, start=(1852, 1)),
        interventions=[Intervention(
            "step", at=1, omega=[0.08], omega_free=[True],
            delta=[0.6], delta_free=[True])],
        ar=[[0.8], [-0.1, -0.1]], boxlam=0.0, mu=0.0, estimate_mu=True,
    )


def _ripc1():
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


# ── Timing helpers ────────────────────────────────────────────────────────────

def _time_c(model_fn, reps=3):
    from fue._engine import estimate as est_c
    times = []
    for _ in range(reps):
        m = model_fn()
        t0 = time.perf_counter()
        est_c(m)
        times.append(time.perf_counter() - t0)
    return min(times) * 1000   # ms


def _time_py(model_fn, reps=3):
    times = []
    for _ in range(reps):
        m = model_fn()
        t0 = time.perf_counter()
        estimate_py(m)
        times.append(time.perf_counter() - t0)
    return min(times) * 1000   # ms


# ── Reference values (fue-1.13.1) ────────────────────────────────────────────

_REF = {
    "ar1":   {"loglik": -23.1683049163, "sigma2": 0.2482607765},
    "sfny2": {"loglik":  13.9573576937, "sigma2": 0.0370593261},
    "ripc1": {"loglik": -100.9274828448, "sigma2": 0.9662469111},
}

# Tolerances match test_estimation.py (reference: fue-1.13.1).
# RIPC.1 has looser bounds because the 14-param optimizer converges to a
# slightly different point depending on initial values and stopping criteria.
_TOL_C = {
    "ar1":   {"loglik": 1e-4, "sigma2": 1e-6},
    "sfny2": {"loglik": 1e-4, "sigma2": 1e-6},
    "ripc1": {"loglik": 1e-3, "sigma2": 1e-4},
}
_TOL_PY = {
    "ar1":   {"loglik": 1e-3, "sigma2": 5e-4},
    "sfny2": {"loglik": 1e-3, "sigma2": 1e-5},
    "ripc1": {"loglik": 5e-3, "sigma2": 5e-4},
}

_HAS_C = pytest.importorskip if True else None
try:
    from fue._fue_engine import ffi as _ffi  # noqa: F401
    _C_AVAILABLE = True
except ImportError:
    _C_AVAILABLE = False

requires_c = pytest.mark.skipif(not _C_AVAILABLE,
                                reason="C extension not compiled")


# ══════════════════════════════════════════════════════════════════════════════
# C engine tests
# ══════════════════════════════════════════════════════════════════════════════

class TestCEngine:
    """Numerical correctness and timing of the C estimator (fue_api.c + GSL)."""

    @requires_c
    def test_ar1_loglik(self):
        from fue._engine import estimate
        r = estimate(_ar1())
        assert abs(r["loglik"] - _REF["ar1"]["loglik"]) < _TOL_C["ar1"]["loglik"]

    @requires_c
    def test_ar1_sigma2(self):
        from fue._engine import estimate
        r = estimate(_ar1())
        assert abs(r["sigma2"] - _REF["ar1"]["sigma2"]) < _TOL_C["ar1"]["sigma2"]

    @requires_c
    def test_sfny2_loglik(self):
        from fue._engine import estimate
        r = estimate(_sfny2())
        assert abs(r["loglik"] - _REF["sfny2"]["loglik"]) < _TOL_C["sfny2"]["loglik"]

    @requires_c
    def test_sfny2_sigma2(self):
        from fue._engine import estimate
        r = estimate(_sfny2())
        assert abs(r["sigma2"] - _REF["sfny2"]["sigma2"]) < _TOL_C["sfny2"]["sigma2"]

    @requires_c
    def test_ripc1_loglik(self):
        from fue._engine import estimate
        r = estimate(_ripc1())
        assert abs(r["loglik"] - _REF["ripc1"]["loglik"]) < _TOL_C["ripc1"]["loglik"]

    @requires_c
    def test_ripc1_sigma2(self):
        from fue._engine import estimate
        r = estimate(_ripc1())
        assert abs(r["sigma2"] - _REF["ripc1"]["sigma2"]) < _TOL_C["ripc1"]["sigma2"]

    @requires_c
    def test_timing_ar1(self, capsys):
        t = _time_c(_ar1)
        with capsys.disabled():
            print(f"\n  C  AR(1)   {t:7.1f} ms")
        assert t < 500, f"C AR(1) unexpectedly slow: {t:.1f} ms"

    @requires_c
    def test_timing_sfny2(self, capsys):
        t = _time_c(_sfny2)
        with capsys.disabled():
            print(f"\n  C  SFNY.2  {t:7.1f} ms")
        assert t < 500, f"C SFNY.2 unexpectedly slow: {t:.1f} ms"

    @requires_c
    def test_timing_ripc1(self, capsys):
        t = _time_c(_ripc1)
        with capsys.disabled():
            print(f"\n  C  RIPC.1  {t:7.1f} ms")
        assert t < 2000, f"C RIPC.1 unexpectedly slow: {t:.1f} ms"


# ══════════════════════════════════════════════════════════════════════════════
# Python engine tests
# ══════════════════════════════════════════════════════════════════════════════

class TestPyEngine:
    """Numerical correctness and timing of the pure-Python estimator."""

    def test_ar1_loglik(self):
        r = estimate_py(_ar1())
        assert abs(r["loglik"] - _REF["ar1"]["loglik"]) < _TOL_PY["ar1"]["loglik"]

    def test_ar1_sigma2(self):
        r = estimate_py(_ar1())
        assert abs(r["sigma2"] - _REF["ar1"]["sigma2"]) < _TOL_PY["ar1"]["sigma2"]

    def test_sfny2_loglik(self):
        r = estimate_py(_sfny2())
        assert abs(r["loglik"] - _REF["sfny2"]["loglik"]) < _TOL_PY["sfny2"]["loglik"]

    def test_sfny2_sigma2(self):
        r = estimate_py(_sfny2())
        assert abs(r["sigma2"] - _REF["sfny2"]["sigma2"]) < _TOL_PY["sfny2"]["sigma2"]

    def test_ripc1_loglik(self):
        r = estimate_py(_ripc1())
        assert abs(r["loglik"] - _REF["ripc1"]["loglik"]) < _TOL_PY["ripc1"]["loglik"]

    def test_ripc1_sigma2(self):
        r = estimate_py(_ripc1())
        assert abs(r["sigma2"] - _REF["ripc1"]["sigma2"]) < _TOL_PY["ripc1"]["sigma2"]

    def test_timing_ar1(self, capsys):
        t = _time_py(_ar1)
        with capsys.disabled():
            print(f"\n  Py AR(1)   {t:7.1f} ms")
        assert t < 500, f"Py AR(1) unexpectedly slow: {t:.1f} ms"

    def test_timing_sfny2(self, capsys):
        t = _time_py(_sfny2)
        with capsys.disabled():
            print(f"\n  Py SFNY.2  {t:7.1f} ms")
        assert t < 5000, f"Py SFNY.2 unexpectedly slow: {t:.1f} ms"

    def test_timing_ripc1(self, capsys):
        t = _time_py(_ripc1, reps=1)
        with capsys.disabled():
            print(f"\n  Py RIPC.1  {t:7.1f} ms")
        assert t < 30000, f"Py RIPC.1 unexpectedly slow: {t:.1f} ms"


# ══════════════════════════════════════════════════════════════════════════════
# Comparison: C vs Python
# ══════════════════════════════════════════════════════════════════════════════

@requires_c
class TestCvsPython:
    """Structural invariants that must hold between both engines."""

    def test_ar1_py_slower_than_c(self):
        t_c  = _time_c(_ar1)
        t_py = _time_py(_ar1)
        assert t_py > t_c, (
            f"Python ({t_py:.1f} ms) should be slower than C ({t_c:.1f} ms)"
        )

    def test_sfny2_py_slower_than_c(self):
        t_c  = _time_c(_sfny2)
        t_py = _time_py(_sfny2)
        assert t_py > t_c

    def test_ripc1_py_slower_than_c(self):
        t_c  = _time_c(_ripc1, reps=1)
        t_py = _time_py(_ripc1, reps=1)
        assert t_py > t_c

    def test_ar1_loglik_agree(self):
        from fue._engine import estimate as est_c
        r_c  = est_c(_ar1())
        r_py = estimate_py(_ar1())
        assert abs(r_c["loglik"] - r_py["loglik"]) < 1e-2

    def test_sfny2_loglik_agree(self):
        from fue._engine import estimate as est_c
        r_c  = est_c(_sfny2())
        r_py = estimate_py(_sfny2())
        assert abs(r_c["loglik"] - r_py["loglik"]) < 1e-2

    def test_ripc1_loglik_agree(self):
        from fue._engine import estimate as est_c
        r_c  = est_c(_ripc1())
        r_py = estimate_py(_ripc1())
        assert abs(r_c["loglik"] - r_py["loglik"]) < 1e-2

    def test_summary_table(self, capsys):
        """Print C vs Python timing table (always runs, never fails)."""
        rows = [
            ("AR(1)",  _ar1,   3,  3),
            ("SFNY.2", _sfny2, 3,  3),
            ("RIPC.1", _ripc1, 1,  1),
        ]
        header = f"\n  {'Case':<8}  {'C (ms)':>8}  {'Py (ms)':>9}  {'factor':>8}"
        lines  = [header, "  " + "-" * 40]
        for label, fn, rc, rp in rows:
            t_c  = _time_c(fn, rc)
            t_py = _time_py(fn, rp)
            lines.append(
                f"  {label:<8}  {t_c:8.1f}  {t_py:9.1f}  {t_py/t_c:7.0f}x"
            )
        with capsys.disabled():
            print("\n".join(lines))
