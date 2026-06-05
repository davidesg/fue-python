"""
Performance tests: fue C engine vs fue Python engine.

Run all:         pytest tests/test_performance.py -v -s
C-only:          pytest tests/test_performance.py -v -s -k "c_engine"
Python-only:     pytest tests/test_performance.py -v -s -k "py_engine"
Summary table:   pytest tests/test_performance.py -v -s -k summary

These tests verify:
  1. Numerical equivalence — C and Python results match reference within tolerance.
  2. Timing ratio — Python is slower than C (structural invariant).
  3. Absolute ceilings — no runaway regressions (loose, environment-independent).

Timing assertions are intentionally loose because wall-clock time is
environment-dependent.  The printed summary table is the real deliverable.

Reference values: fue-1.13.1, Linux x86-64, CPython 3.12, gcc -O2, 2026-06-05.
Reproduce:        pytest tests/test_performance.py::TestCvsPython::test_summary -v -s
"""

import os
import time
import pytest
import numpy as np

from fue import TimeSeries, Model
from fue.intervention import Intervention
from fue.cast_us import estimate_py

# ── C availability ────────────────────────────────────────────────────────────

try:
    from fue._fue_engine import ffi as _ffi   # noqa: F401
    _C_AVAILABLE = True
except ImportError:
    _C_AVAILABLE = False

requires_c = pytest.mark.skipif(not _C_AVAILABLE,
                                reason="C extension not compiled")

# ── Shared data ───────────────────────────────────────────────────────────────

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

_REAL = os.path.join(os.path.dirname(__file__), "real_cases")


# ── Model factories ───────────────────────────────────────────────────────────

def _ar1():
    return Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ar=[[0.5]])

def _ima11():
    return Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)), ma=[[0.3]], d=1)

def _arma11():
    return Model(TimeSeries(_SFNY30, freq=1, start=(1852, 1)),
                 ar=[[0.5]], ma=[[0.3]])

def _sfny2():
    return Model(
        TimeSeries(_SFNY62, freq=1, start=(1852, 1)),
        interventions=[Intervention(
            "step", at=1, omega=[0.08], omega_free=[True],
            delta=[0.6], delta_free=[True])],
        ar=[[0.8], [-0.1, -0.1]], boxlam=0.0, mu=0.0, estimate_mu=True,
    )

def _from_inp(rel):
    import fue
    _, m = fue.load(os.path.join(_REAL, rel + ".inp"))
    return m

def _gdp_r1():
    return _from_inp("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.1")

def _gdp_r2():
    return _from_inp("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.2")

def _ipct_r3():
    return _from_inp("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.3")

def _ipct_r5():
    return _from_inp("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/old/R.5")

def _ripc0():
    return _from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.0")

def _ripc1():
    return _from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1")

def _ripc3():
    return _from_inp("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.3")


# ── Case catalogue ────────────────────────────────────────────────────────────
#
# (id, factory, freq, nobs, npar, ref_loglik, ref_sigma2, tol_c_ll, tol_c_s2,
#                                                          tol_py_ll, tol_py_s2)
#
# Tolerances match test_estimation.py conventions.  RIPC cases with many params
# use looser bounds because the 14-17D optimizer may stop at a slightly
# different point; loglik differences are still tiny (< 5e-3).

_CASES = [
    # id          factory   fr   n   p    ref_loglik          ref_sigma2       tc_ll  tc_s2   tp_ll  tp_s2
    ("AR(1)",      _ar1,     1,  30,  1, -23.1683049163,  0.2482607765,     1e-4, 1e-6, 1e-3, 5e-4),
    ("IMA(1,1)",   _ima11,   1,  30,  1, -18.3455244692,  0.2060862799,     1e-3, 1e-4, 2e-3, 5e-4),
    ("ARMA(1,1)",  _arma11,  1,  30,  2, -20.5630815929,  0.2054046560,     1e-3, 1e-5, 5e-3, 5e-4),
    ("SFNY.2",     _sfny2,   1,  62,  6,  13.9573576937,  0.0370593261,     1e-4, 1e-6, 1e-3, 1e-5),
    ("GDP/R.1",    _gdp_r1,  4,  68,  1,-160.6143098964,  7.0750884543,     1e-3, 1e-4, 2e-3, 5e-3),
    ("GDP/R.2",    _gdp_r2,  4,  68,  3,-141.8868690241,  4.0452844123,     1e-3, 1e-4, 5e-3, 5e-2),
    ("IPC-T/R.3",  _ipct_r3, 4,  68,  4, 228.2017654123,  0.0000658767,     1e-3, 1e-7, 5e-3, 1e-6),
    ("IPC-T/R.5",  _ipct_r5, 4,  68,  7, 182.8429916793,  0.0002658274,     1e-3, 1e-7, 5e-3, 5e-6),
    ("RIPC.0",     _ripc0,  12,  78, 13, -58.6021700483,  0.2486850071,     1e-3, 1e-4, 5e-3, 5e-4),
    ("RIPC.1",     _ripc1,  12,  72, 14,-100.9274828448,  0.9662469111,     1e-3, 1e-4, 5e-3, 5e-4),
    ("RIPC.3",     _ripc3,  12,  78, 17, -39.1764600458,  0.1551235917,     1e-3, 1e-4, 5e-3, 5e-4),
]

# Parametrize IDs (used by pytest -k)
_IDS = [c[0] for c in _CASES]


# ── Timing helpers ────────────────────────────────────────────────────────────

def _time_engine(fn, factory, reps):
    times = []
    for _ in range(reps):
        m = factory()
        t0 = time.perf_counter()
        fn(m)
        times.append(time.perf_counter() - t0)
    return min(times) * 1000   # ms


def _reps(npar):
    """Fewer repetitions for expensive Python runs."""
    return 1 if npar >= 10 else 3


# ══════════════════════════════════════════════════════════════════════════════
# C engine — correctness
# ══════════════════════════════════════════════════════════════════════════════

@requires_c
@pytest.mark.parametrize(
    "case_id,factory,freq,nobs,npar,ref_ll,ref_s2,tol_c_ll,tol_c_s2,tol_py_ll,tol_py_s2",
    _CASES, ids=_IDS,
)
def test_c_loglik(case_id, factory, freq, nobs, npar,
                  ref_ll, ref_s2, tol_c_ll, tol_c_s2, tol_py_ll, tol_py_s2):
    from fue._engine import estimate
    r = estimate(factory())
    assert abs(r["loglik"] - ref_ll) < tol_c_ll, (
        f"{case_id}: loglik {r['loglik']:.8f} vs {ref_ll:.8f} "
        f"(diff {abs(r['loglik']-ref_ll):.2e} > tol {tol_c_ll:.0e})"
    )


@requires_c
@pytest.mark.parametrize(
    "case_id,factory,freq,nobs,npar,ref_ll,ref_s2,tol_c_ll,tol_c_s2,tol_py_ll,tol_py_s2",
    _CASES, ids=_IDS,
)
def test_c_sigma2(case_id, factory, freq, nobs, npar,
                  ref_ll, ref_s2, tol_c_ll, tol_c_s2, tol_py_ll, tol_py_s2):
    from fue._engine import estimate
    r = estimate(factory())
    assert abs(r["sigma2"] - ref_s2) < tol_c_s2, (
        f"{case_id}: sigma2 {r['sigma2']:.8f} vs {ref_s2:.8f} "
        f"(diff {abs(r['sigma2']-ref_s2):.2e} > tol {tol_c_s2:.0e})"
    )


# ══════════════════════════════════════════════════════════════════════════════
# Python engine — correctness
# ══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize(
    "case_id,factory,freq,nobs,npar,ref_ll,ref_s2,tol_c_ll,tol_c_s2,tol_py_ll,tol_py_s2",
    _CASES, ids=_IDS,
)
def test_py_loglik(case_id, factory, freq, nobs, npar,
                   ref_ll, ref_s2, tol_c_ll, tol_c_s2, tol_py_ll, tol_py_s2):
    r = estimate_py(factory())
    assert abs(r["loglik"] - ref_ll) < tol_py_ll, (
        f"{case_id}: loglik {r['loglik']:.8f} vs {ref_ll:.8f} "
        f"(diff {abs(r['loglik']-ref_ll):.2e} > tol {tol_py_ll:.0e})"
    )


@pytest.mark.parametrize(
    "case_id,factory,freq,nobs,npar,ref_ll,ref_s2,tol_c_ll,tol_c_s2,tol_py_ll,tol_py_s2",
    _CASES, ids=_IDS,
)
def test_py_sigma2(case_id, factory, freq, nobs, npar,
                   ref_ll, ref_s2, tol_c_ll, tol_c_s2, tol_py_ll, tol_py_s2):
    r = estimate_py(factory())
    assert abs(r["sigma2"] - ref_s2) < tol_py_s2, (
        f"{case_id}: sigma2 {r['sigma2']:.8f} vs {ref_s2:.8f} "
        f"(diff {abs(r['sigma2']-ref_s2):.2e} > tol {tol_py_s2:.0e})"
    )


# ══════════════════════════════════════════════════════════════════════════════
# C engine — timing ceilings
# ══════════════════════════════════════════════════════════════════════════════

@requires_c
@pytest.mark.parametrize(
    "case_id,factory,freq,nobs,npar,ref_ll,ref_s2,tol_c_ll,tol_c_s2,tol_py_ll,tol_py_s2",
    _CASES, ids=_IDS,
)
def test_c_timing(case_id, factory, freq, nobs, npar,
                  ref_ll, ref_s2, tol_c_ll, tol_c_s2, tol_py_ll, tol_py_s2,
                  capsys):
    ceiling = 500 if npar < 14 else 2000
    t = _time_engine(lambda m: m, lambda: factory(), 1)
    from fue._engine import estimate
    t = _time_engine(estimate, factory, 3)
    with capsys.disabled():
        print(f"\n  C  {case_id:<12} {t:7.1f} ms")
    assert t < ceiling, f"C {case_id}: {t:.1f} ms > ceiling {ceiling} ms"


# ══════════════════════════════════════════════════════════════════════════════
# Python engine — timing ceilings
# ══════════════════════════════════════════════════════════════════════════════

# Ceiling: 500× the reference C time, rounded up generously.
# These are regression guards, not accuracy checks.
_PY_CEILINGS = {
    "AR(1)":     500,
    "IMA(1,1)":  500,
    "ARMA(1,1)": 500,
    "SFNY.2":    5_000,
    "GDP/R.1":   2_000,
    "GDP/R.2":   5_000,
    "IPC-T/R.3": 10_000,
    "IPC-T/R.5": 30_000,
    "RIPC.0":    60_000,
    "RIPC.1":    60_000,
    "RIPC.3":    60_000,
}


@pytest.mark.parametrize(
    "case_id,factory,freq,nobs,npar,ref_ll,ref_s2,tol_c_ll,tol_c_s2,tol_py_ll,tol_py_s2",
    _CASES, ids=_IDS,
)
def test_py_timing(case_id, factory, freq, nobs, npar,
                   ref_ll, ref_s2, tol_c_ll, tol_c_s2, tol_py_ll, tol_py_s2,
                   capsys):
    ceiling = _PY_CEILINGS[case_id]
    t = _time_engine(estimate_py, factory, _reps(npar))
    with capsys.disabled():
        print(f"\n  Py {case_id:<12} {t:8.1f} ms")
    assert t < ceiling, f"Py {case_id}: {t:.1f} ms > ceiling {ceiling} ms"


# ══════════════════════════════════════════════════════════════════════════════
# C vs Python — structural invariants + summary table
# ══════════════════════════════════════════════════════════════════════════════

@requires_c
@pytest.mark.parametrize(
    "case_id,factory,freq,nobs,npar,ref_ll,ref_s2,tol_c_ll,tol_c_s2,tol_py_ll,tol_py_s2",
    _CASES, ids=_IDS,
)
def test_py_slower_than_c(case_id, factory, freq, nobs, npar,
                           ref_ll, ref_s2, tol_c_ll, tol_c_s2,
                           tol_py_ll, tol_py_s2):
    from fue._engine import estimate as est_c
    t_c  = _time_engine(est_c,        factory, 3)
    t_py = _time_engine(estimate_py,  factory, _reps(npar))
    assert t_py > t_c, (
        f"{case_id}: Python ({t_py:.1f} ms) should be slower than C ({t_c:.1f} ms)"
    )


@requires_c
@pytest.mark.parametrize(
    "case_id,factory,freq,nobs,npar,ref_ll,ref_s2,tol_c_ll,tol_c_s2,tol_py_ll,tol_py_s2",
    _CASES, ids=_IDS,
)
def test_loglik_agree(case_id, factory, freq, nobs, npar,
                      ref_ll, ref_s2, tol_c_ll, tol_c_s2, tol_py_ll, tol_py_s2):
    from fue._engine import estimate as est_c
    r_c  = est_c(factory())
    r_py = estimate_py(factory())
    assert abs(r_c["loglik"] - r_py["loglik"]) < 1e-1, (
        f"{case_id}: C loglik {r_c['loglik']:.6f} vs Py {r_py['loglik']:.6f}"
    )


@requires_c
def test_summary(capsys):
    """Print C vs Python timing table for all cases (never fails on timing)."""
    from fue._engine import estimate as est_c

    header = (f"\n  {'Case':<14} {'fr':>3} {'n':>4} {'p':>3}"
              f"  {'C (ms)':>8}  {'Py (ms)':>9}  {'factor':>8}")
    sep    = "  " + "-" * 60
    lines  = [header, sep]

    for (case_id, factory, freq, nobs, npar,
         ref_ll, ref_s2, *_) in _CASES:
        t_c  = _time_engine(est_c,       factory, 3)
        t_py = _time_engine(estimate_py, factory, _reps(npar))
        lines.append(
            f"  {case_id:<14} {freq:3d} {nobs:4d} {npar:3d}"
            f"  {t_c:8.1f}  {t_py:9.1f}  {t_py/t_c:7.0f}x"
        )

    with capsys.disabled():
        print("\n".join(lines))
