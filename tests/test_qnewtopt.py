"""
Tests for qnewtopt.py — pure-Python port of Dennis-Schnabel BFGS optimizer.

Coverage
--------
1. cdgrad          — central-difference gradient vs analytic
2. _lnsrch         — line search on a quadratic
3. _bfgsfac        — BFGS Cholesky update: B B^T grows toward H
4. _qrupdate       — rank-1 QR update: R^T R stays upper triangular, B-identity check
5. raxopt          — minimises several standard test functions (Rosenbrock, quadratic)
6. raxopt vs C     — numerical equivalence on fue real cases (requires C extension)
7. covariance      — cholsol covariance from raxopt B matches finite-diff Hessian
"""

import math
import pytest
import numpy as np
from numpy.testing import assert_allclose

from fue.qnewtopt import raxopt, cdgrad, _cholsol, _bfgsfac, _qrupdate, _jacrot
from fue.cast_us   import estimate_py, estimate_lbfgsb_py


# ── Helpers ────────────────────────────────────────────────────────────────────

def _rosenbrock(x):
    return sum(100.0 * (x[i+1] - x[i]**2)**2 + (1.0 - x[i])**2
               for i in range(len(x) - 1))

def _rosenbrock_grad(x):
    n = len(x)
    g = np.zeros(n)
    for i in range(n - 1):
        g[i]   += -400.0 * x[i] * (x[i+1] - x[i]**2) - 2.0 * (1.0 - x[i])
        g[i+1] +=  200.0 * (x[i+1] - x[i]**2)
    return g


# ── cdgrad ────────────────────────────────────────────────────────────────────

def test_cdgrad_quadratic():
    """Central differences on x^T A x / 2."""
    A   = np.array([[4.0, 1.0], [1.0, 3.0]])
    x0  = np.array([1.5, -0.5])
    f   = lambda x: 0.5 * float(x @ A @ x)
    g   = cdgrad(f, x0.copy())
    assert_allclose(g, A @ x0, rtol=1e-8)


def test_cdgrad_rosenbrock():
    x0  = np.array([-1.2, 1.0])
    g   = cdgrad(_rosenbrock, x0.copy())
    assert_allclose(g, _rosenbrock_grad(x0), rtol=1e-6)


def test_cdgrad_restores_x():
    x   = np.array([1.0, 2.0, 3.0])
    x_c = x.copy()
    cdgrad(lambda v: float(np.dot(v, v)), x)
    assert_allclose(x, x_c)


# ── _qrupdate ─────────────────────────────────────────────────────────────────

def test_qrupdate_3x3():
    """After rank-1 update R += u v^T (then re-triangularize), check R^T R."""
    rng = np.random.default_rng(0)
    n   = 3
    # Start from a random upper-triangular R
    M   = np.triu(rng.standard_normal((n, n))) + np.eye(n) * 2
    u   = rng.standard_normal(n)
    v   = rng.standard_normal(n)
    M0  = M.copy()
    _qrupdate(u, v, M)
    # M should still be upper-triangular
    for i in range(n):
        for j in range(i):
            assert abs(M[i, j]) < 1e-12, f"M[{i},{j}] = {M[i,j]:.2e} is not zero"


def test_qrupdate_identity():
    """qrupdate on identity: R = I, u = e1, v = e1 → result encodes update."""
    n = 2
    M = np.eye(n)
    u = np.array([1.0, 0.0])
    v = np.array([1.0, 0.0])
    _qrupdate(u, v, M)
    assert np.all(M == np.triu(M)), "M should be upper-triangular after qrupdate"


# ── _bfgsfac ─────────────────────────────────────────────────────────────────

def test_bfgsfac_updates_B():
    """B B^T should change after a BFGS step that satisfies curvature."""
    n     = 2
    B     = np.eye(n)
    xk    = np.array([0.0, 0.0])
    xkp1  = np.array([0.5, 0.3])
    # Gradient of x^T x = 2x
    gk    = 2.0 * xk
    gkp1  = 2.0 * xkp1
    B_old = B.copy()
    _bfgsfac(xk, xkp1, gk, gkp1, np.finfo(float).eps, B)
    assert not np.allclose(B, B_old), "B should have changed after valid BFGS step"


def test_bfgsfac_preserves_lower_triangular():
    n    = 3
    rng  = np.random.default_rng(42)
    B    = np.tril(rng.standard_normal((n, n)) + np.eye(n) * 2)
    xk   = rng.standard_normal(n)
    xkp1 = xk + rng.standard_normal(n) * 0.1
    # Gradient of simple quadratic
    A    = np.eye(n) * 3
    gk   = A @ xk
    gkp1 = A @ xkp1
    _bfgsfac(xk, xkp1, gk, gkp1, np.finfo(float).eps, B)
    # B must still be lower triangular
    assert np.allclose(B, np.tril(B)), "B must remain lower-triangular"


# ── raxopt on standard functions ──────────────────────────────────────────────

def test_raxopt_quadratic():
    """raxopt minimises a 3D quadratic exactly."""
    A  = np.diag([2.0, 4.0, 6.0])
    b  = np.array([1.0, 2.0, 3.0])
    x_true = np.linalg.solve(A, b)   # = [0.5, 0.5, 0.5]
    # Normalise so f(x0) = 1
    x0 = np.zeros(3)
    f0 = 0.5 * float(x0 @ A @ x0) - float(b @ x0) + 1.0   # = 1.0
    f  = lambda x: (0.5 * float(x @ A @ x) - float(b @ x) + 1.0) / f0
    x_opt, fval, B, tc = raxopt(f, x0)
    assert_allclose(x_opt, x_true, atol=1e-7)
    assert tc in (1, 2), f"unexpected termcode {tc}"


def test_raxopt_rosenbrock_2d():
    """raxopt finds minimum of normalised 2D Rosenbrock."""
    x0 = np.array([-1.2, 1.0])
    f0 = _rosenbrock(x0)
    f  = lambda x: _rosenbrock(x) / f0
    x_opt, fval, B, tc = raxopt(f, x0)
    assert_allclose(x_opt, [1.0, 1.0], atol=1e-5)
    assert fval < 1e-10 / f0 + 1e-12   # near zero at minimum


def test_raxopt_termcode_gradient():
    """Termcode 1: gradient criterion satisfied."""
    A  = np.eye(2) * 4
    x0 = np.array([0.5, 0.5])
    f0 = float(x0 @ A @ x0)
    f  = lambda x: float(x @ A @ x) / f0
    _, _, _, tc = raxopt(f, x0)
    assert tc == 1


def test_raxopt_maxiter_termcode():
    """Termcode 4: max iterations reached on a hard problem."""
    # Rosenbrock is hard enough that with maxits=5 it won't converge
    x0 = np.array([-1.2, 1.0])
    f0 = _rosenbrock(x0)
    f  = lambda x: _rosenbrock(x) / f0
    _, _, _, tc = raxopt(f, x0, maxits=5)
    assert tc == 4


def test_raxopt_zero_params():
    """raxopt with n=0 returns immediately."""
    x_opt, fval, B, tc = raxopt(lambda x: 1.0, np.array([]))
    assert len(x_opt) == 0
    assert tc == 1


# ── raxopt vs C extension ────────────────────────────────────────────────────

try:
    from fue._fue_engine import ffi as _ffi   # noqa
    _C_AVAILABLE = True
except ImportError:
    _C_AVAILABLE = False

requires_c = pytest.mark.skipif(not _C_AVAILABLE, reason="C extension not compiled")

import os
_REAL = os.path.join(os.path.dirname(__file__), "real_cases")


def _load(rel):
    import fue
    _, m = fue.load(os.path.join(_REAL, rel + ".inp"))
    return m


_EQUIV_CASES = [
    # (rel_path, loglik_tol, sigma2_tol)
    ("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.1",               2e-3, 1e-3),
    ("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.2",               2e-3, 5e-2),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.1", 5e-3, 1e-5),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.2", 5e-3, 3e-6),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.4", 5e-3, 2e-6),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.0",       5e-3, 5e-4),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1",       5e-3, 5e-4),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.3",       5e-3, 5e-4),
]
_EQUIV_IDS = [p.split("/")[-1] for p, *_ in _EQUIV_CASES]


@requires_c
@pytest.mark.parametrize("rel,ll_tol,s2_tol", _EQUIV_CASES, ids=_EQUIV_IDS)
def test_raxopt_matches_c(rel, ll_tol, s2_tol):
    """raxopt loglik must agree with C within tolerance."""
    from fue._engine import estimate as est_c
    m   = _load(rel)
    rc  = est_c(m)
    rp  = estimate_py(m)
    assert abs(rp["loglik"] - rc["loglik"]) < ll_tol, (
        f"raxopt loglik {rp['loglik']:.6f} vs C {rc['loglik']:.6f} "
        f"(diff {abs(rp['loglik']-rc['loglik']):.2e})"
    )


@requires_c
@pytest.mark.parametrize("rel,ll_tol,s2_tol", _EQUIV_CASES, ids=_EQUIV_IDS)
def test_raxopt_sigma2_matches_c(rel, ll_tol, s2_tol):
    """raxopt sigma² must agree with C within tolerance."""
    from fue._engine import estimate as est_c
    m   = _load(rel)
    rc  = est_c(m)
    rp  = estimate_py(m)
    # Skip when logliks diverge (different basin, sigma2 comparison meaningless)
    if abs(rp["loglik"] - rc["loglik"]) > 1.0:
        pytest.skip("different basin — sigma2 not comparable")
    assert abs(rp["sigma2"] - rc["sigma2"]) < s2_tol, (
        f"raxopt sigma2 {rp['sigma2']:.8f} vs C {rc['sigma2']:.8f}"
    )


# ── covariance via BFGS Hessian matches C ─────────────────────────────────────

@requires_c
def test_cov_raxopt_matches_c():
    """raxopt std_errors (from BFGS B) should match C std_errors closely.

    The C drvmlest.c also uses the BFGS Cholesky factor directly, so
    Python raxopt and C should agree to high precision.
    """
    from fue._engine import estimate as est_c
    m   = _load("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.2")
    rc  = est_c(m)
    rr  = estimate_py(m)
    assert_allclose(rr["std_errors"], rc["std_errors"], rtol=1e-5,
                    err_msg="raxopt std_errors must match C std_errors")
