"""
Pure-Python port of qnewtopt.c: Factorized BFGS Quasi-Newton optimizer.

Source: Dennis & Schnabel (1983), Numerical Methods for Unconstrained
        Optimization and Nonlinear Equations, Algorithm A9.4.1.
Original C: Copyright (C) José Alberto Mauricio, 1995.
Python port: fue project.

Public API
----------
raxopt(func, x0, ...) -> (x, f, B, termcode)
    Minimize func starting from x0.  Returns optimal x, function value f,
    lower-triangular Cholesky factor B of the approximate Hessian at x, and
    a convergence code.

cdgrad(func, x, eta) -> g
    Central-difference gradient of func at x.

_cholsol(L, b) -> x
    Solve L @ L.T @ x = b  (L lower-triangular Cholesky factor).

Termination codes
-----------------
1  Gradient criterion satisfied (||scaled_grad||_inf <= gradtol)
2  Step-size criterion satisfied (||scaled_step||_inf <= steptol)
3  Line-search failed to find a lower point
4  Iteration limit reached
5  Five consecutive maximum-length steps taken
"""

import math
import numpy as np
from scipy.linalg import solve_triangular

_MACHEPS = np.finfo(float).eps
_GRTOL   = _MACHEPS ** (1.1 / 3.0)   # ≈ 1.82e-6 — matches fue_defaults() grtol
_SPTOL   = _MACHEPS ** (2.0 / 3.0)   # ≈ 3.67e-11 — matches fue_defaults() sptol
_MAXCMAX = 5


# ── Central-difference gradient ───────────────────────────────────────────────

def cdgrad(func, x, eta=None):
    """Central-difference gradient of func at x.

    x is modified in place during computation and restored on exit.
    Mirrors cdgrad() in qnewtopt.c exactly.
    """
    if eta is None:
        eta = _MACHEPS
    n     = len(x)
    g     = np.empty(n)
    third = eta ** (1.0 / 3.0)
    for i in range(n):
        xi    = x[i]
        stepi = third * (max(xi, 1.0) if xi >= 0.0 else min(xi, -1.0))
        x[i]  = xi + stepi
        stepi = x[i] - xi          # reduces FP errors (same trick as C)
        fpls  = func(x)
        x[i]  = xi - stepi
        fmns  = func(x)
        g[i]  = (fpls - fmns) / (2.0 * stepi)
        x[i]  = xi
    return g


# ── Convergence criteria ──────────────────────────────────────────────────────

def _umstop0(x, f, g, gradtol, maxits):
    """Initial stopping criterion (before first iteration)."""
    if maxits == 0:
        return 4
    fabsf = abs(f) + 1.0
    max1  = 0.0
    for i in range(len(x)):
        tmp = abs(g[i]) * (abs(x[i]) + 1.0) / fabsf
        if tmp > max1:
            max1 = tmp
    return 1 if max1 <= gradtol else 0


def _umstop(xk, xkp1, fkp1, gkp1, retcode,
            gradtol, steptol, k, maxits, maxcmax, maxtaken, consecmax):
    """Convergence check after each iteration.

    Returns (termcode, updated_consecmax).
    """
    fabsf = abs(fkp1) + 1.0
    n     = len(xkp1)
    max1  = 0.0
    for i in range(n):
        tmp = abs(gkp1[i]) * (abs(xkp1[i]) + 1.0) / fabsf
        if tmp > max1:
            max1 = tmp
    max2 = 0.0
    for i in range(n):
        tmp = abs(xkp1[i] - xk[i]) / (abs(xkp1[i]) + 1.0)
        if tmp > max2:
            max2 = tmp
    if max1 <= gradtol:
        return 1, consecmax
    if max2 <= steptol:
        return 2, consecmax
    if retcode == 1:
        return 3, consecmax
    if k >= maxits:
        return 4, consecmax
    if maxtaken:
        consecmax += 1
        return (5 if consecmax == maxcmax else 0), consecmax
    return 0, 0


# ── Line search ───────────────────────────────────────────────────────────────

def _lnsrch(func, xk, fk, gk, dk, maxstep, steptol):
    """Backtracking line search with cubic interpolation.

    Mirrors lnsrch() in qnewtopt.c (Dennis-Schnabel Algorithm A6.3.1).
    Returns (xkp1, fkp1, retcode, maxtaken, lambda).
    retcode: 0=success, 1=step too small, 2=not converged (never returned).
    """
    maxtaken = False
    retcode  = 2
    alpha    = 1.0e-4

    newtlen = math.sqrt(float(np.dot(dk, dk)))
    if newtlen > maxstep:
        dk      = dk * (maxstep / newtlen)
        newtlen = maxstep

    initslp = float(np.dot(gk, dk))

    rellen = 0.0
    for i in range(len(xk)):
        tmp = abs(dk[i]) / max(abs(xk[i]), 1.0)
        if tmp > rellen:
            rellen = tmp
    minlam = steptol / rellen

    lam    = 1.0
    prelam = 0.0
    pfkp1  = 0.0
    xkp1   = np.empty_like(xk)

    while retcode == 2:
        xkp1[:] = xk + lam * dk
        fkp1    = func(xkp1)

        if fkp1 <= fk + alpha * lam * initslp:
            retcode = 0
            if lam == 1.0 and newtlen > 0.99 * maxstep:
                maxtaken = True
        elif lam < minlam:
            retcode  = 1
            xkp1[:]  = xk
            fkp1     = fk
        else:
            if lam == 1.0:
                tlambda = -initslp / (2.0 * (fkp1 - fk - initslp))
            else:
                t1 = fkp1 - fk - lam   * initslp
                t2 = pfkp1 - fk - prelam * initslp
                t3 = 1.0 / (lam - prelam)
                a  = (t1 / (lam*lam) - t2 / (prelam*prelam)) * t3
                b  = (t2 * lam / (prelam*prelam) -
                      t1 * prelam / (lam*lam)) * t3
                if a == 0.0:
                    tlambda = -initslp / (2.0 * b)
                else:
                    disc = b*b - 3.0*a*initslp
                    if disc < 0.0:
                        raise RuntimeError("ROUNDOFF PROBLEM IN LINE SEARCH")
                    tlambda = (-b + math.sqrt(disc)) / (3.0 * a)
                if tlambda > 0.5 * lam:
                    tlambda = 0.5 * lam
            prelam = lam
            pfkp1  = fkp1
            lam    = tlambda if tlambda > 0.1 * lam else 0.1 * lam

    return xkp1, fkp1, retcode, maxtaken, lam


# ── Cholesky solve ────────────────────────────────────────────────────────────

def _cholsol(L, b):
    """Solve L @ L.T @ x = b  (L lower-triangular). Returns x."""
    y = solve_triangular(L, b, lower=True)
    return solve_triangular(L, y, lower=True, trans='T')


# ── Givens rotation ───────────────────────────────────────────────────────────

def _jacrot(i, a, b, M):
    """Apply Givens rotation to rows i and i+1 of M, columns i..n-1 (0-based).

    Mirrors jacrot(n, i, a, b, M) in qnewtopt.c.  C uses 1-based rows;
    here i is 0-based, mapping i_C → i_Py = i_C - 1.
    """
    if a == 0.0:
        c, s = 0.0, (1.0 if b >= 0.0 else -1.0)
    else:
        den  = math.sqrt(a*a + b*b)
        c, s = a / den, b / den
    n         = M.shape[1]
    y         = M[i,     i:n].copy()
    w         = M[i + 1, i:n].copy()
    M[i,     i:n] = c*y - s*w
    M[i + 1, i:n] = s*y + c*w


# ── Rank-1 QR update ─────────────────────────────────────────────────────────

def _qrupdate(u, v, M):
    """Rank-1 QR update of upper-triangular M (in place).

    Mirrors qrupdate(n, u, v, M) in qnewtopt.c.  All arrays are 0-based.
    After the call M is still upper-triangular and encodes the updated R factor.
    """
    n = len(u)
    u = u.copy()                    # C modifies u internally; don't clobber caller

    # Find 0-based last non-zero of u
    k = n - 1
    while k > 0 and u[k] == 0.0:
        k -= 1

    # Apply Givens rotations downward to concentrate u into u[0]
    # C 1-based: i from k-1 down to 1 → 0-based: i from k-1 down to 0
    for i in range(k - 1, -1, -1):
        _jacrot(i, u[i], -u[i + 1], M)
        u[i] = (abs(u[i + 1]) if u[i] == 0.0
                else math.sqrt(u[i]*u[i] + u[i + 1]*u[i + 1]))

    # Rank-1 update to row 0
    M[0, :] += u[0] * v

    # Restore upper-triangularity with Givens rotations upward
    # C 1-based: i from 1 to k-1 → 0-based: i from 0 to k-1
    for i in range(k):
        _jacrot(i, M[i, i], -M[i + 1, i], M)


# ── BFGS Cholesky-factor update ───────────────────────────────────────────────

def _bfgsfac(xk, xkp1, gk, gkp1, eta, B):
    """BFGS update of lower-triangular Cholesky factor B (modified in place).

    B @ B.T approximates the Hessian.  Mirrors bfgsfac() in qnewtopt.c.
    """
    s  = xkp1 - xk
    y  = gkp1 - gk
    ys = float(np.dot(y, s))
    ss = float(np.dot(s, s))
    yy = float(np.dot(y, y))

    if ys <= math.sqrt(_MACHEPS * ss * yy):
        return                                  # curvature condition not satisfied

    # t = B^T s  (B lower → B^T upper triangular)
    t  = B.T @ s
    tt = float(np.dot(t, t))                    # ‖B^T s‖²

    alpha = math.sqrt(ys / tt)
    tol   = math.sqrt(eta)

    Bt = B @ t                                  # = B @ (B^T s)

    # Check whether update can be skipped
    ref  = np.maximum(np.abs(gk), np.abs(gkp1))
    diff = np.abs(y - Bt)
    if np.all(diff < tol * ref):
        return

    u       = y - alpha * Bt
    t_sc    = t / math.sqrt(ys * tt)

    # qrupdate operates on upper-triangular R = B^T
    R    = B.T.copy()
    _qrupdate(t_sc, u, R)
    B[:] = np.tril(R.T)


# ── Main optimizer ────────────────────────────────────────────────────────────

def raxopt(func, x0, maxits=500, nrits=10, gradtol=None, steptol=None):
    """Factorized BFGS Quasi-Newton optimizer (Dennis & Schnabel 1983).

    Minimises func starting from x0.  The objective must satisfy
    func(x0) == 1.0 (normalised) — this matches the fue convention where
    objcfunc is normalised by initial sumsq and det-factor.

    Parameters
    ----------
    func    : callable  f(x) -> float, x modified in-place is fine.
    x0      : array_like  initial parameters (copied internally).
    maxits  : int   maximum iterations (default 500, matches fue_defaults).
    nrits   : int   report interval (accepted but ignored in Python).
    gradtol : float gradient stopping tolerance (default macheps^(1.1/3)).
    steptol : float step-size stopping tolerance (default macheps^(2/3)).

    Returns
    -------
    x        : ndarray  optimal parameters.
    f        : float    objective value at x.
    B        : ndarray (n×n) lower-triangular Cholesky factor of the
               approximate Hessian at x.  Use _cholsol(B, e_i) to get
               columns of H^{-1} for covariance computation.
    termcode : int  see module docstring for codes.
    niter    : int  number of iterations performed.
    gnorm    : float  gradient norm at termination.
    """
    if gradtol is None:
        gradtol = _GRTOL
    if steptol is None:
        steptol = _SPTOL

    x0 = np.asarray(x0, dtype=float)
    n  = len(x0)
    if n == 0:
        return x0.copy(), 1.0, np.eye(0), 1, 0, 0.0

    xk = x0.copy()

    # Dennis-Schnabel parameters
    maxstep = max(float(np.linalg.norm(xk)), 1.0) * 1.0e3
    eta     = _MACHEPS

    # Initialise: fk=1.0 is correct because objective is normalised to 1 at x0
    k          = 0
    fk         = 1.0
    gk         = cdgrad(func, xk, eta)
    B          = np.eye(n)
    consecmax  = 0

    termcode = _umstop0(xk, fk, gk, gradtol, maxits)

    while termcode == 0:
        # Search direction: solve B B^T dk = -gk
        dk = _cholsol(B, -gk)

        # Line search
        xkp1, fkp1, retcode, maxtaken, _lam = _lnsrch(
            func, xk, fk, gk, dk, maxstep, steptol)

        # Gradient at new point
        gkp1 = cdgrad(func, xkp1, eta)

        # Convergence check
        termcode, consecmax = _umstop(
            xk, xkp1, fkp1, gkp1, retcode,
            gradtol, steptol, k + 1, maxits, _MAXCMAX, maxtaken, consecmax)

        # BFGS Cholesky-factor update
        _bfgsfac(xk, xkp1, gk, gkp1, eta, B)

        k   += 1
        xk   = xkp1.copy()
        gk   = gkp1
        fk   = fkp1

    return xk, fk, B, termcode, k, float(np.linalg.norm(gk))
