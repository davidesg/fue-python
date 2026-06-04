"""
cast_us_py — pure-Python estimator core.

Mirrors the C chain:  populate_globals → cast_us → est(cast_us, …) → elf.

Public entry points
-------------------
build_est_spec(model)        Pre-compute fixed data (DataMat, rnsop, …).
cast_us_py(x, est_spec)      Map parameter vector x → (p, q, phi, theta, mu, w, ifault).
estimate_py(model)           Full pure-Python ML estimator; returns the same dict
                             as _engine.estimate().
"""

import math
import numpy as np
from scipy.optimize import minimize

from .elfvarma  import flikam_scalar, elf_scalar
from .forecast  import (
    _reconstruct_params, _unscramble, _nonsop_coefs,
    _ifadf_ornsop, _boxcox,
)

_LOG2PI   = 1.837877066
_SQRT_EPS = math.sqrt(2.2e-16)   # sqrt(machine epsilon)


# ── calcnu_py ─────────────────────────────────────────────────────────────────

def calcnu_py(omega, delta, lags):
    """Impulse response of ω(B)/δ(B) up to lag *lags*.

    Mirrors calcnu() in fue_api.c (extracted from fue.c:4489).

    Parameters
    ----------
    omega : sequence  [ω₀, ω₁, …, ωₛ]  (length Nomega+1)
    delta : sequence  [δ₁, δ₂, …, δᵣ]  (length Ndelta, 0-indexed)
    lags  : int

    Returns
    -------
    nu : ndarray (lags+1,)
    """
    s  = len(omega) - 1        # Nomega
    r  = len(delta)            # Ndelta
    nu = np.zeros(lags + 1)
    nu[0] = omega[0]
    for j in range(1, lags + 1):
        s1 = 0.0
        if r > 0:
            for i in range(1, min(j, r) + 1):
                s1 += delta[i - 1] * nu[j - i]
        s2 = float(omega[j]) if j <= s else 0.0
        nu[j] = s1 - s2
    return nu


# ── indicator series (DataMat[i]) ─────────────────────────────────────────────

def _build_indicator(itv, nobs, freq, begtime):
    """Return 1-indexed indicator array (shape nobs+1; index 0 unused).

    Mirrors the DataMat[i] construction in populate_globals() / fue_api.c.
    begtime = ts.start[1] (the period within the year of observation 1).
    """
    ind = np.zeros(nobs + 1)
    obs = itv.at + 1          # 0-based → 1-based
    t   = itv.type

    if t == "pulse":
        if 1 <= obs <= nobs:
            ind[obs] = 1.0
    elif t == "step":
        if 1 <= obs <= nobs:
            ind[obs:nobs + 1] = 1.0
    elif t == "ramp":
        if 1 <= obs <= nobs:
            for j in range(obs, nobs + 1):
                ind[j] = float(j - obs + 1)
    elif t == "seasonal":
        for j in range(1, nobs + 1):
            if ((j - begtime) % freq) + 1 == obs:
                ind[j] = 1.0
    elif t == "cos":
        k = itv.harmonic
        for j in range(1, nobs + 1):
            ind[j] = math.cos(2.0 * math.pi * k / freq * j)
    elif t == "sin":
        k = itv.harmonic
        for j in range(1, nobs + 1):
            ind[j] = math.sin(2.0 * math.pi * k / freq * j)
    elif t == "alter":
        for j in range(1, nobs + 1):
            ind[j] = 1.0 if j % 2 == 0 else -1.0
    elif t == "custom" and itv.data is not None:
        ind[1:nobs + 1] = itv.data[:nobs]
    return ind


# ── EstSpec ───────────────────────────────────────────────────────────────────

class EstSpec:
    """Pre-computed estimation state (equiv. to C globals Tm + Ts + DataMat).

    Fields
    ------
    model     : Model
    nobs      : int
    sper      : int   seasonal period (≥1)
    eml       : bool  exact ML flag
    chkma     : bool  check MA invertibility
    data0     : ndarray (nobs+1,)  BoxCox-transformed series, 1-indexed
    ind_data  : list[ndarray]  indicator series per intervention, 1-indexed
    rnsop     : ndarray  positive non-stat operator coefs [r1, r2, …]
    ornsop    : int  len(rnsop)
    """
    __slots__ = ("model", "nobs", "sper", "eml", "chkma",
                 "data0", "ind_data", "rnsop", "ornsop")


def build_est_spec(model):
    """Build EstSpec from a Model — called once before optimization."""
    ts      = model.series
    nobs    = ts.nobs
    freq    = ts.freq if ts.freq > 0 else 1
    begtime = ts.start[1] if ts.start else 1

    spec         = EstSpec()
    spec.model   = model
    spec.nobs    = nobs
    spec.sper    = freq
    spec.eml     = model.eml
    spec.chkma   = model.chkma

    # DataMat[0]: BoxCox-transformed series (1-indexed)
    raw   = ts.data
    data0 = np.empty(nobs + 1)
    for i in range(1, nobs + 1):
        data0[i] = _boxcox(raw[i - 1], model.boxlam, model.refactor)
    spec.data0 = data0

    # DataMat[1..k]: intervention indicator series
    spec.ind_data = [
        _build_indicator(itv, nobs, freq, begtime)
        for itv in model.interventions
    ]

    # Non-stationary operator
    ifadf       = model.ifadf or []
    rnsop       = _nonsop_coefs(model.d, model.D, freq, ifadf=ifadf)
    spec.rnsop  = rnsop
    spec.ornsop = len(rnsop)

    return spec


# ── cast_us_py ────────────────────────────────────────────────────────────────

def cast_us_py(x, est_spec):
    """Translate parameter vector x → (p, q, phi, theta, mu, w, ifault).

    Mirrors cast_us() in fue_api.c.  x contains only the *free* parameters
    in the same order as count_npar_build_par():
      1. omega free  per intervention
      2. delta free  per intervention
      3. AR regular  free coefs (factor by factor)
      4. AR seasonal free coefs
      5. MA regular  free coefs  (MA(1) |θ|>1 flip applied)
      6. MA seasonal free coefs
      7. AR fixed-freq free coef2 per factor
      8. MA fixed-freq free coef2 (|θ₂|>1 flip applied)
      9. mu  if estimate_mu

    Returns
    -------
    p, q      : int  AR and MA orders of the expanded polynomial
    phi       : ndarray (p,)
    theta     : ndarray (q,)
    mu        : float
    w         : ndarray (n_eff,)  differenced noise series (0-indexed)
    ifault    : 0 OK, 1 bad fixed-freq coef (coef2 ≥ 0)
    """
    model = est_spec.model
    sper  = est_spec.sper

    # [1] Unpack x → current parameter values (mirrors cast_us step [1])
    itv_omega, itv_delta, ar_est, ar_s_est, ma_est, ma_s_est, \
        ar_f_coefs, ma_f_coefs, mu = _reconstruct_params(model, x)

    # [2] MA(1) invertibility flip: |θ₁| > 1 → invert (mirrors unscramble [6])
    for k, factor in enumerate(ma_est):
        if len(factor) == 1 and abs(factor[0]) > 1.0:
            ma_est[k] = [1.0 / factor[0]]

    # [3] Fixed-freq AR → equivalent 2-lag regular factor [c1, c2]
    #     c1 = 2·cos(2π·freq/sper)·√(−c2),  c2 must be < 0
    ff_ar = []
    for i, ff in enumerate(model.ar_f):
        c2 = ar_f_coefs[i]
        if c2 >= 0.0:
            return 0, 0, np.array([]), np.array([]), 0.0, np.array([]), 1
        c1 = 2.0 * math.cos(2.0 * math.pi * ff.freq / sper) * math.sqrt(-c2)
        ff_ar.append([c1, c2])

    # [4] Fixed-freq MA → [c1, c2] with |θ₂| > 1 flip (mirrors unscramble [7])
    ff_ma = []
    for i, ff in enumerate(model.ma_f):
        c2 = ma_f_coefs[i]
        if c2 < -1.0:
            c2 = 1.0 / c2
        if c2 >= 0.0:
            return 0, 0, np.array([]), np.array([]), 0.0, np.array([]), 1
        c1 = 2.0 * math.cos(2.0 * math.pi * ff.freq / sper) * math.sqrt(-c2)
        ff_ma.append([c1, c2])

    # [5] Expand combined AR/MA polynomials (regular + fixed-freq as extra factors)
    phi   = np.array(_unscramble(list(ar_est) + ff_ar, ar_s_est, sper))
    theta = np.array(_unscramble(list(ma_est) + ff_ma, ma_s_est, sper))
    p, q  = len(phi), len(theta)

    # [6] Build intervention-filtered series z[1..nobs] (1-indexed)
    nobs = est_spec.nobs
    z    = est_spec.data0.copy()   # start from BoxCox data

    for j, itv in enumerate(model.interventions):
        # Lags for impulse response: Nomega if pure FIR, else 40
        lags = len(itv_omega[j]) - 1 if len(itv_delta[j]) == 0 else 40
        nu   = calcnu_py(itv_omega[j], itv_delta[j], lags)
        ind  = est_spec.ind_data[j]
        for t in range(1, nobs + 1):
            acc = 0.0
            for k in range(lags + 1):
                if t - k >= 1:
                    acc += nu[k] * ind[t - k]
            z[t] -= acc

    # [7] Apply non-stationary operator → w[0..n_eff-1]  (0-indexed output)
    #     w[t] = z[ornsop+1+t] − Σ_{i=0}^{ornsop-1} rnsop[i]·z[ornsop+t-i]
    ornsop = est_spec.ornsop
    rnsop  = est_spec.rnsop
    n_eff  = nobs - ornsop
    if n_eff <= 0:
        return p, q, phi, theta, mu, np.array([]), 0

    w = np.empty(n_eff)
    for idx in range(n_eff):
        j   = ornsop + 1 + idx    # 1-indexed position in z
        val = z[j]
        for i in range(ornsop):
            val -= rnsop[i] * z[j - 1 - i]
        w[idx] = val

    return p, q, phi, theta, mu, w, 0


# ── Initial parameter vector ──────────────────────────────────────────────────

def _build_initial_x(model):
    """Build x0 from model initial values — mirrors count_npar_build_par."""
    x = []

    for itv in model.interventions:
        for j, v in enumerate(itv.omega):
            if itv.omega_free[j]:
                x.append(float(v))

    for itv in model.interventions:
        for j, v in enumerate(itv.delta):
            if itv.delta_free[j]:
                x.append(float(v))

    def _add_factors(factors, free_lists):
        for i, factor in enumerate(factors):
            free = free_lists[i] if free_lists is not None else None
            for j, v in enumerate(factor):
                if free is None or free[j]:
                    x.append(float(v))

    _add_factors(model.ar,   model.ar_free)
    _add_factors(model.ar_s, model.ar_s_free)
    _add_factors(model.ma,   model.ma_free)
    _add_factors(model.ma_s, model.ma_s_free)

    for ff in model.ar_f:
        if ff.free:
            x.append(float(ff.coef))
    for ff in model.ma_f:
        if ff.free:
            x.append(float(ff.coef))

    if model.estimate_mu:
        x.append(float(model.mu0))

    return np.array(x, dtype=float)


# ── Finite-difference Hessian (mirrors fdhess in qnewtopt.c) ─────────────────

def _fdhess(f, x, f0, eta):
    """Second-order finite-difference Hessian, symmetric (n×n)."""
    n  = len(x)
    H  = np.zeros((n, n))
    dx = eta * np.maximum(np.abs(x), 1.0)

    # Diagonal
    for i in range(n):
        xi = x[i]
        x[i] = xi + dx[i]; fp = f(x)
        x[i] = xi - dx[i]; fm = f(x)
        x[i] = xi
        H[i, i] = (fp - 2.0 * f0 + fm) / (dx[i] ** 2)

    # Off-diagonal
    for i in range(n):
        for j in range(i + 1, n):
            xi, xj = x[i], x[j]
            x[i] += dx[i]; x[j] += dx[j]; fpp = f(x)
            x[i] -= 2*dx[i];              fmp = f(x)
            x[j] -= 2*dx[j];              fmm = f(x)
            x[i] += 2*dx[i];              fpm = f(x)
            x[i] = xi; x[j] = xj
            H[i, j] = H[j, i] = (fpp - fmp - fpm + fmm) / (4.0 * dx[i] * dx[j])

    return H


# ── logelf matching the C drvmlest.c formula ─────────────────────────────────

def _logelf_c(n, f1, f2):
    """Concentrated exact log-likelihood matching drvmlest.c:

        logelf = −½n·(LOG2π − log n + 1) − ½n·(log f1 + log f2)

    f1 = sumsq from elf_scalar (sigma²=1), f2 = det_factor = exp(log|det|/n).
    """
    if f1 <= 0.0 or f2 <= 0.0:
        return -1e30
    return (-0.5 * n * (_LOG2PI - math.log(n) + 1.0)
            - 0.5 * n * (math.log(f1) + math.log(f2)))


# ── estimate_py ───────────────────────────────────────────────────────────────

def estimate_py(model):
    """Pure-Python ML estimator.  Returns the same dict as _engine.estimate().

    Optimizer : scipy L-BFGS-B minimising (sumsq/sumsq0)·(fact/fact0)
                (normalised like C's objcfunc, using flikam_scalar).
    Final eval : elf_scalar with compute_residuals=True for exact residuals
                and the concentrated logelf matching fue.c.
    """
    ts   = model.series
    sper = ts.freq if ts.freq > 0 else 1

    spec    = build_est_spec(model)
    x0      = _build_initial_x(model)
    npar    = len(x0)
    ornsop  = spec.ornsop
    n_eff   = ts.nobs - ornsop

    _empty = {
        "ifault": -1, "npar": npar, "nresiduals": n_eff,
        "sigma2": 0.0, "loglik": 0.0, "aic": 0.0, "bic": 0.0,
        "params": x0.copy(), "std_errors": np.zeros(npar),
        "cov_matrix": np.zeros((npar, npar)),
        "residuals": np.zeros(n_eff),
    }

    xitol    = 1e-3
    do_chkma = model.chkma

    # [1] Evaluate at initial parameters
    p0, q0, phi0, theta0, mu0, w0, fault0 = cast_us_py(x0, spec)
    if fault0 or len(w0) == 0:
        return {**_empty, "ifault": 6}

    sumsq0, fact0, _, _, iflt0 = flikam_scalar(
        n_eff, p0, q0, phi0, theta0, mu0, w0,
        xitol=xitol, do_chkma=do_chkma,
    )
    if iflt0:
        return {**_empty, "ifault": iflt0}
    if sumsq0 <= 0.0 or fact0 <= 0.0:
        return {**_empty, "ifault": 3}

    # [2] Objective function: normalised like C's objcfunc
    def objective(x):
        p, q, phi, theta, mu, w, fault = cast_us_py(x, spec)
        if fault or len(w) == 0:
            return 1.0
        sumsq, fact, _, _, iflt = flikam_scalar(
            n_eff, p, q, phi, theta, mu, w,
            xitol=xitol, do_chkma=do_chkma,
        )
        if iflt or sumsq <= 0.0 or fact <= 0.0:
            return 1.0
        return (sumsq / sumsq0) * (fact / fact0)

    # [3] Optimize
    if npar > 0:
        opt = minimize(
            objective, x0,
            method="L-BFGS-B",
            options={"maxiter": 500, "ftol": 1e-14, "gtol": 1e-7},
        )
        x_opt    = opt.x
        obj_opt  = opt.fun
        converged = opt.success
    else:
        x_opt    = x0.copy()
        obj_opt  = 1.0
        converged = True

    # [4] Final evaluation with elf_scalar (exact ML)
    p, q, phi, theta, mu, w, fault = cast_us_py(x_opt, spec)
    if fault or len(w) == 0:
        return {**_empty, "ifault": 6}

    _, f1, f2, a_res, ifault_e = elf_scalar(
        n_eff, p, q, phi, theta, w,
        sigma2=1.0, mu=mu, do_chkma=do_chkma, compute_residuals=True,
    )
    if ifault_e:
        return {**_empty, "ifault": ifault_e}

    sigma2_hat = f1 / n_eff if n_eff > 0 else 1.0
    logelf_c   = _logelf_c(n_eff, f1, f2)
    aic        = -2.0 * logelf_c + 2.0 * npar
    bic        = -2.0 * logelf_c + npar * math.log(n_eff) if n_eff > 0 else 0.0

    # [5] Standard errors via finite-difference Hessian of the objective
    #     Mirrors fdhess → cholsol → dev[i] = sqrt(cov[i][i]) in drvmlest.c.
    cov        = np.zeros((npar, npar))
    std_errors = np.zeros(npar)
    if npar > 0:
        H = _fdhess(objective, x_opt.copy(), obj_opt, _SQRT_EPS)
        try:
            # cov[j,i] = 2·obj·(H⁻¹)[j,i] / n_eff
            H_inv = np.linalg.inv(H)
            cov   = 2.0 * obj_opt * H_inv / n_eff
            diag  = np.diag(cov)
            std_errors = np.sqrt(np.maximum(diag, 0.0))
        except np.linalg.LinAlgError:
            converged = False

    ifault_final = 0 if converged else 6

    return {
        "ifault":     ifault_final,
        "npar":       npar,
        "nresiduals": n_eff,
        "sigma2":     sigma2_hat,
        "loglik":     logelf_c,
        "aic":        aic,
        "bic":        bic,
        "params":     x_opt,
        "std_errors": std_errors,
        "cov_matrix": cov,
        "residuals":  a_res,
    }
