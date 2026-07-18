"""
ARMAX forecast engine — pure Python implementation.

Mirrors the usfo.c / fuf.c forecast logic exactly:
  varphi()  → combine stationary AR polynomial with non-stationary operator
  forecast() → recursive level forecasts + psi-weight variances
  point_forecast() → first-difference and seasonal-difference forecasts
"""

import math
import numpy as np
from dataclasses import dataclass


# ── Data class ────────────────────────────────────────────────────────────────

@dataclass
class ForecastResult:
    """Point forecasts and standard errors from Model.forecast()."""
    horizon: int
    level: np.ndarray           # point forecasts, original scale     [L]
    level_std: np.ndarray       # forecast std error, original scale   [L]
    diff1: np.ndarray           # month/quarter change (%)             [L]
    diff1_std: np.ndarray
    seasonal_diff: np.ndarray   # year-on-year change (%)              [L]
    seasonal_diff_std: np.ndarray
    sigma2: float


# ── Box-Cox helpers ───────────────────────────────────────────────────────────

def _boxcox(y, lam, refactor):
    if lam == 0.0:
        return refactor * math.log(y)
    elif lam != 1.0:
        return refactor * (y ** lam - 1.0) / lam
    else:
        return refactor * y

def _inv_boxcox(z, lam, refactor):
    y = z / refactor
    if lam == 0.0:
        return math.exp(y)
    elif lam != 1.0:
        return (y * lam + 1.0) ** (1.0 / lam)
    else:
        return y


# ── Parameter reconstruction ──────────────────────────────────────────────────

def _reconstruct_params(model, params):
    """Return estimated AR/MA/omega/delta/mu from result.params.

    Mirrors the unpack loop in cast_us() (fue_api.c), which matches
    the ordering in count_npar_build_par().

    Returns
    -------
    itv_omega : list[list[float]]  — estimated omega coefs per intervention
    itv_delta : list[list[float]]  — estimated delta coefs per intervention
    ar_est    : list[list[float]]  — estimated regular AR factor coefs
    ar_s_est  : list[list[float]]  — estimated seasonal AR factor coefs
    ma_est    : list[list[float]]  — estimated regular MA factor coefs
    ma_s_est  : list[list[float]]  — estimated seasonal MA factor coefs
    mu        : float
    """
    k = 0

    # Intervention omegas (all interventions, free coefs only)
    itv_omega = []
    for itv in model.interventions:
        om = list(itv.omega)
        for j in range(len(om)):
            if itv.omega_free[j]:
                om[j] = float(params[k]); k += 1
        itv_omega.append(om)

    # Intervention deltas (all interventions, free coefs only)
    itv_delta = []
    for itv in model.interventions:
        dl = list(itv.delta)
        for j in range(len(dl)):
            if itv.delta_free[j]:
                dl[j] = float(params[k]); k += 1
        itv_delta.append(dl)

    def _update_factors(factors, free_lists):
        result = []
        for i, factor in enumerate(factors):
            f = list(factor)
            free = free_lists[i] if free_lists is not None else None
            for j in range(len(f)):
                if free is None or free[j]:
                    f[j] = float(params[k]); globals()['k'] = k + 1
            result.append(f)
        return result

    # Use a mutable container to allow nested assignment of k
    idx = [k]

    def _update(factors, free_lists):
        result = []
        for i, factor in enumerate(factors):
            f = list(factor)
            free = free_lists[i] if free_lists is not None else None
            for j in range(len(f)):
                if free is None or free[j]:
                    f[j] = float(params[idx[0]]); idx[0] += 1
            result.append(f)
        return result

    ar_est   = _update(model.ar,   model.ar_free)
    ar_s_est = _update(model.ar_s, model.ar_s_free)
    ma_est   = _update(model.ma,   model.ma_free)
    ma_s_est = _update(model.ma_s, model.ma_s_free)

    # Fixed-frequency AR coefs (phi2 per factor) — after AR/MA, before mu
    ar_f_coefs = []
    for ff in model.ar_f:
        if ff.free:
            ar_f_coefs.append(float(params[idx[0]])); idx[0] += 1
        else:
            ar_f_coefs.append(ff.coef)

    ma_f_coefs = []
    for ff in model.ma_f:
        if ff.free:
            ma_f_coefs.append(float(params[idx[0]])); idx[0] += 1
        else:
            ma_f_coefs.append(ff.coef)

    mu = float(model.mu0) if model.mu0 else 0.0
    if model.estimate_mu:
        mu = float(params[idx[0]]); idx[0] += 1

    return itv_omega, itv_delta, ar_est, ar_s_est, ma_est, ma_s_est, \
           ar_f_coefs, ma_f_coefs, mu


# ── AR/MA unscramble (replicates unscramble() in fue_api.c) ──────────────────

def _unscramble(factors, s_factors, freq):
    """Convolve regular + seasonal factors into a unified polynomial.

    Returns positive coefficients [c1, c2, ...] for the combined polynomial
    1 - c1*B - c2*B² - ...  (same convention as fue_api.c ArFactor[1..]).

    Regular factors are convolved in B; seasonal in B^freq; then combined.
    """
    # Regular: convolve factor polynomials in B
    # Each factor f = [phi1, phi2, ...] → polynomial [-1, phi1, phi2, ...]
    pr = sum(len(f) for f in factors)
    size_r = max(pr, 1) + 1
    phir = np.zeros(size_r); phir[0] = -1.0
    tmp  = phir.copy()
    p1 = 0
    for factor in factors:
        ar1k = np.empty(len(factor) + 1)
        ar1k[0] = -1.0
        ar1k[1:] = factor
        new_p = np.zeros(size_r); tmp[0] = -1.0
        for i in range(p1 + 1):
            for j in range(len(factor) + 1):
                new_p[j + i] -= ar1k[j] * tmp[i]
        p1 += len(factor)
        for i in range(p1 + 1):
            tmp[i] = new_p[i]
    phir[:p1 + 1] = tmp[:p1 + 1]

    # Seasonal: convolve factor polynomials in B^freq
    pa = sum(len(f) for f in s_factors)
    size_a = max(pa, 1) + 1
    phia = np.zeros(size_a); phia[0] = -1.0
    tmp  = phia.copy()
    p2 = 0
    for factor in s_factors:
        ar2k = np.empty(len(factor) + 1)
        ar2k[0] = -1.0
        ar2k[1:] = factor
        new_p = np.zeros(size_a); tmp[0] = -1.0
        for i in range(p2 + 1):
            for j in range(len(factor) + 1):
                new_p[j + i] -= ar2k[j] * tmp[i]
        p2 += len(factor)
        for i in range(p2 + 1):
            tmp[i] = new_p[i]
    phia[:p2 + 1] = tmp[:p2 + 1]

    # Combine regular × seasonal: ArFactor[j + i*freq] -= phir[j] * phia[i]
    order = pr + freq * pa
    af = np.zeros(order + 1); af[0] = -1.0
    for i in range(p2 + 1):
        for j in range(p1 + 1):
            af[j + i * freq] -= phir[j] * phia[i]

    return af[1:order + 1]   # positive coefs c[1..order]


# ── Non-stationary polynomial ─────────────────────────────────────────────────

# Irreducible factors of (1 - B^s) in standard polynomial form [1, c1, c2, ...].
# Each entry corresponds to ifadf[k] for k = 0, 1, ...
_IFADF_FACTORS = {
    12: [
        [1.0, -1.0],                        # k=0: (1 - B)
        [1.0, -math.sqrt(3.0),  1.0],       # k=1: (1 - √3B + B²)
        [1.0, -1.0,             1.0],       # k=2: (1 - B + B²)
        [1.0,  0.0,             1.0],       # k=3: (1 + B²)
        [1.0,  1.0,             1.0],       # k=4: (1 + B + B²)
        [1.0,  math.sqrt(3.0),  1.0],       # k=5: (1 + √3B + B²)
        [1.0,  1.0],                        # k=6: (1 + B)
    ],
    4: [
        [1.0, -1.0],                        # k=0: (1 - B)
        [1.0,  0.0,  1.0],                  # k=1: (1 + B²)
        [1.0,  1.0],                        # k=2: (1 + B)
    ],
}

# Extra order added to ornsop by each ifadf flag.
_IFADF_ORDERS = {12: [1, 2, 2, 2, 2, 2, 1], 4: [1, 2, 1]}


def _ifadf_ornsop(freq, ifadf):
    """Extra order contributed by ifadf flags."""
    orders = _IFADF_ORDERS.get(freq, [])
    return sum(orders[k] for k in range(min(len(ifadf), len(orders))) if ifadf[k])


def _nonsop_coefs(d, D, freq, ifadf=None):
    """Positive coefficients rnsop[1..ornsop] of the full non-stationary operator.

    Includes (1-B)^d, (1-B^s)^D, and any ifadf individual annual factors.
    Returns 0-indexed array [rnsop1, rnsop2, ...].
    """
    poly = np.array([1.0])
    for _ in range(d):
        poly = np.convolve(poly, [1.0, -1.0])
    if D > 0 and freq > 1:
        seasonal = np.zeros(freq + 1)
        seasonal[0] = 1.0; seasonal[freq] = -1.0
        for _ in range(D):
            poly = np.convolve(poly, seasonal)
    if ifadf and freq in _IFADF_FACTORS:
        fac_list = _IFADF_FACTORS[freq]
        for k, flag in enumerate(ifadf[:len(fac_list)]):
            if flag:
                poly = np.convolve(poly, fac_list[k])
    # poly = [1, -rnsop1, -rnsop2, ...]  →  rnsop[i] = -poly[i]
    return -poly[1:]


# ── varphi (replicates varphi() in usfo.c) ────────────────────────────────────

def _varphi(phi_coefs, rnsop_coefs):
    """Combine stationary AR polynomial with non-stationary operator.

    phi_coefs   : 0-indexed positive AR coefs [phi1, ..., phip]
    rnsop_coefs : 0-indexed positive non-stat coefs [rnsop1, ..., rnsop_ornsop]

    Returns phi0 (0-indexed, length ornsop+p), the combined AR coefficients
    used to make level-series forecasts.  Mirrors varphi() in usfo.c.
    """
    p      = len(phi_coefs)
    ornsop = len(rnsop_coefs)
    order  = ornsop + p
    phi0   = np.zeros(order + 1); phi0[0] = -1.0

    phi_ext   = np.empty(p + 1);      phi_ext[0]   = -1.0; phi_ext[1:]   = phi_coefs
    rnsop_ext = np.empty(ornsop + 1); rnsop_ext[0] = -1.0; rnsop_ext[1:] = rnsop_coefs

    for i in range(ornsop + 1):
        for j in range(p + 1):
            phi0[j + i] -= phi_ext[j] * rnsop_ext[i]

    phi0[0] = 0.0
    return phi0[1:]   # [phi0_1, ..., phi0_{ornsop+p}]


# ── calcnu impulse response (replicates calcnu() in fue_api.c) ────────────────

def _calcnu(omega, delta, lags):
    """Impulse response nu[0..lags] of the ω(B)/δ(B) filter.

    omega : list [ω₀, ω₁, ...] (0-indexed)
    delta : list [δ₁, δ₂, ...] (0-indexed, δ_i = coefficient at lag i)

    Matches calcnu() in fue_api.c: nu[j] = Σ δ_i*nu[j-i] − ω[j].
    """
    s  = len(omega) - 1   # Nomega
    r  = len(delta)        # Ndelta
    nu = np.zeros(lags + 1)
    nu[0] = omega[0]
    for j in range(1, lags + 1):
        sum1 = sum(delta[i - 1] * nu[j - i] for i in range(1, min(j, r) + 1))
        sum2 = omega[j] if j <= s else 0.0
        nu[j] = sum1 - sum2
    return nu


# ── Build indicator and xi ────────────────────────────────────────────────────

def _build_xi(model, nobs, freq, horizon, itv_omega, itv_delta):
    """Total deterministic effect xi[1..nobs+horizon] for the forecast period.

    For each intervention: extend its indicator D to nobs+horizon, compute the
    nu impulse response, and convolve: xi[t] += Σ_k nu[k]*D[t-k].

    Mirrors fue_api.c DataMat[i] generation and fuf.c DataMat extension.
    `itv.at` is 0-based; obs = itv.at + 1 gives the 1-based C obs index.
    """
    T  = nobs + horizon
    xi = np.zeros(T + 1)   # 1-indexed, index 0 unused

    # begtime: starting period within year (1-based), used by seasonal indicator
    begtime = model.series.start[1] if model.series.start else 1

    for idx, itv in enumerate(model.interventions):
        omega  = itv_omega[idx]
        delta  = itv_delta[idx]
        obs    = itv.at + 1    # convert 0-based at → 1-based C obs index

        # NuLag mirrors fue_api.c: Nomega if no delta, else 40
        nu_lag = len(omega) - 1 if not delta else 40
        nu = _calcnu(omega, delta, nu_lag)

        # Build indicator D[1..T]
        D = np.zeros(T + 1)
        itype = itv.type_code

        if itype == 0:    # pulse
            if 1 <= obs <= T:
                D[obs] = 1.0
        elif itype == 1:  # step
            for t in range(max(1, obs), T + 1):
                D[t] = 1.0
        elif itype == 2:  # ramp
            for t in range(max(1, obs), T + 1):
                D[t] = float(t - obs + 1)
        elif itype == 3:  # seasonal: obs is 1-based period within year
            # mirrors fue_api.c: ((j - begtime) % freq) + 1 == obs
            for t in range(1, T + 1):
                if ((t - begtime) % freq) + 1 == obs:
                    D[t] = 1.0
        elif itype == 4:  # cos
            k = itv.harmonic
            for t in range(1, T + 1):
                D[t] = math.cos(2.0 * math.pi * k / freq * t)
        elif itype == 5:  # sin
            k = itv.harmonic
            for t in range(1, T + 1):
                D[t] = math.sin(2.0 * math.pi * k / freq * t)
        elif itype == 6:  # alter: (-1)^t, mirrors fue_api.c j%2==0 → +1, else -1
            for t in range(1, T + 1):
                D[t] = 1.0 if t % 2 == 0 else -1.0
        elif itype == 7:  # custom: use provided data; zero beyond observed range
            if itv.data is not None:
                for t in range(1, min(len(itv.data), T) + 1):
                    D[t] = itv.data[t - 1]

        # Apply filter: xi[t] += Σ_{k=0}^{nu_lag} nu[k] * D[t-k]
        for t in range(1, T + 1):
            for kk in range(nu_lag + 1):
                if t - kk >= 1:
                    xi[t] += nu[kk] * D[t - kk]

    return xi


# ── Main forecast function ────────────────────────────────────────────────────

def forecast(model, result, horizon):
    """Compute L-step-ahead ARMAX forecasts.

    Implements the same algorithm as fuf.c + usfo.c:
      1.  Reconstruct estimated AR/MA/omega/delta from result.params
      2.  Unscramble factor polynomials → phi, theta
      3.  Build combined polynomial phi0 = varphi(phi, rnsop)
      4.  Compute deterministic effects xi[nobs+1..nobs+L]
      5.  Recursive level forecasts using phi0 and the Box-Cox level series
      6.  Psi-weights → forecast error variances (level, 1st diff, seasonal)
      7.  Inverse Box-Cox for original-scale output

    Parameters
    ----------
    model   : Model  (fitted)
    result  : FitResult  (model._result)
    horizon : int

    Returns
    -------
    ForecastResult
    """
    ts      = model.series
    nobs    = ts.nobs
    freq    = ts.freq if ts.freq > 0 else 1
    L       = int(horizon)
    boxlam  = model.boxlam
    refactor = model.refactor
    d       = model.d
    D       = model.D

    params    = np.asarray(result.params)
    sigma2    = result.sigma2
    residuals = np.asarray(result.residuals)   # 0-indexed, length nobs-ornsop

    # [1] Reconstruct estimated coefficients
    itv_omega, itv_delta, ar_est, ar_s_est, ma_est, ma_s_est, \
        ar_f_coefs, ma_f_coefs, mu = _reconstruct_params(model, params)

    # [2a] Expand fixed-freq factors to equivalent [phi1, phi2] regular factors
    #      phi1 = 2·cos(2π·ff.freq/sper)·√(−phi2), phi2 = ar_f_coefs[i]
    sper = freq
    ff_ar = []
    for i, ff in enumerate(model.ar_f):
        phi2 = ar_f_coefs[i]
        phi1 = 2.0 * math.cos(2.0 * math.pi * ff.freq / sper) * math.sqrt(-phi2)
        ff_ar.append([phi1, phi2])

    ff_ma = []
    for i, ff in enumerate(model.ma_f):
        theta2 = ma_f_coefs[i]
        theta1 = 2.0 * math.cos(2.0 * math.pi * ff.freq / sper) * math.sqrt(-theta2)
        ff_ma.append([theta1, theta2])

    # [2b] Unscramble factor polynomials (regular + f-fixed combined, then × seasonal)
    phi_coefs   = _unscramble(list(ar_est) + ff_ar, ar_s_est, freq)   # stationary AR
    theta_coefs = _unscramble(list(ma_est) + ff_ma, ma_s_est, freq)   # MA
    p  = len(phi_coefs)
    q  = len(theta_coefs)

    # [3] Non-stationary polynomial and combined phi0
    ifadf       = getattr(model, 'ifadf', None)
    ornsop      = d + D * freq + _ifadf_ornsop(freq, ifadf or [])
    rnsop_coefs = _nonsop_coefs(d, D, freq, ifadf=ifadf)
    phi0        = _varphi(phi_coefs, rnsop_coefs)   # length ornsop + p
    p0          = len(phi0)

    # [4] Box-Cox level series nt[1..nobs]  (1-indexed)
    raw = ts.data   # 0-indexed numpy array
    nt  = np.empty(nobs + 1)
    for t in range(1, nobs + 1):
        nt[t] = _boxcox(raw[t - 1], boxlam, refactor)

    # [5] Deterministic effects xi[1..nobs+L]
    xi = _build_xi(model, nobs, freq, L, itv_omega, itv_delta)

    # [6] Point forecasts f1[1..L]  (level, Box-Cox space)
    # Mirrors the forecast() loop in usfo.c, but for the univariate case.
    #
    # n = nobs, b = 0  (fuf.c: Fs.b=0)
    # w = varma1.nt = stochastic level = Box-Cox(data) - xi  (NOT raw Box-Cox!)
    # AR part : phi0[i] * (f1[l-i]  if l>i  else  w[nobs-i+l])
    # MA part : theta[j] * a[nobs-j+l-ornsop]  for l≤j  (else 0)
    #
    # BUG-0001: the mean drift is the constant intercept c = μ·φ(1) of the level
    # recursion phi0(B)w = μ·φ(1) + a, added HERE each step — NOT as an accumulated
    # l·μ afterwards (that double-counted the drift, over-shooting the level by
    # μ·φ/(1-φ) in the transient, because the homogeneous recursion is seeded with
    # drift-laden initial conditions). φ(1) = 1 - Σ φ_i, d-independent.
    drift = mu * (1.0 - float(np.sum(phi_coefs))) if mu else 0.0

    f1 = np.zeros(L + 1)   # f1[0] unused

    for l in range(1, L + 1):
        vtmp1 = 0.0
        for i in range(1, p0 + 1):
            if l > i:
                vtmp1 += phi0[i - 1] * f1[l - i]
            else:
                idx_nt = nobs - i + l
                if 1 <= idx_nt <= nobs:
                    vtmp1 += phi0[i - 1] * (nt[idx_nt] - xi[idx_nt])

        vtmp2 = 0.0
        for j in range(1, q + 1):
            if l <= j:
                # C 1-indexed: a[nobs - j + l - ornsop]
                # Python 0-indexed: residuals[nobs - j + l - ornsop - 1]
                py_idx = nobs - j + l - ornsop - 1
                if 0 <= py_idx < len(residuals):
                    vtmp2 += theta_coefs[j - 1] * residuals[py_idx]

        f1[l] = vtmp1 - vtmp2 + drift

    # [7] Add the deterministic effect (interventions/harmonics) at the forecast
    # dates.  The mean drift is NOT added here — it is the intercept in [6].
    if len(model.interventions) > 0:
        for l in range(1, L + 1):
            f1[l] += xi[nobs + l]

    # [8] Psi-weights ψ[0..L]  (impulse response of phi0 and theta)
    psi = np.zeros(L + 1)
    psi[0] = 1.0
    for j in range(1, L + 1):
        for i in range(1, min(j, p0) + 1):
            psi[j] += phi0[i - 1] * psi[j - i]
        if j <= q:
            psi[j] -= theta_coefs[j - 1]

    # [9] Forecast error variances (level v1, 1st-diff v2, seasonal-diff v3)
    psi2 = np.zeros(L + 1)
    psi2[0] = psi[0]
    for j in range(1, L + 1):
        psi2[j] = psi[j] - psi[j - 1]

    psi3 = np.zeros(L + 1)
    for j in range(min(freq, L + 1)):
        psi3[j] = psi[j]
    for j in range(freq, L + 1):
        psi3[j] = psi[j] - psi[j - freq]

    v1 = np.array([sigma2 * sum(psi[j] ** 2 for j in range(l)) for l in range(1, L + 1)])
    v2 = np.array([sigma2 * sum(psi2[j] ** 2 for j in range(l)) for l in range(1, L + 1)])
    v3 = np.array([sigma2 * sum(psi3[j] ** 2 for j in range(l)) for l in range(1, L + 1)])

    # [10] point_forecast: f2[l] = Δf1[l],  f3[l] = Δ_freq f1[l]
    f2 = np.zeros(L + 1)
    f3 = np.zeros(L + 1)
    f2[1] = f1[1] - nt[nobs]
    for l in range(2, L + 1):
        f2[l] = f1[l] - f1[l - 1]
    for l in range(1, min(freq, L) + 1):
        f3[l] = f1[l] - nt[max(1, nobs - freq + l)]
    for l in range(freq + 1, L + 1):
        f3[l] = f1[l] - f1[l - freq]

    # [11] Inverse Box-Cox and output assembly
    level = np.array([_inv_boxcox(f1[l], boxlam, refactor) for l in range(1, L + 1)])
    level_std = np.sqrt(v1) / refactor

    diff1         = 100.0 * np.array([f2[l] for l in range(1, L + 1)]) / refactor
    diff1_std     = 100.0 * np.sqrt(v2) / refactor
    seasonal_diff = 100.0 * np.array([f3[l] for l in range(1, L + 1)]) / refactor
    seasonal_diff_std = 100.0 * np.sqrt(v3) / refactor

    return ForecastResult(
        horizon=L,
        level=level,
        level_std=level_std,
        diff1=diff1,
        diff1_std=diff1_std,
        seasonal_diff=seasonal_diff,
        seasonal_diff_std=seasonal_diff_std,
        sigma2=sigma2,
    )
