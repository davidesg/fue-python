"""
elfvarma.py  —  Log-likelihood for scalar ARMA(p,q) models.

Two algorithms, both specialized to the univariate (m=1) case:

  flikam_scalar   Mélard (1984) AS 197.  Used in the optimizer inner loop.
  elf_scalar      Mauricio (1995) JASA 90, 282-291.  Used for final evaluation.

fue always uses m=1; the general vector-ARMA case is not needed here.

All arrays are 0-indexed NumPy arrays.  phi[0]=φ₁, theta[0]=θ₁, etc.
qq is always 1.0 (fue_api.c normalises by sigma2).

Copyright (C) 1995 José Alberto Mauricio (original C); 2026 Python port.
"""

import numpy as np
from numpy.linalg import solve, eigvals
from scipy.linalg import cholesky, solve_triangular

_LOG2PI = 1.837877066


# ── MA invertibility check ─────────────────────────────────────────────────────

def chekma_scalar(theta):
    """Return 1 if MA companion matrix has an eigenvalue with |λ| ≥ 1.00005.

    The companion matrix eigenvalues are the reciprocals of the roots of
    Θ(z) = 1 − θ₁z − … − θ_q z^q.  Invertibility requires all |λ| < 1.
    """
    q = len(theta)
    if q == 0:
        return 0
    # Companion matrix (same construction as fue's chekma for m=1)
    A = np.zeros((q, q))
    for k in range(q):
        A[k, 0] = theta[k]          # first column = theta values
    for j in range(q - 1):
        A[j, j + 1] = 1.0           # superdiagonal
    return 1 if np.any(np.abs(eigvals(A)) >= 1.00005) else 0


# ── cgamma: autocovariances and cross-covariances for scalar ARMA(p,q) ─────────

def _cgamma_scalar(p, q, phi, theta):
    """Compute initial covariances needed by elf_scalar.

    Returns
    -------
    gamma_arr : ndarray, shape (p,)
        gamma_arr[k] = E[w_t w_{t-k}] for k = 0 … p-1.
    gamwa : dict  {int → float}
        gamwa[k] = E[w_t a_{t-k}] for k = 0, -1, …, -(q-1).
    ifault : int
        1 if the Yule-Walker system is singular (AR near unit-root).
    """
    # ── Cross-covariances gamwa[0 .. -(q-1)] ──────────────────────────────────
    gamwa = {}
    if q > 0:
        gamwa[0] = 1.0                          # = qq = 1.0
    for k in range(1, q):                       # k = 1 … q-1
        s = -theta[k - 1]                       # -θ_k * qq
        for l in range(1, k + 1):
            if l <= p:
                s += phi[l - 1] * gamwa[l - k]
        gamwa[-k] = s

    # ── wzero: RHS for the lag-0 autocovariance equation ──────────────────────
    wzero = 0.0
    for i in range(1, p + 1):
        for j in range(i, q + 1):
            wzero += phi[i - 1] * gamwa[i - j] * theta[j - 1]
    wzero = 1.0 - 2.0 * wzero                  # qq − wzero − wzero
    for j in range(1, q + 1):
        wzero += theta[j - 1] ** 2             # MA variance contribution

    if p == 0:
        return np.array([]), gamwa, 0

    # ── Linear system for gamma[0 … p-1] (size p × p for m=1) ───────────────
    big = p
    mat = np.zeros((big, big))
    rhs = np.zeros(big)

    # Row 0 (C: i=j=1, row=1): lag-0 Yule-Walker row
    s = 0.0
    for r in range(p):
        s -= phi[r] ** 2                        # −Σ φ_r²
    mat[0, 0] = s + 1.0                         # diagonal

    for sv in range(1, p):                      # off-diagonal columns
        s = 0.0
        for r in range(p - sv):
            s -= 2.0 * phi[r] * phi[r + sv]
        mat[0, sv] = s
    rhs[0] = wzero

    # Rows 1 … p-1 (C: s = 1 … p-1)
    for sv in range(1, p):
        row = sv
        mat[row, 0] = -phi[sv - 1]             # first column
        for r in range(1, p):                   # remaining columns
            col = r
            if r + sv <= p:
                mat[row, col] = -phi[r + sv - 1]
            if sv > r:
                mat[row, col] -= phi[sv - r - 1]
        mat[row, row] += 1.0                    # diagonal

        rhs_val = 0.0
        for h in range(sv, q + 1):
            gw = gamwa.get(sv - h, 0.0)
            rhs_val -= gw * theta[h - 1]
        rhs[row] = rhs_val

    try:
        gamma_arr = solve(mat, rhs)
        ifault = 0
    except np.linalg.LinAlgError:
        gamma_arr = np.zeros(big)
        ifault = 1

    return gamma_arr, gamwa, ifault


# ── elf_scalar: exact ML log-likelihood (Mauricio 1995) ───────────────────────

def elf_scalar(n, p, q, phi, theta, w, sigma2=1.0, mu=0.0,
               xitol=1e-3, do_chkma=True, compute_residuals=False):
    """Exact log-likelihood of a scalar ARMA(p,q) model.

    Parameters
    ----------
    n         : number of observations
    p, q      : AR and MA orders
    phi       : ndarray (p,)  AR coefficients φ₁…φ_p  (0-indexed)
    theta     : ndarray (q,)  MA coefficients θ₁…θ_q  (0-indexed)
    w         : ndarray (n,)  mean-subtracted, differenced series
    sigma2    : concentrated noise variance (1.0 during optimisation)
    xitol     : truncation tolerance for Green's function (default 1e-3)
    do_chkma  : check MA invertibility
    compute_residuals : whether to return exact residuals

    Returns
    -------
    logelf : float   exact log-likelihood
    f1     : float   quadratic form  Σa²  (before σ² adjustment)
    f2     : float   determinant factor (see Mauricio 1995 eq. [11])
    a      : ndarray (n,)  residuals (zeros if compute_residuals=False)
    ifault : int    0=OK, 1=Q not PD, 2=AR unit root, 3=non-stationary,
                    4=MA non-invertible, 5=numerical
    """
    phi   = np.asarray(phi,   dtype=float)
    theta = np.asarray(theta, dtype=float)
    w     = np.asarray(w,     dtype=float)
    a     = np.zeros(n)
    g     = max(p, q)
    logelf = f1 = f2 = 0.0

    # [0]: qq = 1.0  →  q1 = 1.0, q1inv = 1.0, detq = 1.0

    # [2]: autocovariances (only needed when p > 0)
    gamma_arr = np.array([])
    gamwa     = {}
    if p > 0:
        gamma_arr, gamwa, ifault = _cgamma_scalar(p, q, phi, theta)
        if ifault:
            return logelf, f1, f2, a, 2

    # [3.1]: mtmp1 = Ω v₁ᵀ  (shape (p+q) × g)
    # For m=1: mtmp1[i-1, j-1] accumulates scalar contributions.
    mtmp1 = np.zeros((p + q, g)) if (p + q) > 0 and g > 0 else np.zeros((1, 1))

    for i in range(1, p + 1):
        for j in range(1, g + 1):
            acc = 0.0
            # AR–AR block
            for k in range(j - i, p - i + 1):
                k_abs = abs(k)
                phi_idx = p - k - i + j - 1    # 0-indexed, range [j-1 .. p-1]
                acc += gamma_arr[k_abs] * phi[phi_idx]
            # AR–MA block
            for k in range(j - i, q - i + 1):
                if p + k <= q:
                    key = -q + p + k
                    theta_idx = q - k - i + j - 1
                    if 0 <= theta_idx < q:
                        acc -= gamwa.get(key, 0.0) * theta[theta_idx]
            mtmp1[i - 1, j - 1] = acc

    for i in range(p + 1, p + q + 1):
        for j in range(1, g + 1):
            acc = 0.0
            # MA–AR block
            for k in range(p + j - i, p + p - i + 1):
                if p - k <= q and 0 <= p + p - k - i + j - 1 < p:
                    key = -q + p - k
                    phi_idx = p + p - k - i + j - 1
                    acc += gamwa.get(key, 0.0) * phi[phi_idx]
            # MA–MA block (pure MA term, using qq=1)
            if p - i + j <= 0:
                theta_idx = q + p - i + j - 1
                if 0 <= theta_idx < q:
                    acc -= theta[theta_idx]         # qq=1 absorbed
            mtmp1[i - 1, j - 1] = acc

    # [3.2]: mtmp0 = v₁ Ω v₁ᵀ  (g × g, upper triangle then symmetrised)
    mtmp0 = np.zeros((g, g))
    for i in range(1, g + 1):
        for j in range(i, g + 1):
            acc = 0.0
            for k in range(0, p - i + 1):              # AR rows of mtmp1
                phi_idx = p - k - 1                    # phi[p-k-1] (0-indexed)
                row1    = k + i - 1                    # 0-indexed row of mtmp1
                acc += phi[phi_idx] * mtmp1[row1, j - 1]
            for k in range(0, q - i + 1):              # MA rows of mtmp1
                theta_idx = q - k - 1
                row1      = k + p + i - 1
                acc -= theta[theta_idx] * mtmp1[row1, j - 1]
            mtmp0[i - 1, j - 1] = acc
    # Symmetrise
    for i in range(g):
        for j in range(i + 1, g):
            mtmp0[j, i] = mtmp0[i, j]

    # [3.3]: Cholesky of mtmp0 → M (lower triangular).
    # C choldcp special-cases a numerically-zero matrix: sets M=0, ifault=0.
    # This happens e.g. for pure MA(q) at theta≈0 (no initial state uncertainty).
    _sqrteps = 1.49e-8   # sqrt(macheps)
    _maxoffl  = max(abs(mtmp0[i, i]) for i in range(g)) if g > 0 else 0.0
    _zero_M   = (_maxoffl <= _sqrteps)

    if _zero_M:
        M = np.zeros((g, g))
    else:
        try:
            M = cholesky(mtmp0, lower=True)
        except np.linalg.LinAlgError:
            return logelf, f1, f2, a, 3

    # [4.1]: check MA invertibility
    if q > 0 and do_chkma:
        if chekma_scalar(theta):
            return logelf, f1, f2, a, 4

    # [4.2]: Green's function rxi[0..nlim]  (scalar, m=1: rxi[k] is a scalar)
    rxi = np.zeros(n)
    rxi[0] = 1.0
    r = 0
    delta = False
    while not delta and r < n - 1:
        r += 1
        for j in range(1, min(q, r) + 1):
            rxi[r] += theta[j - 1] * rxi[r - j]
        if abs(rxi[r]) < xitol:
            nq, delta = 1, True
            while nq <= q and r < n - 1 and delta:
                nq += 1
                r  += 1
                for j in range(1, min(q, r) + 1):
                    rxi[r] += theta[j - 1] * rxi[r - j]
                if abs(rxi[r]) > xitol:
                    delta = False
            if delta:
                r -= nq
    nlim = r
    # Premultiply by q1inv = 1.0: rxi unchanged (m=1 trivial)

    # [5.1]: conditional residuals a[t] = (w[t]-mu) − Σφ_j(w[t-j]-mu) + Σθ_j a[t-j]
    for t in range(n):
        val = w[t] - mu
        for j in range(1, p + 1):
            if t - j >= 0:
                val -= phi[j - 1] * (w[t - j] - mu)
        for j in range(1, q + 1):
            if t - j >= 0:
                val += theta[j - 1] * a[t - j]
        a[t] = val
    # [5.2]: premultiply by q1inv = 1.0: a unchanged

    # [6.1]: h = Σ_{j=1..g} Σ_{i=0..nlim} rxi[i] * a[i+j]
    vechh = np.zeros(g)
    for j in range(1, g + 1):
        acc = 0.0
        for i in range(0, n - j + 1):
            if i <= nlim:
                acc += rxi[i] * a[i + j - 1]
        vechh[j - 1] = acc

    # [6.2]: premultiply by Mᵀ; if M=0, vechh becomes 0
    # C: vechh[i] = Σ_{k>=i} M[k,i] * vechh[k]  ≡  M^T @ vechh  (matrix-vector product)
    if not _zero_M:
        vechh = M.T @ vechh
    else:
        vechh = np.zeros(g)

    # Save M for residuals if needed
    M_save = M.copy() if compute_residuals else None

    # [7]: H'H  (g × g symmetric)
    HtH = np.zeros((g, g))
    for i in range(1, g + 1):
        for k in range(0, n - i + 1):
            if k + i - 1 <= nlim:
                jl = 1 if i == 1 else g      # upper block width
                for jj in range(1, (1 if i == 1 else g) + 1):
                    HtH[i - 1, jj - 1] += rxi[k] * rxi[k + i - 1]
    # Actually redo [7] faithfully:
    HtH = np.zeros((g, g))
    for i in range(1, g + 1):
        for k in range(0, n - i + 1):
            if k + i - 1 <= nlim:
                # first block (i=1): all j from 1..1 (m=1)
                jl = 1
                for jj in range(1, 2):
                    HtH[i - 1, jj - 1] += rxi[k] * rxi[k + i - 1]
    for i in range(2, g + 1):
        for j in range(2, i + 1):
            # Cumulative update: HtH[i-1,j-1] = HtH[i-2,j-2] − rxi[n-i]*rxi[n-j]
            s = 0.0
            if (n - i <= nlim) and (n - j <= nlim):
                s = rxi[n - i] * rxi[n - j]
            HtH[i - 1, j - 1] = HtH[i - 2, j - 2] - s
    # Symmetrise
    for i in range(g):
        for j in range(i + 1, g):
            HtH[i, j] = HtH[j, i]

    # [8.1-8.3]: L = chol(I + Mᵀ H'H M)
    # If M=0 then I + Mᵀ H'H M = I, so L=I and detom=1.
    if _zero_M:
        L     = np.eye(g)
        detom = 1.0
    else:
        MtHtH     = M.T @ HtH
        ImMtHtHM  = np.eye(g) + MtHtH @ M
        try:
            L = cholesky(ImMtHtHM, lower=True)
        except np.linalg.LinAlgError:
            return logelf, f1, f2, a, 5
        detom = float(np.prod(np.diag(L)) ** 2)

    # [9]: forward substitution L λ = vechh  (if M=0, vechh=0 → stays 0)
    if not _zero_M:
        vechh = solve_triangular(L, vechh, lower=True)

    # [10]: quadratic form f1 = Σa² − Σλ²
    f1 = float(np.dot(a, a) - np.dot(vechh, vechh))

    # [11]: f2 = detq * exp(log(detom)/n)  with detq=1.0
    log_detom = np.log(detom) if detom > 0 else 0.0
    f2 = float(np.exp(log_detom / n))

    # [12]: log-likelihood
    logelf = float(-0.5 * (n * (_LOG2PI + np.log(sigma2))
                   + log_detom + f1 / sigma2))

    # [13]: exact residuals via cres
    if compute_residuals:
        a = _cres_scalar(n, g, nlim, rxi, M_save, L, vechh, a)

    return logelf, f1, f2, a, 0


def _cres_scalar(n, g, nlim, rxi, M, L, lambda_vec, a):
    """Compute exact residuals (cres for m=1).

    Steps mirror cres() in elfvarma.c:
      [1] back-solve Lᵀ c = lambda_vec  → overwrite lambda_vec
      [2] d = M c                        → overwrite lambda_vec
      [3] subtract correction from a
    """
    # [1]: cholbak: Lᵀ c = lambda_vec  (L lower → Lᵀ upper)
    c = solve_triangular(L, lambda_vec, trans='T', lower=True)

    # [2]: d = M c  (M lower triangular, stored column-major in C but here as numpy)
    d = np.zeros(g)
    for i in range(g - 1, -1, -1):
        d[i] = np.dot(M[i, :i + 1], c[:i + 1])

    # [3]: a[t] -= Σ_{j=1..g} Σ_{t-j=0..nlim, j<=g} rxi[t-j] * d[j-1]
    for t in range(n):
        for j in range(1, g + 1):
            i_rxi = t - j + 1                   # 0-indexed t+1-j
            if 0 <= i_rxi <= nlim and j <= g:
                a[t] -= rxi[i_rxi] * d[j - 1]
    # q1 = 1.0: premultiplication is trivial
    return a


# ── flikam_scalar: Mélard (1984) approximate ML ───────────────────────────────
#
# Used inside the BFGS optimizer (called at every function evaluation).
# Much faster than elf_scalar; gives exact results as n→∞.

def _twacf_scalar(p, q, phi, theta):
    """Compute ACF and related quantities for ARMA(p,q).

    Direct Python port of twacf() in usmelard.c (m=1 specialisation).

    Returns
    -------
    acf   : ndarray (mxpq+1,)   autocovariance function
    cvli  : ndarray (mxpqp1,)   covariance between w_t and a_{t-k}
    alpha : ndarray (mxpq,)     partial autocorrelation work array
    ifault: int   1 if near unit root in AR, else 0
    """
    mxpq   = max(p, q)
    mxpqp1 = mxpq + 1
    ma     = mxpqp1

    acf   = np.zeros(ma + 1)
    cvli  = np.zeros(mxpqp1 + 1)
    alpha = np.zeros(mxpq + 1)
    epsil2 = 1e-10

    acf[1]  = 1.0
    cvli[1] = 1.0

    if ma == 1:
        return acf, cvli, alpha, 0
    for i in range(2, ma + 1):
        acf[i] = 0.0
    for i in range(2, mxpqp1 + 1):
        cvli[i] = 0.0
    for k in range(1, mxpq + 1):
        alpha[k] = 0.0

    # [1]: MA ACF
    if q > 0:
        for k in range(1, q + 1):
            cvli[k + 1] = -theta[k - 1]
            acf[k + 1]  = -theta[k - 1]
            if q != k:
                for j in range(1, q - k + 1):
                    acf[k + 1] += theta[j - 1] * theta[j + k - 1]
            acf[1] += theta[k - 1] ** 2

    if p == 0:
        return acf, cvli, alpha, 0

    # [2]: initialise cvli with AR coefficients
    for k in range(1, p + 1):
        alpha[k] = phi[k - 1]
        cvli[k]  = phi[k - 1]

    # [3]: T.W.-S ALPHA and DELTA
    for k in range(1, mxpq + 1):
        kc = mxpq - k
        if kc < p:
            divv = 1.0 - alpha[kc + 1] ** 2
            if divv <= epsil2:
                return acf, cvli, alpha, 1
            if kc != 0:
                for j in range(1, kc + 1):
                    alpha[j] = (cvli[j] + alpha[kc + 1] * cvli[kc + 1 - j]) / divv
        if kc < q:
            j1 = max(1, kc + 1 - p)
            for j in range(j1, kc + 1):
                acf[j + 1] += acf[kc + 2] * alpha[kc + 1 - j]
        if kc < p:
            for j in range(1, kc + 1):
                cvli[j] = alpha[j]

    # [4]: NU
    acf[1] *= 0.5
    for k in range(1, mxpq + 1):
        if k <= p:
            divv = 1.0 - alpha[k] ** 2
            for j in range(1, k + 2):
                cvli[j] = (acf[j] + alpha[k] * acf[k + 2 - j]) / divv
            for j in range(1, k + 2):
                acf[j] = cvli[j]

    # [5]: full ACF
    for i in range(1, ma + 1):
        for j in range(1, min(i - 1, p) + 1):
            acf[i] += phi[j - 1] * acf[i - j]
    acf[1] *= 2.0

    # [6]: cvli
    cvli[1] = 1.0
    if q <= 0:
        return acf, cvli, alpha, 0
    for k in range(1, q + 1):
        cvli[k + 1] = -theta[k - 1]
        if p != 0:
            for j in range(1, min(k, p) + 1):
                cvli[k + 1] += phi[j - 1] * cvli[k + 1 - j]

    return acf, cvli, alpha, 0


def flikam_scalar(n, p, q, phi, theta, mu, w, xitol=1e-3, do_chkma=True,
                  compute_residuals=False):
    """Mélard (1984) scalar ARMA log-likelihood.

    Direct Python port of flikam() in usmelard.c (m=1 specialisation).

    Parameters
    ----------
    n, p, q   : obs count, AR/MA orders
    phi       : ndarray (p,)
    theta     : ndarray (q,)
    mu        : float  mean (0 after subtract in cast_us)
    w         : ndarray (n,)  mean-subtracted differenced series
    xitol     : switching tolerance (default 1e-3)
    do_chkma  : check MA invertibility
    compute_residuals : return innovations in `at`

    Returns
    -------
    sumsq  : float  Σ(innovation²/variance)
    fact   : float  determinant factor
    loglik : float  log-likelihood (sigma2=1)
    at     : ndarray (n,)  innovations (residuals)
    ifault : int  0=OK, 1=singular, 2=unit root, 3=non-stationary, 4=non-invertible
    """
    phi   = np.asarray(phi,   dtype=float)
    theta = np.asarray(theta, dtype=float)
    w     = np.asarray(w,     dtype=float)
    at    = np.zeros(n)

    sumsq = fact = loglik = 0.0
    epsil1 = 1e-10

    mxpq   = max(p, q)
    mxpqp1 = mxpq + 1
    mqp1   = q + 1
    mpp1   = p + 1
    mr     = max(p, q + 1)
    mrp1   = mr + 1

    # [1]: MA invertibility check
    if q > 0 and do_chkma:
        if chekma_scalar(theta):
            return sumsq, fact, loglik, at, 4

    # [2]: ACF and innovation covariance sequence
    acf, cvli, alpha, iflt = _twacf_scalar(p, q, phi, theta)
    if iflt:
        return sumsq, fact, loglik, at, 2

    # [3]: first column of P matrix
    vk = np.zeros(mrp1 + 1)
    vk[1] = acf[1]
    for k in range(2, mr + 1):
        s = 0.0
        if k <= p:
            for j in range(k, p + 1):
                s += phi[j - 1] * acf[j + 2 - k]
        if k <= mqp1:
            for j in range(k, mqp1 + 1):
                s -= theta[j - 2] * cvli[j + 1 - k]   # q[j-1][1][1] = theta[j-1]
        vk[k] = s

    # [4]: initial L and K vectors
    vw = np.zeros(mrp1 + 1)
    vl = np.zeros(mrp1 + 1)
    r_var = vk[1]
    vl[mr] = 0.0
    for j in range(1, mr + 1):
        vw[j] = 0.0
        vl[j] = vk[j + 1] if j != mr else 0.0
        if j <= p:
            vl[j] += phi[j - 1] * r_var
        vk[j] = vl[j]

    # [5]: main time loop
    last  = mpp1 - q
    loop  = p
    jfrom = mpp1
    vw[mpp1] = 0.0
    vl[mxpqp1] = 0.0
    nexti = n + 1

    e = np.zeros(n + 1)
    detman = 1.0
    detcar = 0.0

    for i in range(1, n + 1):
        if i == last:
            loop  = min(p, q)
            jfrom = loop + 1
            if q <= 0:
                nexti = i
                break

        if r_var <= epsil1:
            return sumsq, fact, loglik, at, 3

        if abs(r_var - 1.0) < xitol and i > mxpq:
            nexti = i
            break

        # Update determinant
        detman *= r_var
        while abs(detman) >= 1.0:
            detman *= 0.0625; detcar += 4.0
        while abs(detman) < 0.0625:
            detman *= 16.0;  detcar -= 4.0

        vw1 = vw[1]
        a_t = w[i - 1] - mu - vw1
        at[i - 1] = a_t
        e[i]      = a_t / np.sqrt(r_var)
        aor       = a_t / r_var
        sumsq    += a_t * aor
        vl1       = vl[1]
        alf       = vl1 / r_var
        r_var    -= alf * vl1

        # Update vectors
        if loop != 0:
            for j in range(1, loop + 1):
                flj   = vl[j + 1] + phi[j - 1] * vl1
                vw[j] = vw[j + 1] + phi[j - 1] * vw1 + aor * vk[j]
                vl[j] = flj - alf * vk[j]
                vk[j] -= alf * flj
        if jfrom <= q:
            for j in range(jfrom, q + 1):
                vw[j]  = vw[j + 1] + aor * vk[j]
                vl[j]  = vl[j + 1] - alf * vk[j]
                vk[j] -= alf * vl[j + 1]
        if jfrom <= p:
            for j in range(jfrom, p + 1):
                vw[j] = vw[j + 1] + phi[j - 1] * (w[i - 1] - mu)

    # [6]: quick recursions for remaining observations
    if nexti <= n:
        for i in range(nexti, n + 1):
            e[i] = w[i - 1] - mu
        if p != 0:
            for i in range(nexti, n + 1):
                for j in range(1, p + 1):
                    e[i] -= phi[j - 1] * (w[i - j - 1] - mu) if i - j >= 1 else 0.0
        if q != 0:
            for i in range(nexti, n + 1):
                for j in range(1, q + 1):
                    if i - j >= 1:
                        e[i] += theta[j - 1] * e[i - j]
        for i in range(nexti, n + 1):
            sumsq += e[i] ** 2
        if compute_residuals:
            for i in range(nexti, n + 1):
                at[i - 1] = e[i]

    # [7]: determinant factor and log-likelihood
    det_val = detman * (2.0 ** detcar)
    fact    = float(np.exp(np.log(det_val) / n))
    loglik  = float(-0.5 * (n * _LOG2PI + np.log(fact) + sumsq))

    return sumsq, fact, loglik, at, 0
