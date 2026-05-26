"""Diagnostic statistics: ACF, PACF, Jarque-Bera, Ljung-Box."""

import numpy as np


def acf(data, lags=24):
    """
    Sample autocorrelation function.

    Returns array of length *lags* with r[k] = Corr(x_t, x_{t-k}).
    """
    x   = np.asarray(data, dtype=float)
    n   = len(x)
    mu  = x.mean()
    c0  = np.dot(x - mu, x - mu) / n
    out = np.empty(lags)
    for k in range(1, lags + 1):
        out[k - 1] = np.dot(x[k:] - mu, x[:-k] - mu) / (n * c0)
    return out


def pacf(data, lags=24):
    """
    Partial autocorrelation function via Durbin-Levinson recursion.

    Returns array of length *lags*.
    """
    r   = acf(data, lags=lags)
    phi = np.zeros((lags + 1, lags + 1))
    p   = np.empty(lags)
    phi[1, 1] = r[0]
    p[0] = r[0]
    for k in range(2, lags + 1):
        num = r[k - 1] - sum(phi[k - 1, j] * r[k - 1 - j] for j in range(1, k))
        den = 1.0  - sum(phi[k - 1, j] * r[j - 1]         for j in range(1, k))
        phi[k, k] = num / den if abs(den) > 1e-14 else 0.0
        for j in range(1, k):
            phi[k, j] = phi[k - 1, j] - phi[k, k] * phi[k - 1, k - j]
        p[k - 1] = phi[k, k]
    return p


def jarque_bera(data):
    """
    Jarque-Bera normality test.

    Returns (statistic, p-value).
    """
    from scipy import stats as st
    return st.jarque_bera(data)


def ljung_box(data, lags=None, df_correction=0):
    """
    Ljung-Box portmanteau test.

    Parameters
    ----------
    data : array-like
        Residuals.
    lags : int or list of int
        Lag(s) at which to compute the test.  Defaults to min(10, nobs//5).
    df_correction : int
        Number of estimated ARMA parameters (subtracted from degrees of freedom).

    Returns
    -------
    dict with keys 'statistic', 'pvalue', 'lags'.
    """
    from scipy import stats as st
    x   = np.asarray(data, dtype=float)
    n   = len(x)
    r   = acf(x, lags=(max(lags) if hasattr(lags, '__iter__') else (lags or min(10, n // 5))))
    ks  = ([lags] if isinstance(lags, int) else
           list(lags) if lags is not None else [min(10, n // 5)])
    stats_, pvals = [], []
    for k in ks:
        s = n * (n + 2) * sum(r[j] ** 2 / (n - j - 1) for j in range(k))
        df = max(1, k - df_correction)
        stats_.append(s)
        pvals.append(1 - st.chi2.cdf(s, df))
    return {"statistic": stats_, "pvalue": pvals, "lags": ks}
