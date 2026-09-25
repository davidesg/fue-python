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


# ── What the residuals of a MODEL need to be read correctly ──────────────────
#
# BUG-0023. Three things decide whether a residual figure or a Q statistic is
# right, and none of them is in the residuals themselves: how many ARMA
# parameters were ESTIMATED, how many observations the differencing consumed,
# and how many lags to look at. They lived in three places (report.py, plots.py
# and art) and the three disagreed. They live here now, once.

def default_lags(nobs, freq):
    """Default number of acf/pacf lags — the rule of fug C, exactly.

    fug-1.14-proto/src/diagnose.c, default_lags(): «same rule for the .out
    file, the plots and the GUI». At most nobs − 2 lags, at least 1.
    """
    nobs, freq = int(nobs), max(int(freq), 1)
    if nobs < 3 * (freq + 1):
        lags = nobs - freq // 2
    elif freq == 1 and nobs > 200:
        lags = 45
    elif freq == 1:
        lags = 9
    else:
        lags = 3 * (freq + 1)
    if lags > nobs - 2:
        lags = nobs - 2
    return max(lags, 1)


def free_arma_count(model):
    """Number of ESTIMATED ARMA parameters — the Ljung-Box df correction.

    Counts the free coefficients of the regular and seasonal AR/MA factors and
    of the fixed-frequency AR(2)/MA(2) factors. Deterministic regressors
    (harmonics, interventions, the mean) do NOT count: they are not estimated
    from the autocorrelation of the residuals. Fixed coefficients do not count
    either: an AR(1) fixed at 0 estimates nothing.
    """
    n = 0
    for factors, free_lists in (
        (getattr(model, "ar", None),   getattr(model, "ar_free", None)),
        (getattr(model, "ar_s", None), getattr(model, "ar_s_free", None)),
        (getattr(model, "ma", None),   getattr(model, "ma_free", None)),
        (getattr(model, "ma_s", None), getattr(model, "ma_s_free", None)),
    ):
        for i, factor in enumerate(factors or []):
            fl = (free_lists[i]
                  if free_lists is not None and i < len(free_lists)
                  and free_lists[i] is not None
                  else [True] * len(factor))
            n += sum(1 for f in fl if f)
    n += sum(1 for ff in (getattr(model, "ar_f", None) or []) if ff.free)
    n += sum(1 for ff in (getattr(model, "ma_f", None) or []) if ff.free)
    return n


def differencing_offset(model):
    """How many observations the differencing consumes before the 1st residual.

        d          each regular difference consumes 1
        D · s      each seasonal difference consumes s
        ifadf[f]   each stochastic seasonal ROOT consumes the degree of its
                   factor: 2 at an interior frequency (1 − 2cos(ω)B + B²),
                   1 at the Nyquist frequency (1 + B)

    The third term is the one that was missing wherever the count was written
    by hand as d + D·s (art BUG-0172, BUG-0185).
    """
    series = getattr(model, "series", None)
    s = int(getattr(series, "freq", 1) or 1)
    n = int(getattr(model, "d", 0) or 0) + int(getattr(model, "D", 0) or 0) * s
    for f, v in enumerate(list(getattr(model, "ifadf", None) or [])):
        if v == 1:
            n += 1 if (s >= 2 and f == s // 2) else 2
    return n


def residuals_start(model):
    """(year, period) of the FIRST residual: the series start plus the offset."""
    y0, p0 = model.series.start
    freq = max(int(getattr(model.series, "freq", 1) or 1), 1)
    off = (int(p0) - 1) + differencing_offset(model)
    return (int(y0) + off // freq, off % freq + 1)
