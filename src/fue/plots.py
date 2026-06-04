"""matplotlib-based visualization for FUE models and series."""

import numpy as np


def plot_series(series, title=None, ax=None):
    import matplotlib.pyplot as plt
    fig, ax = _get_ax(ax)
    xs = _date_axis(series)
    ax.plot(xs, series.data)
    ax.set_title(title or series.name)
    if series.freq > 1:
        fig.autofmt_xdate()
    _show(fig, ax)


def _date_axis(series):
    """Return a list of x-axis tick labels (decimal year or observation index)."""
    freq = series.freq if series.freq > 0 else 1
    year, period = series.start
    result = []
    y, p = year, period
    for _ in range(series.nobs):
        if freq == 1:
            result.append(y)
        else:
            result.append(y + (p - 1) / freq)
        p += 1
        if p > freq:
            p = 1
            y += 1
    return result


def plot_acf(data, lags=24, title="ACF", confidence=0.95, ax=None):
    from .diagnostics import acf as _acf
    import matplotlib.pyplot as plt
    r   = _acf(data, lags=lags)
    n   = len(data)
    bound = _ci_bound(n, confidence)
    fig, ax = _get_ax(ax)
    _stem_plot(ax, range(1, lags + 1), r, bound, title)
    _show(fig, ax)


def plot_pacf(data, lags=24, title="PACF", confidence=0.95, ax=None):
    from .diagnostics import pacf as _pacf
    import matplotlib.pyplot as plt
    p   = _pacf(data, lags=lags)
    n   = len(data)
    bound = _ci_bound(n, confidence)
    fig, ax = _get_ax(ax)
    _stem_plot(ax, range(1, lags + 1), p, bound, title)
    _show(fig, ax)


def plot_residual_diagnostics(residuals, lags=24, title=""):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f"Residual diagnostics — {title}")

    r = np.asarray(residuals)

    # Time plot
    axes[0, 0].plot(r)
    axes[0, 0].axhline(0, color="k", linewidth=0.8)
    axes[0, 0].set_title("Residuals")

    # Histogram
    axes[0, 1].hist(r, bins="auto", density=True, alpha=0.7)
    from scipy.stats import norm
    xs = np.linspace(r.min(), r.max(), 200)
    axes[0, 1].plot(xs, norm.pdf(xs, r.mean(), r.std()), "r-", lw=1.5)
    axes[0, 1].set_title("Histogram")

    # ACF
    from .diagnostics import acf as _acf
    rc = _acf(r, lags=lags)
    bound = _ci_bound(len(r), 0.95)
    _stem_plot(axes[1, 0], range(1, lags + 1), rc, bound, "ACF of residuals")

    # PACF
    from .diagnostics import pacf as _pacf
    pc = _pacf(r, lags=lags)
    _stem_plot(axes[1, 1], range(1, lags + 1), pc, bound, "PACF of residuals")

    plt.tight_layout()
    plt.show()


# ── Internal helpers ──────────────────────────────────────────────────────

def _ci_bound(n, confidence):
    from scipy.stats import norm
    return norm.ppf((1 + confidence) / 2) / np.sqrt(n)


def _stem_plot(ax, lags, values, bound, title):
    ax.bar(lags, values, width=0.3, color="steelblue")
    ax.axhline( bound, color="r", linestyle="--", linewidth=0.8)
    ax.axhline(-bound, color="r", linestyle="--", linewidth=0.8)
    ax.axhline(0,      color="k", linewidth=0.5)
    ax.set_title(title)


def _get_ax(ax):
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots()
        return fig, ax
    return ax.get_figure(), ax


def _show(fig, ax):
    import matplotlib.pyplot as plt
    if ax.get_figure() == fig:
        plt.show()
