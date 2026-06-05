"""matplotlib-based visualization for FUE models — Treadway-Jenkins design."""

import numpy as np


# ── Public plot functions ──────────────────────────────────────────────────────

def plot_series(series, title=None, ax=None):
    """Raw series plot."""
    import matplotlib.pyplot as plt
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.get_figure()
    xs = _obs_to_decimal_year(series.nobs, *series.start, series.freq)
    _tj_spines(ax)
    ax.plot(xs, series.data, color='k', linewidth=0.9)
    ax.set_title(title or series.name, fontweight='bold', fontsize=11)
    ax.tick_params(direction='out', labelsize=9)
    if standalone:
        fig.tight_layout()
        plt.show()
    return fig


def plot_acf(data, lags=24, title="ACF", confidence=0.95, ax=None):
    """Single ACF panel (legacy helper)."""
    from .diagnostics import acf as _acf
    import matplotlib.pyplot as plt
    r = np.asarray(data, dtype=float)
    rc = _acf(r, lags=lags)
    band = _ci_bound(len(r), confidence)
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6, 3))
    else:
        fig = ax.get_figure()
    _stem_plot(ax, range(1, lags + 1), rc, band, title)
    if standalone:
        fig.tight_layout()
        plt.show()
    return fig


def plot_pacf(data, lags=24, title="PACF", confidence=0.95, ax=None):
    """Single PACF panel (legacy helper)."""
    from .diagnostics import pacf as _pacf
    import matplotlib.pyplot as plt
    r = np.asarray(data, dtype=float)
    pc = _pacf(r, lags=lags)
    band = _ci_bound(len(r), confidence)
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6, 3))
    else:
        fig = ax.get_figure()
    _stem_plot(ax, range(1, lags + 1), pc, band, title)
    if standalone:
        fig.tight_layout()
        plt.show()
    return fig


def plot_residuals_ts(residuals, model=None, title="", ax=None):
    """Standardized residuals time series — Treadway-Jenkins style.

    Draws linespoints against a decimal-year x-axis, ±2σ dashed reference
    lines, and seasonal vertical dividers (for freq > 1).  Xlabel shows the
    sample mean, its standard error, and the sample std dev.
    """
    import matplotlib.pyplot as plt

    r = np.asarray(residuals, dtype=float)
    n = len(r)
    rmean = r.mean()
    rstd  = r.std(ddof=0)
    z = (r - rmean) / rstd if rstd > 0 else r.copy()
    abs_max = max(float(np.abs(z).max()), 2.5)

    freq = 1
    start_year, start_period = 1, 1
    refactor = 100.0
    if model is not None:
        freq = model.series.freq
        start_year, start_period = model.series.start
        refactor = float(model.refactor) if model.refactor != 0 else 100.0

    xs = _obs_to_decimal_year(n, start_year, start_period, freq)

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.get_figure()

    _tj_spines(ax)
    ax.plot(xs, z, color='k', linewidth=0.9,
            marker='o', markersize=3, markerfacecolor='k',
            markeredgewidth=0, zorder=3)
    ax.axhline(0,  color='k',   lw=0.8, linestyle='-',  zorder=2)
    ax.axhline( 2, color='0.3', lw=1.2, linestyle='--', zorder=2)
    ax.axhline(-2, color='0.3', lw=1.2, linestyle='--', zorder=2)

    # seasonal dividers at year boundaries (every 2 years for long series)
    if freq > 1:
        x0, x1 = xs[0], xs[-1]
        step = 2 if (x1 - x0) > 5 else 1
        first_yr = int(np.ceil(x0 - 1e-9))
        for yr in range(first_yr, int(x1) + 2, step):
            if x0 < yr <= x1 + 1.0 / freq:
                ax.axvline(yr, color='k', lw=0.5, zorder=1)

    y_max = int(np.ceil(abs_max / 2.0)) * 2
    ax.set_ylim(-y_max - 0.15, y_max + 0.15)
    ax.set_yticks(range(-y_max, y_max + 1, 2))
    ax.tick_params(axis='y', direction='out', labelsize=9)
    ax.tick_params(axis='x', direction='out', labelsize=8)
    x_pad = 0.3 / freq
    ax.set_xlim(xs[0] - x_pad, xs[-1] + x_pad)

    se = rstd / np.sqrt(n)
    scale = refactor
    ax.set_xlabel(
        f"$\\bar{{w}}$ = {rmean * scale:.2f}  ({se * scale:.2f})"
        f"    $\\hat{{\\sigma}}_w$ = {rstd * scale:.2f}",
        fontsize=9,
    )
    ax.set_title(title, fontweight='bold', fontsize=12)

    if standalone:
        fig.tight_layout()
        plt.show()
    return fig


def plot_acf_pacf(residuals, npar=0, freq=1, lags=None, title="",
                  ax_acf=None, ax_pacf=None):
    """Stacked ACF and PACF panels — impulse style, Ljung-Box Q in ACF xlabel.

    Both panels share the same y-axis range.  Confidence bands are drawn at
    ±2/√n.  Seasonal grid lines mark multiples of freq.
    """
    import matplotlib.pyplot as plt
    from .diagnostics import acf as _acf, pacf as _pacf, ljung_box as _lb

    r = np.asarray(residuals, dtype=float)
    n = len(r)

    if lags is None:
        lags = _default_lags(n, freq)
    lags = max(int(lags), 1)

    rc = _acf(r, lags=lags)
    pc = _pacf(r, lags=lags)

    band = 2.0 / np.sqrt(n)
    raw_max = max(float(np.abs(rc).max()), float(np.abs(pc).max()))
    cmax = _round_cmax(max(raw_max, band * 1.2))

    standalone = ax_acf is None
    if standalone:
        fig, (ax_acf, ax_pacf) = plt.subplots(
            2, 1, figsize=(5, 6),
            gridspec_kw={'hspace': 0.12},
            layout='constrained',
        )
    else:
        fig = ax_acf.get_figure()

    lag_x = np.arange(1, lags + 1)

    _draw_acf_panel(ax_acf,  lag_x, rc, band, cmax, freq, lags, 'acf')
    lb = _lb(r, lags=lags, df_correction=npar)
    lb_stat = lb["statistic"][0]
    ax_acf.set_xlabel(f"Q({lags - npar}) = {lb_stat:.1f}", fontsize=8)

    _draw_acf_panel(ax_pacf, lag_x, pc, band, cmax, freq, lags, 'pacf')
    ax_pacf.set_xlabel('')

    if title:
        ax_acf.set_title(f"acf — {title}", loc='left', fontsize=9, pad=2)

    if standalone:
        plt.show()
    return fig


def plot_histogram(residuals, title="", ax=None):
    """Histogram of standardized residuals with standard-normal overlay.

    Matches gnuplot design: 15% fill, bandwidth=0.5, ±4σ range (or ±8σ for
    extreme outliers), y-axis in %, S/K/JB stats in xlabel.
    """
    import matplotlib.pyplot as plt
    from scipy.stats import norm, skew as _skew, kurtosis as _kurt
    from .diagnostics import jarque_bera as _jb

    r = np.asarray(residuals, dtype=float)
    n = len(r)
    rmean = r.mean()
    rstd  = r.std(ddof=0)
    z = (r - rmean) / rstd if rstd > 0 else r.copy()

    sk = float(_skew(z, bias=False))
    ku = float(_kurt(z, fisher=True, bias=False))   # excess kurtosis
    jb_stat = float(_jb(r)[0])

    abs_max = float(np.abs(z).max())
    xmax = 4.0 if abs_max <= 4.0 else 8.0

    # bins: width=0.5, range [-xmax, xmax]
    bin_edges = np.arange(-xmax, xmax + 1e-9, 0.5)
    counts, _ = np.histogram(z, bins=bin_edges)
    centers   = bin_edges[:-1] + 0.25           # bin center
    prob      = 100.0 * 2.0 * counts / n        # % density (×2 because bw=0.5)

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 5))
    else:
        fig = ax.get_figure()

    _tj_spines(ax)
    ax.bar(centers, prob, width=0.5,
           color='k', alpha=0.15, edgecolor='k', linewidth=1.0, zorder=3)

    xs = np.linspace(-xmax, xmax, 500)
    ax.plot(xs, 100.0 * norm.pdf(xs, 0, 1), 'k-', lw=1.5, zorder=4)

    ax.set_xlim(-xmax, xmax)
    ax.set_xticks(np.arange(-xmax, xmax + 0.001, 2.0))
    ax.set_ylabel('%', fontsize=10)
    ax.yaxis.set_tick_params(direction='out', labelsize=9)
    ax.xaxis.set_tick_params(direction='out', labelsize=9)
    ax.grid(True, which='major', lw=0.4, alpha=0.6, zorder=0)
    ax.set_xlabel(f"S = {sk:.1f}     K = {ku:.1f}     JB = {jb_stat:.1f}", fontsize=10)
    ax.set_title(title, fontweight='bold', fontsize=11)

    if standalone:
        fig.tight_layout()
        plt.show()
    return fig


def plot_model_diagnostics(model, lags=None, save_prefix=None):
    """Full Treadway-Jenkins diagnostic panel.

    Produces two figures:
      1. Residuals time series (left) + stacked ACF/PACF (right)
      2. Histogram with normal overlay

    Parameters
    ----------
    model : Model  (must be fitted)
    lags : int, optional
    save_prefix : str, optional
        If given, saves <save_prefix>_diag.png and <save_prefix>_hist.png

    Returns
    -------
    (fig_diag, fig_hist)
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    model._require_fit()
    r     = model._result.residuals
    npar  = model._result.npar
    freq  = model.series.freq
    title = model.series.name

    if lags is None:
        lags = _default_lags(len(r), freq)

    # Figure 1: residuals series (left) + stacked ACF/PACF (right)
    fig1 = plt.figure(figsize=(13, 5), layout='constrained')
    gs = gridspec.GridSpec(2, 2, figure=fig1,
                           width_ratios=[1.65, 1.0],
                           hspace=0.12, wspace=0.3)
    ax_ser  = fig1.add_subplot(gs[:, 0])
    ax_acf  = fig1.add_subplot(gs[0, 1])
    ax_pacf = fig1.add_subplot(gs[1, 1])

    plot_residuals_ts(r, model=model, title=title, ax=ax_ser)
    plot_acf_pacf(r, npar=npar, freq=freq, lags=lags,
                  ax_acf=ax_acf, ax_pacf=ax_pacf)

    # Figure 2: histogram
    fig2, ax_h = plt.subplots(figsize=(5, 5))
    plot_histogram(r, title=title, ax=ax_h)
    fig2.tight_layout()

    if save_prefix:
        fig1.savefig(f"{save_prefix}_diag.png", dpi=150, bbox_inches='tight')
        fig2.savefig(f"{save_prefix}_hist.png", dpi=150, bbox_inches='tight')

    plt.show()
    return fig1, fig2


# keep old name as alias for backward compat (model.py still calls it)
def plot_residual_diagnostics(residuals, lags=24, title=""):
    """Legacy wrapper — prefer plot_model_diagnostics(model)."""
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    r = np.asarray(residuals, dtype=float)
    freq = 1
    lags = int(lags)

    fig1 = plt.figure(figsize=(13, 5), layout='constrained')
    gs = gridspec.GridSpec(2, 2, figure=fig1,
                           width_ratios=[1.65, 1.0],
                           hspace=0.12, wspace=0.3)
    ax_ser  = fig1.add_subplot(gs[:, 0])
    ax_acf  = fig1.add_subplot(gs[0, 1])
    ax_pacf = fig1.add_subplot(gs[1, 1])

    plot_residuals_ts(r, title=title, ax=ax_ser)
    plot_acf_pacf(r, freq=freq, lags=lags, ax_acf=ax_acf, ax_pacf=ax_pacf)

    fig2, ax_h = plt.subplots(figsize=(5, 5))
    plot_histogram(r, title=title, ax=ax_h)
    fig2.tight_layout()

    plt.show()


# ── Internal helpers ───────────────────────────────────────────────────────────

def _tj_spines(ax, sides=('left', 'bottom')):
    """Remove top/right spines — Treadway-Jenkins minimal border."""
    for sp in ('top', 'right', 'left', 'bottom'):
        ax.spines[sp].set_visible(sp in sides)


def _obs_to_decimal_year(n, start_year, start_period, freq):
    """Decimal-year x-coordinates for n observations."""
    result = []
    y, p = int(start_year), int(start_period)
    for _ in range(n):
        result.append(y + (p - 1) / freq)
        p += 1
        if p > freq:
            p = 1
            y += 1
    return result


def _default_lags(n, freq):
    if n < 3 * (freq + 1):
        return max(n - freq // 2, 1)
    if freq == 1:
        return 9
    return 3 * (freq + 1)


def _round_cmax(val):
    """Round up to nearest 0.1 for a clean ACF y-axis range."""
    return max(round(float(np.ceil(val * 10)) / 10.0, 1), 0.3)


def _draw_acf_panel(ax, lag_x, vals, band, cmax, freq, lags, label):
    """One ACF or PACF panel: impulse style, confidence bands, seasonal grid."""
    _tj_spines(ax, sides=('left',))
    ax.vlines(lag_x, 0, vals, colors='k', linewidth=2.5, zorder=3)
    ax.axhline( band, color='k', lw=1.0, linestyle='--', zorder=2)
    ax.axhline(-band, color='k', lw=1.0, linestyle='--', zorder=2)
    ax.axhline(0,     color='k', lw=1.5, zorder=2)

    n_lags = int(lag_x[-1])
    if freq > 1:
        grid_lags = [freq * m for m in range(1, 4) if freq * m <= n_lags]
    elif n_lags > 9:
        gap = n_lags // 3
        grid_lags = [gap * m for m in range(1, 4) if gap * m <= n_lags]
    else:
        grid_lags = [x for x in (3, 6, 9) if x <= n_lags]

    for xv in grid_lags:
        ax.axvline(xv, color='0.5', lw=0.8, zorder=1)

    ax.set_ylim(-cmax, cmax)
    ax.set_yticks([-cmax, 0.0, cmax])
    ax.yaxis.set_tick_params(direction='out', labelsize=8)
    ax.set_xticks(grid_lags)
    ax.set_xticklabels([str(x) for x in grid_lags], fontsize=7)
    ax.tick_params(axis='x', direction='out', length=3)
    ax.set_xlim(0.5, n_lags + 0.5)
    ax.set_title(label, loc='left', fontsize=9, pad=2)


def _ci_bound(n, confidence):
    from scipy.stats import norm
    return norm.ppf((1 + confidence) / 2) / np.sqrt(n)


def _stem_plot(ax, lags, values, bound, title):
    ax.bar(lags, values, width=0.3, color='steelblue')
    ax.axhline( bound, color='r', linestyle='--', linewidth=0.8)
    ax.axhline(-bound, color='r', linestyle='--', linewidth=0.8)
    ax.axhline(0,      color='k', linewidth=0.5)
    ax.set_title(title)
