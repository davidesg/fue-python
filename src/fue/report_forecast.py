"""
HTML forecast report generator — SPS (Sistema de Previsión y Seguimiento).

Layout
------
LEFT column  : table (hist + first-year forecast + H=L row)  + model details
RIGHT column : single SVG with two panels (mirrors C gnuplot multiplot)
                  top    — annual rate of change (history + full horizon)
                  bottom — ERR (x-axis truncated at forecast origin,
                               aligned with historical half of top panel)

Both panels share the same figure → same column widths → perfect x-alignment.
ERR spine ends at the forecast origin; impulses stop there naturally.

Public API
----------
write_forecast_report(model, fr, path,
                      title=None, source=None, sps_name=None,
                      narrative=None, pdf=False)
"""

import io
import math
from datetime import date as _date

import numpy as np


# ── HTML/CSS template (Jinja2) ────────────────────────────────────────────────

_HTML = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{{ page_title }}</title>
  <style>
    :root {
      --font: system-ui, -apple-system, "Segoe UI", Helvetica, sans-serif;
      --mono: ui-monospace, "Cascadia Code", "Fira Code", monospace;
      --fg:    #111827;
      --muted: #6b7280;
      --border:#e5e7eb;
      --fore-bg:#eff6ff;
      --accent: #1d4ed8;
      --max-w:  1200px;
    }
    *  { box-sizing: border-box; margin: 0; padding: 0 }
    body {
      font-family: var(--font);
      font-size: 14px;
      color: var(--fg);
      background: #fff;
      padding: 2rem 1.5rem;
    }
    .container { max-width: var(--max-w); margin: 0 auto }

    /* ── Header ── */
    header {
      border-bottom: 2.5px solid var(--fg);
      padding-bottom: .9rem;
      margin-bottom: 1.8rem;
    }
    .sps-label {
      font-size: .72rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: .1em;
      color: var(--accent);
      margin-bottom: .3rem;
    }
    header h1 { font-size: 1.55rem; font-weight: 700; line-height: 1.2 }
    .meta {
      display: flex;
      flex-wrap: wrap;
      gap: .25rem 1.6rem;
      margin-top: .55rem;
      font-size: .81rem;
      color: var(--muted);
    }
    .meta strong { color: var(--fg) }

    /* ── Two-column body ─────────────────────────────────────────────────
       LEFT:  table + model details
       RIGHT: forecast chart (top) + ERR chart (bottom) — single SVG
    ─────────────────────────────────────────────────────────────────── */
    .report-body {
      display: grid;
      grid-template-columns: minmax(0, 1fr) minmax(0, 1.2fr);
      gap: 2rem 2.5rem;
      align-items: start;
    }
    .col-left  {}
    .col-right {}
    .col-right figure { margin: 0 }
    .col-right figure svg { width: 100%; height: auto; display: block }
    .col-right figcaption {
      font-size: .73rem;
      color: var(--muted);
      margin-top: .3rem;
      text-align: center;
      font-style: italic;
    }

    /* ── Data table ── */
    .table-wrap { overflow-x: auto }
    .data-table {
      width: 100%;
      border-collapse: collapse;
      font-size: .77rem;
    }
    .data-table caption {
      font-size: .72rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: .07em;
      color: var(--muted);
      text-align: left;
      padding-bottom: .5rem;
    }
    .data-table thead tr:first-child th {
      padding: .3rem .45rem;
      border-bottom: 1px solid var(--border);
      font-size: .65rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: .06em;
      color: var(--muted);
      text-align: center;
    }
    .data-table thead tr:first-child th:first-child { text-align: left }
    .data-table thead tr:last-child th {
      padding: .28rem .45rem;
      border-bottom: 2px solid var(--fg);
      font-size: .68rem;
      font-weight: 600;
      text-align: right;
      white-space: nowrap;
    }
    .data-table thead tr:last-child th:first-child { text-align: left }
    .data-table td {
      padding: .22rem .45rem;
      border-bottom: 1px solid var(--border);
      text-align: right;
      font-family: var(--mono);
      font-size: .74rem;
      font-variant-numeric: tabular-nums;
      white-space: nowrap;
    }
    .data-table td:first-child { text-align: left; font-family: var(--font) }
    .data-table tr.fore td     { background: var(--fore-bg) }
    .data-table tr.sep  td     { border-top: 2px solid var(--fg) }
    .data-table tr.blank td    {
      padding-top: .04rem; padding-bottom: .04rem;
      border-bottom: none;
      background: #fff;
      height: 6px;
    }
    .data-table .na            { color: var(--muted) }

    /* ── Model details ── */
    .model-details { margin-top: 1.2rem }
    details > summary {
      cursor: pointer;
      font-size: .72rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: .08em;
      color: var(--muted);
      padding: .4rem 0;
      border-top: 1px solid var(--border);
      list-style: none;
      user-select: none;
    }
    details > summary::after       { content: " ▸"; font-size: .7rem }
    details[open] > summary::after { content: " ▾" }
    .model-grid {
      display: flex; flex-wrap: wrap;
      gap: .3rem 2rem; padding: .7rem 0; font-size: .79rem;
    }
    .model-grid span   { color: var(--muted) }
    .model-grid strong { font-family: var(--mono) }

    /* ── Narrative ── */
    .narrative {
      margin-top: 2rem; border-top: 1px solid var(--border);
      padding-top: 1.2rem; line-height: 1.7;
    }
    .narrative h2 {
      font-size: .72rem; font-weight: 600; text-transform: uppercase;
      letter-spacing: .08em; color: var(--muted); margin-bottom: .8rem;
    }
    .narrative p { margin-bottom: .7rem }

    /* ── Footer ── */
    footer {
      margin-top: 2.5rem; padding-top: .7rem;
      border-top: 1px solid var(--border);
      font-size: .71rem; color: var(--muted);
    }

    @media print {
      body { padding: 0; font-size: 10pt }
      .report-body { grid-template-columns: 1fr 1.2fr }
      .data-table  { font-size: 7.5pt }
      details > summary { display: none }
      details { display: none }
    }
    @media (max-width: 860px) {
      .report-body { grid-template-columns: 1fr }
    }
  </style>
</head>
<body>
<div class="container">

<header>
  <div class="sps-label">SPS{% if sps_name %}: {{ sps_name }}{% endif %}</div>
  <h1>{{ page_title }}</h1>
  <div class="meta">
    {% if source %}<span>Data Source&#160;<strong>{{ source }}</strong></span>{% endif %}
    <span>Forecast Origin&#160;<strong>{{ origin }}</strong></span>
    <span>Horizon&#160;<strong>{{ horizon }}&#160;{{ freq_label }}</strong></span>
    <span>Generated&#160;<strong>{{ generated }}</strong></span>
  </div>
</header>

<div class="report-body">

  <!-- LEFT: table + model details -->
  <div class="col-left">
    <div class="table-wrap">
    <table class="data-table">
      <caption>Forecast</caption>
      <thead>
        <tr>
          <th rowspan="2">Date</th>
          <th colspan="2">Level</th>
          <th colspan="2">{{ diff1_label }}</th>
          <th colspan="2">Annual (%)</th>
          <th rowspan="2">ERR<br>(%)</th>
        </tr>
        <tr>
          <th>Value</th><th>Std&#160;(%)</th>
          <th>(%)</th><th>Std&#160;(%)</th>
          <th>(%)</th><th>Std&#160;(%)</th>
        </tr>
      </thead>
      <tbody>
        {% for row in hist_rows %}
        <tr>
          <td>{{ row.date }}</td>
          <td>{{ row.level }}</td>
          <td class="na">—</td>
          <td>{{ row.diff1 }}</td>
          <td class="na">—</td>
          <td>{{ row.annual }}</td>
          <td class="na">—</td>
          <td>{{ row.err }}</td>
        </tr>
        {% endfor %}
        {% for row in fore_rows_main %}
        <tr class="fore{% if loop.first %} sep{% endif %}">
          <td>{{ row.date }}</td>
          <td>{{ row.level }}</td>
          <td>{{ row.level_std }}</td>
          <td>{{ row.diff1 }}</td>
          <td>{{ row.diff1_std }}</td>
          <td>{{ row.annual }}</td>
          <td>{{ row.annual_std }}</td>
          <td class="na">—</td>
        </tr>
        {% endfor %}
        {% if fore_row_end %}
        <tr class="fore blank"><td colspan="8"></td></tr>
        <tr class="fore">
          <td>{{ fore_row_end.date }}</td>
          <td>{{ fore_row_end.level }}</td>
          <td>{{ fore_row_end.level_std }}</td>
          <td>{{ fore_row_end.diff1 }}</td>
          <td>{{ fore_row_end.diff1_std }}</td>
          <td>{{ fore_row_end.annual }}</td>
          <td>{{ fore_row_end.annual_std }}</td>
          <td class="na">—</td>
        </tr>
        {% endif %}
      </tbody>
    </table>
    </div>

    <div class="model-details">
      <details>
        <summary>Model details</summary>
        <div class="model-grid">
          <div><span>Model&#160;</span><strong>{{ stem }}</strong></div>
          <div><span>npar&#160;</span><strong>{{ npar }}</strong></div>
          <div><span>σ²&#160;</span><strong>{{ sigma2 }}</strong></div>
          <div><span>AIC&#160;</span><strong>{{ aic }}</strong></div>
          <div><span>BIC&#160;</span><strong>{{ bic }}</strong></div>
          <div><span>N&#160;</span><strong>{{ nobs }}</strong></div>
          <div><span>Sample&#160;</span><strong>{{ sample }}</strong></div>
        </div>
      </details>
    </div>
  </div><!-- .col-left -->

  <!-- RIGHT: forecast + ERR charts (single SVG, panels share x-axis) -->
  <div class="col-right">
    <figure>
      {{ charts_svg | safe }}
      <figcaption>Forecast bands ±1σ &nbsp;·&nbsp; ERR bands ±2σ</figcaption>
    </figure>
  </div>

</div><!-- .report-body -->

{% if narrative %}
<div class="narrative">
  <h2>Analysis</h2>
  {{ narrative | safe }}
</div>
{% endif %}

<footer>
  FUE Python {{ version }} · {{ generated }}
</footer>

</div>
</body>
</html>
"""


# ── Public API ────────────────────────────────────────────────────────────────

def write_forecast_report(model, fr, path,
                          title=None, source=None, sps_name=None,
                          narrative=None, pdf=False):
    """Write a self-contained HTML forecast report."""
    try:
        from jinja2 import Environment
    except ImportError:
        raise ImportError("write_forecast_report requires jinja2 — pip install jinja2")

    if model._result is None:
        raise RuntimeError(
            "model._result not set — call forecast_fuf() first"
        )

    hist_rows, fore_rows_main, fore_row_end, meta = _table_data(model, fr)
    charts_svg = _make_charts_svg(model, fr)

    stem       = meta["stem"]
    page_title = title or f"A.{stem}"

    r   = model._result
    ts  = model.series
    env = Environment(autoescape=False)
    html = env.from_string(_HTML).render(
        page_title     = page_title,
        sps_name       = sps_name or "",
        stem           = stem,
        source         = source or "",
        origin         = meta["origin"],
        horizon        = fr.horizon,
        freq_label     = meta["freq_label"],
        diff1_label    = meta["diff1_label"],
        sigma2         = f"{fr.sigma2:.6f}",
        aic            = f"{r.aic:.2f}"  if hasattr(r, "aic")  else "—",
        bic            = f"{r.bic:.2f}"  if hasattr(r, "bic")  else "—",
        npar           = r.npar          if hasattr(r, "npar")  else "—",
        nobs           = ts.nobs,
        sample         = meta["sample"],
        hist_rows      = hist_rows,
        fore_rows_main = fore_rows_main,
        fore_row_end   = fore_row_end,
        charts_svg     = charts_svg,
        narrative      = narrative,
        version        = "1.13 (Python port)",
        generated      = _date.today().isoformat(),
    )

    path = str(path)
    if not path.endswith(".html"):
        path += ".html"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)

    if pdf:
        _write_pdf(html, path.replace(".html", ".pdf"))

    return path


# ── Chart — combined multiplot ────────────────────────────────────────────────

def _make_charts_svg(model, fr) -> str:
    """Two-panel figure (mirrors C gnuplot multiplot).

    Top    : annual rate of change — history (dots+line) + forecast (line)
             + dashed ±1σ bands + vertical separator at forecast origin.
    Bottom : ERR residuals as impulses + ±2σ dashed bands.
             x-axis ends at the forecast origin (spine truncated); the
             impulse data naturally stops there.

    Both panels share the same figure and same x-range so their x-axes
    are pixel-perfectly aligned.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    from matplotlib.gridspec import GridSpec
    from .forecast import _boxcox
    from .plots    import _tj_spines

    ts       = model.series
    nobs     = ts.nobs
    freq     = ts.freq if ts.freq > 0 else 1
    L        = fr.horizon
    refactor = model.refactor
    boxlam   = model.boxlam
    raw      = ts.data

    # ── Shared x-setup ────────────────────────────────────────────────────────
    xtick_pos, xtick_lbl = _year_ticks(ts, nobs, L, freq)
    xlim = (-0.5, 2 * L - 0.5)

    # ── Figure with two panels ────────────────────────────────────────────────
    fig = plt.figure(figsize=(7, 7.5))
    gs  = GridSpec(2, 1, figure=fig,
                   height_ratios=[2.2, 1],
                   hspace=0.18,
                   left=0.10, right=0.97, top=0.96, bottom=0.07)
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1])

    # ── Top panel: annual rate of change ─────────────────────────────────────
    hist_annual = np.array([
        100.0 * (_boxcox(raw[nobs - L + i],         boxlam, refactor)
                 - _boxcox(raw[nobs - L + i - freq], boxlam, refactor)) / refactor
        for i in range(L)
    ])
    y_all  = np.concatenate([hist_annual, fr.seasonal_diff])
    x_all  = np.arange(2 * L)
    x_hist = np.arange(L)
    x_fore = np.arange(L, 2 * L)

    _tj_spines(ax_top, ('left', 'bottom'))
    # col1: linespoints for all 2L — dashed line + circles (gnuplot ls 6 lt 5, linespoints)
    ax_top.plot(x_all, y_all, color='k', ls='--', marker='o', lw=1.2, ms=4.0, zorder=3)
    # col4: larger circles for historical only (gnuplot: with points ls 1)
    ax_top.plot(x_hist, hist_annual, 'ko', ms=6.0, zorder=4)
    # bands (matches gnuplot ls 3: dashed, lw 1.5)
    ax_top.plot(x_fore, fr.seasonal_diff + fr.seasonal_diff_std, 'k--', lw=1.5, zorder=2)
    ax_top.plot(x_fore, fr.seasonal_diff - fr.seasonal_diff_std, 'k--', lw=1.5, zorder=2)
    ax_top.axvline(L - 0.5, color='0.55', lw=0.9, zorder=1)
    ax_top.axhline(0,        color='k',    lw=0.7, zorder=1)
    ax_top.set_xlim(*xlim)
    ax_top.set_xticks(xtick_pos)
    ax_top.set_xticklabels(xtick_lbl, fontsize=9)
    ax_top.tick_params(direction='out', labelsize=9)
    # grid at xtic positions (matches gnuplot: set grid xtics lt 1)
    ax_top.set_axisbelow(True)
    ax_top.grid(axis='x', color='0.75', lw=0.5, ls='-', zorder=0)
    ax_top.set_title('Annual rate of change (%)', loc='left', fontsize=9, pad=4)

    # ── Bottom panel: ERR ─────────────────────────────────────────────────────
    residuals  = np.asarray(model._result.residuals)
    err_L      = min(L, len(residuals))
    err_vals   = 100.0 * residuals[-err_L:] / refactor
    x_err      = np.arange(err_L)
    sigma_plot = math.sqrt(fr.sigma2)
    prevcmax   = _prevcmax(err_vals, sigma_plot)

    # x-ticks: only historical years (< err_L)
    hist_tick_pos = [p for p in xtick_pos if p < err_L]
    hist_tick_lbl = xtick_lbl[:len(hist_tick_pos)]

    x_end = err_L - 0.5   # data x-coordinate of the forecast origin

    _tj_spines(ax_bot, ('left', 'bottom'))
    ax_bot.vlines(x_err, 0, err_vals, colors='k', lw=1.6, zorder=3)
    # bands and zero-line end at the forecast origin, not at the plot edge
    ax_bot.hlines( 2 * sigma_plot, -0.5, x_end, colors='k', lw=1.0, ls='--', zorder=2)
    ax_bot.hlines(-2 * sigma_plot, -0.5, x_end, colors='k', lw=1.0, ls='--', zorder=2)
    ax_bot.hlines(0,               -0.5, x_end, colors='k', lw=1.2,          zorder=2)
    margin = 0.1 * sigma_plot
    ax_bot.set_ylim(-(prevcmax + margin), prevcmax + margin)
    yt = np.arange(0, prevcmax + 0.05 * sigma_plot, 2 * sigma_plot)
    ax_bot.set_yticks(np.concatenate([-yt[1:][::-1], yt]))
    ax_bot.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))

    # same x-range as top panel → pixel-perfect alignment
    ax_bot.set_xlim(*xlim)
    ax_bot.set_xticks(hist_tick_pos)
    ax_bot.set_xticklabels(hist_tick_lbl, fontsize=9)
    ax_bot.tick_params(direction='out', labelsize=9)
    # grid at historical xtic positions only (matches gnuplot: set grid xtics lt 1)
    ax_bot.set_axisbelow(True)
    ax_bot.grid(axis='x', color='0.75', lw=0.5, ls='-', zorder=0)

    # truncate the bottom spine at the last ERR point (forecast origin)
    ax_bot.spines['bottom'].set_bounds(-0.5, err_L - 0.5)

    ax_bot.set_title('ERR', loc='left', fontsize=9, pad=4)

    return _fig_to_svg(fig)


# ── Table data ────────────────────────────────────────────────────────────────

def _table_data(model, fr):
    """Return (hist_rows, fore_rows_main, fore_row_end, meta_dict).

    hist_rows      : freq+1 rows ending at the forecast origin
    fore_rows_main : first freq forecast rows (one year ahead)
    fore_row_end   : the H=horizon row; None if horizon <= freq
    """
    from .forecast import _boxcox
    from .report   import _fuf_obs_to_date_str

    ts       = model.series
    nobs     = ts.nobs
    freq     = ts.freq if ts.freq > 0 else 1
    L        = fr.horizon
    refactor = model.refactor
    boxlam   = model.boxlam
    raw      = ts.data
    begyear, begtime = ts.start

    residuals = np.asarray(model._result.residuals)
    ornsop    = nobs - len(residuals)
    n_hist    = freq + 1

    # ── Historical rows ───────────────────────────────────────────────────────
    hist_rows = []
    for k in range(nobs - n_hist + 1, nobs + 1):
        date  = _fuf_obs_to_date_str(k, begyear, begtime, freq)
        level = raw[k - 1]
        diff1 = (100.0 * (_boxcox(raw[k-1], boxlam, refactor)
                          - _boxcox(raw[k-2], boxlam, refactor)) / refactor
                 if k > 1 else 0.0)
        annual = (100.0 * (_boxcox(raw[k-1], boxlam, refactor)
                           - _boxcox(raw[k-1-freq], boxlam, refactor)) / refactor
                  if k > freq else
                  100.0 * _boxcox(raw[k-1], boxlam, refactor) / refactor)
        res_idx = k - ornsop - 1
        err_str = _f2(residuals[res_idx]) if 0 <= res_idx < len(residuals) else "—"
        hist_rows.append({
            "date":   date,
            "level":  _f2(level),
            "diff1":  _f2(diff1),
            "annual": _f2(annual),
            "err":    err_str,
        })

    # ── All forecast rows ─────────────────────────────────────────────────────
    all_fore = []
    for h in range(L):
        k    = nobs + h + 1
        date = _fuf_obs_to_date_str(k, begyear, begtime, freq)
        lstd = fr.level_std[h] * refactor
        all_fore.append({
            "date":       date,
            "level":      _f2(fr.level[h]),
            "level_std":  _f2(lstd),
            "diff1":      _f2(fr.diff1[h]),
            "diff1_std":  _f2(fr.diff1_std[h]),
            "annual":     _f2(fr.seasonal_diff[h]),
            "annual_std": _f2(fr.seasonal_diff_std[h]),
        })

    fore_rows_main = all_fore[:freq]
    fore_row_end   = all_fore[-1] if L > freq else None

    # ── Meta ──────────────────────────────────────────────────────────────────
    sample = (f"{_fuf_obs_to_date_str(1,    begyear, begtime, freq)} – "
              f"{_fuf_obs_to_date_str(nobs, begyear, begtime, freq)}")
    stem   = model._inp_stem or ts.name
    origin = _fuf_obs_to_date_str(nobs, begyear, begtime, freq)

    if freq == 12:
        freq_label, diff1_label = "months", "Monthly (%)"
    elif freq == 4:
        freq_label, diff1_label = "quarters", "Quarterly (%)"
    else:
        freq_label, diff1_label = "periods", "Period (%)"

    return hist_rows, fore_rows_main, fore_row_end, {
        "stem":        stem,
        "origin":      origin,
        "freq_label":  freq_label,
        "diff1_label": diff1_label,
        "sample":      sample,
    }


# ── Helpers ───────────────────────────────────────────────────────────────────

def _f2(v) -> str:
    return f"{v:.2f}"


def _fig_to_svg(fig) -> str:
    import matplotlib.pyplot as plt
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches=None)
    plt.close(fig)
    buf.seek(0)
    svg = buf.read()
    idx = svg.find("<svg")
    return svg[idx:] if idx >= 0 else svg


def _year_ticks(ts, nobs, L, freq):
    begyear, begtime = ts.start

    def obs_yr_per(obs1):
        total = begyear * freq + (begtime - 1) + (obs1 - 1)
        return int(total // freq), int(total % freq + 1)

    prevby = previndex = None
    for i in range(freq):
        yr, per = obs_yr_per(nobs - L + 1 + i)
        if per == 1:
            prevby, previndex = yr, i
            break
    if prevby is None:
        prevby, _ = obs_yr_per(nobs - L + 1)
        previndex = 0

    yr_step = 1 if freq == 12 else (2 if freq == 4 else 10)
    x_step  = freq * yr_step if freq > 1 else 10

    pos, lbl = [], []
    cur_yr, cur_x = prevby, previndex
    while cur_x < 2 * L:
        pos.append(cur_x)
        lbl.append(str(cur_yr))
        cur_yr += yr_step
        cur_x  += x_step
    return pos, lbl


def _prevcmax(err_vals, sigma_plot):
    prevcmax = 4.0 * sigma_plot
    for v in np.abs(err_vals):
        if v >= prevcmax:
            prevcmax = float(v)
    if   4.0 * sigma_plot < prevcmax <= 6.0 * sigma_plot:  prevcmax = 6.0 * sigma_plot
    elif 6.0 * sigma_plot < prevcmax <= 7.0 * sigma_plot:  prevcmax = 7.0 * sigma_plot
    elif prevcmax > 7.0 * sigma_plot:                       prevcmax = 10.0 * sigma_plot
    return prevcmax


def _write_pdf(html: str, path: str):
    try:
        from weasyprint import HTML as _W
    except ImportError:
        raise ImportError("PDF output requires weasyprint — pip install weasyprint")
    _W(string=html).write_pdf(path)
