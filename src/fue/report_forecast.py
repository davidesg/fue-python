"""
HTML forecast report generator — SPS (Sistema de Previsión y Seguimiento).

Produces a self-contained .html file with:
  - Descriptive title (series name, not model identifier)
  - Two-column layout: charts LEFT, data table RIGHT
  - Embedded SVG charts (forecast panel + ERR panel)
  - Optional narrative section for LLM-generated analysis
  - Optional PDF via weasyprint

Public API
----------
write_forecast_report(model, fr, path,
                      title=None, source=None, narrative=None, pdf=False)
"""

import io
import math
from datetime import date as _date

import numpy as np


# ── HTML/CSS template (Jinja2) ────────────────────────────────────────────────

_HTML = """\
<!DOCTYPE html>
<html lang="es">
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
      --fore-bg:#f0f6ff;
      --accent: #1d4ed8;
      --max-w:  1160px;
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

    /* ── Two-column body: charts left, table right ── */
    .report-body {
      display: grid;
      grid-template-columns: minmax(0, 1.05fr) minmax(0, 0.95fr);
      gap: 2.5rem;
      align-items: start;
    }

    /* ── Charts column ── */
    .col-charts {}
    figure { margin: 0 0 1.6rem }
    figure:last-child { margin-bottom: 0 }
    figure svg { width: 100%; height: auto; display: block }
    figcaption {
      font-size: .73rem;
      color: var(--muted);
      margin-top: .35rem;
      text-align: center;
      font-style: italic;
    }

    /* ── Table column ── */
    .col-table {}
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
    .data-table th {
      padding: .35rem .45rem;
      border-bottom: 2px solid var(--fg);
      font-size: .7rem;
      font-weight: 600;
      text-align: right;
      letter-spacing: .03em;
      white-space: nowrap;
    }
    .data-table th:first-child { text-align: left }
    .data-table td {
      padding: .25rem .45rem;
      border-bottom: 1px solid var(--border);
      text-align: right;
      font-family: var(--mono);
      font-size: .74rem;
      font-variant-numeric: tabular-nums;
      white-space: nowrap;
    }
    .data-table td:first-child { text-align: left; font-family: var(--font) }
    .data-table tr.fore td   { background: var(--fore-bg) }
    .data-table tr.sep  td   { border-top: 2px solid var(--fg) }
    .data-table .na          { color: var(--muted) }

    /* ── Detalles del modelo (collapsible) ── */
    .model-details {
      margin-top: 2rem;
      font-size: .8rem;
    }
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
    details > summary::after { content: " ▸"; font-size: .7rem }
    details[open] > summary::after { content: " ▾" }
    .model-grid {
      display: flex;
      flex-wrap: wrap;
      gap: .3rem 2rem;
      padding: .7rem 0;
      font-size: .79rem;
    }
    .model-grid span { color: var(--muted) }
    .model-grid strong { font-family: var(--mono) }

    /* ── Narrative ── */
    .narrative {
      margin-top: 2rem;
      border-top: 1px solid var(--border);
      padding-top: 1.2rem;
      line-height: 1.7;
    }
    .narrative h2 {
      font-size: .72rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: .08em;
      color: var(--muted);
      margin-bottom: .8rem;
    }
    .narrative p { margin-bottom: .7rem }

    /* ── Footer ── */
    footer {
      margin-top: 2.5rem;
      padding-top: .7rem;
      border-top: 1px solid var(--border);
      font-size: .71rem;
      color: var(--muted);
    }

    /* ── Print / PDF ── */
    @media print {
      body { padding: 0; font-size: 10pt }
      .report-body { grid-template-columns: 1fr 1fr }
      .data-table  { font-size: 7.5pt }
      details > summary { display: none }
      details { display: none }
    }

    /* ── Responsive: stack on narrow screens ── */
    @media (max-width: 860px) {
      .report-body { grid-template-columns: 1fr }
      .col-table   { order: -1 }
    }
  </style>
</head>
<body>
<div class="container">

<header>
  <div class="sps-label">SPS{% if sps_name %}: {{ sps_name }}{% endif %}</div>
  <h1>{{ page_title }}</h1>
  <div class="meta">
    {% if source %}<span>Fuente&#160;<strong>{{ source }}</strong></span>{% endif %}
    <span>Origen de previsión&#160;<strong>{{ origin }}</strong></span>
    <span>Horizonte&#160;<strong>{{ horizon }}&#160;{{ freq_label }}</strong></span>
    <span>Generado&#160;<strong>{{ generated }}</strong></span>
  </div>
</header>

<div class="report-body">

  <!-- LEFT: charts -->
  <div class="col-charts">
    <figure>
      {{ forecast_svg | safe }}
      <figcaption>Variación anual (%) — histórico y previsión con banda ±1σ</figcaption>
    </figure>
    <figure>
      {{ err_svg | safe }}
      <figcaption>Residuos estandarizados (ERR) — bandas ±2σ</figcaption>
    </figure>
  </div>

  <!-- RIGHT: table -->
  <div class="col-table">
    <div class="table-wrap">
    <table class="data-table">
      <caption>Previsión</caption>
      <thead>
        <tr>
          <th rowspan="2">Fecha</th>
          <th colspan="2">Nivel</th>
          <th colspan="2">Δ {{ diff1_label }}</th>
          <th colspan="2">Δ anual</th>
          <th rowspan="2">ERR</th>
        </tr>
        <tr>
          <th>Valor</th><th>σ</th>
          <th>%</th><th>σ</th>
          <th>%</th><th>σ</th>
        </tr>
      </thead>
      <tbody>
        {% for row in hist_rows %}
        <tr{% if loop.last %} class="sep"{% endif %}>
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
        {% for row in fore_rows %}
        <tr class="fore">
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
      </tbody>
    </table>
    </div>

    <!-- Model details (collapsible, for internal use) -->
    <div class="model-details">
      <details>
        <summary>Detalles del modelo</summary>
        <div class="model-grid">
          <div><span>Modelo&#160;</span><strong>{{ stem }}</strong></div>
          <div><span>npar&#160;</span><strong>{{ npar }}</strong></div>
          <div><span>σ²&#160;</span><strong>{{ sigma2 }}</strong></div>
          <div><span>AIC&#160;</span><strong>{{ aic }}</strong></div>
          <div><span>BIC&#160;</span><strong>{{ bic }}</strong></div>
          <div><span>N&#160;</span><strong>{{ nobs }}</strong></div>
          <div><span>Muestra&#160;</span><strong>{{ sample }}</strong></div>
        </div>
      </details>
    </div>
  </div>

</div><!-- .report-body -->

{% if narrative %}
<div class="narrative">
  <h2>Análisis</h2>
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
    """Write a self-contained HTML forecast report.

    Parameters
    ----------
    model     : Model  (model._result must be set — call forecast_fuf() first)
    fr        : ForecastResult
    path      : str — output .html path
    title     : str, optional — descriptive title, e.g. "Inflación España — IPC"
                If None, falls back to "A.<stem>"
    source    : str, optional — data source, e.g. "INE" or "Eurostat"
    sps_name  : str, optional — SPS label, e.g. "Inflación Eurozona"
    narrative : str, optional — HTML snippet injected in the Análisis section
    pdf       : bool — also write .pdf (requires weasyprint)
    """
    try:
        from jinja2 import Environment
    except ImportError:
        raise ImportError("write_forecast_report requires jinja2 — pip install jinja2")

    if model._result is None:
        raise RuntimeError(
            "write_forecast_report: model._result not set — call forecast_fuf() first"
        )

    hist_rows, fore_rows, meta = _table_data(model, fr)
    forecast_svg = _make_forecast_svg(model, fr)
    err_svg      = _make_err_svg(model, fr)

    stem       = meta["stem"]
    page_title = title or f"A.{stem}"

    r   = model._result
    ts  = model.series
    env = Environment(autoescape=False)
    html = env.from_string(_HTML).render(
        page_title  = page_title,
        sps_name    = sps_name or "",
        stem        = stem,
        source      = source or "",
        origin      = meta["origin"],
        horizon     = fr.horizon,
        freq_label  = meta["freq_label"],
        diff1_label = meta["diff1_label"],
        sigma2      = f"{fr.sigma2:.6f}",
        aic         = f"{r.aic:.2f}"  if hasattr(r, "aic")  else "—",
        bic         = f"{r.bic:.2f}"  if hasattr(r, "bic")  else "—",
        npar        = r.npar          if hasattr(r, "npar")  else "—",
        nobs        = ts.nobs,
        sample      = meta["sample"],
        hist_rows   = hist_rows,
        fore_rows   = fore_rows,
        forecast_svg = forecast_svg,
        err_svg      = err_svg,
        narrative   = narrative,
        version     = "1.13 (Python port)",
        generated   = _date.today().isoformat(),
    )

    path = str(path)
    if not path.endswith(".html"):
        path += ".html"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)

    if pdf:
        _write_pdf(html, path.replace(".html", ".pdf"))

    return path


# ── Chart helpers ─────────────────────────────────────────────────────────────

def _make_forecast_svg(model, fr) -> str:
    """Top-panel only: history + forecast line + dashed confidence bands."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .forecast import _boxcox
    from .plots import _tj_spines

    ts       = model.series
    nobs     = ts.nobs
    freq     = ts.freq if ts.freq > 0 else 1
    L        = fr.horizon
    refactor = model.refactor
    boxlam   = model.boxlam
    raw      = ts.data

    hist_annual = np.array([
        100.0 * (_boxcox(raw[nobs - L + i],        boxlam, refactor)
                 - _boxcox(raw[nobs - L + i - freq], boxlam, refactor)) / refactor
        for i in range(L)
    ])
    y_all  = np.concatenate([hist_annual, fr.seasonal_diff])
    x_all  = np.arange(2 * L)
    x_hist = np.arange(L)
    x_fore = np.arange(L, 2 * L)

    xtick_pos, xtick_lbl = _year_ticks(ts, nobs, L, freq)

    fig, ax = plt.subplots(figsize=(7, 3.5))
    _tj_spines(ax, ('left', 'bottom'))
    ax.plot(x_all,  y_all,       'k-',  lw=1.8, zorder=3)
    ax.plot(x_hist, hist_annual,  'ko', ms=3.5,  zorder=4)
    ax.plot(x_fore, fr.seasonal_diff + fr.seasonal_diff_std, 'k--', lw=1.2, zorder=2)
    ax.plot(x_fore, fr.seasonal_diff - fr.seasonal_diff_std, 'k--', lw=1.2, zorder=2)
    ax.axvline(L - 0.5, color='0.6', lw=0.9, zorder=1)
    ax.axhline(0,        color='k',   lw=0.7, zorder=1)
    if xtick_pos:
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(xtick_lbl, fontsize=9)
    ax.set_xlim(-0.5, 2 * L - 0.5)
    ax.tick_params(direction='out', labelsize=9)
    ax.set_title('Variación anual (%)', loc='left', fontsize=9, pad=4)
    fig.tight_layout(pad=0.6)
    return _fig_to_svg(fig)


def _make_err_svg(model, fr) -> str:
    """ERR panel only: residual impulses with ±2σ bands."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    from .plots import _tj_spines

    ts       = model.series
    freq     = ts.freq if ts.freq > 0 else 1
    L        = fr.horizon
    refactor = model.refactor

    residuals  = np.asarray(model._result.residuals)
    err_L      = min(L, len(residuals))
    err_vals   = 100.0 * residuals[-err_L:] / refactor
    x_err      = np.arange(err_L)
    sigma_plot = math.sqrt(fr.sigma2)
    prevcmax   = _prevcmax(err_vals, sigma_plot)

    xtick_pos, xtick_lbl = _year_ticks(ts, ts.nobs, L, freq)
    err_tick_pos = [p for p in xtick_pos if p < L]
    err_tick_lbl = xtick_lbl[:len(err_tick_pos)]

    fig, ax = plt.subplots(figsize=(7, 2.8))
    _tj_spines(ax, ('left', 'bottom'))
    ax.vlines(x_err, 0, err_vals, colors='k', lw=1.6, zorder=3)
    ax.axhline( 2 * sigma_plot, color='k', lw=1.0, ls='--', zorder=2)
    ax.axhline(-2 * sigma_plot, color='k', lw=1.0, ls='--', zorder=2)
    ax.axhline(0,               color='k', lw=1.2, zorder=2)
    margin = 0.1 * sigma_plot
    ax.set_ylim(-(prevcmax + margin), prevcmax + margin)
    yt = np.arange(0, prevcmax + 0.05 * sigma_plot, 2 * sigma_plot)
    ax.set_yticks(np.concatenate([-yt[1:][::-1], yt]))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
    ax.tick_params(direction='out', labelsize=9)
    if err_tick_pos:
        ax.set_xticks(err_tick_pos)
        ax.set_xticklabels(err_tick_lbl, fontsize=9)
    ax.set_xlim(-0.5, L - 0.5)
    ax.set_title('Residuos (ERR)', loc='left', fontsize=9, pad=4)
    fig.tight_layout(pad=0.6)
    return _fig_to_svg(fig)


# ── Table data ────────────────────────────────────────────────────────────────

def _table_data(model, fr):
    """Return (hist_rows, fore_rows, meta_dict)."""
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
    n_hist    = freq   # one full seasonal cycle of history

    # ── Historical rows ───────────────────────────────────────────────────────
    hist_rows = []
    for k in range(nobs - n_hist + 1, nobs + 1):
        date = _fuf_obs_to_date_str(k, begyear, begtime, freq)
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
            "date":  date,
            "level": _f2(level),
            "diff1": _f2(diff1),
            "annual": _f2(annual),
            "err":   err_str,
        })

    # ── Forecast rows ─────────────────────────────────────────────────────────
    fore_rows = []
    for h in range(L):
        k    = nobs + h + 1
        date = _fuf_obs_to_date_str(k, begyear, begtime, freq)
        lstd = fr.level_std[h] * refactor
        fore_rows.append({
            "date":       date,
            "level":      _f2(fr.level[h]),
            "level_std":  _f2(lstd),
            "diff1":      _f2(fr.diff1[h]),
            "diff1_std":  _f2(fr.diff1_std[h]),
            "annual":     _f2(fr.seasonal_diff[h]),
            "annual_std": _f2(fr.seasonal_diff_std[h]),
        })

    # ── Sample label ─────────────────────────────────────────────────────────
    from .report import _fuf_obs_to_date_str as _d
    sample = f"{_d(1, begyear, begtime, freq)} – {_d(nobs, begyear, begtime, freq)}"

    # ── Meta ──────────────────────────────────────────────────────────────────
    stem   = model._inp_stem or ts.name
    origin = _fuf_obs_to_date_str(nobs, begyear, begtime, freq)
    if freq == 12:
        freq_label, diff1_label = "meses", "mensual"
    elif freq == 4:
        freq_label, diff1_label = "trimestres", "trimestral"
    else:
        freq_label, diff1_label = "periodos", "periodo"

    return hist_rows, fore_rows, {
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
    fig.savefig(buf, format="svg", bbox_inches="tight")
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
