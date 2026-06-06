"""
SPS:Inflation — Generate HTML forecast reports for all countries.

Usage
-----
    python sps_generate.py [--output-dir DIR]

Default output: each report is written next to its source .inp file.
The index page (sps_index.html) is written to the Analisis/ directory.

Dependencies: fue Python package  (pip install -e path/to/fue)
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import date

# ── SPS configuration ─────────────────────────────────────────────────────────

BASE = Path(__file__).parent

SPS_ENTRIES = [
    {
        "id":       "canada",
        "country":  "Canada",
        "title":    "Canada CPI Inflation",
        "source":   "Statistics Canada",
        "sps_name": "G7 Inflation Monitor",
        "fuf":      BASE / "Canada/forecast/CA.1.inp",
        "out_dir":  BASE / "Canada/forecast",
        "md_doc":   BASE / "Canada/forecast/Canada_CA1.md",
    },
    {
        "id":       "emu",
        "country":  "Euro Area",
        "title":    "Euro Area HICP Inflation",
        "source":   "Eurostat",
        "sps_name": "G7 Inflation Monitor",
        "fuf":      BASE / "EMU/forecast/EU.2.inp",
        "out_dir":  BASE / "EMU/forecast",
        "md_doc":   BASE / "EMU/forecast/EMU_EU2.md",
    },
    {
        "id":       "germany",
        "country":  "Germany",
        "title":    "Germany CPI Inflation",
        "source":   "Destatis",
        "sps_name": "G7 Inflation Monitor",
        "fuf":      BASE / "Germany/forecast/G.2.inp",
        "out_dir":  BASE / "Germany/forecast",
        "md_doc":   BASE / "Germany/forecast/Germany_G2.md",
    },
    {
        "id":       "spain",
        "country":  "Spain",
        "title":    "Spain CPI Inflation",
        "source":   "INE",
        "sps_name": "G7 Inflation Monitor",
        "fuf":      BASE / "Spain/forecast_b2025/forecast_S.2.inp",
        "out_dir":  BASE / "Spain/forecast_b2025",
        "md_doc":   BASE / "Spain/forecast_b2025/Spain_S2.md",
    },
    {
        "id":       "uk",
        "country":  "UK",
        "title":    "UK CPI Inflation",
        "source":   "ONS",
        "sps_name": "G7 Inflation Monitor",
        "fuf":      BASE / "UK/Forecast/UK.3.inp",
        "out_dir":  BASE / "UK/Forecast",
        "md_doc":   BASE / "UK/Forecast/UK_UK3.md",
    },
    {
        "id":       "usa",
        "country":  "USA",
        "title":    "USA CPI Inflation",
        "source":   "BLS",
        "sps_name": "G7 Inflation Monitor",
        "fuf":      BASE / "USA/forecast/US.3.inp",
        "out_dir":  BASE / "USA/forecast",
        "md_doc":   BASE / "USA/forecast/USA_US3.md",
    },
]


# ── Index page template ───────────────────────────────────────────────────────

_INDEX_HTML = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>SPS: G7 Inflation Monitor</title>
  <style>
    :root {{
      --font: system-ui, -apple-system, "Segoe UI", sans-serif;
      --fg: #111827; --muted: #6b7280; --border: #e5e7eb;
      --accent: #1d4ed8; --fore-bg: #eff6ff;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0 }}
    body {{ font-family: var(--font); color: var(--fg); background: #fff;
            padding: 2rem 1.5rem; font-size: 14px }}
    .container {{ max-width: 900px; margin: 0 auto }}
    header {{ border-bottom: 2.5px solid var(--fg); padding-bottom: .9rem; margin-bottom: 1.8rem }}
    .sps-label {{ font-size: .72rem; font-weight: 600; text-transform: uppercase;
                  letter-spacing: .1em; color: var(--accent); margin-bottom: .3rem }}
    h1 {{ font-size: 1.55rem; font-weight: 700 }}
    .meta {{ font-size: .81rem; color: var(--muted); margin-top: .5rem }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 1.5rem }}
    th {{ padding: .4rem .7rem; border-bottom: 2px solid var(--fg);
          font-size: .72rem; text-transform: uppercase; letter-spacing: .06em;
          text-align: left; font-weight: 600; color: var(--muted) }}
    td {{ padding: .35rem .7rem; border-bottom: 1px solid var(--border);
          font-size: .82rem; vertical-align: middle }}
    tr:hover td {{ background: #f9fafb }}
    a {{ color: var(--accent); text-decoration: none; font-weight: 500 }}
    a:hover {{ text-decoration: underline }}
    .flag {{ font-size: 1.1rem }}
    .stale {{ color: #b45309; font-size: .75rem; margin-left: .4rem }}
    footer {{ margin-top: 2rem; padding-top: .7rem; border-top: 1px solid var(--border);
              font-size: .71rem; color: var(--muted) }}
  </style>
</head>
<body>
<div class="container">
<header>
  <div class="sps-label">SPS: G7 Inflation Monitor</div>
  <h1>Inflation Forecast Dashboard</h1>
  <div class="meta">Generated {generated} · FUE Python 1.13</div>
</header>

<table>
  <thead>
    <tr>
      <th>Country</th>
      <th>Forecast Origin</th>
      <th>Horizon</th>
      <th>σ²</th>
      <th>Source</th>
      <th>Report</th>
    </tr>
  </thead>
  <tbody>
    {rows}
  </tbody>
</table>

<footer>FUE Python 1.13 (Python port) · {generated}</footer>
</div>
</body>
</html>
"""

_INDEX_ROW = """\
    <tr>
      <td><strong>{country}</strong></td>
      <td>{origin}{stale}</td>
      <td>{horizon} months</td>
      <td>{sigma2}</td>
      <td>{source}</td>
      <td><a href="{rel_path}">View report →</a></td>
    </tr>"""


# ── Generator ─────────────────────────────────────────────────────────────────

def generate_all(output_dir=None, verbose=True):
    """Generate HTML reports for all SPS entries and write the index page."""
    # lazy import — keep script runnable without installing fue in advance
    sys.path.insert(0, str(Path(__file__).parent.parent.parent /
                           "SRC/atws/fue/fue/src"))
    try:
        from fue.inp import load_fuf
        from fue.report_forecast import write_forecast_report
        from fue.report import _fuf_obs_to_date_str
    except ImportError:
        from fue.inp import load_fuf
        from fue.report_forecast import write_forecast_report
        from fue.report import _fuf_obs_to_date_str

    index_rows = []
    results = []

    for entry in SPS_ENTRIES:
        country = entry["country"]
        fuf_path = Path(entry["fuf"])
        out_dir  = Path(output_dir or entry["out_dir"])
        out_dir.mkdir(parents=True, exist_ok=True)
        out_html = out_dir / f"{entry['id']}_forecast.html"

        if not fuf_path.exists():
            print(f"  SKIP  {country}: {fuf_path} not found")
            continue

        try:
            ts, model = load_fuf(str(fuf_path))
            fr = model.forecast_fuf()
        except Exception as exc:
            print(f"  ERROR {country}: {exc}")
            continue

        # Generate HTML report
        write_forecast_report(
            model, fr,
            path=str(out_html),
            title=entry["title"],
            source=entry["source"],
            sps_name=entry["sps_name"],
        )

        # Compute origin date
        by, bt = ts.start
        origin = _fuf_obs_to_date_str(ts.nobs, by, bt, ts.freq)

        # Flag stale data (origin > 6 months ago)
        try:
            o_month, o_year = int(origin.split("/")[0]), int(origin.split("/")[1])
            today = date.today()
            months_ago = (today.year - o_year) * 12 + (today.month - o_month)
            stale_tag = '<span class="stale">⚠ stale</span>' if months_ago > 6 else ""
        except Exception:
            stale_tag = ""

        rel_path = os.path.relpath(out_html, BASE)
        index_rows.append(_INDEX_ROW.format(
            country  = country,
            origin   = origin,
            stale    = stale_tag,
            horizon  = fr.horizon,
            sigma2   = f"{fr.sigma2:.6f}",
            source   = entry["source"],
            rel_path = rel_path,
        ))
        results.append((country, str(out_html)))
        if verbose:
            print(f"  OK    {country:12s}  origin={origin}  → {out_html}")

    # Write index
    index_path = BASE / "sps_index.html"
    index_html = _INDEX_HTML.format(
        generated = date.today().isoformat(),
        rows      = "\n".join(index_rows),
    )
    index_path.write_text(index_html, encoding="utf-8")
    if verbose:
        print(f"\n  INDEX → {index_path}")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SPS inflation reports")
    parser.add_argument("--output-dir", default=None,
                        help="Override output directory for all reports")
    args = parser.parse_args()
    print("SPS: G7 Inflation Monitor — generating reports\n")
    generate_all(output_dir=args.output_dir)
