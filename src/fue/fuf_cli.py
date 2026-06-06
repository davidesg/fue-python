"""
fuf command-line interface.

Mirrors the C fuf binary interface:

    fuf forecast_model

where «forecast_model» is the base name of the forecast .inp file written
by «fue model -f [horizon]».  Extension optional.

Reads:
    forecast_model.inp  — fuf input file with fixed params + horizon/sigma2

Produces:
    forecast_model.out  — forecast report (ASCII)
"""

import os
import sys


_VERSION = "1.13 (Python port)"
_BANNER = f"FUF {_VERSION}: Copyright (C) 2026 A.B. Treadway & D.E. Guerrero"


def _parse_args(argv):
    for a in argv:
        if a in ("-h", "--help"):
            print(_BANNER)
            print()
            print("Usage: fuf forecast_input")
            print()
            print("  forecast_input   fuf input file (omit .inp extension)")
            sys.exit(0)
        if a == "--version":
            print(_BANNER)
            sys.exit(0)

    if not argv:
        return None

    return argv[0]


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        print(_BANNER)
        print()
        print("Usage: fuf forecast_input")
        print("       fuf --help  for details")
        sys.exit(1)

    input_arg = _parse_args(argv)

    base = input_arg
    if base.endswith(".inp"):
        base = base[:-4]
    inp_path = base + ".inp"
    out_path = base + ".out"
    namepart = os.path.basename(base)

    print()
    print(_BANNER)
    print()
    print(f"Input file             : {inp_path}")
    print(f"Output file            : {out_path}")
    print()

    try:
        from .inp import load_fuf
    except ImportError:
        from fue.inp import load_fuf

    if not os.path.isfile(inp_path):
        print(f"fuf: cannot open input file: {inp_path}", file=sys.stderr)
        sys.exit(1)

    ts, model = load_fuf(inp_path)

    horizon = getattr(model, "_fuf_horizon", None)
    sigma2  = getattr(model, "_fuf_sigma2",  None)
    print(f"Forecast horizon       : {horizon}")
    print(f"sigma2                 : {sigma2:.10f}" if sigma2 else "sigma2                 : (estimated)")
    print()

    try:
        fr = model.forecast_fuf()
    except Exception as exc:
        print(f"fuf: forecast failed — {exc}", file=sys.stderr)
        sys.exit(1)

    inp_name = os.path.basename(inp_path)
    out_name = os.path.basename(out_path)
    model.write_fuf_out(fr, path=out_path, inp_name=inp_name, out_name=out_name)
    print(f"Written: {out_path}")
