"""
fue command-line interface.

Mirrors the C fue binary interface:

    fue model [eml|aml] [chk|nochk] [-f [horizon]]

where «model» is the base name of the .inp file (extension optional).

Produces:
    model.out   — estimation report (ASCII)
    model.pre   — pre-estimation file with estimated parameters as initials
    model.tex   — LaTeX report (not yet implemented; skipped with a warning)

With -f [N]:
    forecast_model.inp  — forecast .inp file (same format as .pre, horizon N)
"""

import os
import sys


_VERSION = "1.13 (Python port)"
_BANNER = f"FUE {_VERSION}: Copyright (C) 2026 A.B. Treadway & D.E. Guerrero"


def _parse_args(argv):
    """Manual parser — mirrors fue.c argument conventions exactly."""
    for a in argv:
        if a in ("-h", "--help"):
            print(_BANNER)
            print()
            print("Usage: fue input [eml|aml] [chk|nochk] [-f [horizon]]")
            print()
            print("  input        model-data file (omit .inp extension)")
            print("  eml | aml    exact | approximate maximum likelihood  (default: eml)")
            print("  chk | nochk  check | skip MA invertibility           (default: chk)")
            print("  -f [N]       write forecast .inp file, horizon N     (default: 24)")
            sys.exit(0)
        if a == "--version":
            print(_BANNER)
            sys.exit(0)

    if not argv:
        return None, True, True, False, 24

    input_arg = argv[0]
    eml = True
    chkma = True
    forecast_flag = False
    forecast_horizon = 24

    i = 1
    while i < len(argv):
        a = argv[i]
        if a == "eml":
            eml = True
        elif a == "aml":
            eml = False
        elif a == "chk":
            chkma = True
        elif a == "nochk":
            chkma = False
        elif a == "-f":
            forecast_flag = True
            if i + 1 < len(argv) and argv[i + 1].lstrip("-").isdigit():
                i += 1
                h = int(argv[i])
                forecast_horizon = h if h > 0 else 24
        else:
            print(f"fue: unknown option '{a}' — ignored.", file=sys.stderr)
        i += 1

    return input_arg, eml, chkma, forecast_flag, forecast_horizon


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        print(_BANNER)
        print()
        print("Usage: fue input [eml|aml] [chk|nochk] [-f [horizon]]")
        print("       fue --help  for details")
        sys.exit(1)

    input_arg, eml, chkma, forecast_flag, forecast_horizon = _parse_args(argv)

    # Derive base name and file paths
    base = input_arg
    if base.endswith(".inp"):
        base = base[:-4]
    inp_path = base + ".inp"
    out_path = base + ".out"
    dirpart  = os.path.dirname(base) or "."
    namepart = os.path.basename(base)
    pre_path = (base + ".pre") if not forecast_flag else os.path.join(dirpart, f"forecast_{namepart}.inp")

    # Header
    print()
    print(_BANNER)
    print()
    print(f"Input file             : {inp_path}")
    print(f"Output file            : {out_path}")
    print(f"Estimation method      : {'exact' if eml else 'approximate'} maximum likelihood")
    print(f"Check for invertibility: {'constrained' if chkma else 'unconstrained'} search")
    if forecast_flag:
        print(f"Forecast file          : {pre_path}  (horizon {forecast_horizon})")
    print()

    # Load model
    try:
        from .inp import load
    except ImportError:
        from fue.inp import load

    if not os.path.isfile(inp_path):
        print(f"fue: cannot open input file: {inp_path}", file=sys.stderr)
        sys.exit(1)

    ts, model = load(inp_path)

    # Apply CLI overrides
    model.eml   = eml
    model.chkma = chkma

    # Fit
    print("Estimating parameters...")
    try:
        model.fit()
    except RuntimeError as exc:
        print(f"fue: estimation failed — {exc}", file=sys.stderr)
        sys.exit(1)

    r = model._result
    print(f"  npar    = {r.npar}")
    print(f"  loglik  = {r.loglik:.10f}")
    print(f"  sigma²  = {r.sigma2:.10f}")
    print(f"  AIC     = {r.aic:.4f}")
    print(f"  BIC     = {r.bic:.4f}")
    print()

    # Write .out
    model.write_out(out_path)
    print(f"Written: {out_path}")

    # Write .pre / forecast .inp
    model.write_pre(pre_path)
    print(f"Written: {pre_path}")

    # LaTeX (.tex) — not yet implemented
    tex_path = base + ".tex"
    print(f"Skipped: {tex_path}  (LaTeX output not yet implemented)")
