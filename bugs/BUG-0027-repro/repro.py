#!/usr/bin/env python3
"""BUG-0027: the forecast table of the .out is not fuf C's.

fue writes the table of `write_fuf_out` (report._fuf_section_forecast_table)
and atsw-gui's fuf writes it in its .out. fuf fixed its table on 15-nov-2026
(its BUG-0001, BUG-0002, BUG-0003); the port must say the same, byte for byte.

Two cases, each a fuf input and the table fuf writes for it:

- CPI_USA_model_AR: in levels (lambda = 1). The port multiplied the variations
  by 100 and wrote the level's std as sqrt(v1) without the scale.
- MINK_m01: an annual series (Jenkins and Alavi's mink). The port filled the
  annual columns, and the first observed rows' "annual change" was the level.

Run with the fue to check (`python repro.py`). It prints the first line that
differs and exits 1 while the table is not fuf's; 0 when both are.
"""
import pathlib
import sys

import fue

AQUI = pathlib.Path(__file__).resolve().parent


def tabla(texto):
    lineas = texto.splitlines()
    i = next(k for k, l in enumerate(lineas) if l.startswith("FORECAST REPORT"))
    out, bordes = [], 0
    for l in lineas[i:]:
        out.append(l)
        bordes += l.startswith(" +--")
        if bordes == 3:
            break
    return out


def main():
    abierto = False
    for caso in ("CPI_USA_model_AR", "MINK_m01"):
        ts, m = fue.load_fuf(str(AQUI / f"{caso}.fuf.inp"))
        got = tabla(m.write_fuf_out(m.forecast_fuf()))
        ref = (AQUI / f"{caso}.tabla").read_text(encoding="latin-1").splitlines()
        if got == ref:
            print(f"{caso}: the table is fuf's")
            continue
        abierto = True
        k = next((i for i, (a, b) in enumerate(zip(got, ref)) if a != b),
                 min(len(got), len(ref)))
        print(f"{caso}: the table is not fuf's (line {k + 1})")
        print(f"  fuf: {ref[k] if k < len(ref) else '<none>'}")
        print(f"  fue: {got[k] if k < len(got) else '<none>'}")
    return 1 if abierto else 0


if __name__ == "__main__":
    sys.exit(main())
