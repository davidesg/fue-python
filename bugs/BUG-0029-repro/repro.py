#!/usr/bin/env python3
"""BUG-0029: on macOS Intel the Python port (raxopt) and the C engine stop at
different points on a file outside the optimiser's band.

tests/test_qnewtopt.py compares fue's two engines ON THE SAME MACHINE. In the
atsw 2.1.3 preflight (2026-11-17, fue 0.1.22 from PyPI, macos-15-intel,
Python 3.12) two comparisons fail there and nowhere else (Linux, macOS arm64,
Windows):

  test_raxopt_sigma2_matches_c[R.4]: raxopt sigma2 0.00009672 vs C 0.00009453
  test_cov_raxopt_matches_c (GDP R.2): the standard errors, rtol 1e-5

R.4 (IPC trimestral, Coint) is lambda = 1 with refactor = 1: sigma2 ~ 1e-4, a
spread of ~0.01, below the band [0.1, 1000] of fue.scale -- where the stopping
point loses precision (atsw-gui docs/ESTUDIO-escala.md).

Runs both comparisons. Exits 1 if they disagree here; 0 if they agree on
macOS x86_64 (fixed); 2 if they agree on another platform (it does not
reproduce here: run it on a Mac Intel).
"""
import os
import platform
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REAL = os.path.join(os.path.dirname(os.path.dirname(HERE)), "tests", "real_cases")

import fue                                              # noqa: E402
from fue._engine import estimate as est_c               # noqa: E402
from fue.cast_us import estimate_py                     # noqa: E402


def load(rel):
    return fue.load(os.path.join(REAL, rel + ".inp"))[1]


def main():
    from fue._engine import engine_backend
    if engine_backend() != "c":
        print(f"fue's C engine is not loaded here ({engine_backend()}): the comparison "
              "would set the port against itself. Install fue's binary wheel.")
        return 2
    bad = []
    m = load("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.4")
    rc, rp = est_c(m), estimate_py(m)
    print(f"R.4 (lambda {m.boxlam:g}, refactor {m.refactor:g}): sigma2 C {rc['sigma2']:.8f} "
          f"Python {rp['sigma2']:.8f}; loglik {rc['loglik']:.6f} / {rp['loglik']:.6f}")
    if abs(rp["loglik"] - rc["loglik"]) <= 1.0 and abs(rp["sigma2"] - rc["sigma2"]) >= 2e-6:
        bad.append("R.4 sigma2")
    m = load("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.2")
    rc, rp = est_c(m), estimate_py(m)
    se_c, se_p = np.asarray(rc["std_errors"]), np.asarray(rp["std_errors"])
    print(f"GDP R.2: s.e. C {np.round(se_c, 6)}  Python {np.round(se_p, 6)}")
    if not np.allclose(se_p, se_c, rtol=1e-5, atol=0):
        bad.append("R.2 std errors")
    mac_intel = platform.system() == "Darwin" and platform.machine() == "x86_64"
    if bad:
        print("OPEN: the two engines disagree here:", ", ".join(bad))
        return 1
    if mac_intel:
        print("FIXED on macOS x86_64")
        return 0
    print(f"they agree on {platform.system()} {platform.machine()}: the defect shows on "
          "macOS x86_64 only -- run this there")
    return 2


if __name__ == "__main__":
    sys.exit(main())
