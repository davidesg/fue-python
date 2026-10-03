"""BUG-0005 — the same fits on every platform.

The cross-platform half of BUG-0005: in July 2026 the Windows wheel of fue 0.1.7
sent US_CPI's AR(2)×(2,0,0)_12 to a spurious optimum that Linux did not reach.
art BUG-0118 (2026-09-08) reran it on Windows with the corrected seed and it
matched Linux. This script checks it again with the current optimizer (with the
restart) and adds R.4, the fit whose stopping point depended on the last bit.

It prints one line per fit, with every number to many digits, so the Linux and
Windows logs of the `platform-check` workflow can be compared line by line, and
exits 1 if a fit misses its reference optimum.

    python scripts/platform_check_bug0005.py
"""
import copy
import os
import platform
import sys
import warnings

import numpy as np

import fue

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
US_CPI = os.path.join(ROOT, "tests/data/bug_0005/US_CPI.pre")
COINT = os.path.join(ROOT, "tests/real_cases/PRICES/IPC/Trimestral/"
                     "Sample_1.2003_4.2019/Mod/Coint")

#: Reference optima (Linux). US_CPI: the correct basin, μ̂≈0.0021, σ̂ₐ≈0.0026
#: in proportion; the July spurious one had μ̂≈−0.144.
REF = {"US_CPI": None, "R.4": 251.682958, "R.1": None}


def fit(name, path, ar=None, ar_s=None, mu=None):
    _, m = fue.load(path)
    if ar is not None:
        m.ar = copy.deepcopy(ar)
    if ar_s is not None:
        m.ar_s = copy.deepcopy(ar_s)
    if mu is not None:
        m.mu0 = mu
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit()
    r = m._result
    p = np.asarray(r.params, float)
    print(f"{name:24} loglik {r.loglik:.6f}  sigma {np.sqrt(r.sigma2):.8f}  "
          f"conv {r.converged}  ifault {r.ifault}  "
          f"term {getattr(r, 'termcode', None)}  "
          f"restarts {getattr(r, 'restarts', None)}")
    print(f"{'':24} params {np.array2string(p, precision=6, max_line_width=200)}")
    return r


def main():
    print(f"fue {fue.__version__}  backend {fue.engine_backend()}  "
          f"{platform.system()} {platform.machine()}  "
          f"python {platform.python_version()}  numpy {np.__version__}\n")
    bad = []
    if fue.engine_backend() != "c":
        # The divergence under test is the C engine's (MSVC vs GCC); the Python
        # port would pass for the wrong reason.
        bad.append(f"no C engine: {fue.engine_load_error()}")

    # US_CPI from three starts: its optimum, the July wrong-sign seasonal seed,
    # and the identified seed. All three must land on the same optimum.
    us = [fit("US_CPI from .pre", US_CPI),
          fit("US_CPI wrong-sign seed", US_CPI, ar=[[0.60, -0.17]],
              ar_s=[[0.04, 0.08]], mu=0.0),
          fit("US_CPI identified seed", US_CPI, ar=[[0.60, -0.17]],
              ar_s=[[-0.11, -0.09]])]
    for r in us:
        mu = float(np.asarray(r.params, float)[-1])
        if abs(mu - 0.0021) > 5e-4 or abs(r.loglik - us[0].loglik) > 1e-4:
            bad.append(f"US_CPI: mu {mu:+.6f}, loglik {r.loglik:.6f}")

    r4 = fit("R.4", os.path.join(COINT, "R.4.inp"))
    if abs(r4.loglik - REF["R.4"]) > 1e-4:
        bad.append(f"R.4: loglik {r4.loglik:.6f} != {REF['R.4']}")
    r1 = fit("R.1", os.path.join(COINT, "R.1.inp"))
    if not r1.converged:
        bad.append("R.1 did not converge")

    print()
    if bad:
        print("FAIL  " + "\n      ".join(bad))
        sys.exit(1)
    print("OK  every fit reaches its reference optimum")


if __name__ == "__main__":
    main()
