"""BUG-0025: the line search never returned when the objective was NaN.

`lnsrch` (Dennis & Schnabel A6.3.1) decides with comparisons, and with a NaN
every comparison is false: in the C the step was neither accepted nor
abandoned, lambda became NaN and the loop spun for ever. The Python port
survived only through the direction of one comparison, and then crashed on the
NaN gradient. Both now treat a non-finite point as inadmissible, as drvarma,
drtran and the atsw-gui engines do. Each run is a subprocess with a timeout: a
regression must fail this test, not hang the suite.
"""
import os
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPRO = os.path.join(ROOT, "bugs", "BUG-0025-repro")


def test_python_raxopt_returns_on_a_nan_objective():
    r = subprocess.run([sys.executable, "-W", "error",
                        os.path.join(REPRO, "repro_python.py")],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr[-500:]
    # it stops against the edge of the admissible region, x0 -> 1, x1 = -1/3
    assert "termcode 2" in r.stdout
    assert "0.9999" in r.stdout


def test_python_lnsrch_shrinks_past_a_nan_without_interpolating():
    import numpy as np
    from fue.qnewtopt import _lnsrch

    def f(x):
        return float("nan") if x[0] > 1.0 else (x[0] - 3.0) ** 2 / 9.0

    xkp1, fkp1, retcode, _maxtaken, lam = _lnsrch(
        f, np.array([0.0]), 1.0, np.array([-6.0 / 9.0]), np.array([6.0]),
        100.0, 1e-7)
    assert retcode == 0 and lam == pytest.approx(0.1)
    assert xkp1[0] == pytest.approx(0.6)


@pytest.mark.skipif(shutil.which("gcc") is None or shutil.which("pkg-config") is None,
                    reason="needs gcc and GSL to build the C harness")
def test_c_lnsrch_returns_on_a_nan_objective(tmp_path):
    gsl = subprocess.run(["pkg-config", "--cflags", "--libs", "gsl"],
                         capture_output=True, text=True)
    if gsl.returncode:
        pytest.skip("GSL not found by pkg-config")
    exe = tmp_path / "repro_c"
    internal = os.path.join(ROOT, "csrc", "internal")
    subprocess.run(["gcc", "-I" + internal, os.path.join(REPRO, "repro_c.c"),
                    os.path.join(internal, "qnewtopt.c"),
                    os.path.join(internal, "nlatools.c"), *gsl.stdout.split(),
                    "-lm", "-o", str(exe)], check=True)
    r = subprocess.run([str(exe)], capture_output=True, text=True, timeout=20)
    assert "retcode 0" in r.stdout and "x 0.6" in r.stdout, r.stdout
