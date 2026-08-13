"""The five graded examples must run, and must still say what they claim.

An example that stops running is worse than none: it is the first thing a reader
executes. So each one is run end to end here, and the numbers it prints are
asserted — not the prose, but the results the prose is about.

    01  ARIMA(1,1,0) on the CPI            the minimal flow
    02  airline (0,1,1)(0,1,1)_12          the canonical seasonal model
    03  ARMAX with step and impulse        the "X" half of the thesis
    04  harmonics against ∇_12             the seasonality half
    05  a mixed MEG                        the class only fue can specify

Examples 3-5 are simulated with a fixed seed precisely so they can be checked:
the parameters are known, and what the example prints must be near them. Three
specification errors were caught this way while they were being written — a
1-based intervention index, a harmonic phase off by one period, and a missing
mean on an undifferenced series (which reported sigma = 100 for a true 1).
"""
import os
import re
import subprocess
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_EX = os.path.join(_ROOT, "examples")

_EXAMPLES = [
    "01_minimal_arima.py",
    "02_airline.py",
    "03_armax_interventions.py",
    "04_harmonics_vs_annual_difference.py",
    "05_mixed_meg.py",
]


def _run(name):
    p = os.path.join(_EX, name)
    if not os.path.exists(p):
        pytest.skip(f"{name} not present")
    r = subprocess.run([sys.executable, p], capture_output=True, text=True,
                       cwd=_ROOT, timeout=600)
    assert r.returncode == 0, f"{name} failed:\n{r.stderr[-1500:]}"
    return r.stdout


def _num(texto, patron):
    m = re.search(patron, texto)
    assert m, f"pattern not found in the output: {patron}"
    return float(m.group(1))


@pytest.mark.parametrize("name", _EXAMPLES)
def test_the_example_runs(name):
    salida = _run(name)
    assert salida.strip(), f"{name} printed nothing"


def test_01_reports_a_converged_fit():
    salida = _run("01_minimal_arima.py")
    assert "converged : True" in salida
    phi = _num(salida, r"phi_1\s+=\s+([\d.-]+)")
    assert 0.0 < phi < 1.0, "the AR(1) is outside the stationarity region"
    assert "criterio del gradiente satisfecho" in salida, (
        "the minimal example no longer stops on the gradient — read it before "
        "changing this assertion")


def test_02_airline_is_invertible_and_seasonal():
    salida = _run("02_airline.py")
    Theta = _num(salida, r"Theta\s+=\s+([\d.-]+)")
    assert 0.0 < Theta < 1.0
    assert "59 usable" in salida, "the arithmetic of the annual difference moved"


def test_03_recovers_the_intervention_parameters():
    """The example's whole point: with n=600 the estimates land on the truth."""
    salida = _run("03_armax_interventions.py")
    omega_step = _num(salida, r"omega step\s+12\.000\s+([\d.-]+)")
    omega_pulse = _num(salida, r"omega pulse\s+-8\.000\s+([\d.-]+)")
    phi = _num(salida, r"phi\s+0\.700\s+([\d.-]+)")

    assert omega_step == pytest.approx(12.0, abs=1.0)
    assert omega_pulse == pytest.approx(-8.0, abs=2.5)
    assert phi == pytest.approx(0.7, abs=0.06)
    assert "absorbed by the AR" in salida


def test_04_the_ifadf_equivalence_holds():
    """ifadf all set IS the annual difference — the identity the example shows."""
    salida = _run("04_harmonics_vs_annual_difference.py")
    dif = _num(salida, r"difference ([\deE.+-]+)\s+← the same operator")
    assert dif < 1e-8, f"ifadf-all and D=1 differ by {dif}, they must not"

    cos1 = _num(salida, r"cos\(f=1\)\s+=\s+([\d.-]+)")
    sin1 = _num(salida, r"sin\(f=1\)\s+=\s+([\d.-]+)")
    assert cos1 == pytest.approx(3.0, abs=0.4), "the harmonic phase convention moved"
    assert sin1 == pytest.approx(1.5, abs=0.4)
    assert "NOT comparable" in salida, (
        "the warning against comparing likelihoods across differencing is gone")


def test_05_the_mixed_meg_separates_the_two_frequencies():
    salida = _run("05_mixed_meg.py")
    cos2 = _num(salida, r"cos\(f=2\)\s+=\s+([\d.-]+)")
    sigma = _num(salida, r"sigma\s+=\s+([\d.-]+)")
    radius = _num(salida, r"radius r = ([\d.-]+)")

    assert cos2 == pytest.approx(2.5, abs=0.4), "the deterministic harmonic moved"
    assert sigma == pytest.approx(1.0, abs=0.5), (
        "sigma is far from the truth — check that mu is being estimated: with "
        "ifadf slot 0 off, nothing else removes the level")
    assert radius < 0.95, (
        "the MA_f witness is at the boundary, which would say the integration "
        "at f=1 is not needed — but the series is integrated there by "
        "construction")
    assert "[0, 1, 0, 0, 0, 0, 0]" in salida
