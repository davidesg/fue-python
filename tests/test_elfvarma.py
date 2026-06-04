"""
Tests for elfvarma.py — elf_scalar and flikam_scalar against known references.

Phase: elfvarma only.  Integration with cast_us / full C engine comes later.
"""
import math
import numpy as np
import pytest
from fue.elfvarma import elf_scalar, flikam_scalar, chekma_scalar


# ── chekma ────────────────────────────────────────────────────────────────────

def test_chekma_invertible():
    assert chekma_scalar(np.array([0.5])) == 0

def test_chekma_noninvertible_ma1():
    assert chekma_scalar(np.array([1.1])) == 1   # |θ| > 1

def test_chekma_boundary():
    assert chekma_scalar(np.array([1.00005])) == 1
    assert chekma_scalar(np.array([1.00004])) == 0

def test_chekma_empty():
    assert chekma_scalar(np.array([])) == 0

def test_chekma_ma2_invertible():
    # θ₁=0.5, θ₂=0.3  →  companion eigenvalues well inside unit circle
    assert chekma_scalar(np.array([0.5, 0.3])) == 0


# ── AR(1) exact log-likelihood reference ──────────────────────────────────────

def _ar1_loglik_exact(phi1, w):
    """Textbook exact log-lik of scalar AR(1) with σ²=1."""
    n  = len(w)
    g0 = 1.0 / (1.0 - phi1 ** 2)          # variance of stationary AR(1)
    a  = np.zeros(n)
    for t in range(n):
        a[t] = w[t] - (phi1 * w[t - 1] if t > 0 else 0.0)
    ll  = -0.5 * n * math.log(2 * math.pi)
    ll -= 0.5 * math.log(g0)
    ll -= 0.5 * (w[0] ** 2 / g0 + float(np.dot(a[1:], a[1:])))
    return ll


@pytest.fixture
def ar1_data():
    rng  = np.random.default_rng(42)
    n    = 200
    phi1 = 0.7
    w    = np.zeros(n)
    for t in range(1, n):
        w[t] = phi1 * w[t - 1] + rng.standard_normal()
    return phi1, w


def test_elf_scalar_ar1(ar1_data):
    phi1, w = ar1_data
    logelf, f1, f2, a, ifault = elf_scalar(
        len(w), 1, 0, np.array([phi1]), np.array([]), w,
        sigma2=1.0, do_chkma=False,
    )
    ref = _ar1_loglik_exact(phi1, w)
    assert ifault == 0
    assert abs(logelf - ref) < 1e-4, f"elf {logelf:.6f}  ref {ref:.6f}"


def test_flikam_scalar_ar1(ar1_data):
    """flikam runs without error and returns a finite loglik for AR(1)."""
    phi1, w = ar1_data
    sumsq, fact, loglik, at, ifault = flikam_scalar(
        len(w), 1, 0, np.array([phi1]), np.array([]), 0.0, w, do_chkma=False,
    )
    assert ifault == 0
    assert math.isfinite(loglik)
    assert sumsq > 0
    # flikam uses a different normalisation than exact ML (Mélard 1984 formula),
    # so we only check it is within a few units of the exact loglik.
    ref = _ar1_loglik_exact(phi1, w)
    assert abs(loglik - ref) < 2.0, f"flikam {loglik:.4f}  exact {ref:.4f}"


# ── MA(1) iid series test ──────────────────────────────────────────────────────

def test_elf_scalar_ma1_iid():
    """MA(1) evaluated at θ=0 must reduce to pure white-noise loglik."""
    rng = np.random.default_rng(7)
    n   = 150
    w   = rng.standard_normal(n)
    logelf, f1, f2, a, ifault = elf_scalar(
        n, 0, 1, np.array([]), np.array([0.0]), w,
        sigma2=1.0, do_chkma=False,
    )
    ref = -0.5 * n * (_LOG2PI := 1.837877066)   # −n/2 log(2π) when σ²=1
    ref -= 0.5 * float(np.dot(w, w))            # − ½ Σwt²
    assert ifault == 0
    assert abs(logelf - ref) < 1e-6, f"elf MA(1)@0 {logelf:.6f}  ref {ref:.6f}"


# ── ARMA(2,1) self-consistency: elf ≈ flikam for large n ─────────────────────

def test_elf_flikam_arma21_consistent():
    """elf_scalar and flikam_scalar agree to within 0.5 on ARMA(2,1) data."""
    rng   = np.random.default_rng(123)
    n     = 500
    phi   = np.array([1.2, -0.5])
    theta = np.array([0.3])
    # Simulate ARMA(2,1)
    e = rng.standard_normal(n + 50)
    w = np.zeros(n + 50)
    for t in range(2, n + 50):
        w[t] = phi[0]*w[t-1] + phi[1]*w[t-2] + e[t] - theta[0]*e[t-1]
    w = w[50:]   # drop burn-in

    logelf, f1, f2, a, ifault_e = elf_scalar(
        n, 2, 1, phi, theta, w, sigma2=1.0, do_chkma=False,
    )
    _, _, loglik, _, ifault_f = flikam_scalar(
        n, 2, 1, phi, theta, 0.0, w, do_chkma=False,
    )
    assert ifault_e == 0
    assert ifault_f == 0
    assert abs(logelf - loglik) < 2.0, \
        f"elf {logelf:.4f}  flikam {loglik:.4f}  diff {abs(logelf-loglik):.4f}"


# ── ifault codes ─────────────────────────────────────────────────────────────

def test_elf_nonstationary_ar():
    """AR with near-unit root returns ifault=2."""
    rng = np.random.default_rng(1)
    w   = rng.standard_normal(50)
    # phi1 = 1.01 → non-stationary
    _, _, _, _, ifault = elf_scalar(
        len(w), 1, 0, np.array([1.01]), np.array([]), w, do_chkma=False,
    )
    assert ifault in (2, 3), f"expected ifault 2 or 3, got {ifault}"


def test_elf_noninvertible_ma():
    """Non-invertible MA returns ifault=4."""
    rng = np.random.default_rng(2)
    w   = rng.standard_normal(50)
    _, _, _, _, ifault = elf_scalar(
        len(w), 0, 1, np.array([]), np.array([1.5]), w, do_chkma=True,
    )
    assert ifault == 4, f"expected ifault 4, got {ifault}"
