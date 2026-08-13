"""AS 311 verified the only way it can be: by its published identities.

`literature/Mauricio-Algorithm311Exact-1997.pdf` (Applied Statistics 46,
157-171) publishes the algorithm but **not** the listing — unlike AS 197, whose
FORTRAN is printed in full and is executed in
`tests/test_as197_published_fortran.py`. What it does publish, on page 158, is
the exact log-likelihood in closed form:

    l(Φ,Θ,μ,σ²,Q|w) = −½{ nm·log(2πσ²) + n·log|Q| + log|ΛᵀΛ| + S(·)/σ² }   (2)

    S(Φ,Θ,μ,Q|w) = ηᵀη − (Mᵀh̃)ᵀ(I + MᵀHᵀHM)⁻¹(Mᵀh̃)                       (3)

    |ΛᵀΛ| = |I + MᵀHᵀHM|                                                   (4)

and, in the working paper this article summarises (`literature/9316.pdf`,
pp. 7-8), the procedure as **ten numbered steps** with equations (2.15)-(2.22).

So the verification here is of three kinds, and none of them needs the paper to
have shipped code:

1. **Equation (2) on our own outputs.** `elf_scalar` returns `f1` = S and
   `f2` = |ΛᵀΛ|^(1/n); the published formula must reproduce `logelf` for any σ².

2. **Equations (3) and (4) against a different author.** For a scalar model,
   Melard's `SUMSQ` *is* Mauricio's S and Melard's `FACT` *is* |ΛᵀΛ|^(1/n).
   Two published algorithms, two implementations that share nothing, computing
   the same two quantities. Measured 13 August 2026, over six Box-Jenkins
   specifications: **S agrees to 1e-14 and |ΛᵀΛ|^(1/n) to 4.4e-16**.

3. **The ten steps, traceable.** Mauricio's C numbers its blocks `[1]`…`[9]`
   in the order the working paper numbers the steps, and the Python port marks
   `(a)`…`(k)` with the AS 311 equation numbers. That correspondence is what
   makes the claim "this code is that algorithm" checkable by a reader, so it
   is asserted here rather than left to good intentions.

What this cannot do is check the paper's *numerical example* (WP 9316, Tables 4
and 5): the estimates are published to two decimals but the series is not, so
there is nothing to run.
"""
import math
import os
import shutil
import subprocess

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_BJ = os.path.expanduser(
    "~/Dropbox/SRC/drvus-source/1.2.01/drvus/src/Box_y_Jenkings")

_CASES = ["SeriesA/a2.inp", "SeriesB/b.inp", "SeriesC/c2.inp",
          "SeriesA/a1.inp", "SeriesC/c.inp", "SeriesE/e1.inp"]


def _fit(rel):
    import fue

    ts, m = fue.load(os.path.join(_BJ, rel))
    y = np.array(ts.data)
    d = int(m.d)
    if m.estimate_mu:
        m.mu0 = float((np.diff(y, d) if d else y).mean())
    m.fit()
    w = np.diff(y, d) if d else y
    if m.estimate_mu:
        w = w - m._result.params[-1]
    phi = np.array([c for f in m.ar for c in f])
    theta = np.array([c for f in m.ma for c in f])
    return np.asarray(w, float), phi, theta


def _elf(w, phi, theta, sigma2=1.0):
    """Exact evaluation: xitol NEGATIVE, which is the AS 197 TOLER contract
    that fue reuses. With a positive xitol the Ξ sequence is truncated and the
    quantities below stop being the exact ones — measured, it moves S by 1e-5
    relative, which is enough to fail every assertion here."""
    from fue.elfvarma import elf_scalar

    return elf_scalar(len(w), len(phi), len(theta), phi, theta, w,
                      sigma2=sigma2, xitol=-1e-3)


needs_bj = pytest.mark.skipif(not os.path.isdir(_BJ),
                              reason="DRVUS Box-Jenkins fixtures not present")


# ── 1. equation (2), on the engine's own outputs ───────────────────────────

@needs_bj
@pytest.mark.parametrize("rel", _CASES)
@pytest.mark.parametrize("sigma2", [1.0, 0.5, 2.5])
def test_equation_2_reproduces_the_log_likelihood(rel, sigma2):
    """l = −½{n·log(2πσ²) + n·log|Q| + log|ΛᵀΛ| + S/σ²}, with m=1 and |Q|=1.

    The tolerance is 1e-6 and not machine epsilon for a reason worth knowing:
    the engine carries log|ΛᵀΛ| internally but exports its n-th root, so the
    round trip through `f2 = exp(log|ΛᵀΛ|/n)` and back costs about 2e-8 on a
    series of this length. It is the export that loses it, not the algorithm.
    """
    w, phi, theta = _fit(rel)
    logelf, f1, f2, _a, ifault = _elf(w, phi, theta, sigma2)
    assert ifault == 0

    n = len(w)
    eq2 = -0.5 * (n * math.log(2 * math.pi * sigma2)
                  + n * math.log(f2) + f1 / sigma2)
    assert logelf == pytest.approx(eq2, abs=1e-6)


@needs_bj
@pytest.mark.parametrize("rel", _CASES)
def test_the_concentrated_variance_is_S_over_n(rel):
    """From (2), ∂l/∂σ² = 0 gives σ̂² = S/(nm) — the relation the estimation
    driver relies on to concentrate σ² out of the search."""
    w, phi, theta = _fit(rel)
    _ll, f1, _f2, _a, _if = _elf(w, phi, theta)
    n = len(w)

    grid = [f1 / n * k for k in (0.9, 0.99, 1.0, 1.01, 1.1)]
    lls = [_elf(w, phi, theta, s2)[0] for s2 in grid]
    assert lls.index(max(lls)) == 2, "σ̂² = S/n is not the maximiser"


@needs_bj
def test_minimising_S_times_the_determinant_is_maximising_the_likelihood():
    """AS 311, *Additional Comments*: maximising the exact likelihood is
    equivalent to minimising S^m·|Q|^m·|ΛᵀΛ|^(1/n) — for m=1 and |Q|=1, S·f2.

    This is the identity the optimiser's objective function rests on, so it is
    checked along a grid and not only at the optimum.
    """
    w, _phi, theta0 = _fit("SeriesA/a2.inp")
    n = len(w)
    objective, concentrated = [], []
    for th in np.linspace(0.40, 0.95, 23):
        _ll, f1, f2, _a, ifault = _elf(w, np.array([]), np.array([th]))
        assert ifault == 0
        objective.append(f1 * f2)
        concentrated.append(-0.5 * n * (math.log(2 * math.pi) + 1.0
                                        + math.log(f1 / n) + math.log(f2)))

    assert int(np.argmin(objective)) == int(np.argmax(concentrated))
    # and the maximiser is where fue's own estimate sits
    grid = np.linspace(0.40, 0.95, 23)
    assert grid[int(np.argmax(concentrated))] == pytest.approx(theta0[0],
                                                              abs=0.03)


# ── 2. equations (3) and (4), against Melard's published FORTRAN ───────────

@needs_bj
@pytest.mark.skipif(shutil.which("gfortran") is None,
                    reason="gfortran not available")
@pytest.mark.parametrize("rel", _CASES)
def test_S_and_the_determinant_match_algorithm_as197(rel, tmp_path):
    """The strongest check available for a paper that published no code.

    Melard (1984) and Mauricio (1997) derive the exact ARMA likelihood by
    different routes — Kalman-type recursions against the innovations form with
    Cholesky factorisations — and both express it through a sum of squares and
    a determinant factor. If `SUMSQ == S` and `FACT == |ΛᵀΛ|^(1/n)` to machine
    precision, the two are the same function, computed twice.
    """
    exe = tmp_path / "as197d"
    r = subprocess.run(
        ["gfortran", "-std=legacy", "-fdefault-real-8", "-O2", "-o", str(exe),
         os.path.join(_HERE, "fortran", "as197.f"),
         os.path.join(_HERE, "fortran", "as197_driver.f")],
        capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip("AS 197 did not compile here")

    w, phi, theta = _fit(rel)
    stdin = f"{len(phi)} {len(theta)} {len(w)} -0.001\n"
    if len(phi):
        stdin += " ".join(repr(float(v)) for v in phi) + "\n"
    if len(theta):
        stdin += " ".join(repr(float(v)) for v in theta) + "\n"
    stdin += "\n".join(repr(float(v)) for v in w) + "\n"
    out = subprocess.run([str(exe)], input=stdin, capture_output=True,
                         text=True).stdout
    mel = {k: float(v) for k, v in
           (line.split(None, 1) for line in out.strip().split("\n"))}

    _ll, f1, f2, _a, _if = _elf(w, phi, theta)
    assert mel["SUMSQ"] == pytest.approx(f1, rel=1e-10), "eq. (3) disagrees"
    assert mel["FACT"] == pytest.approx(f2, rel=1e-12), "eq. (4) disagrees"


# ── 3. the ten steps, traceable in both engines ────────────────────────────

def test_the_c_numbers_the_steps_the_working_paper_numbers():
    """`elfvarma.c` blocks [1]…[9] against the ten steps of WP 9316, pp. 7-8.

    Step (10) — the quadratic form S = ηᵀη − λᵀλ — has no block of its own in
    the C: it is computed where λ is, which is why only nine are asserted.
    """
    src = open(os.path.join(_ROOT, "csrc", "internal", "elfvarma.c"),
               encoding="utf-8").read()
    esperado = {
        1: "Cholesky factor of qq",
        2: "autocovariances",
        3: "Cholesky factor of v1 * omega * v1",
        4: "xi(k)",
        5: "vector eta",
        6: "M'h",
        7: "H'H",
        8: "I+M'H'HM",
        9: "lambda",
    }
    for k, marca in esperado.items():
        bloque = f"/* [{k}]:"
        assert bloque in src, f"step [{k}] is no longer marked in elfvarma.c"
        i = src.index(bloque)
        assert marca in src[i:i + 200], (
            f"step [{k}] no longer says '{marca}' — the numbering and the "
            f"working paper's steps have drifted apart")


def test_the_python_port_cites_the_equations_it_implements():
    """(a)…(k) with AS 311 equation numbers. A translation that stops saying
    which equation it implements is a translation nobody can check."""
    src = open(os.path.join(_ROOT, "src", "fue", "elfvarma.py"),
               encoding="utf-8").read()
    for paso in "abcdefghij":
        assert f"({paso})" in src or f"[{paso}]" in src, (
            f"step ({paso}) is not marked in elfvarma.py")
    for eq in ("AS311 eq.2", "AS311 eq.3", "9316 eq.2.15", "9316 eq.2.16"):
        assert eq in src, f"the citation of {eq} has been lost"
