"""fue's likelihood against Melard's own FORTRAN, as printed in AS 197.

This is the check `docs/PROVENANCE.md` §6 listed as missing, and it is a
different kind of witness from everything else in the suite. `usmelard.c` is
Mauricio's C, the Python port is a translation of it, and the reference `.out`
files were produced by it: they all descend from one implementation. The
listing in

    G. Melard (1984) "Algorithm AS 197: A Fast Algorithm for the Exact
    Likelihood of Autoregressive-Moving Average Models",
    Applied Statistics 33(1), 104-114

is the algorithm as its author published it, printed in full on pages 110-114
and transcribed in `tests/fortran/as197.f`. Agreement with it tests the
**algorithm**, not the port.

Measured 13 August 2026 on the nine Box-Jenkins specifications, log-likelihood
of AS 197's FORTRAN against fue's engine:

    a1  ARMA(1,1)     -50.745091555   -50.745091514   -4.0e-08
    a2  IMA(0,1,1)    -53.508690319   -53.508690279   -4.0e-08
    b   IMA(0,1,1)  -1249.974932650 -1249.974932574   -7.5e-08
    c   ARI(1,1,0)    131.668146734   131.668146780   -4.6e-08
    c2  IMA(0,2,2)    123.399306445   123.399306491   -4.6e-08
    d   AR(1)         -67.751684568   -67.751684505   -6.3e-08
    d1  IMA(0,1,1)    -76.691867132   -76.691867069   -6.3e-08
    e1  AR(2)        -414.617408735  -414.617408715   -2.1e-08
    e2  AR(3)        -412.494817402  -412.494817381   -2.1e-08

FLIKAM returns `FACT` and `SUMSQ`; the log-likelihood is built here, from the
paper's own definition — likelihood ∝ FACT·SUMSQ with FACT = (Π hₜ²)^(1/n) and
SUMSQ = Σ(âₜ/hₜ)² — concentrating σ̂² = SUMSQ/n:

    log L = −(n/2)·[ log 2π + 1 + log(SUMSQ/n) + log FACT ]

**A note for whoever maintains this.** The first transcription put label `170`
one line too early, on `IF (MQ .LE. 0) GOTO 300` instead of on
`IF (R .LE. EPSIL1) GOTO 400`. Everything with an MA part still agreed to 1e-8;
every model with an AR part was wrong by 2 to 7 in log-likelihood, because the
routine switched to quick recursions at t=1 instead of at t=p+1. The numbers
caught it — a structural read of the listing would not have. That is the
argument for running the paper's code rather than admiring it.
"""
import math
import os
import shutil
import subprocess

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_FORTRAN = os.path.join(_HERE, "fortran")
_BJ = os.path.expanduser(
    "~/Dropbox/SRC/drvus-source/1.2.01/drvus/src/Box_y_Jenkings")

#: case -> (relative .inp, expected IFAULT from FLIKAM)
#:
#: IFAULT is not an error here. The paper documents `-m` as "not a failure:
#: indicates that quick recursions took place from t = m". A pure AR(p) needs
#: the exact treatment only for its first p observations, so switching at
#: t = p+1 is the exact likelihood, and the codes below say exactly that:
#: AR(1) -> -2, AR(2) -> -3, AR(3) -> -4. Models with an MA part run the
#: exact recursion throughout and return 0.
_CASES = {
    "a1": ("SeriesA/a1.inp", 0),
    "a2": ("SeriesA/a2.inp", 0),
    "b":  ("SeriesB/b.inp", 0),
    "c":  ("SeriesC/c.inp", -2),
    "c2": ("SeriesC/c2.inp", 0),
    "d":  ("SeriesD/d.inp", -2),
    "d1": ("SeriesD/d1.inp", 0),
    "e1": ("SeriesE/e1.inp", -3),
    "e2": ("SeriesE/e2.inp", -4),
}

pytestmark = [
    pytest.mark.skipif(shutil.which("gfortran") is None,
                       reason="gfortran not available"),
    pytest.mark.skipif(not os.path.isdir(_BJ),
                       reason="DRVUS Box-Jenkins fixtures not present"),
]


@pytest.fixture(scope="module")
def as197(tmp_path_factory):
    """Compile the transcribed listing.

    `-fdefault-real-8` is the paper's own Precision note (page 108) applied
    without editing the source: "all the real variables should be replaced by
    double precision variables". `-std=legacy` is for the fixed form and the
    computed GOTOs.
    """
    out = tmp_path_factory.mktemp("as197") / "as197d"
    cmd = ["gfortran", "-std=legacy", "-fdefault-real-8", "-O2", "-o", str(out),
           os.path.join(_FORTRAN, "as197.f"),
           os.path.join(_FORTRAN, "as197_driver.f")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip(f"AS 197 did not compile here:\n{r.stderr[:400]}")
    return str(out)


def _flikam(exe, w, phi, theta, toler=-1e-3):
    """Run FLIKAM. `toler` negative asks for the exact likelihood — the
    paper's contract for TOLER, and the reason fue signs `xitol` the way it
    does (`fue_api.c:951-956`)."""
    stdin = f"{len(phi)} {len(theta)} {len(w)} {toler}\n"
    if phi:
        stdin += " ".join(repr(v) for v in phi) + "\n"
    if theta:
        stdin += " ".join(repr(v) for v in theta) + "\n"
    stdin += "\n".join(repr(v) for v in w) + "\n"
    out = subprocess.run([exe], input=stdin, capture_output=True,
                         text=True).stdout
    return {k: float(v) for k, v in
            (line.split(None, 1) for line in out.strip().split("\n"))}


def _loglik(res, n):
    return -0.5 * n * (math.log(2 * math.pi) + 1.0
                       + math.log(res["SUMSQ"] / n) + math.log(res["FACT"]))


def _fue_fit(rel):
    """Fit with fue, seeding mu at the mean of the differenced variable.

    The seed matters — see bugs/BUG-0012 — and `a1.inp` ships a stale one.
    """
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
    phi = [c for f in m.ar for c in f]
    theta = [c for f in m.ma for c in f]
    return m, list(w), phi, theta


# ── the algorithm, against its own publication ─────────────────────────────

@pytest.mark.parametrize("case", sorted(_CASES))
def test_fue_agrees_with_melards_published_fortran(case, as197):
    rel, _ifault = _CASES[case]
    m, w, phi, theta = _fue_fit(rel)
    res = _flikam(as197, w, phi, theta)
    assert res["IFAULT"] <= 0, f"{case}: FLIKAM reported fault {res['IFAULT']:.0f}"

    ll = _loglik(res, len(w))
    assert ll == pytest.approx(m.loglik, abs=1e-6), (
        f"{case}: AS 197 {ll:.9f} vs fue {m.loglik:.9f}")


@pytest.mark.parametrize("case", sorted(_CASES))
def test_the_switching_point_is_the_one_the_paper_documents(case, as197):
    """`IFAULT = -m` means quick recursions took place from t = m.

    Asserted because it is the property that caught a transcription error: with
    the label one line out of place every AR model switched at t=1, which is
    the conditional likelihood, and the log-likelihoods moved by 2 to 7.
    """
    rel, want = _CASES[case]
    _m, w, phi, theta = _fue_fit(rel)
    res = _flikam(as197, w, phi, theta)
    assert int(res["IFAULT"]) == want


def test_a_positive_toler_gives_the_approximate_likelihood(as197):
    """The other half of the TOLER contract, and it must be visible.

    AS 197: TOLER "should be negative if the exact likelihood is desired.
    Otherwise, switching to approximate recursions occurs when h²ₜ < 1 + δ."
    With δ > 0 the routine switches early and the likelihood changes — which is
    why `fue` signs `xitol` negative for `eml=True` and positive otherwise.
    """
    _m, w, phi, theta = _fue_fit("SeriesA/a2.inp")
    exact = _flikam(as197, w, phi, theta, toler=-1e-3)
    approx = _flikam(as197, w, phi, theta, toler=1e-1)

    assert int(exact["IFAULT"]) == 0
    assert int(approx["IFAULT"]) < 0, "a positive TOLER did not switch at all"
    assert _loglik(exact, len(w)) != pytest.approx(_loglik(approx, len(w)),
                                                   abs=1e-9)


def test_the_transcription_is_the_paper_and_not_our_code():
    """A guard on the file itself, so it cannot quietly become our code.

    Three markers from the printed listing: the header line the paper carries,
    the FORTRAN 66 spellings that gfortran does not even know (they are shimmed
    in the driver rather than edited here), and the label placement that the
    first transcription got wrong.
    """
    src = open(os.path.join(_FORTRAN, "as197.f")).read()
    assert "ALGORITHM AS 197  APPL. STATIST. (1984) VOL.33, NO.1" in src
    assert "MAXO(" in src and "MINO(" in src, (
        "MAXO/MINO were replaced by MAX0/MIN0 — that is editing the paper")
    assert "      IF (MQ .LE. 0) GOTO 300\n  170 IF (R .LE. EPSIL1) GOTO 400" in src, (
        "label 170 has moved; page 111 puts it on the EPSIL1 line, and the "
        "difference is the exact likelihood of every AR model")
