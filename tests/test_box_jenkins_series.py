"""The Box & Jenkins series, against Mauricio's own C.

`drvus-source/1.2.01/drvus/src/Box_y_Jenkings/` carries Series A–E of Box &
Jenkins (1976) as `.inp` files with the `.out` that DRVUS produced around 2001.
It is the canonical benchmark of the field, fitted by the author of the engine,
preserved for twenty-five years — and until 12 August 2026 none of it loaded,
because of `bugs/BUG-0010` (Latin-1) and `bugs/BUG-0011` (the bands section that
did not exist yet).

Measured once both were fixed — log-likelihood, `fue` against DRVUS:

    a2   Series A, IMA(0,1,1)     -53.508690   1.5e-11
    b    Series B, IMA(0,1,1)   -1249.974933   4.0e-11
    c    Series C, ARI(1,1,0)     131.668147   4.7e-12
    c2   Series C, IMA(0,2,1)     123.399306   3.9e-11
    d    Series D, AR(1)          -67.751685   2.6e-10
    d1   Series D, IMA(0,1,1)     -76.691867   2.8e-11
    e1   Series E, AR(1)         -414.617409   3.4e-12
    e2   Series E, AR(1)         -412.494817   8.2e-12

and one that does **not** agree, `a1`, which is `bugs/BUG-0012` and is asserted
here as a known divergence rather than hidden by a wide tolerance.

The files live outside the repository, so these tests skip when DRVUS is absent —
the same pattern as the thesis and precipitation fixtures elsewhere in the suite.
"""
import os
import re

import pytest

_BJ = os.path.expanduser(
    "~/Dropbox/SRC/drvus-source/1.2.01/drvus/src/Box_y_Jenkings")

pytestmark = pytest.mark.skipif(not os.path.isdir(_BJ),
                                reason="DRVUS Box-Jenkins fixtures not present")

#: case -> (relative .inp, DRVUS's logelf from the .out beside it)
_CASES = {
    "a2": ("SeriesA/a2.inp", -53.5086902793),
    "b":  ("SeriesB/b.inp", -1249.974933),
    "c":  ("SeriesC/c.inp", 131.668147),
    "c2": ("SeriesC/c2.inp", 123.399306),
    "d":  ("SeriesD/d.inp", -67.751685),
    "d1": ("SeriesD/d1.inp", -76.691867),
    "e1": ("SeriesE/e1.inp", -414.617409),
    "e2": ("SeriesE/e2.inp", -412.494817),
}


def _drvus_loglik(inp_path):
    """`logelf:` from the .out that ships beside the .inp."""
    out = inp_path[:-4] + ".out"
    if not os.path.exists(out):
        return None
    txt = open(out, encoding="latin-1").read().replace("\r\n", "\n")
    m = re.search(r"logelf:\s*([-\d.]+)", txt)
    return float(m.group(1)) if m else None


def _fit(rel):
    import fue

    ts, m = fue.load(os.path.join(_BJ, rel))
    m.fit()
    return ts, m


# ── the two barriers, as regressions ───────────────────────────────────────

def test_the_files_load_untouched():
    """BUG-0010 and BUG-0011: Latin-1 and the missing bands section. These files
    are read as they were written in 1996 — no conversion, no hand-editing."""
    import fue

    for rel, _ in _CASES.values():
        ts, m = fue.load(os.path.join(_BJ, rel))
        assert ts.nobs > 0


def test_the_annual_header_with_two_fields_is_read():
    """Series E is annual and its header is `100 1770` — two fields, the DRVUS
    form. The parser knew the three- and four-field forms only."""
    ts, _m = _fit("SeriesE/e1.inp")
    assert ts.nobs == 100
    assert ts.start[0] == 1770
    assert ts.freq == 1


# ── the benchmark ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("case", sorted(_CASES))
def test_fue_reproduces_drvus(case):
    """The port against Mauricio's own C, on the canonical series."""
    rel, expected = _CASES[case]
    _ts, m = _fit(rel)
    assert m.loglik == pytest.approx(expected, abs=1e-4), (
        f"{case}: fue {m.loglik:.6f} vs DRVUS {expected:.6f}")


@pytest.mark.parametrize("case", sorted(_CASES))
def test_the_reference_really_comes_from_the_out_file(case):
    """The expected values above are transcribed. This checks the transcription
    against the file, so a typo here fails rather than passing quietly."""
    rel, expected = _CASES[case]
    ref = _drvus_loglik(os.path.join(_BJ, rel))
    if ref is None:
        pytest.skip(f"{case}: no logelf in the .out")
    assert ref == pytest.approx(expected, abs=1e-4)


@pytest.mark.parametrize("case", sorted(_CASES))
def test_the_termination_criterion_matches_drvus(case):
    """Not just the destination — how each run ended.

    DRVUS recorded its criterion, iteration count and gradient norm in the
    `.out`, and since BUG-0012 the port returns the same three. The criterion
    must agree on all eight; the iteration count is deliberately NOT asserted,
    because it drifts by a few (50 against 47 on `d`, 2 against 3 on `c`)
    between a 2001 binary and this one. That drift is itself the evidence for
    what happens on `a1`: two slightly different paths, which matters only where
    the valley is flat.
    """
    rel, _expected = _CASES[case]
    out = os.path.join(_BJ, rel)[:-4] + ".out"
    txt = open(out, encoding="latin-1").read()
    assert "GRADIENT STOPPING CRITERIUM" in txt, (
        f"{case}: DRVUS did not stop on the gradient — reread the .out")

    _ts, m = _fit(rel)
    assert m._result.termcode == 1, (
        f"{case}: DRVUS stopped on the gradient and fue stopped by "
        f"termcode {m._result.termcode} (|g|={m._result.gnorm:.3g})")
    assert m._result.gnorm == pytest.approx(0.0, abs=1e-4)


# ── the one that disagrees, asserted rather than hidden ────────────────────

def test_series_a_arma11_still_stops_on_the_boundary():
    """BUG-0012, resolved — and kept as a two-sided marker.

    On `a1` fue stops at φ→1 by the STEP criterion, 6.86 in log-likelihood
    below the `.out` DRVUS left in 2001. **So does DRVUS itself when compiled
    today**: same source, `gcc -O2`, same 23 iterations, same |g|=0.01009, same
    −57.60386. Nothing in the port is implicated.

    The 2001 binary was 32-bit x86, whose FPU carried intermediates in 80-bit
    registers; x86-64 uses SSE2 and works in 64. Rebuilt today with
    `-O0 -mfpmath=387` the original C reproduces its own 2001 trace iteration by
    iteration and lands on −50.745092 again. Eleven extra bits of mantissa
    decide, in a valley this flat, whether steptol fires at iteration 23.

    The marker stays because the *behaviour* must not drift unnoticed in either
    direction — and because the day this platform changes its floating point,
    this test is what will say so.
    """
    import warnings

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _ts, m = _fit("SeriesA/a1.inp")
    drvus = _drvus_loglik(os.path.join(_BJ, "SeriesA/a1.inp"))
    assert drvus is not None

    assert m.ar[0][0] > 0.999, "fue used to stop on the AR boundary — did it escape?"
    # And it must say so: the stop is the STEP criterion with a live gradient,
    # not a maximum. Silence here is the half of BUG-0012 that is now fixed.
    assert m._result.termcode == 2
    assert m._result.converged is False
    assert any(issubclass(w.category, RuntimeWarning) for w in caught)
    assert m.loglik < drvus - 1.0, (
        f"fue {m.loglik:.6f} is no longer well below DRVUS {drvus:.6f} — "
        f"BUG-0012 may be fixed; check and update this test")


def test_the_optimum_drvus_reached_in_2001_is_the_real_one():
    """And fue reaches it — the half that matters for the estimate.

    Seed mu at the sample mean instead of the stale 2.5 the file carries, and
    fue lands on Mauricio's 2001 answer. `statsmodels` — different authors,
    different algorithm, no shared code with any of this — lands on the same
    point, which is what makes it the optimum rather than a third opinion:

        statsmodels   phi=0.908685  theta=0.575841  logL=-50.745092
        DRVUS 2001    phi=0.908683  theta=0.575839  logL=-50.745092
        fue, mu0=17   phi=0.908685  theta=0.575841  logL=-50.745092
    """
    sm = pytest.importorskip("statsmodels.tsa.arima.model")
    import numpy as np

    ts, m = fue_load_a1 = _fit_with_mu("SeriesA/a1.inp", 17.0)
    assert m._result.termcode == 1
    assert m.loglik == pytest.approx(-50.745092, abs=1e-5)

    r = sm.ARIMA(np.array(ts.data), order=(1, 0, 1),
                 trend="c").fit(method="innovations_mle")
    assert m.ar[0][0] == pytest.approx(r.arparams[0], abs=1e-5)
    assert m.ma[0][0] == pytest.approx(-r.maparams[0], abs=1e-5)
    assert m.loglik == pytest.approx(r.llf, abs=1e-4)


def _fit_with_mu(rel, mu0):
    import fue

    ts, m = fue.load(os.path.join(_BJ, rel))
    m.mu0 = mu0
    m.fit()
    return ts, m
