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


# ── the one that disagrees, asserted rather than hidden ────────────────────

def test_series_a_arma11_still_stops_on_the_boundary():
    """BUG-0012, as a two-sided marker.

    On `a1` — ARMA(1,1) on the level of Series A — fue stops at φ→1, where the
    model degenerates into the IMA(1,1) of `a2`, while DRVUS escapes to
    φ=0.9087, θ=0.5758 from the same starting values: 6.86 better in
    log-likelihood, and the published estimate.

    If this ever starts agreeing, the assertion below fails and that is the
    point: a known defect must not become a silent one in either direction.
    """
    _ts, m = _fit("SeriesA/a1.inp")
    drvus = _drvus_loglik(os.path.join(_BJ, "SeriesA/a1.inp"))
    assert drvus is not None

    assert m.ar[0][0] > 0.999, "fue used to stop on the AR boundary — did it escape?"
    assert m.loglik < drvus - 1.0, (
        f"fue {m.loglik:.6f} is no longer well below DRVUS {drvus:.6f} — "
        f"BUG-0012 may be fixed; check and update this test")
