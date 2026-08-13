"""The C likelihood core must remain José Alberto Mauricio's, verbatim.

`csrc/internal/{elfvarma,usmelard,drvmlest}.c` are Mauricio's DRVUS sources
(`SRC/drvus/src/`), carrying the implementations of

    Mauricio (1997) Algorithm AS 311   — exact VARMA likelihood
    Mauricio (1995) JASA 90, 282-291   — the estimation method
    Melard  (1984) Algorithm AS 197    — fast ARMA likelihood

with a GPL header, an encoding change and two mechanical edits on top. That is
the central claim of `docs/PROVENANCE.md`, and until now it was a claim: a
property nobody checked, which is exactly how the synchronisation of the C copies
broke in `drvarma` (`drvarma/bugs/BUG-0002` — `qnewtopt.c` drifted 17 lines and
nothing noticed).

**This test does not count differing lines.** A count would pass if twenty-three
benign lines were replaced by twenty-three malicious ones. Every differing line
must match a DECLARED exception; anything else fails and is printed.

The one functional change in the whole core has its own exception, written out in
full, so that changing it further breaks this test rather than sliding through:

    elfvarma.c:513   eigenqr(...)  →  if (n>1) {gsl_eigenqr(...);}

`nlatools.c` is deliberately out of scope: it was rewritten around GSL (764 lines
against Mauricio's 1355) and is not a verbatim copy. See `test_nlatools_is_not_claimed_verbatim`.
"""
import difflib
import os
import re

import pytest

_DRVUS = os.path.expanduser("~/Dropbox/SRC/drvus/src")
_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_INTERNAL = os.path.join(_HERE, "csrc", "internal")

#: The modules claimed to be verbatim. `nlatools.c` is NOT among them.
_VERBATIM = ["elfvarma", "usmelard", "drvmlest"]

pytestmark = pytest.mark.skipif(
    not os.path.isdir(_DRVUS),
    reason="DRVUS sources not present (expected at ~/Dropbox/SRC/drvus/src)")


# ── the declared exceptions ────────────────────────────────────────────────
#
# Each entry is (name, regex). A differing line is accepted only if it matches
# one of them. Keep the names meaningful: they are what a failure prints.

_LICENCE = re.compile(
    r"^/\*\s*("
    r"This program is free software|"
    r"it under the terms of the GNU|"
    r"the Free Software Foundation; either version|"
    r"\(at your option\) any later version|"
    r"This program is distributed in the hope|"
    r"but WITHOUT ANY WARRANTY|"
    r"MERCHANTABILITY or FITNESS|"
    r"GNU General Public License for more details|"
    r"You should have received a copy|"
    r"along with this program; if not|"
    r"Inc\., 51 Franklin Street"
    r")")

_EXCEPTIONS = [
    # The GPL header added to Mauricio's files when the suite was licensed.
    ("gpl-header", _LICENCE.match),
    ("gpl-header-blank", lambda l: l.strip() in ("/*", "/*  */") or
                                   re.fullmatch(r"/\*\s*\*/", l.strip()) is not None),
    # The copyright line itself: same text, different encoding of "José".
    ("copyright-encoding",
     lambda l: "Copyright (C)" in l and "Alberto Mauricio" in l),
    # drvus.h split into fue.h + nlatools.h.
    ("include-rename",
     lambda l: re.match(r'^#include "(drvus|fue|nlatools)\.h"', l) is not None),
    # A blank line dropped next to the includes, and the closing comment banner.
    ("blank-line", lambda l: l.strip() == ""),
    ("banner-comment", lambda l: re.fullmatch(r"/\*+/", l.strip()) is not None),
    # Another encoding artefact: "nº of parameters" in drvmlest.c.
    ("text-encoding",
     lambda l: "of parameters to estimate" in l),
    # ── THE ONE FUNCTIONAL CHANGE, written out in full ────────────────────
    # Numerical Recipes' eigenvalue routine replaced by GSL's, guarded n>1.
    ("gsl-eigenqr-old",
     lambda l: l.strip().startswith("eigenqr( a, n, wr, wi );")),
    ("gsl-eigenqr-new",
     lambda l: l.strip().startswith("if ( n>1 ) {gsl_eigenqr( a, n, wr, wi );}")),
]


def _lines(path, encoding):
    """Normalised lines: CRLF stripped, decoded, no trailing whitespace."""
    with open(path, "rb") as fh:
        raw = fh.read().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return [l.rstrip() for l in raw.decode(encoding, errors="replace").split("\n")]


def _classify(line):
    """The exception a differing line falls under, or None."""
    for name, matches in _EXCEPTIONS:
        try:
            if matches(line):
                return name
        except Exception:                       # a bad regex must not pass a line
            continue
    return None


def _differing(module):
    """Every line that differs, as (side, text)."""
    a = _lines(os.path.join(_DRVUS, f"{module}.c"), "latin-1")
    b = _lines(os.path.join(_INTERNAL, f"{module}.c"), "utf-8")
    out = []
    for line in difflib.unified_diff(a, b, n=0, lineterm=""):
        if line.startswith(("---", "+++", "@@")):
            continue
        if line[:1] in "-+":
            out.append((line[0], line[1:]))
    return out


# ── the invariant ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("module", _VERBATIM)
def test_every_difference_is_a_declared_exception(module):
    """The claim of docs/PROVENANCE.md §2, as an invariant."""
    undeclared = [(side, text) for side, text in _differing(module)
                  if _classify(text) is None]
    assert not undeclared, (
        f"\n{module}.c has changed against Mauricio's DRVUS source in ways this "
        f"test does not know about.\n"
        f"Either the change is legitimate — and belongs in _EXCEPTIONS with a "
        f"name and a reason, and in docs/PROVENANCE.md §2 — or the core is no "
        f"longer verbatim and the documentation is now false.\n\n"
        + "\n".join(f"  {s} {t}" for s, t in undeclared[:20]))


@pytest.mark.parametrize("module", _VERBATIM)
def test_the_source_files_are_present_on_both_sides(module):
    assert os.path.exists(os.path.join(_DRVUS, f"{module}.c"))
    assert os.path.exists(os.path.join(_INTERNAL, f"{module}.c"))


def test_the_single_functional_change_is_still_the_only_one():
    """Two of the three modules must differ ONLY in licence and encoding.

    If a functional edit ever appears in `usmelard.c` (AS 197) or `drvmlest.c`,
    it will show up here even if someone adds it to _EXCEPTIONS for elfvarma.
    """
    functional = {"gsl-eigenqr-old", "gsl-eigenqr-new"}
    for module in ("usmelard", "drvmlest"):
        kinds = {_classify(t) for _s, t in _differing(module)}
        assert not (kinds & functional), (
            f"{module}.c now carries the eigenqr change, which belonged only to "
            f"elfvarma.c")

    kinds = {_classify(t) for _s, t in _differing("elfvarma")}
    assert functional <= kinds, (
        "elfvarma.c no longer shows the GSL eigenvalue substitution. If it was "
        "reverted or rewritten, docs/PROVENANCE.md §2 needs updating.")


def test_nlatools_is_not_claimed_verbatim():
    """The exclusion is deliberate and is recorded here so it stays deliberate.

    `nlatools.c` was rewritten around GSL — 764 lines against Mauricio's 1355 —
    so it is not a copy and this test suite makes no claim about it.
    """
    assert "nlatools" not in _VERBATIM
    a = _lines(os.path.join(_DRVUS, "nlatools.c"), "latin-1")
    b = _lines(os.path.join(_INTERNAL, "nlatools.c"), "utf-8")
    assert len(b) < len(a), "nlatools.c is no longer the shorter rewrite"


# ── the optimizer, against its own ancestor ────────────────────────────────
#
# `qnewtopt.c` is raxopt: Mauricio (1995) JASA 90, 282-291. It was outside the
# invariant above until BUG-0012, which is precisely the file that then had to
# be edited — and `drvarma/bugs/BUG-0002` is the record of what happens when
# this copy drifts unwatched (17 lines, nobody noticed). Its ancestor is
# `fue-1.13.1/src/qnewtopt.c`, not DRVUS: the fue line already wrote its
# progress to `outputv` rather than stdout.

_FUE_1131 = os.path.expanduser("~/Dropbox/SRC/atws/fue/fue-1.13.1/src")

_OPT_EXCEPTIONS = [
    # "José" — the same text in two encodings.
    ("copyright-encoding",
     lambda l: "Copyright (C)" in l and "Alberto Mauricio" in l),
    # Progress lines redirected: the binding has no stdout to write to.
    ("printf-to-outputv",
     lambda l: '"%4d F: %0.10f' in l),
    # BUG-0012: recording what raxopt already computed. No criterion changes.
    ("termcode-recording",
     lambda l: re.match(r"^\s*((int|double)\s+)?qn_last_(termcode|nit|gnorm)\s*=",
                        l) is not None),
]


def _code_only(path, encoding):
    """The file with comments and blank lines removed.

    Comparing code rather than prose is deliberate: a comment cannot change what
    raxopt does, and the note explaining the BUG-0012 recording is long. What
    must not move is the executable text.
    """
    src = "\n".join(_lines(path, encoding))
    src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    src = re.sub(r"//[^\n]*", "", src)
    return [l.rstrip() for l in src.split("\n") if l.strip()]


def _classify_opt(line):
    for name, matches in _OPT_EXCEPTIONS:
        try:
            if matches(line):
                return name
        except Exception:
            continue
    return None


@pytest.mark.skipif(not os.path.isdir(_FUE_1131),
                    reason="fue-1.13.1 sources not present")
def test_the_optimizer_is_still_mauricios():
    """raxopt must keep its own rules: BUG-0012 records, it does not decide."""
    a = _code_only(os.path.join(_FUE_1131, "qnewtopt.c"), "latin-1")
    b = _code_only(os.path.join(_INTERNAL, "qnewtopt.c"), "utf-8")
    undeclared = []
    for line in difflib.unified_diff(a, b, n=0, lineterm=""):
        if line.startswith(("---", "+++", "@@")) or line[:1] not in "-+":
            continue
        if _classify_opt(line[1:]) is None:
            undeclared.append(line)
    assert not undeclared, (
        "\nqnewtopt.c has changed against fue-1.13.1 in ways this test does not "
        "know about. raxopt is published work (Mauricio 1995, JASA 90, 282-291):"
        " a change to its stopping criteria is a STUDY, not a bug fix.\n\n"
        + "\n".join(f"  {l}" for l in undeclared[:20]))


@pytest.mark.skipif(not os.path.isdir(_FUE_1131),
                    reason="fue-1.13.1 sources not present")
def test_the_stopping_criteria_themselves_are_untouched():
    """The two routines that DECIDE when to stop, byte for byte.

    `umstop0` (the first-iterate test) and `umstop` (the rest) are where the
    gradient and step criteria live. BUG-0012 is about reporting their verdict,
    so their text must be identical — and if a future study does change them,
    it will have to change this test and say so.
    """
    def routine(path, enc, name):
        src = "\n".join(_lines(path, enc))
        m = re.search(r"^\s*int\s+" + name + r"\s*\(", src, re.M)
        assert m, f"{name} not found in {path}"
        i, depth, started = m.start(), 0, False
        for j in range(i, len(src)):
            if src[j] == "{":
                depth += 1
                started = True
            elif src[j] == "}":
                depth -= 1
                if started and depth == 0:
                    return src[i:j + 1]
        raise AssertionError(f"{name}: unbalanced braces")

    for name in ("umstop0", "umstop"):
        old = routine(os.path.join(_FUE_1131, "qnewtopt.c"), "latin-1", name)
        new = routine(os.path.join(_INTERNAL, "qnewtopt.c"), "utf-8", name)
        assert old == new, f"{name}() differs — a stopping criterion changed"


def test_the_provenance_document_says_the_same_thing():
    """A document that drifts from the test it describes is worse than none."""
    doc = os.path.join(_HERE, "docs", "PROVENANCE.md")
    if not os.path.exists(doc):
        pytest.skip("docs/PROVENANCE.md not present")
    text = open(doc, encoding="utf-8").read()
    for module in _VERBATIM:
        assert f"`{module}.c`" in text, f"{module}.c is not listed in PROVENANCE §2"
    assert "gsl_eigenqr" in text, "the one functional change is not documented"
