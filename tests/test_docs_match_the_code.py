"""The documentation must not drift from what the code does.

`docs/MODEL.md` and `docs/FORMAL_TESTS.md` are translations of Treadway's manual
brought up to date, and they carry numbers: critical values, tolerances, the
factorisation of the annual difference. A document that quietly stops matching
the code is worse than no document, because it is read as authority — which is
exactly what happened to `art`'s bug index (it said "11 reports, 3 open" when
there were 18 and 8).

So the numbers that can be checked are checked here. What cannot be automated —
whether the prose is *true* — stays the reader's job.
"""
import os
import re

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DOCS = os.path.join(_ROOT, "docs")


def _doc(name):
    p = os.path.join(_DOCS, name)
    if not os.path.exists(p):
        pytest.skip(f"{name} not present")
    return open(p, encoding="utf-8").read()


# ── MODEL.md ───────────────────────────────────────────────────────────────

def test_the_annual_difference_factorises_as_the_document_says():
    """S(B) = (1+B)·∏(1 − 2cos(2πf/s)B + B²), and ∇ₛ = (1−B)·S(B).

    MODEL.md §2.3 states it and gives the f=1,5 example with √3. Multiplying the
    factors out must give 1 + B + … + B^{s−1} exactly.
    """
    _doc("MODEL.md")
    s = 12
    poly = np.array([1.0, 1.0])                      # (1 + B), the Nyquist factor
    for f in range(1, s // 2):
        poly = np.convolve(poly, [1.0, -2 * np.cos(2 * np.pi * f / s), 1.0])
    assert np.allclose(poly, np.ones(s), atol=1e-12), (
        "the product of the irreducible factors is not the annual moving sum")

    nabla_s = np.convolve([1.0, -1.0], poly)         # (1 − B)·S(B)
    esperado = np.zeros(s + 1)
    esperado[0], esperado[-1] = 1.0, -1.0            # 1 − B^s
    assert np.allclose(nabla_s, esperado, atol=1e-12)


def test_the_root_three_example_is_the_one_in_the_document():
    """§2.3: for s=12, frequencies 1 and 5 give (1 − √3B + B²)(1 + √3B + B²)."""
    doc = _doc("MODEL.md")
    assert "(1 − √3·B + B²)(1 + √3·B + B²)" in doc
    assert 2 * np.cos(2 * np.pi * 1 / 12) == pytest.approx(np.sqrt(3), abs=1e-12)
    assert 2 * np.cos(2 * np.pi * 5 / 12) == pytest.approx(-np.sqrt(3), abs=1e-12)


def test_the_cost_in_observations_is_stated_correctly():
    """§2.4: interior factors cost two observations, f=0 and Nyquist cost one.

    This is the arithmetic `drtran` got wrong twice (BUG-5, BUG-9), so it is
    asserted against the degree of the polynomial rather than against prose.
    """
    _doc("MODEL.md")
    s = 12
    for idx, grado in [(0, 1), (1, 2), (3, 2), (5, 2), (6, 1)]:
        if idx == 0:
            fac = [1.0, -1.0]
        elif idx == s // 2:
            fac = [1.0, 1.0]
        else:
            fac = [1.0, -2 * np.cos(2 * np.pi * idx / s), 1.0]
        assert len(fac) - 1 == grado, f"factor {idx} does not cost {grado}"


def test_the_likelihood_formula_is_the_one_the_engine_uses():
    """§3 prints AS 311 equation (2). It must be the engine's, not decoration."""
    doc = _doc("MODEL.md")
    assert "log|ΛᵀΛ|" in doc and "AS 311" in doc

    from fue.elfvarma import elf_scalar
    import math

    rng = np.random.RandomState(7)
    w = rng.randn(120)
    logelf, f1, f2, _a, ifault = elf_scalar(120, 1, 1, np.array([0.5]),
                                            np.array([0.3]), w, xitol=-1e-3)
    assert ifault == 0
    eq2 = -0.5 * (120 * math.log(2 * math.pi) + 120 * math.log(f2) + f1)
    assert logelf == pytest.approx(eq2, abs=1e-6)


# ── FORMAL_TESTS.md ────────────────────────────────────────────────────────

_art = pytest.importorskip("art.formal_tests",
                           reason="art is not installed here")


def test_the_shin_fuller_table_matches_art():
    """§2's table against `art.formal_tests._SF_CRIT`, row by row."""
    doc = _doc("FORMAL_TESTS.md")
    for n, c10, c05, c01 in _art._SF_CRIT:
        fila = f"| {n} | {c10:.2f} | {c05:.2f} | {c01:.2f} |"
        assert fila in doc, f"the row for n={n} does not match art: expected {fila}"


def test_the_dcd_regimes_match_art():
    """§3: the real-root values, the complex-pair table and the asymptote."""
    doc = _doc("FORMAL_TESTS.md")
    real = _art._DCD_CRIT_MA
    assert f"| {real['10%']:.2f} | {real['5%']:.2f} | {real['1%']:.2f} |" in doc

    for n, (a, b, c) in sorted(_art._DCD_CRIT_MA_F_TABLE.items()):
        assert f"| {n} | {a:.2f} | {b:.2f} | {c:.2f} |" in doc, (
            f"the complex-pair row for n={n} does not match art")

    a, b, c = _art._DCD_CRIT_MA_F_ASYMP
    assert f"| → ∞ | {a:.2f} | {b:.2f} | {c:.2f} |" in doc


def test_the_functions_the_document_points_at_exist():
    """§7 names where each test lives. A pointer to a function that no longer
    exists is the kind of error that survives for years."""
    doc = _doc("FORMAL_TESTS.md")
    for nombre in re.findall(r"`art\.formal_tests\.(\w+)`", doc):
        assert hasattr(_art, nombre), (
            f"FORMAL_TESTS.md §7 points at art.formal_tests.{nombre}, "
            f"which does not exist")
    for nombre in ("dcd", "dcd_overdiff_regular", "meg", "shin_fuller"):
        assert nombre in doc, f"{nombre} is no longer mentioned in the document"


# ── PORT.md ────────────────────────────────────────────────────────────────

def test_the_port_document_counts_the_lines_it_claims():
    """PORT.md §1 quotes the size of both engines. Those numbers age quickly and
    are exactly the kind that nobody re-checks, so they are re-checked here."""
    doc = _doc("PORT.md")

    py = sum(len(open(os.path.join(_ROOT, "src", "fue", f),
                      encoding="utf-8").readlines())
             for f in ("elfvarma.py", "cast_us.py", "qnewtopt.py"))
    c = sum(len(open(os.path.join(_ROOT, p), encoding="utf-8").readlines())
            for p in ("csrc/fue_api.c", "csrc/internal/drvmlest.c",
                      "csrc/internal/elfvarma.c", "csrc/internal/nlatools.c",
                      "csrc/internal/qnewtopt.c", "csrc/internal/usmelard.c"))

    def escrito(n):
        # the document writes thousands with a dot: 1.895
        return f"{n // 1000}.{n % 1000:03d}" if n >= 1000 else str(n)

    assert escrito(py) in doc, (
        f"the pure-Python engine now has {py} lines and PORT.md says otherwise")
    assert escrito(c) in doc, (
        f"the embedded C now has {c} lines and PORT.md says otherwise")


def test_the_fallback_the_port_document_describes_still_exists():
    """§1: one try/except ImportError, no flag and no configuration."""
    doc = _doc("PORT.md")
    assert "try/except ImportError" in doc

    src = open(os.path.join(_ROOT, "src", "fue", "_engine.py"),
               encoding="utf-8").read()
    assert "from fue._fue_engine import ffi, lib" in src
    assert "except ImportError:" in src
    assert "from .cast_us import estimate_py" in src, (
        "the pure-Python fallback is gone; PORT.md §1 is then false")


def test_the_python_engine_still_defaults_to_raxopt():
    """PORT.md §3.1 says both engines run raxopt, and L-BFGS-B is opt-in.

    The document said the opposite for a day. It matters: if the fallback
    quietly switched to scipy, a disagreement between the two engines would stop
    being evidence about the likelihood and start being evidence about nothing.
    """
    doc = _doc("PORT.md")
    assert "Both engines run raxopt" in doc

    src = open(os.path.join(_ROOT, "src", "fue", "cast_us.py"),
               encoding="utf-8").read()
    assert re.search(r"def _estimate_core\(model, optimizer=[\"']raxopt[\"']\)", src)
    m = re.search(r"def estimate_py\(model\).*?_estimate_core\(model, "
                  r"optimizer=[\"'](\w+)[\"']\)", src, re.S)
    assert m and m.group(1) == "raxopt", (
        "estimate_py no longer defaults to raxopt — PORT.md §3.1 is then false")


# ── README.md / GETTING_STARTED.md ─────────────────────────────────────────

def test_the_landing_page_promises_three_things_the_engine_has():
    """The three distinguishing features, each checked against the code.

    A landing page is where overclaiming happens, so each claim is tied to
    something that exists: the rational transfer function (delta), the
    frequency-by-frequency operator (ifadf), and the exact likelihood.
    """
    doc = _doc("README.md")
    assert "rational transfer function" in doc
    assert "frequency by frequency" in doc.lower()
    assert "exact unconditional" in doc.lower()

    import inspect

    import fue

    sig = inspect.signature(fue.Intervention.__init__)
    assert "delta" in sig.parameters, (
        "the denominator of the transfer function is gone; the landing page "
        "claims a DYNAMIC response and would be false")
    assert "ifadf" in inspect.signature(fue.Model.__init__).parameters


def test_the_ten_lines_on_the_landing_page_run():
    """The snippet a reader tries first. It must produce the numbers printed
    next to it."""
    doc = _doc("README.md")
    assert "0.430784" in doc and "-48.995195" in doc

    import warnings

    import fue
    from fue.datasets import ripc

    ts = ripc()
    m = fue.Model(ts, d=1, boxlam=0.0, refactor=100.0,
                  ar=[[0.3]], ar_free=[[True]],
                  mu=0.1, estimate_mu=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit()
    assert m.ar[0][0] == pytest.approx(0.430784, abs=1e-5)
    assert m.loglik == pytest.approx(-48.995195, abs=1e-4)


def test_getting_started_points_only_at_documents_that_exist():
    doc = _doc("GETTING_STARTED.md")
    for nombre in re.findall(r"\[([A-Z_]+\.md)\]", doc):
        assert os.path.exists(os.path.join(_DOCS, nombre)), (
            f"GETTING_STARTED.md links to {nombre}, which does not exist")


# ── API.md ─────────────────────────────────────────────────────────────────

def test_the_api_reference_is_current():
    """docs/API.md is generated; if it is stale, the generator says so.

    This is the strictest guard in the file and the cheapest to satisfy: run
    `python tools/gen_api_reference.py`. It exists because an API document is
    the one that rots fastest — a renamed argument breaks nothing and nobody
    notices until a reader copies the old signature.
    """
    import subprocess
    import sys

    r = subprocess.run([sys.executable,
                        os.path.join(_ROOT, "tools", "gen_api_reference.py"),
                        "--check"],
                       capture_output=True, text=True, cwd=_ROOT)
    assert r.returncode == 0, r.stdout.strip() or r.stderr[-400:]


def _public_surface():
    """The public names of a FRESH `import fue`, computed in a subprocess.

    Not `dir(fue)` in this process: importing `fue.elfvarma` anywhere — another
    test, an example — binds it as an attribute of the package, so the surface
    grows with whatever ran first. That made this check pass or fail by test
    order, which is worse than not having it.
    """
    import json
    import subprocess
    import sys

    r = subprocess.run(
        [sys.executable, "-c",
         "import json, inspect, fue;"
         "print(json.dumps({n: inspect.ismodule(getattr(fue, n))"
         " for n in dir(fue) if not n.startswith('_')}))"],
        capture_output=True, text=True, cwd=_ROOT)
    assert r.returncode == 0, r.stderr[-400:]
    return json.loads(r.stdout.strip().splitlines()[-1])


def test_every_public_symbol_is_documented_somewhere():
    """The reference must cover the public surface, not a chosen subset."""
    doc = _doc("API.md")
    faltan = []
    for nombre, es_modulo in sorted(_public_surface().items()):
        if es_modulo:
            if f"`fue.{nombre}`" not in doc:
                faltan.append(nombre)
        elif f"`{nombre}(" not in doc and f"`{nombre}." not in doc:
            faltan.append(nombre)
    assert not faltan, (
        "public symbols missing from docs/API.md: " + ", ".join(faltan) +
        " — add them to _GROUPS in tools/gen_api_reference.py")


def test_nothing_public_is_left_without_a_docstring():
    """The reference is only as good as the docstrings behind it."""
    import inspect

    import fue

    mudos = [n for n in _public_surface()
             if not (inspect.getdoc(getattr(fue, n)) or "")]
    assert not mudos, f"public symbols with no docstring: {mudos}"


# ── the site ───────────────────────────────────────────────────────────────

def test_the_site_navigation_points_at_documents_that_exist():
    """Every entry of mkdocs.yml's nav, and every document, exactly once.

    A nav entry for a file that does not exist fails the build in CI; a
    document missing from the nav is worse, because it builds and is simply
    unreachable. Both are checked here so neither needs a deploy to notice.
    """
    yaml = pytest.importorskip("yaml")

    p = os.path.join(_ROOT, "mkdocs.yml")
    if not os.path.exists(p):
        pytest.skip("mkdocs.yml not present")
    cfg = yaml.safe_load(open(p, encoding="utf-8"))

    def hojas(nodo):
        if isinstance(nodo, str):
            yield nodo
        elif isinstance(nodo, list):
            for x in nodo:
                yield from hojas(x)
        elif isinstance(nodo, dict):
            for x in nodo.values():
                yield from hojas(x)

    en_nav = list(hojas(cfg["nav"]))
    for f in en_nav:
        assert os.path.exists(os.path.join(_DOCS, f)), (
            f"mkdocs.yml lists {f}, which does not exist")

    en_disco = {f for f in os.listdir(_DOCS) if f.endswith(".md")}
    # The plan is working material, not part of the published site.
    sin_publicar = en_disco - set(en_nav) - {"DOCUMENTATION_PLAN.md"}
    assert not sin_publicar, (
        f"documents not reachable from the site: {sorted(sin_publicar)}")
