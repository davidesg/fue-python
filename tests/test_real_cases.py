"""
Regression tests against real .inp files from dolarization/Analisis/RATIOS.

Reference values were produced by fue-1.13 (installed system-wide).
Tests are skipped if the .inp file is not present.

FUG-format files (Idem/*.inp, RIPC.inp) are explicitly skipped: they use
a different program format that fue-1.13 itself cannot execute (segfault).
"""

import os
import re
import pytest

import fue

try:
    from fue._fue_engine import ffi as _ffi  # noqa
    _C_AVAILABLE = True
except ImportError:
    _C_AVAILABLE = False

requires_c = pytest.mark.skipif(not _C_AVAILABLE, reason="C extension not compiled")

_CASES_DIR = os.path.join(os.path.dirname(__file__), "real_cases")

# FUG-format files: different program, not supported by fue or this parser.
_FUG_FILES = {
    "PRICES/GDP/Sample_1.2003_4.2019/Idem/R.1.inp",
    "PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.inp",
    "PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Idem/R.1.inp",
    "PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Idem/R1.inp",
    "PRICES/PCE/Sample_1.2003_4.2019/Idem/R.1.inp",
}

# Reference values from fue-1.13 output files.
# (rel_path, ref_sigma2, ref_loglik, sigma2_tol, loglik_tol)
_CASES = [
    ("PRICES/GDP/Sample_1.2003_4.2019/Mod/PE.1.inp",                        0.4595942631,  -69.0256059911, 1e-4, 1e-3),
    ("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.1.inp",                         7.0750884543, -160.6143098964, 1e-3, 1e-2),
    ("PRICES/GDP/Sample_1.2003_4.2019/Mod/R.2.inp",                         4.0452844123, -141.8868690241, 1e-3, 1e-2),
    ("PRICES/GDP/Sample_1.2003_4.2019/Mod/SF/R.2.inp",                      0.6231368945,  -80.3087042225, 1e-4, 1e-3),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.0.inp",                 0.2486850071,  -58.6021700483, 1e-4, 1e-3),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.1.inp",                 0.9662469111, -100.9274828448, 1e-4, 1e-3),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.2.inp",                 0.1566201392,  -39.5332743432, 1e-4, 1e-3),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.3.1.inp",               0.1563311728,  -39.8199920374, 1e-4, 1e-3),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.3.inp",                 0.1551235917,  -39.1764600458, 1e-4, 1e-3),
    ("PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.4.inp",                 0.1587564448,  -39.4822974885, 1e-4, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/PE.1.inp",             0.4595942631,  -69.0256059911, 1e-4, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/R.1.inp",              0.6553326752,  -80.9113706253, 1e-4, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/R.2.inp",              0.5046684635,  -72.2974038918, 1e-4, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.1.inp",        0.0008762276,  142.8682603432, 1e-6, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.2.inp",        0.0002154545,  189.8809776965, 1e-6, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.3.inp",        6.58767e-05,   228.2017654123, 1e-7, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.4.inp",        9.45335e-05,   211.2147955479, 1e-7, 1e-2),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/SF/R.2.inp",           0.6231368945,  -80.3087042225, 1e-4, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/old/R.2.inp",          0.011381552,    54.8691279423, 1e-5, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/old/R.3.inp",          0.0010939222,  135.3236947743, 1e-6, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/old/R.4.inp",          0.0003433078,  174.1661689579, 1e-6, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/old/R.5.inp",          0.0002658274,  182.8429916793, 1e-6, 1e-3),
    ("PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/old/R.6.inp",          0.0003485012,  170.713068238,  1e-6, 1e-3),
    ("PRICES/PCE/Sample_1.2003_4.2019/Mod/PE.1.inp",                        0.4595942631,  -69.0256059911, 1e-4, 1e-3),
    ("PRICES/PCE/Sample_1.2003_4.2019/Mod/R.1.inp",                         1.4774986393, -108.2340392827, 1e-3, 1e-2),
    ("PRICES/PCE/Sample_1.2003_4.2019/Mod/R.2.inp",                         0.7296537323,  -84.5787453474, 1e-4, 1e-3),
    ("PRICES/PCE/Sample_1.2003_4.2019/Mod/SF/R.1.inp",                      1.931265089,  -118.202709099,  1e-3, 1e-2),
    ("PRICES/PCE/Sample_1.2003_4.2019/Mod/SF/R.2.inp",                      0.6231368945,  -80.3087042225, 1e-4, 1e-3),
]


def _inp_path(rel):
    return os.path.join(_CASES_DIR, rel)


def _skip_if_missing(rel):
    if not os.path.exists(_inp_path(rel)):
        pytest.skip(f"real-case file not present: {rel}")


@pytest.mark.parametrize("rel,ref_sigma2,ref_loglik,stol,ltol", _CASES,
                         ids=[c[0].split("/")[-1].replace(".inp", "") + "_" + c[0].split("/")[-3]
                              for c in _CASES])
def test_real_case_sigma2(rel, ref_sigma2, ref_loglik, stol, ltol):
    _skip_if_missing(rel)
    _, m = fue.load(_inp_path(rel))
    m.fit()
    r = m._result
    assert abs(r.sigma2 - ref_sigma2) < stol, (
        f"{rel}: sigma2 {r.sigma2:.8f} != ref {ref_sigma2:.8f} (tol {stol})"
    )


@pytest.mark.parametrize("rel,ref_sigma2,ref_loglik,stol,ltol", _CASES,
                         ids=[c[0].split("/")[-1].replace(".inp", "") + "_" + c[0].split("/")[-3]
                              for c in _CASES])
def test_real_case_loglik(rel, ref_sigma2, ref_loglik, stol, ltol):
    _skip_if_missing(rel)
    _, m = fue.load(_inp_path(rel))
    m.fit()
    r = m._result
    assert abs(r.loglik - ref_loglik) < ltol, (
        f"{rel}: loglik {r.loglik:.8f} != ref {ref_loglik:.8f} (tol {ltol})"
    )


# ── write_out / write_pre regression tests (RIPC.1) ──────────────────────────

_RIPC1_DIR = os.path.join(_CASES_DIR,
                          "PRICES/IPC/Mensual/sample_1.2002_12.2007")


@requires_c
def test_write_out_ripc1():
    """write_out() for RIPC.1 must match the reference .out file exactly.

    Requires C because the reference was produced with the C engine; the
    pure-Python raxopt gives parameter values that differ by ~1e-4, which
    propagates through all formatted numbers in the output.
    """
    inp = os.path.join(_RIPC1_DIR, "RIPC.1.inp")
    ref = os.path.join(_RIPC1_DIR, "RIPC.1.out")
    if not os.path.exists(inp) or not os.path.exists(ref):
        pytest.skip("RIPC.1.inp or RIPC.1.out not present")
    _, m = fue.load(inp)
    m.fit()
    assert (m.write_out(inp_name="RIPC.1.inp", out_name="RIPC.1.out").rstrip('\n')
            == open(ref).read().rstrip('\n'))


@requires_c
def test_write_pre_ripc1(tmp_path):
    """write_pre() for RIPC.1 must match the reference .pre file exactly.

    Same caveat as test_write_out_ripc1 regarding the C reference.
    """
    inp = os.path.join(_RIPC1_DIR, "RIPC.1.inp")
    ref = os.path.join(_RIPC1_DIR, "RIPC.1.pre")
    if not os.path.exists(inp) or not os.path.exists(ref):
        pytest.skip("RIPC.1.inp or RIPC.1.pre not present")
    _, m = fue.load(inp)
    m.fit()
    pre_path = tmp_path / "RIPC.1.pre"
    m.write_pre(str(pre_path))
    assert pre_path.read_text() == open(ref).read()
