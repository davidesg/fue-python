"""BUG-0027: the forecast table of the .out says what fuf C says.

The cases and where they come from: tests/data/tablas_de_fuf/README.md.

- On Linux x86_64 with glibc and the C engine, the table is fuf's byte for
  byte, except in the six corpus files whose forecasts explode (a rescaling
  factor of 1 with parameters estimated at 100).
- Everywhere, and for every case, the same rows and fields, the numbers to
  1e-4 relative and a zero's sign ignored.
"""
import math
import pathlib

import pytest

import fue

from tests._plataforma import exact_trajectory

DATOS = pathlib.Path(__file__).parent / "data" / "tablas_de_fuf"
CASOS = sorted(p.stem for p in DATOS.glob("*.tabla"))
EXPLOTAN = {"CPI_USA_model.fuf", "CPI_USA_model_AR.fuf", "IPC_DE_model.fuf",
            "IPC_DE_model_AR.fuf", "IPC_ES_m01_ar1.fuf", "IPC_FR_model.fuf"}


def _tabla(texto):
    lineas = texto.splitlines()
    i = next(k for k, l in enumerate(lineas) if l.startswith("FORECAST REPORT"))
    out, bordes = [], 0
    for l in lineas[i:]:
        out.append(l)
        bordes += l.startswith(" +--")
        if bordes == 3:
            break
    return out


def _del_puerto(caso):
    ts, m = fue.load_fuf(str(DATOS / f"{caso}.inp"))
    return _tabla(m.write_fuf_out(m.forecast_fuf()))


def _de_fuf(caso):
    return (DATOS / f"{caso}.tabla").read_text(encoding="latin-1").splitlines()


def test_hay_casos():
    assert len(CASOS) == 17


@exact_trajectory
@pytest.mark.parametrize("caso", [c for c in CASOS if c not in EXPLOTAN])
def test_byte_a_byte(caso):
    assert _del_puerto(caso) == _de_fuf(caso)


@pytest.mark.parametrize("caso", CASOS)
def test_mismas_filas_y_campos(caso):
    got, ref = _del_puerto(caso), _de_fuf(caso)
    assert len(got) == len(ref)
    for g, r in zip(got, ref):
        if g == r:
            continue
        gc, rc = g.split(), r.split()
        assert len(gc) == len(rc), f"{g!r} vs {r!r}"
        for x, y in zip(gc, rc):
            if x == y:
                continue
            fx, fy = float(x), float(y)          # a word that differs fails here
            assert math.isclose(fx, fy, rel_tol=1e-4, abs_tol=1e-4), f"{g!r} vs {r!r}"


def test_la_escala_en_niveles():
    """lambda = 1: the variation is the change of the level, not 100 times it,
    and the first forecast's std is sigma (fuf BUG-0001)."""
    ts, m = fue.load_fuf(str(DATOS / "CPI_USA_model_AR_refactor100.inp"))
    fr = m.forecast_fuf()
    filas = [l.split() for l in _tabla(m.write_fuf_out(fr)) if l[:9].strip()[:1].isdigit()]
    obs = [f for f in filas if f[-1] != "-"]
    a, b = obs[-2], obs[-1]
    assert float(b[3]) == pytest.approx(float(b[1]) - float(a[1]), abs=2e-4)
    primera = next(f for f in filas if f[-1] == "-" and f[2] != "-")
    assert float(primera[2]) == pytest.approx(math.sqrt(fr.sigma2) / m.refactor, abs=1e-4)


def test_una_serie_anual():
    """One period per year: the date is the year, the annual columns are "-"
    (fuf BUG-0003)."""
    t = _del_puerto("MINK_m01_anual")
    assert "FORECAST ORIGIN : 1911" in t
    filas = [l.split() for l in t if l.strip()[:4].isdigit()]
    assert filas and all(f[0].isdigit() and len(f[0]) == 4 for f in filas)
    assert all(f[5] == "-" and f[6] == "-" for f in filas)
