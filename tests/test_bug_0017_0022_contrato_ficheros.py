"""BUG-0017…0022 — el contrato de ficheros: lo que se lee y se escribe sin mentir.

Cada prueba es la reproducción del informe reducida a un fichero sintético, para
no depender del corpus de conformidad (`atws/conformidad/corpus`), que vive
fuera de este repositorio.
"""
import warnings

import numpy as np
import pytest

import fue
from fue.report import write_pre


def _inp(freq_line, obs_line, det_block, arma_ar, mu_line, bc_line, bands,
         rows, ifadf=" 0"):
    return "\n".join([
        "** Frequency of time series: either 1(A), 4(Q) or 12(M):",
        freq_line,
        "** Number of observations and starting date of time series:",
        obs_line,
        "** Number of deterministic variables (including seasonal components):",
        *det_block,
        "**Number and orders of regular AR operators:",
        *arma_ar,
        "** Number and orders of annual AR operators:",
        "0",
        "** Number and orders of regular MA operators:",
        "0",
        "** Number and orders of anual MA operators:",
        "0",
        "** Number and frequencies of regular AR(2) operators with fixed frequency:",
        "0",
        "** Number and frequencies of regular MA(2) operators with fixed frequency:",
        "0",
        "** Mean parameter (mu):",
        mu_line,
        "** Box-Cox lambda, regular differences and complete annual differences:",
        bc_line,
        "** Individual factors of the annual difference (from freq 0.0):",
        ifadf,
        "** ACF/PACF bands (0 Automatic) and reescaling factor:",
        bands,
        "** Time series (stochastic and non-standard deterministic variables):",
        *rows,
    ]) + "\n"


def _serie(n=60, seed=3):
    rng = np.random.default_rng(seed)
    y = np.empty(n)
    y[0] = 0.0
    for t in range(1, n):
        y[t] = 0.6 * y[t - 1] + rng.standard_normal()
    return y + 10.0


# ── BUG-0017 ────────────────────────────────────────────────────────────────

def _con_regresor(tmp_path, columnas):
    y = _serie()
    x = np.linspace(-1.0, 1.0, len(y))
    rows = [f"{a!r} {b!r}" if columnas else f"{a!r}" for a, b in zip(y, x)]
    p = tmp_path / "reg.inp"
    p.write_text(_inp(" 1", f" {len(y)} 1 1950 REG",
                      ["1", "**", "custom", "**", "0", "**", "0.0 1", "**", "0"],
                      ["1 1", "**", "0.5 1"], "0", "1.00 0 0", " 0 1.00", rows))
    return p, x


def test_0017_la_columna_del_regresor_se_lee(tmp_path):
    p, x = _con_regresor(tmp_path, columnas=True)
    _, m = fue.load(str(p))
    assert np.array_equal(m.interventions[0].data, x)


def test_0017_sin_la_columna_el_lector_falla_en_vez_de_inventar_ceros(tmp_path):
    p, _ = _con_regresor(tmp_path, columnas=False)
    with pytest.raises(ValueError, match="BUG-0017"):
        fue.load(str(p))


# ── BUG-0018 ────────────────────────────────────────────────────────────────

def _simple(tmp_path, freq_line=" 1", bands=" 0 1.00", obs=None):
    y = _serie()
    p = tmp_path / "s.inp"
    p.write_text(_inp(freq_line, obs or f" {len(y)} 1 1950 S", [" 0"],
                      ["1 1", "**", "0.5 1"], "0", "1.00 0 0", bands,
                      [f"{v!r}" for v in y]))
    return p


def test_0018_number_y_cbands_se_conservan_y_vuelven(tmp_path):
    p = _simple(tmp_path, freq_line=" number", bands=" 2.50 50.00")
    ts, m = fue.load(str(p))
    assert ts.freq == 1 and ts.numbering
    assert m.cbands == 2.5 and m.refactor == 50.0
    m.fit()
    out = tmp_path / "s.pre"
    write_pre(m, str(out))
    ts2, m2 = fue.load(str(out))
    assert ts2.numbering and m2.cbands == 2.5
    assert " number\n" in out.read_text()


def test_0018_una_serie_anual_fechada_no_se_marca_como_numerada(tmp_path):
    ts, m = fue.load(str(_simple(tmp_path)))
    assert not ts.numbering and m.cbands == 0.0


# ── BUG-0019 ────────────────────────────────────────────────────────────────

def test_0019_phi2_menos_cero_se_acepta_como_en_el_motor_y_avisa():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ff = fue.FixedFreqFactor(freq=2, coef=-0.0, free=True)
    assert ff.coef == 0.0
    assert any("BUG-0019" in str(x.message) for x in w)


def test_0019_phi2_positivo_se_sigue_rechazando():
    with pytest.raises(ValueError):
        fue.FixedFreqFactor(freq=2, coef=1e-9)


# ── BUG-0020 / BUG-0021 ─────────────────────────────────────────────────────

def _con_delta2_y_fijos(tmp_path):
    y = _serie(80)
    y[40:] += 3.0
    p = tmp_path / "d.inp"
    p.write_text(_inp(" 1", f" {len(y)} 1 1950 D",
                      ["1", "**", "step 1990", "**", "0", "**", "0.0 1",
                       "**", "2", "**", "0.3 1", "-0.2 1"],
                      ["1 1", "**", "0.941176"],        # AR FIJO (sin bandera)
                      "-7.123456789 0",                  # μ FIJA no nula
                      "0.333333333333 0 0", " 0 1.00",   # λ = 1/3
                      [f"{v!r}" for v in y]))
    return p


def test_0020_dos_deltas_salen_en_lineas_distintas_y_el_pre_se_relee(tmp_path):
    _, m = fue.load(str(_con_delta2_y_fijos(tmp_path)))
    m.fit()
    out = tmp_path / "d.pre"
    write_pre(m, str(out))
    _, m2 = fue.load(str(out))                     # antes: ValueError '1-0.8631'
    assert len(m2.interventions[0].delta) == 2


def test_0021_lo_fijo_vuelve_identico(tmp_path):
    _, m = fue.load(str(_con_delta2_y_fijos(tmp_path)))
    assert m.ar_free == [[False]]
    m.fit()
    out = tmp_path / "d.pre"
    write_pre(m, str(out))
    _, m2 = fue.load(str(out))
    assert m2.ar == [[0.941176]] and m2.ar_free == [[False]]
    assert m2.boxlam == m.boxlam
    assert m2.mu0 == -7.123456789 and not m2.estimate_mu


def test_0021_los_datos_del_pre_releen_identicos(tmp_path):
    ts, m = fue.load(str(_con_delta2_y_fijos(tmp_path)))
    m.fit()
    out = tmp_path / "d.pre"
    write_pre(m, str(out))
    ts2, _ = fue.load(str(out))
    assert np.array_equal(ts.data, ts2.data)


# ── BUG-0022 ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("obs, year, name", [
    (" 60 1 1766 2020", 1766, "2020"),   # el del informe: nombre numérico
    (" 60 1 1766 EN.1", 1766, "EN.1"),
    (" 60 1768 1768 GE", 1768, "GE"),    # la forma con el año repetido
    (" 60 1770 GE", 1770, "GE"),         # DRVUS con nombre
])
def test_0022_la_cabecera_anual(tmp_path, obs, year, name):
    ts, _ = fue.load(str(_simple(tmp_path, obs=obs)))
    assert ts.start == (year, 1) and ts.name == name


def test_0022_sin_nombre_se_avisa(tmp_path):
    with pytest.warns(RuntimeWarning, match="BUG-0022"):
        ts, _ = fue.load(str(_simple(tmp_path, obs=" 60 1770")))
    assert ts.start == (1770, 1) and ts.name == "series"
