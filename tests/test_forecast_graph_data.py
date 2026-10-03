"""forecast_graph_data — the data of fuf's forecast graph, as usfo.c prepares them.

The graph itself is pyfug's (`pyfug.plot_forecast`): pyfug is the one graphics
engine. fue only decides which series, which bands and in which units.
"""
import math

import numpy as np
import pytest

import fue
from fue.forecast import forecast_graph_data


def _modelo(boxlam, n=200, seed=3):
    rng = np.random.default_rng(seed)
    y = np.exp(np.cumsum(rng.normal(0.002, 0.01, n)) + 4.6)
    ts = fue.TimeSeries(data=y.tolist(), freq=12, start=[2005, 3], name="SIM")
    m = fue.Model(ts, d=1, D=0, boxlam=boxlam, ar=[[0.3]], ar_free=[[True]],
                  ma=[], ma_free=[], ar_s=[], ma_s=[], interventions=[],
                  ifadf=[0] * 7, mu=0.0, estimate_mu=False)
    m.fit()
    return m, y


def test_en_logaritmos_es_la_tasa_anual_en_por_ciento():
    m, y = _modelo(0.0)
    fr = m.forecast_fuf(12)
    g = forecast_graph_data(m, fr)
    L = 12
    assert g["title"] == "LRC anual (%)" and g["L"] == L
    assert len(g["y"]) == 2 * L and len(g["err"]) == L
    # observed part: 100·(ln y_t − ln y_{t−12})
    t = len(y) - L
    assert g["y"][0] == pytest.approx(100 * (math.log(y[t]) - math.log(y[t - 12])),
                                      rel=1e-9)
    # forecast part and one-sd bands come from fr.seasonal_diff
    np.testing.assert_allclose(g["y"][L:], fr.seasonal_diff)
    np.testing.assert_allclose(g["band"][L:] - g["y"][L:], fr.seasonal_diff_std)
    assert g["sigma"] == pytest.approx(100 * math.sqrt(fr.sigma2) / m.refactor)


def test_en_niveles_es_el_cambio_anual_en_sus_unidades():
    """usfo.c: multiplicar por cien un cambio que no es una tasa sacaba el
    gráfico de escala."""
    m, y = _modelo(1.0)
    fr = m.forecast_fuf(12)
    g = forecast_graph_data(m, fr)
    t = len(y) - 12
    assert g["title"] == "Annual change"
    assert g["y"][0] == pytest.approx(y[t] - y[t - 12], rel=1e-9)


def test_la_fecha_del_primer_punto_es_la_de_la_observacion_n_mas_1_menos_L():
    m, y = _modelo(0.0)
    g = forecast_graph_data(m, m.forecast_fuf(12))
    # start 2005/3, n = 200: observation 189 (1-based) is 2020/11
    total = 2005 * 12 + 2 + (len(y) - 12)
    assert (g["first_year"], g["first_season"]) == (total // 12, total % 12 + 1)


def test_pyfug_la_dibuja():
    pytest.importorskip("pyfug")
    from pyfug.graphics import plot_forecast
    m, _ = _modelo(0.0)
    fig = plot_forecast(**forecast_graph_data(m, m.forecast_fuf(12)))
    assert len(fig.axes) == 2
