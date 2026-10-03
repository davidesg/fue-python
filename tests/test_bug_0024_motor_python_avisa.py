"""BUG-0024 — when the C engine does not load, fue says so.

It used to fall back to the Python port in silence: same model, another
optimizer, another speed, and on a flat likelihood possibly another optimum
(BUG-0005). Now the first fallback in a process warns with the import error,
and `fue.engine_backend()` tells the engine in use.
"""
import sys
import warnings

import numpy as np
import pytest

import fue
from fue import _engine


def _modelo():
    rng = np.random.default_rng(2)
    y = np.exp(4 + np.cumsum(rng.normal(0.0, 0.01, 120)))
    ts = fue.TimeSeries(data=y.tolist(), freq=12, start=[2000, 1], name="S")
    return fue.Model(ts, d=1, D=0, boxlam=0.0, ar=[[0.2]], ar_free=[[True]],
                     ma=[], ma_free=[], ar_s=[], ma_s=[], interventions=[],
                     ifadf=[0] * 7, mu=0.0, estimate_mu=False)


@pytest.fixture
def sin_motor_c(monkeypatch):
    """The extension's import fails, as with a broken wheel or no GSL."""
    monkeypatch.setitem(sys.modules, "fue._fue_engine", None)
    monkeypatch.setattr(_engine, "_C_ENGINE", None)
    monkeypatch.setattr(_engine, "_C_ERROR", None)
    monkeypatch.setattr(_engine, "_WARNED", False)


def test_avisa_con_la_causa_y_estima_igual(sin_motor_c):
    with pytest.warns(RuntimeWarning, match=r"C engine did not load \((Import|ModuleNotFound)Error"):
        m = _modelo().fit()
    assert m._result is not None
    assert fue.engine_backend() == "python"
    assert "Error" in fue.engine_load_error()


def test_avisa_una_sola_vez_por_proceso(sin_motor_c):
    with pytest.warns(RuntimeWarning):
        _modelo().fit()
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=r".*C engine did not load.*")
        _modelo().fit()


def test_con_el_motor_c_no_hay_aviso():
    pytest.importorskip("fue._fue_engine")
    assert fue.engine_backend() == "c"
    assert fue.engine_load_error() is None
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=r".*C engine did not load.*")
        _modelo().fit()
