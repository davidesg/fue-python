"""BUG-0026 — el Jarque-Bera del `.out` y de `describe()` truncaba n/6.

El puerto copió el C antiguo (`nobs / 6`, los dos int: división entera). Con
n = 215, ⌊215/6⌋ = 35 frente a 35,83: 359,19 en lugar de 367,74 (IPC_US AR(2),
art/BUG-0214). El C de atsw-gui ya usaba `nobs / 6.0`, y scipy, pyfug y art
también: la fórmula es n/6 · (S² + K²/4).
"""
import numpy as np
from scipy import stats

from fue import TimeSeries


def _jb(x):
    x = np.asarray(x, dtype=float)
    n = len(x)
    z = (x - x.mean()) / x.std(ddof=0)
    return n / 6.0 * ((z ** 3).mean() ** 2 + ((z ** 4).mean() - 3.0) ** 2 / 4.0)


def test_describe_usa_n_entre_6_exacto():
    rng = np.random.default_rng(26)
    x = rng.standard_t(4, size=215)              # 215 no es múltiplo de 6
    out = TimeSeries.from_array(x.tolist(), freq=12, start=(2002, 2)).describe()
    assert f"{_jb(x):.6f}" in out
    assert f"{(215 // 6) / (215 / 6) * _jb(x):.6f}" not in out


def test_es_el_de_scipy():
    rng = np.random.default_rng(6)
    x = rng.standard_t(5, size=101)
    assert abs(_jb(x) - stats.jarque_bera(x).statistic) < 1e-9
