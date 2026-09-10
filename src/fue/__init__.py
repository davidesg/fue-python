"""
fue — Python interface to the FUE exact maximum likelihood estimation engine.

Typical usage::

    import fue

    ts = fue.TimeSeries.from_array(data, freq=12, start=(1990, 1))

    m = fue.Model(ts, ar=[[1]], ma=[[1]], d=1, D=1)
    m.fit()
    print(m.summary())

    m.plot_residuals()
    m.residuals.plot_acf()
"""

from .series import TimeSeries
from .intervention import Intervention
from .model import Model, FixedFreqFactor
from .forecast import ForecastResult
from .diagnostics import acf, pacf, jarque_bera, ljung_box
from .inp import load, load_fuf
from .report import write_out, write_fuf, write_fuf_out
from .report_forecast import write_forecast_report
from . import datasets

def _version() -> str:
    """La versión, de UNA sola fuente.

    Estaba escrita a mano aquí y se quedó tres versiones atrás: la 0.1.14
    publicada en PyPI declara `__version__ = "0.1.11"` mientras su metadata dice
    0.1.14. Es la enfermedad de siempre —un dato escrito dos veces y la segunda
    copia rezagada— y ya había mordido una vez: la página de API de la 0.1.10
    decía «fue 0.1.9», y aquello se arregló en el generador (leyendo el
    `pyproject`) sin tocar la raíz.

    El orden importa y es el mismo que usa `art.version_instrumento`:

      1. el `pyproject.toml` del árbol, SI existe — en una instalación editable
         la metadata se escribió al instalar y no se regenera al subir la
         versión, así que sería ella la que mentiría;
      2. la metadata del paquete instalado, que es lo correcto para una rueda;
      3. `"0+desconocida"`, que es honesto: dice que no consta, no un número.

    Nunca levanta: esto corre en cada `import fue`.
    """
    import os
    import re
    raiz = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    try:
        with open(os.path.join(raiz, "pyproject.toml"), encoding="utf-8") as fh:
            m = re.search(r'^version\s*=\s*"([^"]+)"', fh.read(), re.M)
        if m:
            return m.group(1)
    except OSError:
        pass
    try:
        from importlib.metadata import version as _v
        return _v("fue")
    except Exception:
        return "0+desconocida"


__version__ = _version()
__all__ = ["TimeSeries", "Intervention", "Model", "FixedFreqFactor",
           "ForecastResult", "acf", "pacf", "jarque_bera", "ljung_box",
           "load", "load_fuf", "write_out", "write_fuf", "write_fuf_out",
           "write_forecast_report", "datasets"]
