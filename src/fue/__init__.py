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
from .report import write_out, write_fuf

__version__ = "0.1.0"
__all__ = ["TimeSeries", "Intervention", "Model", "FixedFreqFactor",
           "ForecastResult", "acf", "pacf", "jarque_bera", "ljung_box",
           "load", "load_fuf", "write_out", "write_fuf"]
