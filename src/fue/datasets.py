"""
Built-in datasets for fue examples and tests.
"""

import numpy as np
from .series import TimeSeries


def sfny() -> TimeSeries:
    """
    SFNY annual precipitation index, 1852–1913 (62 observations).

    The series is a sunspot-New York precipitation proxy used as a
    standard test case for the FUE estimation engine (Mauricio 1995,
    JASA §4 example SFNY.2).

    The recommended model is an ARMAX with a level shift at 1853:

        log(y_t) = ω/(1 − δB) · S_t  +  AR(1) × AR(2)  +  μ  +  ε_t

    where S_t is a step function starting at t=2 (1853).

    Returns
    -------
    TimeSeries
        Annual series (freq=1), start=(1852, 1), name="SFNY".
    """
    data = np.array([
        3.91505848, 2.02125792, 0.81208771, 0.60807414, 1.21576447,
        1.43763055, 1.78032601, 0.82841058, 0.65433228, 0.74324607,
        0.93394905, 0.60094494, 0.80840161, 0.90899270, 0.40822203,
        0.41975993, 0.50368768, 0.57248427, 0.72970370, 0.90175445,
        0.61763439, 0.63607641, 0.67670827, 0.81812744, 0.78095914,
        0.82024104, 0.86103433, 0.84442843, 0.74566075, 0.63347579,
        0.72637557, 0.81351610, 0.79142754, 0.80305873, 0.83867533,
        0.98678814, 0.80485863, 0.81651553, 0.75960093, 0.84070968,
        0.89480882, 0.89407591, 0.84323646, 0.77215182, 0.82509544,
        0.87384443, 0.81360106, 0.78497496, 0.71323360, 0.70688522,
        0.81090348, 0.94831097, 0.72598922, 0.80337325, 0.84011493,
        0.89247202, 0.89328246, 0.90942424, 0.82871189, 0.88647340,
        0.82251497, 0.94737336,
    ])
    return TimeSeries(data, freq=1, start=(1852, 1), name="SFNY")


def ripc() -> TimeSeries:
    """
    RIPC monthly CPI index, January 2002 – December 2007 (72 observations).

    Standard monthly test case for fue, used to verify seasonal ARMAX
    estimation with Fourier harmonics and alternator interventions.

    The series is the log of the Spanish CPI rescaled by 100:
    y_t = 100 · log(CPI_t).

    Returns
    -------
    TimeSeries
        Monthly series (freq=12), start=(2002, 1), name="RIPC".
    """
    data = np.array([
        0.413459, 0.416226, 0.418544, 0.422442, 0.424508, 0.425892,
        0.425137, 0.425577, 0.427322, 0.429367, 0.432350, 0.434795,
        0.443644, 0.443617, 0.443454, 0.448741, 0.450270, 0.448844,
        0.448505, 0.447079, 0.449154, 0.449650, 0.452374, 0.452679,
        0.452326, 0.452981, 0.453226, 0.454730, 0.449935, 0.447131,
        0.445082, 0.444966, 0.445047, 0.443958, 0.445585, 0.446936,
        0.447090, 0.445735, 0.443441, 0.444167, 0.445403, 0.445490,
        0.442748, 0.439849, 0.437653, 0.438303, 0.442595, 0.445719,
        0.444463, 0.446707, 0.447143, 0.443674, 0.440874, 0.438993,
        0.437829, 0.437908, 0.442587, 0.446552, 0.447956, 0.447148,
        0.447111, 0.445032, 0.441439, 0.438548, 0.436016, 0.436859,
        0.438797, 0.439924, 0.441830, 0.441481, 0.441057, 0.443876,
    ])
    return TimeSeries(data, freq=12, start=(2002, 1), name="RIPC")
