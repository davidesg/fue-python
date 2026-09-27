"""BUG-0025: raxopt's line search never returns when the objective is NaN.

A quadratic whose value is NaN beyond x = 1: the first full step lands there.
Normalised to 1 at x0, as raxopt assumes (fue's objcfunc convention).
"""
import numpy as np
from fue.qnewtopt import raxopt


def f(x):
    x = np.asarray(x, float)
    if x[0] > 1.0:
        return float("nan")
    return float((x[0] - 3.0) ** 2 + (x[1] + 1.0) ** 2) / 10.0     # f(x0) = 1


x, fx, B, termcode, niter, gnorm = raxopt(f, np.array([0.0, 0.0]))
print("returned:", x, fx, "termcode", termcode, "niter", niter)
