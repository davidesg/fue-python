"""The rescaling factor: one rule for the whole suite.

fue estimates on ``refactor * BoxCox_lambda(y)``. The factor is a numerical
device -- advice on precision from Treadway's *DRVUS* manual (2001) -- and not
a unit: with logs and 100 everything happens to read in per cent, and that is
all. This module holds what every program of the suite needs to know about
it, so that none of them decides on its own:

* the transformation itself (`boxcox`, `inv_boxcox`): fue's, which is the C's;
* how a number estimated on that scale is PRESENTED (`display_scale`,
  `cross_scale`): in the units of its series, the rule atsw-gui's GUIs follow;
* whether the optimiser can work on a series as it is (`check_scale`) and the
  factor to write for it (`recommended_refactor`);
* a log-likelihood that does not depend on the factor (`loglik_y`).

What the optimiser needs was measured on fue, drtran and drvarma, with the C
engines and with the Python ones (atsw-gui, docs/ESTUDIO-escala.md). Its
typical parameter size is fixed at 1, and only the mean, the deterministic
coefficients and the coefficients between series depend on the scale. Two
conditions:

* the SPREAD of the stationary series: safe from 0.1 to 1000; trouble starts
  a decade beyond on either side;
* with d = D = 0, the LEVEL-TO-SPREAD ratio: marginal at 1e6, failing at 1e8
  whatever the factor. No factor mends it.

With logs and 100 every model of the suite's corpus is inside (0.29 to 39),
which is why 100 stays the convention and the factor is not chosen from the
data unless 100 leaves the series outside the band.
"""

from __future__ import annotations

import math
from collections import namedtuple

import numpy as np

#: The suite's convention.
CONVENTION = 100.0
#: Spread of the stationary series the optimiser works on comfortably.
SPREAD_BAND = (0.1, 1000.0)
#: Level-to-spread ratio (d = D = 0) beyond which the mean cannot be resolved.
RATIO_MAX = 1.0e6
#: Factors a writer may choose when 100 does not do, by groups in order of
#: preference: those every program of the suite holds (atsw-gui's fue GUI
#: keeps the factor in a control from 1 to 1000); the small ones the file
#: format still writes (the C writers print two decimals: fue C BUG-0003, so
#: nothing under 0.01); the large ones.
ADMISSIBLE = ((10.0, 1000.0, 1.0), (0.1, 0.01), (1.0e4, 1.0e5, 1.0e6))

DisplayScale = namedtuple("DisplayScale", "k unit percent")
ScaleCheck = namedtuple(
    "ScaleCheck", "ok spread ratio refactor recommended problems")


# ── the transformation ───────────────────────────────────────────────────────

def boxcox(y, lam, refactor=1.0):
    """``refactor * BoxCox_lambda(y)``: log for 0, y itself for 1 (not y - 1),
    ``(y**lam - 1)/lam`` otherwise. Scalar or array. This is the C's
    transformation; nothing in the suite should write another."""
    lam = float(lam)
    if np.ndim(y) == 0:
        y = float(y)
        if lam == 0.0:
            return refactor * math.log(y)
        if lam != 1.0:
            return refactor * (y ** lam - 1.0) / lam
        return refactor * y
    y = np.asarray(y, dtype=float)
    if lam == 0.0:
        return refactor * np.log(y)
    if lam != 1.0:
        return refactor * (np.power(y, lam) - 1.0) / lam
    return refactor * y


def inv_boxcox(z, lam, refactor=1.0):
    """The inverse of `boxcox`. Scalar or array."""
    lam = float(lam)
    if np.ndim(z) == 0:
        v = float(z) / refactor
        if lam == 0.0:
            return math.exp(v)
        if lam != 1.0:
            return (v * lam + 1.0) ** (1.0 / lam)
        return v
    v = np.asarray(z, dtype=float) / refactor
    if lam == 0.0:
        return np.exp(v)
    if lam != 1.0:
        return np.power(v * lam + 1.0, 1.0 / lam)
    return v


def stationary(y, lam, d=0, D=0, freq=1, refactor=1.0):
    """``nabla^d nabla_s^D [refactor * BoxCox(y)]``: what the optimiser sees,
    before any deterministic term or mean is taken out."""
    w = np.asarray(boxcox(np.asarray(y, dtype=float), lam, refactor), dtype=float)
    for _ in range(int(d or 0)):
        w = np.diff(w)
    s = int(freq or 1)
    for _ in range(int(D or 0)):
        w = w[s:] - w[:-s]
    return w


# ── presentation ─────────────────────────────────────────────────────────────

def display_scale(lam, refactor):
    """How a quantity on the engine's scale is shown: multiply by ``k``.

    For a residual, a standard deviation or a variation of the transformed
    series: ``k = 100/refactor`` with logs, a percentage; ``1/refactor``
    otherwise, in the units of the transformed series -- which with lambda = 1
    are the series' own. ``unit`` is "%", "series" or "transformed".

    A mean or a deterministic coefficient is divided by the factor whatever
    lambda is (`param_scale`): with logs it is then a log, and 100 times it a
    percentage.
    """
    lam = float(lam)
    r = float(refactor) if refactor else 1.0
    if abs(lam) < 1e-8:
        return DisplayScale(100.0 / r, "%", True)
    return DisplayScale(1.0 / r, "series" if abs(lam - 1.0) < 1e-12 else "transformed",
                        False)


def param_scale(refactor):
    """Multiplier that takes a mean or a deterministic coefficient from the
    engine's scale to the transformed series' units: ``1/refactor``."""
    r = float(refactor) if refactor else 1.0
    return 1.0 / r


def cross_scale(refactor_out, refactor_in):
    """Multiplier for a coefficient of series ``in`` in the equation of series
    ``out`` (a transfer weight, a cross AR/MA term): ``refactor_in /
    refactor_out``. It is 1 when the files share a factor."""
    ro = float(refactor_out) if refactor_out else 1.0
    ri = float(refactor_in) if refactor_in else 1.0
    return ri / ro


def loglik_y(loglik, y, lam, refactor, n=None):
    """The log-likelihood in the units of y:

        l_y = l + n ln(refactor) + (lambda - 1) sum(ln y_t)

    over the last ``n`` observations (those the likelihood is of; all by
    default). It does not depend on the factor and is comparable between
    lambdas -- the Jacobian of ``refactor * BoxCox``. For a joint model, add
    the terms of each series.
    """
    y = np.asarray(y, dtype=float)
    n = len(y) if n is None else int(n)
    r = float(refactor) if refactor else 1.0
    jac = n * math.log(r)
    if float(lam) != 1.0:
        jac += (float(lam) - 1.0) * float(np.sum(np.log(y[len(y) - n:])))
    return float(loglik) + jac


# ── the check and the factor to write ────────────────────────────────────────

def _spread_ratio(y, lam, d, D, freq, refactor):
    w = stationary(y, lam, d, D, freq, refactor)
    w = w[np.isfinite(w)]
    if len(w) < 3:
        return float("nan"), float("nan")
    sd = float(np.std(w))
    if not sd > 0.0:
        return 0.0, float("inf")
    ratio = float(abs(np.mean(w)) / sd) if not (int(d or 0) or int(D or 0)) else 0.0
    return sd, ratio


def recommended_refactor(y, lam, d=0, D=0, freq=1, current=None):
    """The factor to write for this series and this transformation.

    ``current`` (the base file's factor) if it leaves the spread inside the
    band; else 100, the convention, if that does; else, from the first group of
    `ADMISSIBLE` that has one, the power of ten that brings the spread nearest 1. None when no admissible factor puts
    it in the band: the series' units have to change.
    """
    lo, hi = SPREAD_BAND
    sd1, _ = _spread_ratio(y, lam, d, D, freq, 1.0)
    if not (sd1 == sd1) or not sd1 > 0.0:
        return float(current) if current else CONVENTION
    for f in ([float(current)] if current else []) + [CONVENTION]:
        if lo <= f * sd1 <= hi:
            return f
    for group in ADMISSIBLE:
        inside = [f for f in group if lo <= f * sd1 <= hi]
        if inside:
            return min(inside, key=lambda f: abs(math.log10(f * sd1)))
    return None


def check_scale(y, lam, d=0, D=0, freq=1, refactor=CONVENTION):
    """Can the optimiser work on this series as its file has it?

    Returns a `ScaleCheck`: ``ok``; the ``spread`` and the level-to-spread
    ``ratio`` measured (the ratio only with d = D = 0); the ``refactor`` it was
    asked about and the ``recommended`` one (None if no factor does); and
    ``problems``, a tuple of ("spread" | "ratio", text).
    """
    lo, hi = SPREAD_BAND
    r = float(refactor) if refactor else 1.0
    sd, ratio = _spread_ratio(y, lam, d, D, freq, r)
    rec = recommended_refactor(y, lam, d, D, freq, current=r)
    problems = []
    if sd == sd and not (lo <= sd <= hi):
        what = "small" if sd < lo else "large"
        if rec is None:
            how = ("no factor the file format can hold brings it into the band: "
                   "express the series in other units")
        else:
            how = f"refactor={rec:g} brings it to {sd / r * rec:.3g}"
        problems.append(("spread",
                         f"the stationary series is too {what} for the optimiser: its "
                         f"standard deviation is {sd:.3g} with refactor={r:g} (it works "
                         f"between {lo:g} and {hi:g}); {how}"))
    if ratio == ratio and ratio > RATIO_MAX:
        problems.append(("ratio",
                         f"the level of the stationary series is {ratio:.2g} times its "
                         f"standard deviation: the mean cannot be resolved, whatever the "
                         f"factor. Difference the series, or express it in other units "
                         f"(with lambda < 0 the transformation sits at -1/lambda plus a "
                         f"tiny variation)"))
    return ScaleCheck(not problems, sd, ratio, r, rec, tuple(problems))


def format_refactor(refactor):
    """The factor as a file carries it: two decimals when that is exact (so
    100 stays "100.00", byte for byte), the shortest exact form otherwise.
    Two decimals alone write 0.001 as 0.00, which every reader takes for 1."""
    r = float(refactor)
    s = f"{r:.2f}"
    return s if float(s) == r else f"{r:.10g}"
