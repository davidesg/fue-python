"""Intervention: linear transfer function ω(B)/δ(B) applied to an indicator."""

import numpy as np


class Intervention:
    """
    Deterministic component with linear transfer function.

    The effect on the series is  ω(B)/δ(B) · x_t  where x_t is a binary
    indicator determined by *type* and *at*.

    Parameters
    ----------
    type : str
        ``'impulse'``  — isolated impulse at *at*.  ``'pulse'`` is accepted as a
                         deprecated alias and normalised to ``'impulse'``: the
                         school's vocabulary — and the one in the `.pre`/`.inp`
                         format and in fue C — is ``impulse``, and having two
                         names for one thing is how a file ends up with a
                         keyword the other interpreter does not know.
        ``'compimp'``  — COMPENSATED impulse: +1 at *at* and **−1 at *at*+1**
                         (``compimp`` in the .pre/.inp format).  A pulse that is
                         undone the next period: the level returns to where it
                         was, so it does not shift the mean of a differenced
                         series the way a plain pulse does.
        ``'step'``     — permanent level shift starting at *at*
        ``'ramp'``     — linear ramp starting at *at*
        ``'seasonal'`` — periodic seasonal dummy (*at* = 0-based period within
                         year).  **Python-only, and outside the methodology:**
                         deterministic seasonality is parameterised with
                         HARMONICS (``cos``/``sin`` plus the Nyquist
                         ``alter``), not with dummies, which is why fue C has no
                         such regressor and the format has no keyword for it.
                         Writing one raises rather than emit a word fue C would
                         silently take for a non-standard variable.
        ``'easter'``   — Easter-holiday variable; **monthly series only**
                         (freq == 12), *at* unused.  See `_build_indicator`.
        ``'trend'``    — deterministic linear trend, 1, 2, …, n; *at* unused
        ``'cos'``      — cosine component cos(2π·harmonic/freq·j); *at* unused
        ``'sin'``      — sine component   sin(2π·harmonic/freq·j); *at* unused
        ``'alter'``    — alternating sign (-1)^j; *at* unused
        ``'custom'``   — external indicator supplied as *data* array
    at : int
        0-based observation index for pulse/compimp/step/ramp (0 = first
        observation); 0-based period within year for seasonal.  Unused for
        easter/trend/cos/sin/alter/custom.
    omega : list of float
        Numerator polynomial coefficients [ω₀, ω₁, …].  Default ``[1.0]``.
    delta : list of float
        Denominator polynomial coefficients [δ₁, δ₂, …].  Default ``[]``
        (no denominator → pure FIR).
    omega_free : list of bool, optional
        Which omega coefficients to estimate.  Defaults to all True.
    delta_free : list of bool, optional
        Which delta coefficients to estimate.  Defaults to all True.
    data : array-like, optional
        Pre-computed indicator values, length nobs.  Required for type='custom'.
    """

    # The codes are the contract with the C engine (`FUE_ITV_*` in
    # `csrc/fue_api.h`): only ever APPENDED, never renumbered, or an extension
    # compiled earlier would build a different regressor than the one asked for.
    TYPES = {"impulse": 0, "step": 1, "ramp": 2, "seasonal": 3,
             "cos": 4, "sin": 5, "alter": 6, "custom": 7,
             "compimp": 8, "easter": 9, "trend": 10}

    #: Deprecated spellings, normalised on construction. `pulse` was this
    #: package's own name for the school's `impulse`; both are accepted, one is
    #: stored, so there is a single vocabulary shared with fue C.
    ALIASES = {"pulse": "impulse"}

    def __init__(self, type, at=0, omega=None, delta=None,
                 omega_free=None, delta_free=None, harmonic=1.0, data=None):
        type = self.ALIASES.get(type, type)
        if type not in self.TYPES:
            raise ValueError(f"type must be one of {list(self.TYPES)}")
        if type == "custom" and data is None:
            raise ValueError("data must be provided for type='custom'")
        self.type     = type
        self.at       = int(at)
        self.harmonic = float(harmonic)
        self.data     = np.asarray(data, dtype=float) if data is not None else None
        self.omega  = list(omega) if omega is not None else [1.0]
        self.delta  = list(delta) if delta is not None else []
        self.omega_free = (list(omega_free) if omega_free is not None
                           else [True] * len(self.omega))
        self.delta_free = (list(delta_free) if delta_free is not None
                           else [True] * len(self.delta))
        if len(self.omega_free) != len(self.omega):
            raise ValueError("omega_free must have the same length as omega")
        if len(self.delta_free) != len(self.delta):
            raise ValueError("delta_free must have the same length as delta")

    @property
    def type_code(self):
        return self.TYPES[self.type]

    def __repr__(self):
        return (f"Intervention(type={self.type!r}, at={self.at}, "
                f"omega={self.omega}, delta={self.delta})")
