"""Intervention: linear transfer function ω(B)/δ(B) applied to an indicator."""


class Intervention:
    """
    Deterministic component with linear transfer function.

    The effect on the series is  ω(B)/δ(B) · x_t  where x_t is a binary
    indicator determined by *type* and *at*.

    Parameters
    ----------
    type : str
        ``'pulse'``    — isolated impulse at *at*
        ``'step'``     — permanent level shift starting at *at*
        ``'ramp'``     — linear ramp starting at *at*
        ``'seasonal'`` — periodic seasonal dummy (*at* = 0-based period within year)
        ``'cos'``      — cosine component cos(2π·harmonic/freq·j); *at* unused
        ``'sin'``      — sine component   sin(2π·harmonic/freq·j); *at* unused
        ``'alter'``    — alternating sign (-1)^j; *at* unused
    at : int
        0-based observation index for pulse/step/ramp (0 = first observation);
        0-based period within year for seasonal.  Unused for cos/sin/alter.
    omega : list of float
        Numerator polynomial coefficients [ω₀, ω₁, …].  Default ``[1.0]``.
    delta : list of float
        Denominator polynomial coefficients [δ₁, δ₂, …].  Default ``[]``
        (no denominator → pure FIR).
    omega_free : list of bool, optional
        Which omega coefficients to estimate.  Defaults to all True.
    delta_free : list of bool, optional
        Which delta coefficients to estimate.  Defaults to all True.
    """

    TYPES = {"pulse": 0, "step": 1, "ramp": 2, "seasonal": 3,
             "cos": 4, "sin": 5, "alter": 6}

    def __init__(self, type, at=0, omega=None, delta=None,
                 omega_free=None, delta_free=None, harmonic=1.0):
        if type not in self.TYPES:
            raise ValueError(f"type must be one of {list(self.TYPES)}")
        self.type     = type
        self.at       = int(at)
        self.harmonic = float(harmonic)
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
