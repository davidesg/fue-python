"""
Bridge between the Python Model API and the cffi-compiled C extension.

Populated in Phase 1 once cast_us() is extracted and fue_api.c is complete.
"""


def estimate(model):
    """
    Convert a Model instance to a FueModelSpec, call fue_estimate(), and
    return the raw FueResult* cdata object.

    Raises ImportError if the C extension has not been compiled yet.
    """
    try:
        from fue._fue_engine import ffi, lib   # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "The FUE C extension (_fue_engine) has not been compiled. "
            "Run: pip install -e .  (requires GSL and a C compiler)"
        ) from exc

    spec = _build_spec(ffi, lib, model)
    result = lib.fue_estimate(spec)
    return result


def _build_spec(ffi, lib, model):
    """Translate a Model object into a FueModelSpec cdata struct."""
    import numpy as np

    spec = ffi.new("FueModelSpec *")
    lib.fue_defaults(spec)

    ts = model.series
    _data = np.ascontiguousarray(ts.data, dtype=np.float64)
    spec.nobs     = ts.nobs
    spec.data     = ffi.cast("double *", _data.ctypes.data)
    spec.sper     = ts.freq if ts.freq > 0 else 1
    spec.numbering = 1 if ts.numbering else 0
    spec.begyear  = ts.start[0]
    spec.begtime  = ts.start[1]

    spec.boxlam      = model.boxlam
    spec.nrdiff      = model.d
    spec.nadiff      = model.D
    spec.mu0         = model.mu0
    spec.estimate_mu = 1 if model.estimate_mu else 0
    spec.chkma       = 1 if model.chkma else 0
    spec.eml         = 1 if model.eml else 0

    # AR factors (regular)
    spec.nar1 = len(model.ar)
    for i, factor in enumerate(model.ar):
        spec.ar1[i].order = len(factor)
        for j, v in enumerate(factor):
            spec.ar1[i].coefs[j]    = v
            spec.ar1[i].coef_free[j] = 1

    # MA factors (regular)
    spec.nma1 = len(model.ma)
    for i, factor in enumerate(model.ma):
        spec.ma1[i].order = len(factor)
        for j, v in enumerate(factor):
            spec.ma1[i].coefs[j]    = v
            spec.ma1[i].coef_free[j] = 1

    # AR seasonal
    spec.nar2 = len(model.ar_s)
    for i, factor in enumerate(model.ar_s):
        spec.ar2[i].order = len(factor)
        for j, v in enumerate(factor):
            spec.ar2[i].coefs[j]    = v
            spec.ar2[i].coef_free[j] = 1

    # MA seasonal
    spec.nma2 = len(model.ma_s)
    for i, factor in enumerate(model.ma_s):
        spec.ma2[i].order = len(factor)
        for j, v in enumerate(factor):
            spec.ma2[i].coefs[j]    = v
            spec.ma2[i].coef_free[j] = 1

    # Interventions
    spec.ninterventions = len(model.interventions)
    for i, itv in enumerate(model.interventions):
        spec.interventions[i].type      = itv.type_code
        spec.interventions[i].obs_index = itv.at
        spec.interventions[i].nomega    = len(itv.omega)
        for j, v in enumerate(itv.omega):
            spec.interventions[i].omega[j]      = v
            spec.interventions[i].omega_free[j] = 1 if itv.omega_free[j] else 0
        spec.interventions[i].ndelta = len(itv.delta)
        for j, v in enumerate(itv.delta):
            spec.interventions[i].delta[j]      = v
            spec.interventions[i].delta_free[j] = 1 if itv.delta_free[j] else 0

    return spec
