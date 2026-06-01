"""
Bridge between the Python Model API and the cffi-compiled C extension.
"""

import numpy as np


def estimate(model):
    """
    Convert a Model instance to a FueModelSpec, call fue_estimate(), and
    return the raw FueResult* cdata object.

    Raises ImportError if the C extension has not been compiled yet.
    """
    try:
        from fue._fue_engine import ffi, lib
    except ImportError as exc:
        raise ImportError(
            "The FUE C extension (_fue_engine) has not been compiled. "
            "Run: pip install -e .  (requires GSL and a C compiler)"
        ) from exc

    spec = ffi.new("FueModelSpec *")
    lib.fue_defaults(spec)

    ts = model.series

    # Keep the numpy array alive via ffi.from_buffer so spec.data stays valid.
    _data = np.ascontiguousarray(ts.data, dtype=np.float64)
    _data_buf = ffi.from_buffer("double[]", _data)

    spec.nobs      = ts.nobs
    spec.data      = _data_buf
    spec.sper      = ts.freq if ts.freq > 0 else 1
    spec.numbering = 1 if ts.numbering else 0
    spec.begyear   = ts.start[0]
    spec.begtime   = ts.start[1]

    spec.boxlam      = model.boxlam
    spec.nrdiff      = model.d
    spec.nadiff      = model.D
    spec.mu0         = model.mu0
    spec.estimate_mu = 1 if model.estimate_mu else 0
    spec.chkma       = 1 if model.chkma else 0
    spec.eml         = 1 if model.eml else 0

    spec.nar1 = len(model.ar)
    for i, factor in enumerate(model.ar):
        spec.ar1[i].order = len(factor)
        for j, v in enumerate(factor):
            spec.ar1[i].coefs[j]    = v
            spec.ar1[i].coef_free[j] = 1

    spec.nma1 = len(model.ma)
    for i, factor in enumerate(model.ma):
        spec.ma1[i].order = len(factor)
        for j, v in enumerate(factor):
            spec.ma1[i].coefs[j]    = v
            spec.ma1[i].coef_free[j] = 1

    spec.nar2 = len(model.ar_s)
    for i, factor in enumerate(model.ar_s):
        spec.ar2[i].order = len(factor)
        for j, v in enumerate(factor):
            spec.ar2[i].coefs[j]    = v
            spec.ar2[i].coef_free[j] = 1

    spec.nma2 = len(model.ma_s)
    for i, factor in enumerate(model.ma_s):
        spec.ma2[i].order = len(factor)
        for j, v in enumerate(factor):
            spec.ma2[i].coefs[j]    = v
            spec.ma2[i].coef_free[j] = 1

    spec.ninterventions = len(model.interventions)
    for i, itv in enumerate(model.interventions):
        spec.interventions[i].type      = itv.type_code
        spec.interventions[i].obs_index = itv.at
        spec.interventions[i].nomega    = len(itv.omega)
        for j, v in enumerate(itv.omega):
            spec.interventions[i].omega[j]       = v
            spec.interventions[i].omega_free[j]  = 1 if itv.omega_free[j] else 0
        spec.interventions[i].ndelta = len(itv.delta)
        for j, v in enumerate(itv.delta):
            spec.interventions[i].delta[j]       = v
            spec.interventions[i].delta_free[j]  = 1 if itv.delta_free[j] else 0

    # _data and _data_buf remain alive until fue_estimate returns.
    raw = lib.fue_estimate(spec)
    if raw == ffi.NULL:
        return {'ifault': -1, 'npar': 0, 'nresiduals': 0,
                'sigma2': 0.0, 'loglik': 0.0, 'aic': 0.0, 'bic': 0.0,
                'params': np.array([]), 'std_errors': np.array([]),
                'cov_matrix': np.zeros((0, 0)), 'residuals': np.array([])}

    n  = raw.npar
    nr = raw.nresiduals
    try:
        data = {
            'ifault':     raw.ifault,
            'npar':       n,
            'nresiduals': nr,
            'sigma2':     raw.sigma2,
            'loglik':     raw.loglik,
            'aic':        raw.aic,
            'bic':        raw.bic,
            'params':     np.array([raw.params[i]     for i in range(n)],    dtype=float),
            'std_errors': np.array([raw.std_errors[i] for i in range(n)],    dtype=float),
            'cov_matrix': np.array([raw.cov_matrix[i] for i in range(n * n)],
                                   dtype=float).reshape(n, n) if n > 0 else np.zeros((0, 0)),
            'residuals':  np.array([raw.residuals[i]  for i in range(nr)],   dtype=float),
        }
    finally:
        lib.fue_result_free(raw)

    return data
