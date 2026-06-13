"""
Bridge between the Python Model API and the cffi-compiled C extension.

Falls back to the pure-Python estimator (cast_us.estimate_py) when the
C extension (_fue_engine) is not available.
"""

import numpy as np


def estimate(model):
    """
    Estimate model parameters by exact ML.

    Tries the C extension first; falls back to the pure-Python estimator
    (scipy L-BFGS-B + elf_scalar) when the extension is not compiled.

    When all ARMA/intervention parameters are fixed (npar=0), uses
    eval_at_params to evaluate the likelihood at the fixed values without
    optimisation — the C backend crashes in this case.
    """
    from .cast_us import _build_initial_x
    if len(_build_initial_x(model)) == 0:
        from .cast_us import eval_at_params
        return eval_at_params(model)

    try:
        from fue._fue_engine import ffi, lib
    except ImportError:
        from .cast_us import estimate_py
        return estimate_py(model)

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
    spec.refactor    = model.refactor
    spec.nrdiff      = model.d
    spec.nadiff      = model.D
    spec.mu0         = model.mu0
    spec.estimate_mu = 1 if model.estimate_mu else 0
    spec.chkma       = 1 if model.chkma else 0
    spec.eml         = 1 if model.eml else 0

    def _fill_factors(spec_arr, factors, free_lists):
        for i, factor in enumerate(factors):
            spec_arr[i].order = len(factor)
            free = free_lists[i] if free_lists is not None else None
            for j, v in enumerate(factor):
                spec_arr[i].coefs[j]     = v
                spec_arr[i].coef_free[j] = (
                    0 if (free is not None and not free[j]) else 1
                )

    spec.nar1 = len(model.ar);   _fill_factors(spec.ar1, model.ar,   model.ar_free)
    spec.nma1 = len(model.ma);   _fill_factors(spec.ma1, model.ma,   model.ma_free)
    spec.nar2 = len(model.ar_s); _fill_factors(spec.ar2, model.ar_s, model.ar_s_free)
    spec.nma2 = len(model.ma_s); _fill_factors(spec.ma2, model.ma_s, model.ma_s_free)

    for i, v in enumerate(model.ifadf[:8]):
        spec.ifadf[i] = 1 if v else 0

    spec.nar1f = len(model.ar_f)
    for i, ff in enumerate(model.ar_f):
        spec.ar1f_freq[i] = ff.freq
        spec.ar1f_coef[i] = ff.coef
        spec.ar1f_free[i] = 1 if ff.free else 0

    spec.nma1f = len(model.ma_f)
    for i, ff in enumerate(model.ma_f):
        spec.ma1f_freq[i] = ff.freq
        spec.ma1f_coef[i] = ff.coef
        spec.ma1f_free[i] = 1 if ff.free else 0

    spec.ninterventions = len(model.interventions)
    # Custom indicator buffers must stay alive until fue_estimate returns.
    _custom_bufs = []
    for i, itv in enumerate(model.interventions):
        spec.interventions[i].type      = itv.type_code
        spec.interventions[i].obs_index = itv.at
        spec.interventions[i].harmonic  = itv.harmonic
        spec.interventions[i].nomega    = len(itv.omega)
        for j, v in enumerate(itv.omega):
            spec.interventions[i].omega[j]       = v
            spec.interventions[i].omega_free[j]  = 1 if itv.omega_free[j] else 0
        spec.interventions[i].ndelta = len(itv.delta)
        for j, v in enumerate(itv.delta):
            spec.interventions[i].delta[j]       = v
            spec.interventions[i].delta_free[j]  = 1 if itv.delta_free[j] else 0
        if itv.type == "custom" and itv.data is not None:
            _arr = np.ascontiguousarray(itv.data, dtype=np.float64)
            _buf = ffi.from_buffer("double[]", _arr)
            spec.interventions[i].indicator_data = _buf
            _custom_bufs.append((_arr, _buf))
        else:
            spec.interventions[i].indicator_data = ffi.NULL

    # _data, _data_buf, and _custom_bufs remain alive until fue_estimate returns.
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
