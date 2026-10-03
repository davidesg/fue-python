"""
Bridge between the Python Model API and the cffi-compiled C extension.

Falls back to the pure-Python estimator (cast_us.estimate_py) when the
C extension (_fue_engine) is not available, and says so (BUG-0024).
"""

import numpy as np

# Transport-struct capacities — must match the #defines in csrc/fue_api.h and
# the cdef literals in _build_cffi.py.  Used only to turn an over-capacity model
# into a clear ValueError instead of a raw cffi IndexError (see BUG-0002).
_MAX_FACTORS = 32   # FUE_MAX_FACTORS: max AR/MA factors per block
_MAX_POLYORD = 64   # FUE_MAX_POLYORD: max polynomial order per factor


#: The Hessian behind the standard errors, as the C engine codes it
#: (est_se_how in drvmlest.c) and as every report of the family writes it.
SE_METHODS = {
    0: "bfgs",
    1: "fdhess",
    2: "bfgs (fdhess: the optimum is on the boundary of the admissible region)",
    3: "bfgs (fdhess: the Hessian is not positive definite)",
    4: ("none (fdhess: the optimum is on the boundary of the admissible region; "
        "the search did not move, so it built no BFGS Hessian)"),
    5: ("none (fdhess: the Hessian is not positive definite; "
        "the search did not move, so it built no BFGS Hessian)"),
}


# ── Which engine estimates (BUG-0024) ────────────────────────────────────────
#
# If the C extension does not load, fue estimates with the Python port. That
# is a homologated engine, but not the same path: another optimizer, another
# speed, and on a flat likelihood possibly another optimum (BUG-0005). It used
# to happen in silence. Now the first fallback in a process warns, with the
# import error, and `engine_backend()` says which engine is in use.

_C_ENGINE = None          # (ffi, lib) once loaded
_C_ERROR = None           # the ImportError text when it does not load
_WARNED = False


def _load_c():
    """(ffi, lib) of the C extension, or None. Tried once per process."""
    global _C_ENGINE, _C_ERROR
    if _C_ENGINE is None and _C_ERROR is None:
        try:
            from fue._fue_engine import ffi, lib
            _C_ENGINE = (ffi, lib)
        except ImportError as e:
            _C_ERROR = f"{type(e).__name__}: {e}"
    return _C_ENGINE


def engine_backend() -> str:
    """"c" when the compiled engine loads, "python" when fue falls back to
    the Python port. art seals it in the guion with the instrument's version."""
    return "c" if _load_c() is not None else "python"


def engine_load_error():
    """Why the C engine did not load (the ImportError), or None."""
    _load_c()
    return _C_ERROR


def _warn_python_fallback():
    global _WARNED
    if _WARNED:
        return
    _WARNED = True
    import warnings
    warnings.warn(
        "fue: the C engine did not load "
        f"({_C_ERROR}); estimating with the Python port. Same model, another "
        "optimizer and speed, and on a flat likelihood possibly another "
        "optimum (fue BUG-0005). If fue was installed without the extension "
        "on purpose (FUE_SKIP_C=1) this is expected; otherwise reinstall the "
        "wheel. fue.engine_backend() tells the engine in use (BUG-0024).",
        RuntimeWarning, stacklevel=3)


def se_method_label(code, npar):
    """The method line for a fit with `npar` free parameters."""
    if npar == 0:
        return "none (no free parameters)"
    return SE_METHODS.get(int(code), f"unknown ({code})")


def estimate(model):
    """
    Estimate model parameters by exact ML.

    Tries the C extension first; falls back to the pure-Python estimator
    (scipy L-BFGS-B + elf_scalar) when the extension is not compiled.

    When all ARMA/intervention parameters are fixed (npar=0), uses
    eval_at_params to evaluate the likelihood at the fixed values without
    optimisation — the C backend crashes in this case.

    When the specification carries NO ARMA factor at all, the pure-Python
    engine estimates it — the C backend segfaults (BUG-0013).
    """
    from .cast_us import _build_initial_x
    if len(_build_initial_x(model)) == 0:
        from .cast_us import eval_at_params
        return eval_at_params(model)

    # BUG-0013: a model with deterministic inputs and no ARMA factor at all
    # used to kill the interpreter in the C engine, and was diverted here to
    # the Python engine. Fixed at its root in 0.1.17: nlatools.c accepts the
    # 0 x 0 matrices elf() works with when max(p,q) = 0 (the fix fue-1.14 and
    # the atsw-gui CLI carry). The C now gives the same fit as the Python engine
    # and as the same model with one AR factor pinned at zero.

    c = _load_c()
    if c is None:
        _warn_python_fallback()
        from .cast_us import estimate_py
        return estimate_py(model)
    ffi, lib = c

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
    # BUG-0015: fdhess at the optimum unless the model asks for the BFGS
    # Hessian of the search (Model(hessian="bfgs")).
    spec.hessian_bfgs = 1 if getattr(model, "hessian", "fd") == "bfgs" else 0

    def _fill_factors(spec_arr, factors, free_lists):
        # Capacity of the cffi transport struct (FueFactor coefs[FUE_MAX_POLYORD]
        # and ar1/ar2/ma1/ma2[FUE_MAX_FACTORS] in csrc/fue_api.h).  The engine
        # itself is dynamic; these only bound the marshalling buffer.  Raise a
        # clear error instead of a raw cffi IndexError when a model exceeds them.
        if len(factors) > _MAX_FACTORS:
            raise ValueError(
                f"{len(factors)} factors exceed the binding capacity of "
                f"{_MAX_FACTORS} factors per block; rebuild the extension with a "
                f"larger FUE_MAX_FACTORS in csrc/fue_api.h and _build_cffi.py."
            )
        for i, factor in enumerate(factors):
            if len(factor) > _MAX_POLYORD:
                raise ValueError(
                    f"factor of order {len(factor)} exceeds the binding capacity "
                    f"of {_MAX_POLYORD}; rebuild the extension with a larger "
                    f"FUE_MAX_POLYORD in csrc/fue_api.h and _build_cffi.py."
                )
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
                'termcode': 0, 'niter': 0, 'gnorm': 0.0,
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
            # BUG-0012: el veredicto del optimizador, que hasta ahora se perdia.
            # 1=gradiente 2=paso 3=sin mejora 4=limite de iteraciones 5=paso maximo
            'termcode':   raw.termcode,
            'niter':      raw.niter,
            'gnorm':      raw.gnorm,
            'sgrad':      raw.sgrad,
            'se_method':  se_method_label(raw.se_method, n),
            'params':     np.array([raw.params[i]     for i in range(n)],    dtype=float),
            'std_errors': np.array([raw.std_errors[i] for i in range(n)],    dtype=float),
            'cov_matrix': np.array([raw.cov_matrix[i] for i in range(n * n)],
                                   dtype=float).reshape(n, n) if n > 0 else np.zeros((0, 0)),
            'residuals':  np.array([raw.residuals[i]  for i in range(nr)],   dtype=float),
        }
    finally:
        lib.fue_result_free(raw)

    return data
