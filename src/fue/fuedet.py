"""
fuedet — the deterministic component of a fue model, with fue's rule.

    noise_t = data_t − Σ_j ν_j(B) x_{j,t}

* ν_j(B) = ω_j(B)/δ_j(B), Box-Jenkins signs, by ``calcnu``;
* ν_j is truncated at Nomega_j lags when the variable has no δ, and at
  ``NULAG_DELTA`` = 40 when it has one;
* the variables are summed (``tmp2`` per variable, ``tmp1 += tmp2``) and the
  sum is subtracted ONCE.

This is step [5] of ``cast_us`` in Mauricio's fue.c. The C code lives in
``csrc/fuedet.c``, a byte-identical copy of atsw-gui's ``lib/fuedet/fuedet.c``
(``tests/test_fuedet_vendored.py`` checks its sha256): the same routine serves
fue, fuf, drtran, drvarma and drvec in C, and ``cast_us.cast_us_py`` and
``forecast._build_xi`` here.

Both functions call the compiled routine when the extension has it. The
pure-Python fallback does the same operations in the same order, so on Linux
x86_64 the two give the same bits (``tests/test_fuedet_vendored.py``).

Conventions: ``omegas[j] = [ω0, ω1, …, ωs]``, ``deltas[j] = [δ1, …, δr]``
(0-indexed lists, as ``cast_us`` keeps them), ``inds[j]`` and ``data`` are
1-indexed arrays (index 0 unused) of length ≥ n+1. The result is 1-indexed,
length n+1, with index 0 set to 0.
"""

import numpy as np

#: Lags of ν when the variable has a δ (fue.c [5], FUEDET_NULAG_DELTA).
NULAG_DELTA = 40


def nulag(nomega, ndelta):
    """The last lag of ν: Nomega without δ, 40 with one."""
    return nomega if ndelta == 0 else NULAG_DELTA


def nu_weights(omega, delta, lags):
    """calcnu of fue.c: ν[0..lags] of ω(B)/δ(B), in calcnu's order."""
    s = len(omega) - 1
    r = len(delta)
    om = [float(v) for v in omega]
    de = [float(v) for v in delta]
    nu = [0.0] * (lags + 1)
    nu[0] = om[0]
    for j in range(1, lags + 1):
        sum1 = 0.0
        if r > 0:
            for i in range(1, j + 1):
                if i <= r:
                    sum1 = sum1 + de[i - 1] * nu[j - i]
        sum2 = 0.0
        if s > 0:
            if j <= s:
                sum2 = om[j]
        nu[j] = sum1 - sum2
    return np.array(nu)


def _lib():
    """(ffi, lib) when the compiled extension carries fuedet, else None."""
    from ._engine import _load_c
    eng = _load_c()
    if eng is None:
        return None
    ffi, lib = eng
    if not hasattr(lib, "fuedet_noise"):      # an extension built before fuedet
        return None
    return ffi, lib


def _core_c(eng, omegas, deltas, inds, data, n):
    ffi, lib = eng
    nvar = len(omegas)
    keep = []
    om_p = ffi.new("double *[]", nvar + 1)
    de_p = ffi.new("double *[]", nvar + 1)
    x_p = ffi.new("double *[]", nvar + 1)
    nom = ffi.new("int[]", nvar + 1)
    nde = ffi.new("int[]", nvar + 1)
    for j in range(nvar):
        om = ffi.new("double[]", [float(v) for v in omegas[j]])
        de = ffi.new("double[]", [0.0] + [float(v) for v in deltas[j]])
        x = np.ascontiguousarray(inds[j][: n + 1], dtype=np.float64)
        keep += [om, de, x]
        om_p[j + 1] = om
        de_p[j + 1] = de
        x_p[j + 1] = ffi.cast("double *", ffi.from_buffer(x))
        nom[j + 1] = len(omegas[j]) - 1
        nde[j + 1] = len(deltas[j])
    out = np.zeros(n + 1)
    out_p = ffi.cast("double *", ffi.from_buffer(out))
    if data is None:
        rc = lib.fuedet_sum(nvar, om_p, nom, de_p, nde, x_p, n, out_p)
    else:
        d = np.ascontiguousarray(data[: n + 1], dtype=np.float64)
        keep.append(d)
        rc = lib.fuedet_noise(nvar, om_p, nom, de_p, nde, x_p,
                              ffi.cast("double *", ffi.from_buffer(d)), n, out_p)
    if rc != 0:
        raise MemoryError("fuedet: allocation failure in the C routine")
    out[0] = 0.0
    return out


def _core_py(omegas, deltas, inds, data, n):
    """fue.c [5] in Python, operation by operation: the fallback."""
    nvar = len(omegas)
    out = np.zeros(n + 1)
    if nvar == 0:
        if data is not None:
            out[1:] = np.asarray(data[1: n + 1], dtype=np.float64)
        return out
    lags = [nulag(len(omegas[j]) - 1, len(deltas[j])) for j in range(nvar)]
    nus = [nu_weights(omegas[j], deltas[j], lags[j]).tolist() for j in range(nvar)]
    xs = [[float(v) for v in inds[j][: n + 1]] for j in range(nvar)]
    for i in range(1, n + 1):
        tmp1 = 0.0
        for j in range(nvar):
            nu, x = nus[j], xs[j]
            tmp2 = 0.0
            for k in range(lags[j] + 1):
                if i - k >= 1:
                    tmp2 += nu[k] * x[i - k]
            tmp1 += tmp2
        out[i] = (float(data[i]) - tmp1) if data is not None else tmp1
    return out


def det_sum(omegas, deltas, inds, n, use_c=True):
    """det[t] = Σ_j ν_j(B) x_j[t], t = 1..n (fuedet_sum)."""
    eng = _lib() if use_c else None
    if eng is not None and len(omegas) > 0:
        return _core_c(eng, omegas, deltas, inds, None, n)
    return _core_py(omegas, deltas, inds, None, n)


def det_noise(omegas, deltas, inds, data, n, use_c=True):
    """noise[t] = data[t] − Σ_j ν_j(B) x_j[t], t = 1..n, the sum subtracted
    once (fuedet_noise; fue.c [5])."""
    eng = _lib() if use_c else None
    if eng is not None and len(omegas) > 0:
        return _core_c(eng, omegas, deltas, inds, data, n)
    return _core_py(omegas, deltas, inds, data, n)
