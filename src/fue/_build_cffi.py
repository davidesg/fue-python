"""
cffi build script for the FUE estimation engine.

Activated in pyproject.toml via:
    cffi_modules = ["src/fue/_build_cffi.py:ffi"]

This file is run at build time (pip install / python setup.py build_ext).
It compiles csrc/fue_api.c + csrc/internal/*.c against GSL and produces
the _fue_engine extension module imported by _engine.py at runtime.
"""

import os
import sys
import cffi

# Compute paths relative to the project root (CWD during pip/setup.py build).
# os.path.relpath normalises away the ".." components that confuse distutils.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT   = os.path.relpath(os.path.join(_THIS_DIR, "..", "..", "csrc"))
_INTERN = os.path.join(_ROOT, "internal")

# cffi cannot process C preprocessor directives (#define, #ifdef, etc.).
# We provide a macro-expanded cdef string with literal values instead.
# Macro values: FUE_MAX_DETVARS=64, FUE_MAX_FACTORS=32, FUE_MAX_POLYORD=64.
# These literals MUST match the #defines in csrc/fue_api.h exactly, or the
# struct layout the extension compiles will disagree with what cffi marshals.
_CDEF = """
typedef struct {
    int    type;
    int    obs_index;
    double harmonic;

    int    nomega;
    double omega[64];
    int    omega_free[64];

    int    ndelta;
    double delta[64];
    int    delta_free[64];

    double *indicator_data;
} FueIntervention;

typedef struct {
    int    order;
    double coefs[64];
    int    coef_free[64];
} FueFactor;

typedef struct {
    int     nobs;
    double *data;
    int     sper;
    int     numbering;
    int     begyear;
    int     begtime;

    double  boxlam;
    double  refactor;

    int     nrdiff;
    int     nadiff;

    double  mu0;
    int     estimate_mu;

    int            ninterventions;
    FueIntervention interventions[64];

    int       nar1;
    FueFactor ar1[32];

    int       nar2;
    FueFactor ar2[32];

    int       nma1;
    FueFactor ma1[32];

    int       nma2;
    FueFactor ma2[32];

    int    nar1f;
    double ar1f_freq[32];
    double ar1f_coef[32];
    int    ar1f_free[32];

    int    nma1f;
    double ma1f_freq[32];
    double ma1f_coef[32];
    int    ma1f_free[32];

    int    ifadf[8];

    int    maxits;
    double grtol;
    double sptol;
    double xitol;
    int    chkma;
    int    eml;
} FueModelSpec;

typedef struct {
    int     ifault;
    int     npar;
    double *params;
    double *std_errors;
    double *cov_matrix;
    double *residuals;
    int     nresiduals;
    double  sigma2;
    double  loglik;
    double  aic;
    double  bic;
} FueResult;

FueResult *fue_estimate(const FueModelSpec *spec);
void       fue_defaults(FueModelSpec *spec);
void       fue_result_free(FueResult *r);
const char *fue_strerror(int ifault);
"""

_SOURCES = [
    os.path.join(_ROOT, "fue_api.c"),
    os.path.join(_INTERN, "elfvarma.c"),
    os.path.join(_INTERN, "usmelard.c"),
    os.path.join(_INTERN, "qnewtopt.c"),
    os.path.join(_INTERN, "nlatools.c"),
    os.path.join(_INTERN, "drvmlest.c"),
]

if sys.platform == "win32":
    _libraries    = ["gsl", "gslcblas"]   # no -lm on Windows
    _compile_args = ["/O2", "/W2"]
    # conda-build sets LIBRARY_INC / LIBRARY_LIB; prefer those when present.
    _lib_inc = os.environ.get("LIBRARY_INC")
    _lib_lib = os.environ.get("LIBRARY_LIB")
    if _lib_inc and _lib_lib:
        _include_dirs = [_ROOT, _INTERN, _lib_inc]
        _library_dirs = [_lib_lib]
    else:
        # GSL via vcpkg (x64-windows-static-md triplet) — used by cibuildwheel.
        # Set GSL_ROOT to override the default vcpkg install prefix.
        _vcpkg_default = os.path.join(
            os.environ.get("VCPKG_INSTALLATION_ROOT", r"C:\vcpkg"),
            "installed", "x64-windows-static-md",
        )
        _gsl_root     = os.environ.get("GSL_ROOT", _vcpkg_default)
        _include_dirs = [_ROOT, _INTERN, os.path.join(_gsl_root, "include")]
        _library_dirs = [os.path.join(_gsl_root, "lib")]
else:
    _libraries    = ["gsl", "gslcblas", "m"]
    _compile_args = ["-O2", "-std=c99", "-Wall"]
    # conda-build sets PREFIX to the environment where host deps are installed.
    _prefix = os.environ.get("PREFIX")
    if _prefix:
        _include_dirs = [_ROOT, _INTERN, os.path.join(_prefix, "include")]
        _library_dirs = [os.path.join(_prefix, "lib")]
    else:
        _include_dirs = [_ROOT, _INTERN]
        _library_dirs = []

ffi = cffi.FFI()
ffi.cdef(_CDEF)

ffi.set_source(
    "fue._fue_engine",
    r'#include "fue_api.h"',
    sources=_SOURCES,
    include_dirs=_include_dirs,
    library_dirs=_library_dirs,
    libraries=_libraries,
    extra_compile_args=_compile_args,
)

if __name__ == "__main__":
    ffi.compile(verbose=True)
