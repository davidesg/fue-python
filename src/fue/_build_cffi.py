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
# Macro values: FUE_MAX_DETVARS=64, FUE_MAX_FACTORS=8, FUE_MAX_POLYORD=16.
_CDEF = """
typedef struct {
    int    type;
    int    obs_index;
    double harmonic;

    int    nomega;
    double omega[16];
    int    omega_free[16];

    int    ndelta;
    double delta[16];
    int    delta_free[16];
} FueIntervention;

typedef struct {
    int    order;
    double coefs[16];
    int    coef_free[16];
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
    FueFactor ar1[8];

    int       nar2;
    FueFactor ar2[8];

    int       nma1;
    FueFactor ma1[8];

    int       nma2;
    FueFactor ma2[8];

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
    # GSL via vcpkg (x64-windows-static-md triplet).
    # Set GSL_ROOT to override the default vcpkg install prefix.
    _vcpkg_default = os.path.join(
        os.environ.get("VCPKG_INSTALLATION_ROOT", r"C:\vcpkg"),
        "installed", "x64-windows-static-md",
    )
    _gsl_root    = os.environ.get("GSL_ROOT", _vcpkg_default)
    _include_dirs = [_ROOT, _INTERN, os.path.join(_gsl_root, "include")]
    _library_dirs = [os.path.join(_gsl_root, "lib")]
    _libraries    = ["gsl", "gslcblas"]   # no -lm on Windows
    _compile_args = ["/O2", "/W2"]
else:
    _include_dirs = [_ROOT, _INTERN]
    _library_dirs = []
    _libraries    = ["gsl", "gslcblas", "m"]
    _compile_args = ["-O2", "-std=c99", "-Wall"]

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
