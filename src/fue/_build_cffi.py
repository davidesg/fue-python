"""
cffi build script for the FUE estimation engine.

Activated in pyproject.toml via:
    cffi_modules = ["src/fue/_build_cffi.py:ffi"]

This file is run at build time (pip install / python setup.py build_ext).
It compiles csrc/fue_api.c + csrc/internal/*.c against GSL and produces
the _fue_engine extension module imported by _engine.py at runtime.
"""

import os
import cffi

_ROOT   = os.path.join(os.path.dirname(__file__), "..", "..", "csrc")
_INTERN = os.path.join(_ROOT, "internal")

# cffi cannot process C preprocessor directives (#define, #ifdef, etc.).
# We provide a macro-expanded cdef string with literal values instead.
# Macro values: FUE_MAX_DETVARS=64, FUE_MAX_FACTORS=8, FUE_MAX_POLYORD=16.
_CDEF = """
typedef struct {
    int    type;
    int    obs_index;

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

ffi = cffi.FFI()
ffi.cdef(_CDEF)

ffi.set_source(
    "fue._fue_engine",
    r'#include "fue_api.h"',
    sources=[
        os.path.join(_ROOT, "fue_api.c"),
        os.path.join(_INTERN, "elfvarma.c"),
        os.path.join(_INTERN, "usmelard.c"),
        os.path.join(_INTERN, "qnewtopt.c"),
        os.path.join(_INTERN, "nlatools.c"),
        os.path.join(_INTERN, "drvmlest.c"),
    ],
    include_dirs=[_ROOT, _INTERN],
    libraries=["gsl", "gslcblas", "m"],
    extra_compile_args=["-O2", "-std=c99", "-Wall"],
)

if __name__ == "__main__":
    ffi.compile(verbose=True)
