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

_ROOT   = os.path.join(os.path.dirname(__file__), "..", "..", "..", "csrc")
_INTERN = os.path.join(_ROOT, "internal")

# Read the public header — this is what cffi parses to build the bindings.
with open(os.path.join(_ROOT, "fue_api.h")) as f:
    _API_HEADER = f.read()

ffi = cffi.FFI()

# Strip preprocessor directives cffi cannot handle (#include, #ifdef, etc.)
import re
_cleaned = re.sub(r'#[^\n]*', '', _API_HEADER)
ffi.cdef(_cleaned)

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
