import os
import sysconfig

from setuptools import setup

# Set FUE_SKIP_C=1 to install without the C extension (pure-Python estimator
# is used automatically as a fallback).  Useful when GSL is not available.
_cffi_modules = []
if not os.environ.get("FUE_SKIP_C"):
    try:
        import cffi  # noqa: F401
        _cffi_modules = ["src/fue/_build_cffi.py:ffi"]
    except ImportError:
        pass

# abi3 (BUG-0029): cffi's setuptools integration already builds the extension
# for the limited API, so one wheel per platform serves every CPython from 3.10
# on, the ones released after it included. Up to 0.1.22 the wheels were per
# minor version and stopped at cp314: on 3.15 pip took the py3-none-any wheel.
# Not on a free-threaded CPython, which has no limited API. The pure wheel
# (FUE_SKIP_C=1) stays py3-none-any.
_options = {}
if _cffi_modules and not sysconfig.get_config_var("Py_GIL_DISABLED"):
    _options = {"bdist_wheel": {"py_limited_api": "cp310"}}

setup(cffi_modules=_cffi_modules, options=_options)
