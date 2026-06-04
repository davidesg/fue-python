import os
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

setup(cffi_modules=_cffi_modules)
