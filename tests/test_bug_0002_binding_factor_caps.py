"""
Regression test for BUG-0002 — the cffi binding capped AR/MA factors at 8 and
factor order at 16 (fixed cdata arrays FueFactor[8] / coefs[16]), so long-order
and many-factor models crashed with IndexError in fue-Python where fue-C runs.

The engine (Tusmodel) allocates factors dynamically; the caps were only on the
FueModelSpec transport struct.  The fix raises them (FUE_MAX_FACTORS=32,
FUE_MAX_POLYORD=64) and adds a clear ValueError guard past the new cap.

The two real fixtures are the England annual series (n=258), which the C engine
estimates as EN.4 = AR(18) unfactored and as its 9×AR(2) factorisation.  The
reference log-likelihoods below were produced by fue-C 1.13 (/usr/local/bin/fue).

See bugs/BUG-0002-binding-fixed-factor-arrays.md.
"""

import os

import numpy as np
import pytest

from fue import TimeSeries, Model
from fue.inp import load
from fue._engine import _MAX_FACTORS, _MAX_POLYORD

try:
    from fue._fue_engine import ffi  # noqa: F401
    _HAS_C_EXT = True
except ImportError:
    _HAS_C_EXT = False

# The caps and their guard live on the C marshalling path; the pure-Python
# fallback has no such limit, so these tests are meaningful only with the ext.
requires_c_ext = pytest.mark.skipif(
    not _HAS_C_EXT, reason="C extension unavailable; BUG-0002 is a C-binding cap"
)

_HERE = os.path.dirname(__file__)
_EN_AR18 = os.path.join(_HERE, "real_cases", "en4_ar18.inp")
_EN_FAC  = os.path.join(_HERE, "real_cases", "en4_fac9x2.inp")

# fue-C 1.13 references (exact ML), England n=258.
_REF_AR18 = -1590.8128831276   # EN.4, one AR(18) factor
_REF_FAC  = -1597.8503371512   # EN.4 as 9 × AR(2) factors


@requires_c_ext
def test_ar18_unfactored_matches_fue_c():
    # Previously: IndexError: index too large for cdata 'double[16]'.
    ts, m = load(_EN_AR18)
    assert len(m.ar) == 1 and len(m.ar[0]) == 18
    m.fit()
    r = m._result
    assert r.ifault == 0
    assert abs(r.loglik - _REF_AR18) < 1e-6


@requires_c_ext
def test_nine_ar2_factors_matches_fue_c():
    # Previously: IndexError: index too large for cdata 'FueFactor[8]'.
    ts, m = load(_EN_FAC)
    assert len(m.ar) == 9 and all(len(f) == 2 for f in m.ar)
    m.fit()
    r = m._result
    assert r.ifault == 0
    assert abs(r.loglik - _REF_FAC) < 1e-6


@requires_c_ext
def test_order_and_factor_count_at_new_cap_marshal():
    # The maximum order / factor count must still marshal (no IndexError); one
    # past the cap is what the Python guard converts into a ValueError.
    from fue._fue_engine import ffi
    spec = ffi.new("FueModelSpec *")
    spec.ar1[0].coefs[_MAX_POLYORD - 1] = 1.0     # last valid coefficient
    spec.ar1[_MAX_FACTORS - 1].order = 2          # last valid factor slot
    with pytest.raises(IndexError):
        spec.ar1[0].coefs[_MAX_POLYORD] = 1.0     # one past the cap
    with pytest.raises(IndexError):
        spec.ar1[_MAX_FACTORS].order = 2          # one past the cap


@requires_c_ext
def test_guard_raises_valueerror_past_cap():
    # A model past the (new) cap must fail with a clear ValueError naming the
    # limit, not a raw cffi IndexError.
    y = np.cumsum(np.random.default_rng(0).normal(0, 1, 200)) + 100.0
    ts = TimeSeries(y, freq=1, start=(1, 1900))

    too_many = Model(ts, ar=[[0.0, 0.0]] * (_MAX_FACTORS + 1), d=0, boxlam=1.0)
    with pytest.raises(ValueError, match="exceed"):
        too_many.fit()

    too_deep = Model(ts, ar=[[0.0] * (_MAX_POLYORD + 1)], d=0, boxlam=1.0)
    with pytest.raises(ValueError, match="exceeds"):
        too_deep.fit()
