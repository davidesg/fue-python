"""
Phase 1.3 — numerical equivalence tests.

Each test runs the Python API on the same data and model spec as fue-1.13.1
and compares estimated parameters, sigma², and log-likelihood against
hard-coded reference values obtained from the reference binary.

Reference values were captured by running:
    fue-1.13.1/bin/fue <case>
and reading the corresponding .out file.

Tolerances:
    params  : abs 1e-6
    loglik  : abs 1e-4
    sigma2  : abs 1e-6
"""

import numpy as np
import pytest

pytest.importorskip("fue._fue_engine",
                    reason="C extension not compiled — skip estimation tests")

from fue import TimeSeries, Model

# First 30 annual observations from the SFNY.2 dataset (1852-1881)
_SFNY30 = np.array([
    3.91505848, 2.02125792, 0.81208771, 0.60807414, 1.21576447,
    1.43763055, 1.78032601, 0.82841058, 0.65433228, 0.74324607,
    0.93394905, 0.60094494, 0.80840161, 0.90899270, 0.40822203,
    0.41975993, 0.50368768, 0.57248427, 0.72970370, 0.90175445,
    0.61763439, 0.63607641, 0.67670827, 0.81812744, 0.78095914,
    0.82024104, 0.86103433, 0.84442843, 0.74566075, 0.63347579,
])

# Full SFNY.2 dataset (1852-1913, 62 annual observations)
_SFNY62 = np.array([
    3.91505848, 2.02125792, 0.81208771, 0.60807414, 1.21576447,
    1.43763055, 1.78032601, 0.82841058, 0.65433228, 0.74324607,
    0.93394905, 0.60094494, 0.80840161, 0.90899270, 0.40822203,
    0.41975993, 0.50368768, 0.57248427, 0.72970370, 0.90175445,
    0.61763439, 0.63607641, 0.67670827, 0.81812744, 0.78095914,
    0.82024104, 0.86103433, 0.84442843, 0.74566075, 0.63347579,
    0.72637557, 0.81351610, 0.79142754, 0.80305873, 0.83867533,
    0.98678814, 0.80485863, 0.81651553, 0.75960093, 0.84070968,
    0.89480882, 0.89407591, 0.84323646, 0.77215182, 0.82509544,
    0.87384443, 0.81360106, 0.78497496, 0.71323360, 0.70688522,
    0.81090348, 0.94831097, 0.72598922, 0.80337325, 0.84011493,
    0.89247202, 0.89328246, 0.90942424, 0.82871189, 0.88647340,
    0.82251497, 0.94737336,
])


# ── Case 1: AR(1), annual, no interventions ───────────────────────────────────
#
# fue-1.13.1 reference (test_ar1.inp, boxlam=1.0, nrdiff=0):
#   npar=1, phi[1]=0.9747519833, sigma2=0.2482607765, logelf=-23.1683049163

def test_ar1_params():
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m  = Model(ts, ar=[[0.5]])
    m.fit()

    assert m._result.npar == 1
    assert abs(m._result.params[0] - 0.9747519833) < 1e-6, (
        f"phi[1] mismatch: got {m._result.params[0]}"
    )


def test_ar1_loglik():
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m  = Model(ts, ar=[[0.5]])
    m.fit()

    assert abs(m._result.loglik - (-23.1683049163)) < 1e-4, (
        f"loglik mismatch: got {m._result.loglik}"
    )


def test_ar1_sigma2():
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m  = Model(ts, ar=[[0.5]])
    m.fit()

    assert abs(m._result.sigma2 - 0.2482607765) < 1e-6, (
        f"sigma2 mismatch: got {m._result.sigma2}"
    )


# ── Case 2: IMA(1,1), annual, no interventions ────────────────────────────────
#
# fue-1.13.1 reference (test_ima11.inp, boxlam=1.0, nrdiff=1):
#   npar=1, theta[1]=-0.4228241648, sigma2=0.2060862798, logelf=-18.3455244692

def test_ima11_params():
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m  = Model(ts, ma=[[-0.5]], d=1)
    m.fit()

    assert m._result.npar == 1
    assert abs(m._result.params[0] - (-0.4228241648)) < 1e-6, (
        f"theta[1] mismatch: got {m._result.params[0]}"
    )


def test_ima11_loglik():
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m  = Model(ts, ma=[[-0.5]], d=1)
    m.fit()

    assert abs(m._result.loglik - (-18.3455244692)) < 1e-4, (
        f"loglik mismatch: got {m._result.loglik}"
    )


def test_ima11_sigma2():
    ts = TimeSeries(_SFNY30, freq=1, start=(1852, 1), name="SFNY30")
    m  = Model(ts, ma=[[-0.5]], d=1)
    m.fit()

    assert abs(m._result.sigma2 - 0.2060862798) < 1e-6, (
        f"sigma2 mismatch: got {m._result.sigma2}"
    )


# ── Case 3: SFNY.2 full model — step intervention + 2 AR factors + mu ────────
#
# fue-1.13.1 reference (SFNY.2.inp, boxlam=0.0, nrdiff=0):
#   npar=6
#   params[0] = omega       = -0.600361  (step omega)
#   params[1] = delta       =  0.587261  (step delta)
#   params[2] = AR1 fac1 φ1 =  0.131073
#   params[3] = AR1 fac2 φ1 =  0.487958
#   params[4] = AR1 fac2 φ2 = -0.258762
#   params[5] = mu          =  1.209415
#   sigma2 = 0.0370593261
#   logelf = 13.9573576937

_REF_SFNY2_PARAMS = np.array([
    -0.600361,  0.587261,
     0.131073,
     0.487958, -0.258762,
     1.209415,
])


def _sfny2_model():
    from fue import Intervention
    ts = TimeSeries(_SFNY62, freq=1, start=(1852, 1), name="SFNY")
    return Model(ts,
        interventions=[
            Intervention("step", at=1,
                         omega=[0.08], omega_free=[True],
                         delta=[0.6],  delta_free=[True])
        ],
        ar=[[0.8], [-0.1, -0.1]],
        boxlam=0.0,
        mu=0.0,
        estimate_mu=True,
    )


def test_sfny2_npar():
    m = _sfny2_model()
    m.fit()
    assert m._result.npar == 6


def test_sfny2_params():
    m = _sfny2_model()
    m.fit()
    for i, (got, ref) in enumerate(zip(m._result.params, _REF_SFNY2_PARAMS)):
        assert abs(got - ref) < 1e-5, (
            f"params[{i}] mismatch: got {got:.6f}, expected {ref:.6f}"
        )


def test_sfny2_loglik():
    m = _sfny2_model()
    m.fit()
    assert abs(m._result.loglik - 13.9573576937) < 1e-4, (
        f"loglik mismatch: got {m._result.loglik}"
    )


def test_sfny2_sigma2():
    m = _sfny2_model()
    m.fit()
    assert abs(m._result.sigma2 - 0.0370593261) < 1e-6, (
        f"sigma2 mismatch: got {m._result.sigma2}"
    )
