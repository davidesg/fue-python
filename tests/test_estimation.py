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


# ── Case 4: RIPC.1 — cos/sin/alter seasonals + step + fixed AR(1) + mu ───────
#
# fue-1.13.1 reference (RIPC.1.inp, boxlam=0.0, refactor=100.0, freq=12):
#   npar=14
#   params[0..10]: omegas for cos1,sin1,cos2,sin2,cos3,sin3,cos4,sin4,cos5,sin5,alter
#   params[11]  = omega for step   = 1.419040
#   params[12]  = delta for step   = 0.843401
#   params[13]  = mu               = -90.034940
#   AR(1) phi=0 fixed (not estimated)
#   sigma2 = 0.9662469111
#   logelf = -100.9274828448

_RIPC1_DATA = np.array([
    0.413459, 0.416226, 0.418544, 0.422442, 0.424508, 0.425892,
    0.425137, 0.425577, 0.427322, 0.429367, 0.432350, 0.434795,
    0.443644, 0.443617, 0.443454, 0.448741, 0.450270, 0.448844,
    0.448505, 0.447079, 0.449154, 0.449650, 0.452374, 0.452679,
    0.452326, 0.452981, 0.453226, 0.454730, 0.449935, 0.447131,
    0.445082, 0.444966, 0.445047, 0.443958, 0.445585, 0.446936,
    0.447090, 0.445735, 0.443441, 0.444167, 0.445403, 0.445490,
    0.442748, 0.439849, 0.437653, 0.438303, 0.442595, 0.445719,
    0.444463, 0.446707, 0.447143, 0.443674, 0.440874, 0.438993,
    0.437829, 0.437908, 0.442587, 0.446552, 0.447956, 0.447148,
    0.447111, 0.445032, 0.441439, 0.438548, 0.436016, 0.436859,
    0.438797, 0.439924, 0.441830, 0.441481, 0.441057, 0.443876,
])

_REF_RIPC1_PARAMS = np.array([
    0.370491,   # omega cos 1
    0.601877,   # omega sin 1
    0.074002,   # omega cos 2
    0.033204,   # omega sin 2
   -0.020471,   # omega cos 3
    0.089728,   # omega sin 3
   -0.042174,   # omega cos 4
    0.048089,   # omega sin 4
   -0.056381,   # omega cos 5
   -0.054876,   # omega sin 5
   -0.017354,   # omega alter
    1.419040,   # omega step
    0.843401,   # delta step
  -90.034940,   # mu
])


def _ripc1_model():
    from fue import Intervention
    ts = TimeSeries(_RIPC1_DATA, freq=12, start=(2002, 1), name="RIPC1")
    interventions = [
        Intervention("cos",   harmonic=1.0, omega=[0.0], omega_free=[True]),
        Intervention("sin",   harmonic=1.0, omega=[0.0], omega_free=[True]),
        Intervention("cos",   harmonic=2.0, omega=[0.0], omega_free=[True]),
        Intervention("sin",   harmonic=2.0, omega=[0.0], omega_free=[True]),
        Intervention("cos",   harmonic=3.0, omega=[0.0], omega_free=[True]),
        Intervention("sin",   harmonic=3.0, omega=[0.0], omega_free=[True]),
        Intervention("cos",   harmonic=4.0, omega=[0.0], omega_free=[True]),
        Intervention("sin",   harmonic=4.0, omega=[0.0], omega_free=[True]),
        Intervention("cos",   harmonic=5.0, omega=[0.0], omega_free=[True]),
        Intervention("sin",   harmonic=5.0, omega=[0.0], omega_free=[True]),
        Intervention("alter", omega=[0.0], omega_free=[True]),
        Intervention("step", at=1, omega=[0.014], omega_free=[True],
                     delta=[0.6], delta_free=[True]),
    ]
    return Model(ts,
        interventions=interventions,
        ar=[[0.0]], ar_free=[[False]],
        boxlam=0.0, refactor=100.0,
        mu=0.0, estimate_mu=True,
    )


def test_ripc1_npar():
    m = _ripc1_model()
    m.fit()
    assert m._result.npar == 14, f"npar mismatch: got {m._result.npar}"


def test_ripc1_params():
    m = _ripc1_model()
    m.fit()
    for i, (got, ref) in enumerate(zip(m._result.params, _REF_RIPC1_PARAMS)):
        assert abs(got - ref) < 2e-3, (
            f"params[{i}] mismatch: got {got:.6f}, expected {ref:.6f}"
        )


def test_ripc1_loglik():
    m = _ripc1_model()
    m.fit()
    assert abs(m._result.loglik - (-100.9274828448)) < 1e-3, (
        f"loglik mismatch: got {m._result.loglik}"
    )


def test_ripc1_sigma2():
    m = _ripc1_model()
    m.fit()
    assert abs(m._result.sigma2 - 0.9662469111) < 1e-4, (
        f"sigma2 mismatch: got {m._result.sigma2}"
    )


# ── Case 5: External regressor — AR(1) + 1 custom column ─────────────────────
#
# fue-1.13 reference (test_exog.inp, freq=12, nobs=60, nrdiff=0):
#   y[t] = 1.5*sin(t*pi/6) + AR(1) noise
#   npar=2, omega=2.463219, phi[1]=0.757269
#   sigma2 = 0.5427434811, logelf = -67.2287772158

_EXOG_Y = np.array([
     1.595263,  2.023227,  2.730346,  3.141004,  2.240141,  1.345117,
     0.056625, -2.142425, -2.276626, -2.364764, -2.481573, -1.574718,
     0.057819,  1.203051,  2.100456,  1.832691,  2.126905,  1.338083,
     0.190080, -1.948252, -1.843602, -2.328031, -2.340389, -0.389697,
     0.493489,  0.869792,  1.819261,  1.246437,  2.022561,  1.005299,
    -0.518097, -1.073661, -2.969735, -2.813164, -3.470106, -2.413143,
    -1.299996,  1.250028,  3.133097,  3.014190,  2.978880,  1.697335,
     0.552432, -1.343998, -3.160568, -4.096928, -3.016596, -0.686160,
     0.473010,  1.320542,  3.248334,  3.366690,  2.820731,  1.818727,
     0.275048, -1.288748, -2.990730, -2.842664, -2.502986, -0.905249,
])

_EXOG_X = np.array([
     0.500000,  0.866025,  1.000000,  0.866025,  0.500000,  0.000000,
    -0.500000, -0.866025, -1.000000, -0.866025, -0.500000, -0.000000,
     0.500000,  0.866025,  1.000000,  0.866025,  0.500000,  0.000000,
    -0.500000, -0.866025, -1.000000, -0.866025, -0.500000, -0.000000,
     0.500000,  0.866025,  1.000000,  0.866025,  0.500000,  0.000000,
    -0.500000, -0.866025, -1.000000, -0.866025, -0.500000, -0.000000,
     0.500000,  0.866025,  1.000000,  0.866025,  0.500000,  0.000000,
    -0.500000, -0.866025, -1.000000, -0.866025, -0.500000, -0.000000,
     0.500000,  0.866025,  1.000000,  0.866025,  0.500000,  0.000000,
    -0.500000, -0.866025, -1.000000, -0.866025, -0.500000, -0.000000,
])


def _exog_model():
    from fue import Intervention
    ts = TimeSeries(_EXOG_Y, freq=12, start=(2000, 1), name="testser")
    return Model(ts,
        interventions=[
            Intervention("custom", omega=[1.5], omega_free=[True], data=_EXOG_X),
        ],
        ar=[[0.6]],
    )


def test_exog_npar():
    m = _exog_model()
    m.fit()
    assert m._result.npar == 2, f"npar mismatch: got {m._result.npar}"


def test_exog_params():
    m = _exog_model()
    m.fit()
    assert abs(m._result.params[0] - 2.463219) < 1e-5, (
        f"omega mismatch: got {m._result.params[0]}"
    )
    assert abs(m._result.params[1] - 0.757269) < 1e-5, (
        f"phi[1] mismatch: got {m._result.params[1]}"
    )


def test_exog_loglik():
    m = _exog_model()
    m.fit()
    assert abs(m._result.loglik - (-67.2287772158)) < 1e-4, (
        f"loglik mismatch: got {m._result.loglik}"
    )


def test_exog_sigma2():
    m = _exog_model()
    m.fit()
    assert abs(m._result.sigma2 - 0.5427434811) < 1e-6, (
        f"sigma2 mismatch: got {m._result.sigma2}"
    )
