"""BUG-0015: the standard errors came from the path of the optimizer.

The BFGS matrix raxopt accumulates depends on where the search went: two runs of
one model gave SE(mu) 0.073 and 0.028, and a search that starts at the optimum
never builds it. Since 0.1.17 the covariance is taken from Mauricio's `fdhess`
AT the optimum (the call his drvmlest.c left commented out), with the guards the
rest of the family uses. Pins: the exact GLS of ES_CPI_m10 (drtran battery §1d),
reached from the .pre and from a start 20 % away, by the C engine and by the
pure-Python path.
"""
import os

import numpy as np
import pytest

import fue
from fue.cast_us import _fdhess, estimate_py

PRE = os.path.join(os.path.dirname(__file__), "data", "bug_0015", "ES_CPI_m10.pre")

#: exact GLS standard errors: cos1, sin1 (first harmonic), cos2, alter, phi, mu
GLS = {0: 0.068328, 1: 0.068294, 2: 0.027692, 10: 0.006094, 11: 0.062421, 12: 0.028502}


def _model(hessian="fd", perturb=False):
    m = [x for x in fue.load(PRE) if hasattr(x, "fit")][0]
    m.hessian = hessian
    if perturb:
        m.ar = [[c * 1.2 for c in f] for f in m.ar]
        m.mu0 *= 1.2
    return m


def _check_gls(se):
    for i, ref in GLS.items():
        assert se[i] == pytest.approx(ref, rel=5e-3), (i, se[i], ref)


@pytest.mark.parametrize("perturb", [False, True])
def test_c_engine_gives_the_exact_gls(perturb):
    m = _model(perturb=perturb).fit()
    assert m._result.se_method == "fdhess"
    _check_gls(m._result.std_errors)


@pytest.mark.parametrize("perturb", [False, True])
def test_python_path_gives_the_exact_gls(perturb):
    r = estimate_py(_model(perturb=perturb))
    assert r["se_method"] == "fdhess"
    _check_gls(r["std_errors"])


def test_bfgs_on_request_is_the_old_path_dependent_answer():
    m = _model("bfgs").fit()
    assert m._result.se_method == "bfgs"
    assert m._result.std_errors[12] > 2.0 * GLS[12]       # 0.073, "run A"


def test_hessian_must_be_fd_or_bfgs():
    ts, m = fue.load(PRE)
    with pytest.raises(ValueError):
        fue.Model(ts, hessian="exact")


def test_fdhess_is_the_c_step_and_exact_on_a_quadratic():
    A = np.array([[4.0, 1.0], [1.0, 3.0]])
    f = lambda x: 1.0 + 0.5 * x @ A @ x                    # noqa: E731
    x = np.array([0.3, -2.0])
    np.testing.assert_allclose(_fdhess(f, x, f(x)), A, rtol=1e-5)
