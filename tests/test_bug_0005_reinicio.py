"""BUG-0005 — Model.fit restarts the optimizer when it stops without zeroing
the gradient.

On IPC-T/Coint/R.4 the BFGS search stopped on the step test with ‖g‖≈1.2e5,
40 loglik units below the optimum, and where exactly it stopped depended on the
last bit (numpy 2 took two more steps). Restarting from where it stopped
reaches the optimum, the same under any numpy. Fits that converge the first
time are not touched.
"""
import os
import warnings

import pytest

import fue

_REAL = os.path.join(os.path.dirname(__file__), "real_cases")
R4 = os.path.join(_REAL, "PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.4.inp")
R1 = os.path.join(_REAL, "PRICES/IPC/Trimestral/Sample_1.2003_4.2019/Mod/Coint/R.1.inp")


def _fit(path):
    _, m = fue.load(path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit()
    return m._result


def test_r4_llega_al_optimo_reiniciando():
    r = _fit(R4)
    assert r.restarts >= 1
    assert r.loglik == pytest.approx(251.682958, abs=1e-4)


def test_un_ajuste_que_converge_no_se_reinicia():
    r = _fit(R1)
    assert r.converged and r.restarts == 0


def test_el_out_dice_los_reinicios(tmp_path):
    _, m = fue.load(R4)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit()
    out = tmp_path / "r4.out"
    fue.write_out(m, str(out))
    assert "OPTIMIZER RESTARTED" in out.read_text()


# ── (B) un alto por el paso con el gradiente anulado ES un máximo ──────────

D1 = "/home/david/Dropbox/SRC/atsw-gui/engines/fue/tests/corpus/D.1.inp"


def test_el_motor_c_registra_el_gradiente_escalado():
    from fue._engine import estimate
    pytest.importorskip("fue._fue_engine")
    _, m = fue.load(R1)
    r = estimate(m)
    assert r["sgrad"] is not None and 0.0 <= r["sgrad"] <= 1.9e-6, \
        "un alto por gradiente tiene el escalado bajo la tolerancia (≈1.82e-6)"


def test_el_puerto_python_registra_el_mismo_gradiente_escalado():
    from fue.cast_us import estimate_py
    _, m = fue.load(R1)
    r = estimate_py(m)
    assert 0.0 <= r["sgrad"] <= 1.9e-6


def test_un_alto_por_el_paso_en_el_optimo_cuenta_como_convergido():
    """R.4 tras su reinicio para por el criterio del paso con el gradiente
    escalado en 4.2e-5: es un máximo y no se avisa."""
    _, m = fue.load(R4)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        m.fit()
    r = m._result
    assert r.termcode == 2 and r.converged
    assert not any("sin anular el gradiente" in str(x.message) for x in w)


@pytest.mark.skipif(not os.path.exists(D1), reason="corpus de atsw-gui")
def test_un_alto_por_el_paso_en_el_optimo_no_se_reinicia():
    """D.1 ya estaba en el óptimo (escalado 3.9e-6): antes se reiniciaba en
    balde y avisaba de que no había convergido."""
    _, m = fue.load(D1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit()
    r = m._result
    assert r.termcode == 2 and r.converged and r.restarts == 0


def test_el_umbral_separa_lo_legitimo_de_lo_atascado():
    from fue.model import FitResult, SGRAD_CONVERGIDO
    base = dict(ifault=0, npar=1, nresiduals=10, sigma2=1.0, loglik=0.0,
                aic=0.0, bic=0.0, params=[], std_errors=[], cov_matrix=[],
                residuals=[], termcode=2)
    assert FitResult({**base, "sgrad": 4.8e-6}).converged
    assert not FitResult({**base, "sgrad": 3.0e4}).converged
    assert not FitResult({**base, "sgrad": None}).converged
    assert 4.2e-5 < SGRAD_CONVERGIDO < 26.9
