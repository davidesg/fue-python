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
