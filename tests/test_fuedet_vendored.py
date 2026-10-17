"""
The deterministic filter: one routine, fue.c's step [5].

csrc/fuedet.c and csrc/fuedet.h are VENDORED: byte-identical copies of
atsw-gui's lib/fuedet/fuedet.c and lib/fuedet/fuedet.h, the one
implementation of the deterministic filter that fue, fuf, drtran, drvarma and
drvec use (atsw-gui, drvec BUG-40). The source lives THERE. To update: copy the
two files from atsw-gui's lib/fuedet/ into csrc/, and record their new sha256
below in the same commit.

The rest checks what fuedet.py promises: the compiled routine and the
pure-Python fallback give the same bits, the sum is subtracted once, and ν is
truncated at 40 lags when there is a δ.
"""

import hashlib
import os
import platform

import numpy as np
import pytest

from fue import fuedet

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: sha256 of atsw-gui's lib/fuedet/fuedet.[ch] (atsw-gui branch fue-det-unico,
#: 2026-10-16). The vendored copies must match them byte for byte.
VENDORED = {
    "csrc/fuedet.c": "6f5acf3313538ea4b836b5a74cf1a1c8c4c92ae3193805071757ad9a89904223",
    "csrc/fuedet.h": "880d044ac466115c4073c1c9028989b2c3b5da5a5ce4a80bc8285be8801320a6",
}


@pytest.mark.parametrize("rel", sorted(VENDORED))
def test_vendored_copy_is_atsw_gui_lib_fuedet(rel):
    with open(os.path.join(ROOT, rel), "rb") as f:
        got = hashlib.sha256(f.read()).hexdigest()
    assert got == VENDORED[rel], (
        f"{rel} is not the recorded copy of atsw-gui lib/fuedet: the filter "
        "lives in atsw-gui; copy it from there and record the new sha256")


def test_vendored_copy_matches_a_local_atsw_gui_checkout():
    """When an atsw-gui checkout is at hand (ATSW_GUI, or a sibling of
    atws/), compare with it directly."""
    cands = [os.environ.get("ATSW_GUI", "")]
    cands.append(os.path.join(ROOT, "..", "..", "..", "atsw-gui"))
    src = next((c for c in cands if c and
                os.path.isfile(os.path.join(c, "lib", "fuedet", "fuedet.c"))), None)
    if src is None:
        pytest.skip("no atsw-gui checkout at hand")
    for name in ("fuedet.c", "fuedet.h"):
        with open(os.path.join(src, "lib", "fuedet", name), "rb") as f:
            theirs = f.read()
        with open(os.path.join(ROOT, "csrc", name), "rb") as f:
            ours = f.read()
        if theirs != ours:
            pytest.skip(f"the atsw-gui checkout at {src} carries another "
                        f"{name} (another branch?): the sha256 test decides")


# ── a model with three deterministic variables ───────────────────────────────

N = 120


def _case():
    t = np.arange(N + 1, dtype=float)
    data = 100.0 + 0.3 * t + np.sin(0.7 * t) + 0.01 * ((t * 7919) % 101)
    data[0] = 0.0
    x1 = np.cos(2.0 * np.pi * t / 12.0); x1[0] = 0.0
    x2 = (t >= 2).astype(float); x2[0] = 0.0
    x3 = (t == 30).astype(float)
    omegas = [[0.37], [1.419], [0.61, -0.25]]
    deltas = [[], [0.8434], [0.5, -0.3]]
    return omegas, deltas, [x1, x2, x3], data


def _fue_c_step5(omegas, deltas, inds, data, n):
    """fue.c [5], written out once more, as the reference."""
    out = np.zeros(n + 1)
    lags = [len(o) - 1 if not d else 40 for o, d in zip(omegas, deltas)]
    nus = [fuedet.nu_weights(o, d, L).tolist() for o, d, L in zip(omegas, deltas, lags)]
    for i in range(1, n + 1):
        tmp1 = 0.0
        for nu, x, L in zip(nus, inds, lags):
            tmp2 = 0.0
            for k in range(L + 1):
                if i - k >= 1:
                    tmp2 += nu[k] * float(x[i - k])
            tmp1 += tmp2
        out[i] = float(data[i]) - tmp1
    return out


def test_python_fallback_is_fue_c_step5():
    om, de, x, y = _case()
    got = fuedet.det_noise(om, de, x, y, N, use_c=False)
    assert np.array_equal(got[1:], _fue_c_step5(om, de, x, y, N)[1:])


def test_sum_subtracted_once_is_the_noise_bit_for_bit():
    om, de, x, y = _case()
    for use_c in (False, True):
        det = fuedet.det_sum(om, de, x, N, use_c=use_c)
        noise = fuedet.det_noise(om, de, x, y, N, use_c=use_c)
        assert np.array_equal(y[1:] - det[1:], noise[1:])


@pytest.mark.skipif(fuedet._lib() is None, reason="the C extension has no fuedet")
def test_compiled_routine_and_fallback_agree():
    om, de, x, y = _case()
    c = fuedet.det_noise(om, de, x, y, N, use_c=True)
    p = fuedet.det_noise(om, de, x, y, N, use_c=False)
    if platform.machine() in ("x86_64", "AMD64") and platform.system() == "Linux":
        assert np.array_equal(c, p)           # the same bits
    else:
        assert np.allclose(c, p, rtol=0, atol=1e-12)


def test_nu_is_truncated_at_40_lags_with_a_delta():
    """A step from t = 1 with δ = 0.95 stops growing at lag 40 (BUG-40)."""
    s = np.ones(N + 1)
    for use_c in (False, True):
        det = fuedet.det_sum([[1.0]], [[0.95]], [s], N, use_c=use_c)
        flat = sum(0.95 ** k for k in range(41))
        assert np.allclose(det[41:], flat, rtol=0, atol=1e-12)
        exact = (1 - 0.95 ** N) / 0.05
        assert exact - det[N] > 2.0
    assert fuedet.nulag(2, 0) == 2 and fuedet.nulag(0, 1) == 40


def test_without_deterministic_terms_the_noise_is_the_data():
    _, _, _, y = _case()
    assert np.array_equal(fuedet.det_noise([], [], [], y, N)[1:], y[1:])
    assert not fuedet.det_sum([], [], [], N).any()
