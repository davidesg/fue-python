"""BUG-0004 — forecast_fuf must forecast from the FIT, not from the pre-fit seed
attributes.

Model.fit() does not sync the estimate back into self.ar/ar_s/ma/ma_s/mu0 (they keep
the seed values). forecast_fuf -> eval_at_params built x0 from those attributes, so it
forecast from the seeds. When a seed is far from the fit — e.g. a mu0 seed carried in
a rescaled space, as art._make_model does (×100) — the level forecast runs away. The
fix: eval_at_params uses model._result.params (the fit) when present, so the forecast
is invariant to the (unsynced) attributes.
"""
import numpy as np
import fue


def _driftless_ar1_level(n=220, phi=0.4, seed=5):
    """A ~driftless log-level with AR(1) noise; level ~100."""
    rng = np.random.default_rng(seed)
    a = rng.normal(0.0, 0.01, n)
    w = np.zeros(n)
    for i in range(1, n):
        w[i] = phi * w[i - 1] + a[i]
    return np.exp(np.cumsum(w) + 4.6)


def _fit_model():
    y = _driftless_ar1_level()
    ts = fue.TimeSeries(data=y.tolist(), freq=12, start=[2000, 1], name="SIM")
    m = fue.Model(ts, d=1, D=0, boxlam=0.0,
                  ar=[[0.4]], ar_free=[[True]], ma=[], ma_free=[],
                  ar_s=[], ma_s=[], interventions=[],
                  ifadf=[0] * (12 // 2 + 1), mu=0.0, estimate_mu=True)
    m.fit()
    return m, float(y[-1])


def test_forecast_fuf_is_sane_and_invariant_to_stale_attributes():
    m, last = _fit_model()
    fc_fit = list(m.forecast_fuf(6).level[:6])

    # (1) A fitted model forecasts near the last observation (mean-reverting), not away.
    assert all(abs(v - last) < 0.10 * last for v in fc_fit), \
        f"forecast diverged from the last obs: {fc_fit}"

    # (2) The forecast reads the FIT, not the attributes: corrupt the (unsynced)
    # parameter attributes — same structure, absurd values, like a ×100 mu0 seed —
    # and the forecast must be UNCHANGED. Pre-fix, mu0=5.0 made the level explode.
    m.mu0 = 5.0
    m.ar  = [[0.95]]
    fc_corrupt = list(m.forecast_fuf(6).level[:6])
    assert np.allclose(fc_fit, fc_corrupt), \
        f"forecast_fuf used the stale attributes, not the fit: {fc_fit} vs {fc_corrupt}"
