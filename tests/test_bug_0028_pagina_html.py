"""BUG-0028: the HTML forecast page says a percentage only with logs."""
import numpy as np
import pytest

import fue
from fue.report_forecast import _table_data


def _modelo(lam, f, seed=7, n=160):
    rng = np.random.default_rng(seed)
    y = np.zeros(n)
    y[0] = 100.0
    for t in range(1, n):
        y[t] = 100.0 + 0.6 * (y[t - 1] - 100.0) + rng.normal(0, 1.0)
    ts = fue.TimeSeries(list(map(float, y)), freq=12, start=(2000, 1), name="Y")
    w = fue.boxcox(y, lam, f)
    m = fue.Model(ts, d=0, D=0, boxlam=lam, ar=[[0.1]], ar_free=[[True]],
                  mu=float(np.mean(w)), estimate_mu=True, refactor=f)
    m.fit()
    return y, m, m.forecast(24)


@pytest.mark.parametrize("f", [1.0, 100.0])
def test_in_levels_the_variation_is_the_change(f):
    y, m, fr = _modelo(1.0, f)
    hist, fore, _, meta = _table_data(m, fr)
    assert float(hist[-1]["diff1"]) == pytest.approx(y[-1] - y[-2], abs=0.006)
    assert float(hist[-1]["annual"]) == pytest.approx(y[-1] - y[-13], abs=0.006)
    sigma = float(np.sqrt(fr.sigma2)) / f
    assert float(fore[0]["level_std"]) == pytest.approx(sigma, abs=0.006)
    assert float(fore[0]["diff1"]) == pytest.approx(float(fore[0]["level"]) - y[-1], abs=0.011)
    assert meta["percent"] is False and "%" not in meta["diff1_label"] + meta["annual_label"]


@pytest.mark.parametrize("f", [1.0, 100.0])
def test_in_logs_it_is_a_percentage_whatever_the_factor(f):
    y, m, fr = _modelo(0.0, f)
    hist, fore, _, meta = _table_data(m, fr)
    assert float(hist[-1]["diff1"]) == pytest.approx(100 * np.log(y[-1] / y[-2]), abs=0.006)
    assert meta["percent"] is True and meta["diff1_label"] == "Monthly (%)"
    assert 0.5 < float(fore[0]["level_std"]) < 2.0          # about 1 % of 100


def test_the_page_is_written_and_headed(tmp_path):
    pytest.importorskip("jinja2")
    pytest.importorskip("matplotlib")
    for lam, palabra in ((1.0, "Monthly change"), (0.0, "Monthly (%)")):
        _, m, fr = _modelo(lam, 100.0)
        p = tmp_path / f"f{lam}.html"
        fue.write_forecast_report(m, fr, str(p))
        html = p.read_text(encoding="utf-8")
        assert palabra in html
        assert ("Std&#160;(%)" in html) is (lam == 0.0)
