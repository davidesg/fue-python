#!/usr/bin/env python3
"""BUG-0028: the HTML forecast page (report_forecast) still multiplies the
variations by 100 whatever lambda -- the defect BUG-0027 fixed in the .out.

The page's table and its chart take 100 * delta(refactor * BoxCox(y)) / refactor
for the period and annual variations, their stds, ERR and the level's std
(`level_std * 100`). With lambda = 0 that is a percentage. With lambda = 1 it is
the change of the level times 100: a series about 100 that moves by 1 reads as
moving by 100, under a column headed "Monthly (%)".

A synthetic monthly series about 100 in levels (an AR(1), innovations of sd 1),
fitted by fue with lambda = 1 and refactor = 100; the rows the page renders
(`_table_data`, the template's input). Exits 1 while the last observed row's
variation is 100 times the change of the level; 0 when it is the change.
"""
import sys
import warnings

import numpy as np

import fue
from fue.report_forecast import _table_data

warnings.filterwarnings("ignore")


def main():
    rng = np.random.default_rng(7)
    n = 160
    y = np.zeros(n)
    y[0] = 100.0
    for t in range(1, n):
        y[t] = 100.0 + 0.6 * (y[t - 1] - 100.0) + rng.normal(0, 1.0)
    ts = fue.TimeSeries(list(map(float, y)), freq=12, start=(2000, 1), name="NIVELES")
    m = fue.Model(ts, d=0, D=0, boxlam=1.0, ar=[[0.1]], ar_free=[[True]],
                  mu=100.0 * float(y.mean()), estimate_mu=True, refactor=100.0)
    m.fit()
    fr = m.forecast(24)
    hist, fore, _, meta = _table_data(m, fr)
    shown = float(hist[-1]["diff1"])
    change = y[-1] - y[-2]
    print(f"last observation: {y[-2]:.4f} -> {y[-1]:.4f}, change {change:+.4f}")
    print(f"the page's '{meta['diff1_label']}' column: {shown:+.2f}")
    print(f"first forecast: level {fore[0]['level']}, its std {fore[0]['level_std']}, "
          f"variation {fore[0]['diff1']} (std {fore[0]['diff1_std']})")
    abierto = abs(shown) > 10 * abs(change) + 0.01
    print("x100 in levels" if abierto else "the change of the level")
    return 1 if abierto else 0


if __name__ == "__main__":
    sys.exit(main())
