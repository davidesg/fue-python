"""TimeSeries: lightweight wrapper around a numpy array with date metadata."""

import numpy as np


class TimeSeries:
    """
    A univariate time series with frequency and start-date metadata.

    Parameters
    ----------
    data : array-like
        Observations in chronological order.
    freq : int
        Observations per year: 1 (annual), 4 (quarterly), 12 (monthly).
        Pass 0 to use plain observation numbering.
    start : tuple (year, period)
        First observation.  Period is 1-based (Jan = 1 for monthly).
    name : str, optional
        Series label used in summaries and plots.
    """

    def __init__(self, data, freq=12, start=(1900, 1), name="series"):
        self.data  = np.asarray(data, dtype=float)
        self.freq  = int(freq)
        self.start = tuple(start)   # (year, period)
        self.name  = str(name)

    # ── Convenience constructors ──────────────────────────────────────────

    @classmethod
    def from_array(cls, data, freq=12, start=(1900, 1), name="series"):
        return cls(data, freq=freq, start=start, name=name)

    @classmethod
    def from_csv(cls, path, freq=12, start=(1900, 1), name=None,
                 column=0, **read_csv_kwargs):
        import pandas as pd
        df = pd.read_csv(path, **read_csv_kwargs)
        col = df.iloc[:, column] if isinstance(column, int) else df[column]
        label = name or (col.name if hasattr(col, "name") else "series")
        return cls(col.values, freq=freq, start=start, name=label)

    @classmethod
    def from_pandas(cls, series, freq=None, name=None):
        """
        Build a TimeSeries from a ``pandas.Series``.

        If the series has a ``DatetimeIndex`` or ``PeriodIndex``, freq and
        start are inferred automatically (annual / quarterly / monthly).
        Pass *freq* explicitly to override.

        Parameters
        ----------
        series : pandas.Series
        freq   : int, optional
            1, 4, or 12.  Inferred from index if omitted.
        name   : str, optional
            Defaults to ``series.name``.
        """
        import pandas as pd

        label = name or (str(series.name) if series.name is not None else "series")
        data  = series.to_numpy(dtype=float)
        idx   = series.index

        if freq is None:
            if hasattr(idx, "freqstr") and idx.freqstr:
                fs = idx.freqstr.upper()
                if fs.startswith("A") or fs.startswith("Y"):
                    freq = 1
                elif fs.startswith("Q"):
                    freq = 4
                elif fs.startswith("M"):
                    freq = 12
                else:
                    freq = 1
            else:
                freq = 1

        # Derive start from the first index value
        if isinstance(idx, pd.PeriodIndex):
            p0 = idx[0]
            year = p0.year
            if freq == 4:
                period = p0.quarter
            elif freq == 12:
                period = p0.month
            else:
                period = 1
        elif isinstance(idx, pd.DatetimeIndex):
            d0 = idx[0]
            year = d0.year
            if freq == 4:
                period = (d0.month - 1) // 3 + 1
            elif freq == 12:
                period = d0.month
            else:
                period = 1
        else:
            year, period = 1900, 1

        return cls(data, freq=freq, start=(year, period), name=label)

    # ── Properties ────────────────────────────────────────────────────────

    @property
    def nobs(self):
        return len(self.data)

    @property
    def numbering(self):
        """True when freq=0 (plain observation indices, no calendar dates)."""
        return self.freq == 0

    # ── Plotting ──────────────────────────────────────────────────────────

    def plot(self, title=None, ax=None):
        """Time-series line plot with calendar x-axis."""
        from .plots import plot_series
        plot_series(self, title=title, ax=ax)

    def plot_acf(self, lags=24, confidence=0.95, ax=None):
        """Autocorrelation function stem plot."""
        from .plots import plot_acf
        plot_acf(self.data, lags=lags, confidence=confidence,
                 title=f"ACF — {self.name}", ax=ax)

    def plot_pacf(self, lags=24, confidence=0.95, ax=None):
        """Partial autocorrelation function stem plot."""
        from .plots import plot_pacf
        plot_pacf(self.data, lags=lags, confidence=confidence,
                  title=f"PACF — {self.name}", ax=ax)

    # ── Descriptive statistics ────────────────────────────────────────────

    def describe(self):
        """
        Sample statistics matching fue's File_StatSer output.

        Uses population moments (divisor n) to match the C implementation.
        Returns the formatted string (also printed to stdout).
        """
        import math
        x   = self.data
        n   = self.nobs
        mu  = float(x.mean())
        std = float(x.std(ddof=0))         # population std (C uses sum/n)
        se  = std / math.sqrt(n)
        skew = (((x - mu) / std) ** 3).mean() if std > 1e-20 else 0.0
        kurt = (((x - mu) / std) ** 4).mean() - 3.0 if std > 1e-20 else 0.0
        jb   = (n // 6) * (skew ** 2 + kurt ** 2 / 4.0)   # C uses integer n/6

        imin = int(x.argmin())             # 0-based
        imax = int(x.argmax())
        ey, ep = self._obs_to_date(imin + 1)
        ay, ap = self._obs_to_date(imax + 1)

        freq = self.freq if self.freq > 0 else 1
        begyear, begtime = self.start
        endyear, endtime = self._obs_to_date(n)

        if freq > 1:
            span = f"from {begtime}/{begyear} to {endtime}/{endyear}"
            min_at = f"at {ep:2d}/{ey} (observation {imin + 1:3d})"
            max_at = f"at {ap:2d}/{ay} (observation {imax + 1:3d})"
        else:
            span = f"from {begyear} to {endyear}"
            min_at = f"at {ey} (observation {imin + 1:3d})"
            max_at = f"at {ay} (observation {imax + 1:3d})"

        lines = [
            f"{self.name} (seasonal period: {freq})",
            f"{n} observations: {span}",
            "",
            f"{'Mean':>24s}: {mu:18.6f}",
            f"{'Standard error of mean':>24s}: {se:18.6f}",
            f"{'Variance':>24s}: {std**2:18.6f}",
            f"{'Standard deviation':>24s}: {std:18.6f}",
            f"{'Skewness':>24s}: {skew:18.6f}",
            f"{'Kurtosis':>24s}: {kurt:18.6f}",
            f"{'Jarque-Bera':>24s}: {jb:18.6f}",
            f"{'Minimum':>24s}: {x[imin]:18.6f}  {min_at}",
            f"{'Maximum':>24s}: {x[imax]:18.6f}  {max_at}",
        ]
        out = "\n".join(lines)
        print(out)
        return out

    def _obs_to_date(self, obs_1based):
        """Convert 1-based observation number to (year, period)."""
        freq = self.freq if self.freq > 0 else 1
        begyear, begtime = self.start
        total = begyear * freq + (begtime - 1) + (obs_1based - 1)
        return total // freq, total % freq + 1

    def __len__(self):
        return self.nobs

    def __repr__(self):
        return (f"TimeSeries(name={self.name!r}, nobs={self.nobs}, "
                f"freq={self.freq}, start={self.start})")
