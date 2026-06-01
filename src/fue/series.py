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

    def __len__(self):
        return self.nobs

    def __repr__(self):
        return (f"TimeSeries(name={self.name!r}, nobs={self.nobs}, "
                f"freq={self.freq}, start={self.start})")
