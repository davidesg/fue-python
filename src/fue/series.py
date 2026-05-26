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

    # ── Properties ────────────────────────────────────────────────────────

    @property
    def nobs(self):
        return len(self.data)

    @property
    def numbering(self):
        """True when freq=0 (plain observation indices, no calendar dates)."""
        return self.freq == 0

    def __len__(self):
        return self.nobs

    def __repr__(self):
        return (f"TimeSeries(name={self.name!r}, nobs={self.nobs}, "
                f"freq={self.freq}, start={self.start})")
