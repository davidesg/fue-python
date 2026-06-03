"""
Parser for fue .inp model-data files.

Reproduces the reading order of fue-1.13.1/src/fue.c [sections 3.0–3.7].
Sections that are commented-out in fue.c (seasonal AR/MA, annual f-fixed)
are silently consumed and discarded so that both old and new .inp files work.
"""

import numpy as np
from .series import TimeSeries
from .model import Model, FixedFreqFactor
from .intervention import Intervention


# ── Public API ────────────────────────────────────────────────────────────────

def load(path):
    """
    Parse a fue .inp file and return (TimeSeries, Model).

    The returned Model is unfitted; call .fit() to estimate parameters.

    Parameters
    ----------
    path : str or path-like
        Path to the .inp file (with or without the .inp extension).

    Returns
    -------
    ts : TimeSeries
    model : Model  (unfitted)
    """
    path = str(path)
    if not path.endswith(".inp") and not path.endswith(".pre"):
        path += ".inp"
    return _InpParser(path).parse()


# ── Internal parser ───────────────────────────────────────────────────────────

class _InpParser:
    """
    Sequential parser that mirrors fue.c's fscanf / fgets reading order.

    Internally, lines are classified as:
    - 'sep'  : lines starting with '**' (section headers / unnamed separators)
    - 'data' : all other non-blank, non-pure-comment lines
    Lines starting with a single '*' (not '**') are treated as pure header
    comments and completely skipped.
    """

    def __init__(self, path):
        blocks = []
        with open(path) as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith("**"):
                    key = s.lstrip("*").strip().lower()
                    blocks.append(("sep", key))
                elif s.startswith("*"):
                    pass  # pure header comment
                else:
                    blocks.append(("data", s.split()))
        self._b = blocks
        self._p = 0

    # ── Low-level helpers ─────────────────────────────────────────────────────

    def _next_data(self):
        """Return tokens of the next data line, skipping all separators."""
        while self._p < len(self._b):
            kind, val = self._b[self._p]
            self._p += 1
            if kind == "data":
                return val
        raise ValueError("Unexpected end of .inp file")

    def _skip_sep(self):
        """Consume (and return the key of) the next separator line."""
        while self._p < len(self._b):
            kind, val = self._b[self._p]
            self._p += 1
            if kind == "sep":
                return val
        return None

    def _peek_sep(self):
        """Return True if the next item is a separator."""
        return self._p < len(self._b) and self._b[self._p][0] == "sep"

    def _peek_key(self):
        """Return the key of the next separator without consuming it, or ''."""
        for k in range(self._p, len(self._b)):
            if self._b[k][0] == "sep":
                return self._b[k][1]
        return ""

    # ── Section readers ───────────────────────────────────────────────────────

    def _read_factor_block(self, nfactors, orders):
        """
        Read nfactors AR or MA factor blocks.  Each block has:
          **            (separator consumed by caller already for the first block)
          coef1 free1
          coef2 free2   (one pair per coefficient)

        Returns (factors, free_lists) where each element is a list of floats/bools.
        """
        factors = []
        free_lists = []
        for i in range(nfactors):
            self._skip_sep()
            order = orders[i]
            coefs, frees = [], []
            # read 'order' coef/free pairs on potentially multiple lines
            remaining = order
            while remaining > 0:
                toks = self._next_data()
                for j in range(0, len(toks), 2):
                    if remaining <= 0:
                        break
                    coefs.append(float(toks[j]))
                    frees.append(bool(int(toks[j + 1])) if j + 1 < len(toks) else False)
                    remaining -= 1
            factors.append(coefs)
            free_lists.append(frees)
        return factors, free_lists

    def _read_arma_section(self):
        """
        Read the next AR or MA section (count + orders + coefficient blocks).
        Returns (factors, free_lists, orders) or ([], [], []) if count is 0.
        """
        toks = self._next_data()
        count = int(toks[0])
        if count == 0:
            return [], [], []
        orders = [int(toks[i + 1]) for i in range(count)]
        factors, free_lists = self._read_factor_block(count, orders)
        return factors, free_lists, orders

    def _read_ffixed_section(self):
        """
        Read a fixed-frequency AR(2) or MA(2) section.
        Format (active in fue-1.13.1):
          count freq1 freq2 ...
          ** coef1 free1
          ** coef2 free2
          ...
        Returns list of FixedFreqFactor.
        """
        toks = self._next_data()
        count = int(toks[0])
        if count == 0:
            return []
        freqs = [float(toks[i + 1]) for i in range(count)]
        result = []
        for i in range(count):
            self._skip_sep()
            coef_toks = self._next_data()
            coef = float(coef_toks[0])
            free = bool(int(coef_toks[1]))
            result.append(FixedFreqFactor(freq=freqs[i], coef=coef, free=free))
        return result

    def _skip_ffixed_section(self):
        """Skip a f-fixed section (commented-out in fue-1.13.1)."""
        toks = self._next_data()
        count = int(toks[0])
        for _ in range(count):
            self._skip_sep()
            self._next_data()

    # ── Main parser ───────────────────────────────────────────────────────────

    def parse(self):
        # [3.0] Skip first 5 lines → in our representation that's at most
        # 5 separator/data blocks before the Frequency section.
        # We find the frequency section by key.
        while self._p < len(self._b):
            kind, val = self._b[self._p]
            if kind == "sep" and "frequency" in val:
                self._p += 1   # consume the separator
                break
            self._p += 1

        # [3.1] Frequency
        freq_toks = self._next_data()
        freq_str = freq_toks[0].lower()
        if freq_str == "number":
            freq, numbering = 1, True
        else:
            freq, numbering = int(freq_str), False

        # [3.1] Observations + start date + name
        self._skip_sep()
        obs_toks = self._next_data()
        nobs = int(obs_toks[0])
        if freq > 1:
            begtime = int(obs_toks[1])
            begyear = int(obs_toks[2])
            name    = obs_toks[3] if len(obs_toks) > 3 else "series"
        else:
            # annual: field order is nobs, outyear (ignored), begyear, name
            begyear = int(obs_toks[2])
            begtime = 1
            name    = obs_toks[3] if len(obs_toks) > 3 else "series"

        # [3.2] Number of deterministic variables
        self._skip_sep()
        ndet = int(self._next_data()[0])

        interventions = []
        custom_col_order = []

        if ndet > 0:
            # [3.2.0] Read intervention type + date for each detvar
            self._skip_sep()
            raw_types = []        # (type_str, at_1based, harmonic, is_custom)
            for _ in range(ndet):
                toks = self._next_data()
                t = toks[0].lower()
                if t in ("impulse", "pulse", "compimp"):
                    py_type = "pulse"
                    at_1 = _date_to_obs(begyear, begtime, freq, toks[1:])
                    raw_types.append((py_type, at_1, None, False))
                elif t == "step":
                    at_1 = _date_to_obs(begyear, begtime, freq, toks[1:])
                    raw_types.append(("step", at_1, None, False))
                elif t == "ramp":
                    at_1 = _date_to_obs(begyear, begtime, freq, toks[1:])
                    raw_types.append(("ramp", at_1, None, False))
                elif t == "cos":
                    raw_types.append(("cos", None, float(toks[1]), False))
                elif t == "sin":
                    raw_types.append(("sin", None, float(toks[1]), False))
                elif t in ("alter", "easter", "trend"):
                    raw_types.append((t, None, None, False))
                else:
                    # Non-standard (external) deterministic variable — data
                    # comes as extra columns in the time series block.
                    custom_col_order.append(len(raw_types))
                    raw_types.append(("custom", None, None, True))

            # [3.2.2] Omega orders for all detvars (one line)
            self._skip_sep()
            omega_orders = [int(x) for x in self._next_data()[:ndet]]

            # [3.2.2] Omega coefficients (one ** block per detvar)
            all_omega = []
            all_omega_free = []
            for s in omega_orders:
                self._skip_sep()
                coefs, frees = [], []
                remaining = s + 1   # s+1 pairs: omega[0..s]
                while remaining > 0:
                    toks = self._next_data()
                    for j in range(0, len(toks), 2):
                        if remaining <= 0:
                            break
                        coefs.append(float(toks[j]))
                        frees.append(bool(int(toks[j + 1])) if j + 1 < len(toks) else False)
                        remaining -= 1
                all_omega.append(coefs)
                all_omega_free.append(frees)

            # [3.2.3] Delta orders for all detvars
            self._skip_sep()
            delta_orders = [int(x) for x in self._next_data()[:ndet]]

            # [3.2.3] Delta coefficients (only for detvars with Ndelta > 0)
            all_delta = [[] for _ in range(ndet)]
            all_delta_free = [[] for _ in range(ndet)]
            for i, r in enumerate(delta_orders):
                if r > 0:
                    self._skip_sep()
                    coefs, frees = [], []
                    remaining = r
                    while remaining > 0:
                        toks = self._next_data()
                        for j in range(0, len(toks), 2):
                            if remaining <= 0:
                                break
                            coefs.append(float(toks[j]))
                            frees.append(bool(int(toks[j + 1])) if j + 1 < len(toks) else False)
                            remaining -= 1
                    all_delta[i] = coefs
                    all_delta_free[i] = frees

            # Build Intervention objects (data for custom filled in after reading)
            for i, (py_type, at_1, harm, is_custom) in enumerate(raw_types):
                if is_custom:
                    # Placeholder — data array assigned after reading the data block
                    interventions.append(Intervention(
                        "custom", at=0, data=np.zeros(nobs),
                        omega=all_omega[i], omega_free=all_omega_free[i],
                        delta=all_delta[i], delta_free=all_delta_free[i],
                    ))
                elif at_1 is not None:
                    at_0 = at_1 - 1
                    interventions.append(Intervention(
                        py_type, at=at_0,
                        omega=all_omega[i], omega_free=all_omega_free[i],
                        delta=all_delta[i], delta_free=all_delta_free[i],
                    ))
                else:
                    # cos / sin / alter / trend / easter
                    interventions.append(Intervention(
                        py_type, at=0,
                        omega=all_omega[i], omega_free=all_omega_free[i],
                        delta=all_delta[i], delta_free=all_delta_free[i],
                        harmonic=harm if harm is not None else 0.0,
                    ))

        # [3.3.1-3.3.2] Regular AR
        self._skip_sep()
        ar, ar_free, _ = self._read_arma_section()

        # Skip optional seasonal AR section (commented-out in fue-1.13.1)
        if "seasonal ar" in self._peek_key() or "frequencies of seasonal ar" in self._peek_key():
            self._skip_sep()
            self._skip_ffixed_section()

        # [3.3.3-3.3.4] Annual / seasonal AR
        self._skip_sep()
        ar_s, ar_s_free, _ = self._read_arma_section()

        # [3.3.5-3.3.6] Regular MA
        self._skip_sep()
        ma, ma_free, _ = self._read_arma_section()

        # Skip optional seasonal MA section
        if "seasonal ma" in self._peek_key() or "frequencies of seasonal ma" in self._peek_key():
            self._skip_sep()
            self._skip_ffixed_section()

        # [3.3.7-3.3.8] Annual / seasonal MA
        self._skip_sep()
        ma_s, ma_s_free, _ = self._read_arma_section()

        # [3.3.9-3.3.10] f-fixed regular AR
        self._skip_sep()
        ar_f = self._read_ffixed_section()

        # Skip optional f-fixed annual AR (commented-out in fue-1.13.1)
        if "anual ar" in self._peek_key() or "annual ar" in self._peek_key():
            self._skip_sep()
            self._skip_ffixed_section()

        # [3.3.13-3.3.14] f-fixed regular MA
        self._skip_sep()
        ma_f = self._read_ffixed_section()

        # Skip optional f-fixed annual MA
        if "anual ma" in self._peek_key() or "annual ma" in self._peek_key():
            self._skip_sep()
            self._skip_ffixed_section()

        # [3.4] Mean parameter
        self._skip_sep()
        mu_toks = self._next_data()
        mu = float(mu_toks[0])
        estimate_mu = bool(int(mu_toks[1])) if len(mu_toks) > 1 else False

        # [3.5] Box-Cox lambda + differences
        self._skip_sep()
        bc_toks = self._next_data()
        boxlam = float(bc_toks[0])
        d      = int(bc_toks[1])
        D      = int(bc_toks[2])

        # [3.6] Individual annual-difference factors
        self._skip_sep()
        toks = self._next_data()
        if freq > 1:
            n_ifadf = freq // 2 + 1
            ifadf = [int(toks[k]) for k in range(min(n_ifadf, len(toks)))]
            ifadf += [0] * (n_ifadf - len(ifadf))
        else:
            ifadf = []

        # [3.7] cbands + refactor
        self._skip_sep()
        rf_toks = self._next_data()
        refactor = float(rf_toks[1]) if len(rf_toks) > 1 else 1.0
        if refactor == 0.0:
            refactor = 1.0

        # [3.7] Time series data (col 0 = stochastic; extra cols = custom detvars)
        self._skip_sep()
        n_custom = len(custom_col_order)
        custom_data = [[] for _ in range(n_custom)]
        data = []
        while len(data) < nobs:
            toks = self._next_data()
            data.append(float(toks[0]))
            for k in range(n_custom):
                col_idx = k + 1
                custom_data[k].append(float(toks[col_idx]) if col_idx < len(toks) else 0.0)

        # Fill in the custom intervention data arrays
        for k, itv_idx in enumerate(custom_col_order):
            interventions[itv_idx].data = np.array(custom_data[k], dtype=float)

        ts = TimeSeries(
            np.array(data),
            freq=freq,
            start=(begyear, begtime),
            name=name,
        )
        model = Model(
            ts,
            ar=ar,           ar_free=ar_free if ar else None,
            ma=ma,           ma_free=ma_free if ma else None,
            ar_s=ar_s,       ar_s_free=ar_s_free if ar_s else None,
            ma_s=ma_s,       ma_s_free=ma_s_free if ma_s else None,
            ar_f=ar_f,
            ma_f=ma_f,
            d=d, D=D, ifadf=ifadf,
            interventions=interventions,
            mu=mu, estimate_mu=estimate_mu,
            boxlam=boxlam, refactor=refactor,
        )
        return ts, model


# ── Date helpers ──────────────────────────────────────────────────────────────

def _date_to_obs(begyear, begtime, freq, arg_tokens):
    """
    Reproduce fue.c's DateToObs: convert (year [, subperiod]) from the .inp
    token list to a 1-based observation index.

    For annual data (freq==1): arg_tokens = [year]
    For sub-annual data:       arg_tokens = [subperiod, year]
    """
    if freq == 1:
        year = int(arg_tokens[0])
        sub  = 1
    else:
        sub  = int(arg_tokens[0])
        year = int(arg_tokens[1])

    # Exact port of diagnose.c DateToObs:
    srest = freq - begtime + 1
    if sub == freq:
        pcad   = year - begyear
        obs_no = srest + freq * pcad
    else:
        pcad   = year - begyear - 1
        obs_no = srest + freq * pcad + sub
    return obs_no
