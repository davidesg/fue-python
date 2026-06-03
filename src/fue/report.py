"""
Generate .out-style estimation reports matching fue's ASCII output format.
"""

import math
import numpy as np


# ── Public API ─────────────────────────────────────────────────────────────────

def write_out(model, path=None):
    """
    Generate an estimation report in fue .out format.

    Parameters
    ----------
    model : Model  (must be fitted)
    path : str or None
        If given, write to this file path.  If None, return the text.

    Returns
    -------
    str
    """
    if model._result is None:
        raise RuntimeError("Model has not been fitted. Call .fit() first.")
    text = "\n".join(_build_report(model))
    if path is not None:
        with open(path, "w") as fh:
            fh.write(text + "\n")
    return text


# ── Top-level report builder ───────────────────────────────────────────────────

def _build_report(model):
    r      = model._result
    ts     = model.series
    freq   = ts.freq if ts.freq > 0 else 1
    ornsop = model.d + model.D * freq
    n_eff  = ts.nobs - ornsop
    sigma2 = r.sigma2
    sigma  = math.sqrt(sigma2)
    res    = r.residuals          # length n_eff

    fitted = _extract_fitted(model, r)

    lines = []
    _section_params(lines, model, r, fitted)
    _section_boxcox(lines, model, freq)
    _section_arma_ops(lines, model, fitted, freq)
    _section_sigma(lines, r, sigma, sigma2, n_eff)
    _section_hq(lines, model, sigma2, n_eff)
    _section_matrices(lines, r)
    _section_residual_stats(lines, ts, res, freq, ornsop)
    _section_plot(lines, ts, res, freq, ornsop)
    _section_histogram(lines, res)
    _section_corr(lines, model, res, n_eff, freq)
    return lines


# ── Parameter extraction ───────────────────────────────────────────────────────

def _extract_fitted(model, r):
    """
    Walk the free-parameter list in the same order as count_npar_build_par
    and return a dict with fitted arrays for AR/MA/omega/delta/etc.
    """
    idx = [0]   # mutable so inner function can advance it

    def _next():
        v = r.params[idx[0]]
        idx[0] += 1
        return v

    def _next_se():
        return r.std_errors[idx[0] - 1]

    # omegas and deltas: collect (value, se) pairs
    omega_vals = []   # list of lists
    omega_se   = []
    delta_vals = []
    delta_se   = []

    for itv in model.interventions:
        ov, os_ = [], []
        for w, wf in zip(itv.omega, itv.omega_free):
            if wf:
                idx[0] += 1   # advance
                ov.append(r.params[idx[0] - 1])
                os_.append(r.std_errors[idx[0] - 1])
            else:
                ov.append(w)
                os_.append(None)
        omega_vals.append(ov)
        omega_se.append(os_)

    for itv in model.interventions:
        dv, ds = [], []
        for d, df in zip(itv.delta, itv.delta_free):
            if df:
                idx[0] += 1
                dv.append(r.params[idx[0] - 1])
                ds.append(r.std_errors[idx[0] - 1])
            else:
                dv.append(d)
                ds.append(None)
        delta_vals.append(dv)
        delta_se.append(ds)

    def _fit_factor_list(factors, free_lists):
        result, result_se = [], []
        for i, factor in enumerate(factors):
            fl = (free_lists[i]
                  if free_lists is not None else [True] * len(factor))
            fv, fs = [], []
            for phi, pf in zip(factor, fl):
                if pf:
                    idx[0] += 1
                    fv.append(r.params[idx[0] - 1])
                    fs.append(r.std_errors[idx[0] - 1])
                else:
                    fv.append(phi)
                    fs.append(None)
            result.append(fv)
            result_se.append(fs)
        return result, result_se

    ar,    ar_se    = _fit_factor_list(model.ar,   model.ar_free)
    ar_s,  ar_s_se  = _fit_factor_list(model.ar_s, model.ar_s_free)
    ma,    ma_se    = _fit_factor_list(model.ma,   model.ma_free)
    ma_s,  ma_s_se  = _fit_factor_list(model.ma_s, model.ma_s_free)

    # f-fixed AR/MA (only phi2/theta2 is free; phi1 is derived)
    ar_f_phi2, ar_f_se = [], []
    for ff in model.ar_f:
        if ff.free:
            idx[0] += 1
            ar_f_phi2.append(r.params[idx[0] - 1])
            ar_f_se.append(r.std_errors[idx[0] - 1])
        else:
            ar_f_phi2.append(ff.coef)
            ar_f_se.append(None)

    ma_f_th2, ma_f_se = [], []
    for ff in model.ma_f:
        if ff.free:
            idx[0] += 1
            ma_f_th2.append(r.params[idx[0] - 1])
            ma_f_se.append(r.std_errors[idx[0] - 1])
        else:
            ma_f_th2.append(ff.coef)
            ma_f_se.append(None)

    # mu
    if model.estimate_mu:
        idx[0] += 1
        mu_val = r.params[idx[0] - 1]
        mu_se  = r.std_errors[idx[0] - 1]
    else:
        mu_val = model.mu0
        mu_se  = None

    return {
        'omega_vals': omega_vals, 'omega_se': omega_se,
        'delta_vals': delta_vals, 'delta_se': delta_se,
        'ar': ar,   'ar_se': ar_se,
        'ar_s': ar_s, 'ar_s_se': ar_s_se,
        'ma': ma,   'ma_se': ma_se,
        'ma_s': ma_s, 'ma_s_se': ma_s_se,
        'ar_f_phi2': ar_f_phi2, 'ar_f_se': ar_f_se,
        'ma_f_th2':  ma_f_th2,  'ma_f_se': ma_f_se,
        'mu_val': mu_val, 'mu_se': mu_se,
    }


# ── Section helpers ────────────────────────────────────────────────────────────

def _emit(lines, val, se, par_num):
    if se is not None:
        lines.append(f"{val:14.6f}  ({se:8.6f}) [{par_num:2d}]")
    else:
        lines.append(f"{val:14.6f}")


def _section_params(lines, model, r, fitted):
    par_num = 0

    def _emit_par(val, se):
        nonlocal par_num
        if se is not None:
            par_num += 1
            lines.append(f"{val:14.6f}  ({se:8.6f}) [{par_num:2d}]")
        else:
            lines.append(f"{val:14.6f}")

    # Omegas
    for i, itv in enumerate(model.interventions):
        lines.append(f"Omegas for deterministic variable {i + 1}:")
        for val, se in zip(fitted['omega_vals'][i], fitted['omega_se'][i]):
            _emit_par(val, se)

    # Deltas (only for interventions that have them)
    for i, itv in enumerate(model.interventions):
        if not itv.delta:
            continue
        lines.append(f"Deltas for deterministic variable {i + 1}:")
        for val, se in zip(fitted['delta_vals'][i], fitted['delta_se'][i]):
            _emit_par(val, se)

    # Regular AR
    for i, factor in enumerate(model.ar):
        lines.append(f"Coefficients for regular AR factor {i + 1}:")
        for val, se in zip(fitted['ar'][i], fitted['ar_se'][i]):
            _emit_par(val, se)

    # Annual/seasonal AR
    for i, factor in enumerate(model.ar_s):
        lines.append(f"Coefficients for annual AR factor {i + 1}:")
        for val, se in zip(fitted['ar_s'][i], fitted['ar_s_se'][i]):
            _emit_par(val, se)

    # Regular MA
    for i, factor in enumerate(model.ma):
        lines.append(f"Coefficients for regular MA factor {i + 1}:")
        for val, se in zip(fitted['ma'][i], fitted['ma_se'][i]):
            _emit_par(val, se)

    # Annual/seasonal MA
    for i, factor in enumerate(model.ma_s):
        lines.append(f"Coefficients for annual MA factor {i + 1}:")
        for val, se in zip(fitted['ma_s'][i], fitted['ma_s_se'][i]):
            _emit_par(val, se)

    # f-fixed AR: phi1 (derived, no SE), phi2 (free or fixed)
    sper = model.series.freq if model.series.freq > 0 else 1
    for i, ff in enumerate(model.ar_f):
        lines.append(
            f"Coefficients for regular f-fixed AR factor {i + 1}"
            f" [f = {ff.freq:.1f}]:"
        )
        phi2 = fitted['ar_f_phi2'][i]
        phi1 = 2.0 * math.cos(2.0 * math.pi * ff.freq / sper) * math.sqrt(-phi2)
        lines.append(f"{phi1:14.6f}")
        _emit_par(phi2, fitted['ar_f_se'][i])

    # f-fixed MA
    for i, ff in enumerate(model.ma_f):
        lines.append(
            f"Coefficients for regular f-fixed MA factor {i + 1}"
            f" [f = {ff.freq:.1f}]:"
        )
        th2 = fitted['ma_f_th2'][i]
        th1 = 2.0 * math.cos(2.0 * math.pi * ff.freq / sper) * math.sqrt(-th2)
        lines.append(f"{th1:14.6f}")
        _emit_par(th2, fitted['ma_f_se'][i])

    # Mean parameter (always printed)
    lines.append("Mean parameter (mu):")
    _emit_par(fitted['mu_val'], fitted['mu_se'])
    lines.append("")


def _section_boxcox(lines, model, freq):
    lines.append(f"Box-Cox lambda     : {model.boxlam:4.1f}")
    lines.append(f"Seasonal period    : {freq:2d}")
    lines.append(f"Regular differences: {model.d:2d}")
    lines.append(f"Annual differences : {model.D:2d}")
    lines.append("")


def _section_arma_ops(lines, model, fitted, freq):
    sper = freq

    # Combined AR polynomial
    ar_poly = [1.0]
    for factor in fitted['ar']:
        ar_poly = _poly_mul(ar_poly, [1.0] + [-c for c in factor])
    for factor in fitted['ar_s']:
        ar_poly = _poly_mul(ar_poly, _seasonal_poly(factor, sper))
    for i, ff in enumerate(model.ar_f):
        phi2 = fitted['ar_f_phi2'][i]
        phi1 = 2.0 * math.cos(2.0 * math.pi * ff.freq / sper) * math.sqrt(-phi2)
        ar_poly = _poly_mul(ar_poly, [1.0, -phi1, -phi2])

    lines.append("Coefficients of the Autoregressive operator: ")
    for k in range(1, len(ar_poly)):
        lines.append(f"  phi[{k:2d}]   = {-ar_poly[k]:15.10f}")

    # Combined MA polynomial
    ma_poly = [1.0]
    for factor in fitted['ma']:
        ma_poly = _poly_mul(ma_poly, [1.0] + [-c for c in factor])
    for factor in fitted['ma_s']:
        ma_poly = _poly_mul(ma_poly, _seasonal_poly(factor, sper))
    for i, ff in enumerate(model.ma_f):
        th2 = fitted['ma_f_th2'][i]
        th1 = 2.0 * math.cos(2.0 * math.pi * ff.freq / sper) * math.sqrt(-th2)
        ma_poly = _poly_mul(ma_poly, [1.0, -th1, -th2])

    lines.append("Coefficients of the Moving-Average operator: ")
    for k in range(1, len(ma_poly)):
        lines.append(f"  theta[{k:2d}] = {-ma_poly[k]:15.10f}")
    lines.append("")


def _section_sigma(lines, r, sigma, sigma2, n_eff):
    se_sigma2 = sigma2 * math.sqrt(2.0 / n_eff)
    lines.append(f"sigma2: {sigma2:14.10f} ({se_sigma2:15.10f})")
    lines.append(f"sigma : {sigma:14.10f}")
    lines.append(f"logelf: {r.loglik:14.10f}")
    lines.append("")


def _section_hq(lines, model, sigma2, n_eff):
    nparma = _count_nparma(model)
    ln_s2  = math.log(sigma2)
    factor = 2.0 * (1 + nparma) / n_eff
    hq     = ln_s2 + factor * math.log(math.log(n_eff))
    schwarz = ln_s2 + factor * math.log(n_eff)
    lines.append("Selection Model Criterium:")
    lines.append(f"Hannan-Quinn = {hq:6.2f}")
    lines.append(f"Schwarz   = {schwarz:6.2f}")
    lines.append("")


def _section_matrices(lines, r):
    npar = r.npar
    cov  = r.cov_matrix   # (npar, npar)
    se   = r.std_errors

    lines.append("Estimated covariance matrix:")
    lines.append("")
    for i in range(npar):
        row = f"x[{i + 1:2d}] ->"
        for j in range(i + 1):
            row += f"{cov[i, j]:13.9f}"
        lines.append(row)
    lines.append("")

    lines.append("")
    lines.append("Estimated correlation matrix:")
    lines.append("")
    for i in range(npar):
        row = f"x[{i + 1:2d}] ->"
        for j in range(i + 1):
            corr_ij = cov[i, j] / (se[i] * se[j])
            row += f"{corr_ij:6.2f}"
        lines.append(row)
    lines.append("")

    lines.append("Correlations greater than or equal to 0.7 in absolute value:")
    lines.append("")
    for i in range(npar):
        for j in range(i):
            c = cov[i, j] / (se[i] * se[j])
            if abs(c) >= 0.7:
                lines.append(f"corr[{i + 1:2d}][{j + 1:2d}] ={c:6.2f}")
    lines.append("")


def _section_residual_stats(lines, ts, res, freq, ornsop):
    """Residual statistics block (mirrors TimeSeries.describe())."""
    n   = len(res)
    mu  = float(res.mean())
    std = float(res.std(ddof=0))
    se  = std / math.sqrt(n)
    skew = float((((res - mu) / std) ** 3).mean()) if std > 1e-20 else 0.0
    kurt = float((((res - mu) / std) ** 4).mean() - 3.0) if std > 1e-20 else 0.0
    jb   = (n // 6) * (skew ** 2 + kurt ** 2 / 4.0)   # C uses integer n/6

    imin = int(res.argmin())
    imax = int(res.argmax())

    # Residuals start at observation ornsop+1 of the original series
    res_begyear, res_begtime = ts._obs_to_date(ornsop + 1)
    res_endyear, res_endtime = ts._obs_to_date(ornsop + n)

    def _obs_to_date_res(k_1based):
        """1-based index within residuals → (year, period)."""
        total = res_begyear * freq + (res_begtime - 1) + (k_1based - 1)
        return total // freq, total % freq + 1

    ey, ep = _obs_to_date_res(imin + 1)
    ay, ap = _obs_to_date_res(imax + 1)

    if freq > 1:
        span   = f"from {res_begtime}/{res_begyear} to {res_endtime}/{res_endyear}"
        min_at = f"at {ep:2d}/{ey} (observation {imin + 1:3d})"
        max_at = f"at {ap:2d}/{ay} (observation {imax + 1:3d})"
    else:
        span   = f"from {res_begyear} to {res_endyear}"
        min_at = f"at {ey} (observation {imin + 1:3d})"
        max_at = f"at {ay} (observation {imax + 1:3d})"

    lines.append(f"Unconditional residuals (seasonal period: {freq})")
    lines.append(f"{n} observations: {span}")
    lines.append("")
    lines.append(f"{'Mean':>24s}: {mu:18.6f}")
    lines.append(f"{'Standard error of mean':>24s}: {se:18.6f}")
    lines.append(f"{'Variance':>24s}: {std**2:18.6f}")
    lines.append(f"{'Standard deviation':>24s}: {std:18.6f}")
    lines.append(f"{'Skewness':>24s}: {skew:18.6f}")
    lines.append(f"{'Kurtosis':>24s}: {kurt:18.6f}")
    lines.append(f"{'Jarque-Bera':>24s}: {jb:18.6f}")
    lines.append(f"{'Minimum':>24s}: {res[imin]:18.6f}  {min_at}")
    lines.append(f"{'Maximum':>24s}: {res[imax]:18.6f}  {max_at}")
    lines.append("")


def _section_plot(lines, ts, res, freq, ornsop):
    """Standardized time series plot (File_PlotSer port)."""
    n   = len(res)
    mu  = float(res.mean())
    std = float(res.std(ddof=0))

    if std < 1e-20:
        return

    # AbsMax: largest |z| among all residuals (forced >= 3.0)
    z_all  = (res - mu) / std
    absmax = float(np.abs(z_all).max())
    if absmax <= 2.0:
        absmax = 3.0
    if absmax > 8.0:
        lines.append("Warning: at least one observation above 8 sigmas")
        return

    horinc   = 25.0 / absmax
    bandpos1 = horinc
    bandpos2 = 2.0 * horinc

    # Build Guions and Marcas strings (length 77)
    guions = list("-------------+-------------------------+-------------------------+--------------")
    marcas = list("                                       0                          ")

    for sigma_i in range(1, 9):
        if absmax >= sigma_i:
            pos = _iround(sigma_i * horinc)
            guions[39 - pos] = '+'
            guions[39 + pos] = '+'
            digit = str(sigma_i)
            marcas[39 - pos]     = digit
            marcas[39 - pos - 1] = '-'
            marcas[39 + pos]     = digit
            marcas[39 + pos - 1] = '+'

    guions_s = "".join(guions)
    marcas_s = "".join(marcas)

    lines.append("Standardized time series plot "
                 "(original values on right-side column):")
    lines.append("")
    lines.append(marcas_s)
    lines.append(guions_s)

    res_begyear, res_begtime = ts._obs_to_date(ornsop + 1)

    for i in range(1, n + 1):
        val = res[i - 1]
        z   = (val - mu) / std

        # Date of this residual observation
        total = res_begyear * freq + (res_begtime - 1) + (i - 1)
        aper  = total // freq
        asub  = total % freq + 1

        if freq == 1:
            prefix = f"{i:4d}{aper:7d} "
        else:
            prefix = f"{i:4d}{asub:3d}/{aper:4d}"

        # Build 55-char Tmpstr (indices 0-54)
        bar = [' '] * 55
        is_outlier = abs(z) >= 2.0

        # Period boundary marker
        if freq != 1 and asub == freq:
            bar[1]  = '+'
            bar[53] = '+'
        else:
            bar[1]  = '|'
            bar[53] = '|'

        if is_outlier:
            bar[0]  = '@'
            bar[54] = '@'

        # Star position
        idx = 27 + _iround(z * horinc)
        if 0 <= idx <= 54:
            bar[idx] = '*'

        # Center marker (only if position 27 is still space)
        if bar[27] == ' ':
            bar[27] = '|'

        # Band markers
        for bp in (bandpos1, bandpos2):
            p_plus  = 27 + _iround(bp)
            p_minus = 27 - _iround(bp)
            if 0 <= p_plus  <= 54 and bar[p_plus]  == ' ':
                bar[p_plus]  = ':'
            if 0 <= p_minus <= 54 and bar[p_minus] == ' ':
                bar[p_minus] = ':'

        bar_s = "".join(bar)
        lines.append(f"{prefix}{bar_s}{val:13.10f}")

    lines.append(guions_s)
    lines.append(marcas_s)
    lines.append("")

    # Outlier table
    lines.append("                 +------------------------------------------+")
    lines.append("                 |       Table of standardized values       |")
    lines.append("                 |       greater than or equal to 2.0       |")
    lines.append("                 +------------------------------------------+")
    lines.append("                 |                                          |")
    lines.append("                 | Observation    Date   Standardized value |")
    lines.append("                 |                                          |")

    for i in range(1, n + 1):
        val = res[i - 1]
        z   = (val - mu) / std
        if abs(z) >= 2.0:
            total = res_begyear * freq + (res_begtime - 1) + (i - 1)
            aper  = total // freq
            asub  = total % freq + 1
            if freq == 1:
                lines.append(
                    f"                 |{i:7d}{aper:13d} {z:13.2f}        |"
                )
            else:
                lines.append(
                    f"                 |{i:7d}{asub:9d}/{aper:4d} {z:13.2f}        |"
                )

    lines.append("                 +------------------------------------------+")
    lines.append("")


def _section_histogram(lines, res):
    """Standardized histogram (File_HistSer port)."""
    n   = len(res)
    mu  = float(res.mean())
    std = float(res.std(ddof=0))

    if std < 1e-20:
        return

    z_all  = (res - mu) / std
    xmax_raw = float(np.abs(z_all).max())

    if xmax_raw > 8.0:
        lines.append("Warning: at least one observation above 8 sigmas")
        return
    xmax = 4.0 if xmax_raw <= 4.0 else 8.0

    if xmax == 4.0:
        nphor = 4
        no_s  = "    "
        yes_s = "...."
        base1 = "        +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+"
        base2 = "       -4      -3      -2      -1       0      +1      +2      +3      +4"
    else:
        nphor = 2
        no_s  = "  "
        yes_s = ".."
        base1 = "        +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
        base2 = "       -8  -7  -6  -5  -4  -3  -2  -1   0  +1  +2  +3  +4  +5  +6  +7  +8"

    # Breakpoints: from -xmax+0.5 to xmax in steps of 0.5
    bk = []
    b  = -xmax + 0.5
    while b <= xmax + 1e-10:
        bk.append(b)
        b += 0.5
    num_cat = len(bk)

    # Frequencies
    freqs = [0] * (num_cat + 1)
    atip1 = 0
    atip2 = 0
    for zi in z_all:
        if zi <= bk[0]:
            freqs[0] += 1
        else:
            for j in range(1, num_cat):
                if bk[j - 1] < zi <= bk[j]:
                    freqs[j] += 1
                    break
        if abs(zi) >= 1.0:
            atip1 += 1
        if abs(zi) >= 2.0:
            atip2 += 1

    # Trim to num_cat entries
    freqs = freqs[:num_cat]
    fmax  = max(freqs) if freqs else 1
    obs_per_fil = fmax / 16.0

    NUM_FIL = 17
    NUM_COL = 64

    # Build shist rows (index 0 = top/numbers, 1..16 = dot rows)
    shist = [None] * NUM_FIL

    # Rows 2..17 (j in 2..17) → shist[j-1] (indices 1..16)
    for j in range(2, NUM_FIL + 1):
        row = []
        for i in range(num_cat):
            if freqs[i] > obs_per_fil * (NUM_FIL - j):
                row.append(yes_s)
            else:
                row.append(no_s)
        s = "".join(row)
        shist[j - 1] = s[: NUM_COL - 1] + "|"

    # Build aux[0..15]: frequency label strings, one slot per j (j=2..17 → aux[0..15]).
    # Label for category i goes into aux[j_first-2] where j_first is the first j
    # (lowest threshold) where freqs[i] > ObsPerFil*(NumFil-j_first).
    # aux[0] → shist[0] (numbers row, above the tallest bars)
    # aux[k] → merged into shist[k] for k=1..15

    def _freq_label(f):
        s = str(f)
        if nphor == 4:
            if len(s) == 1:
                return "  " + s + " "
            elif len(s) == 2:
                return " " + s + " "
            elif len(s) == 3:
                return s + " "
            else:
                return s[:4]
        else:
            return (s + " ")[:2]

    aux  = [list(no_s * num_cat) for _ in range(NUM_FIL - 1)]   # 16 rows of blanks
    chk  = [False] * num_cat
    for j in range(2, NUM_FIL + 1):
        aux_row = aux[j - 2]  # 0-indexed
        col_pos = 0
        for i in range(num_cat):
            if freqs[i] > obs_per_fil * (NUM_FIL - j) and not chk[i]:
                chk[i] = True
                label = _freq_label(freqs[i])
                # Write label into aux row at col_pos
                for ci, ch in enumerate(label):
                    if col_pos + ci < len(aux_row):
                        aux_row[col_pos + ci] = ch
            col_pos += nphor

    # Collapse aux rows to strings
    aux_strs = ["".join(row) for row in aux]

    # shist[0] = aux[0] with '|' at end (the numbers row above tallest bar)
    shist[0] = (aux_strs[0][: NUM_COL - 1] + "|")

    # Merge aux[k] into shist[k] for k=1..15 (only non-space chars replace spaces)
    for k in range(1, NUM_FIL - 1):
        row_list = list(shist[k])
        for ci, ch in enumerate(aux_strs[k]):
            if ch != ' ' and ci < len(row_list) and row_list[ci] == ' ':
                row_list[ci] = ch
        shist[k] = "".join(row_list)

    lines.append("Standardized time series histogram:")
    lines.append("")
    lines.append(base2)
    lines.append(base1)
    for row in shist:
        lines.append("        |" + row)
    lines.append(base1)
    lines.append(base2)
    lines.append("")

    lines.append(
        f"{atip1:16d} values outside (-1,+1):"
        f" {atip1 * 100.0 / n:5.2f} % (31.74 % expected)"
    )
    lines.append(
        f"{atip2:16d} values outside (-2,+2):"
        f" {atip2 * 100.0 / n:5.2f} % ( 4.56 % expected)"
    )
    lines.append("")


def _section_corr(lines, model, res, n_eff, freq):
    """ACF and PACF bar plots with Ljung-Box Q (PlotCor port)."""
    from .diagnostics import acf as _acf, pacf as _pacf

    std = float(np.std(res))
    if std < 1e-20:
        return

    n = n_eff

    if n < 3 * (freq + 1):
        lags = n - freq // 2
    elif freq == 1 and n > 200:
        lags = 45
    elif freq == 1:
        lags = 9
    else:
        lags = 3 * (freq + 1)

    lags = max(1, lags)
    nparma = _count_nparma(model)

    acf_vals  = _acf(res, lags=lags)
    pacf_vals = _pacf(res, lags=lags)

    _plot_corr(lines, acf_vals, lags, n, freq, nparma, is_acf=True)
    _plot_corr(lines, pacf_vals, lags, n, freq, nparma, is_acf=False)


def _plot_corr(lines, corr, lags, nobs, freq, nparma, is_acf):
    """ASCII ACF/PACF bar plot (PlotCor port)."""
    HORINC = 25.0
    band   = 2.0 / math.sqrt(nobs)

    GUIONS = "-------------+-------------------------+-------------------------+--------------"
    if is_acf:
        MARCAS = "            -1                         0                         1  L-B Q  DF"
        lines.append(f"Autocorrelation function (acf bands = +- {band:5.3f}):")
    else:
        MARCAS = "            -1                         0                         1"
        lines.append(f"Partial autocorrelation function (pacf bands = +- {band:5.3f}):")
    lines.append("")
    lines.append(MARCAS)
    lines.append(GUIONS)

    band_pos = _iround(band * HORINC)

    for i in range(1, lags + 1):
        r_i = corr[i - 1]

        # Period boundary marker
        if freq != 1 and i % freq == 0:
            prefix = f"{i:4d} {r_i:7.3f} +"
            end_sym = '+'
        else:
            prefix = f"{i:4d} {r_i:7.3f} |"
            end_sym = '*'

        # Build 53-char TmpStr (indices 0-52)
        bar = [' '] * 53
        bar[51] = end_sym if (freq != 1 and i % freq == 0) else '|'

        pos  = r_i * HORINC
        posi = abs(_iround(pos))
        sym  = '+' if (freq != 1 and i % freq == 0) else '*'
        if pos <= 0.0:
            for j in range(25 - posi, 26):
                bar[j] = sym
        else:
            for j in range(25, 26 + posi):
                bar[j] = sym
        bar[25] = '|'

        if bar[25 + band_pos] == ' ':
            bar[25 + band_pos] = ':'
        if bar[25 - band_pos] == ' ':
            bar[25 - band_pos] = ':'

        bar_s = "".join(bar)

        # Ljung-Box Q
        lb_str = ""
        if is_acf and (i - nparma) >= 1:
            show_lb = (
                (freq != 1 and i % freq == 0) or
                (freq != 1 and i == lags and i % freq != 0) or
                (freq == 1 and i == lags)
            )
            if show_lb:
                q_stat = _chi_test(corr, i, nobs)
                df     = i - nparma
                lb_str = f"{q_stat:6.2f} {df:3d}"

        lines.append(f"{prefix}{bar_s}{lb_str}")

    lines.append(GUIONS)
    lines.append(MARCAS)
    lines.append("")


# ── Utility functions ──────────────────────────────────────────────────────────

def _iround(x):
    """C-style round-half-away-from-zero."""
    return int(x + 0.5) if x >= 0.0 else -int(-x + 0.5)


def _poly_mul(p, q):
    """Multiply two polynomials [1, c1, c2, ...]."""
    result = [0.0] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            result[i + j] += pi * qj
    return result


def _seasonal_poly(coefs, sper):
    """
    Convert seasonal factor [Phi1, Phi2, ...] → polynomial array.
    Represents 1 - Phi1*B^s - Phi2*B^(2s) - ...
    """
    poly = [0.0] * (len(coefs) * sper + 1)
    poly[0] = 1.0
    for k, c in enumerate(coefs, 1):
        poly[k * sper] = -c
    return poly


def _count_nparma(model):
    """Count free ARMA parameters (AR1+AR2+MA1+MA2+AR1f+MA1f), excluding omega/delta/mu."""
    n = 0
    for factors, free_lists in [
        (model.ar,   model.ar_free),
        (model.ar_s, model.ar_s_free),
        (model.ma,   model.ma_free),
        (model.ma_s, model.ma_s_free),
    ]:
        for i, factor in enumerate(factors):
            fl = (free_lists[i] if free_lists is not None
                  else [True] * len(factor))
            n += sum(1 for f in fl if f)
    n += sum(1 for ff in model.ar_f if ff.free)
    n += sum(1 for ff in model.ma_f if ff.free)
    return n


def _chi_test(corr, lags, nobs):
    """Ljung-Box Q statistic (ChiTest port)."""
    q = 0.0
    for i in range(1, lags + 1):
        q += corr[i - 1] ** 2 / (nobs - i)
    return q * nobs * (nobs + 2)
