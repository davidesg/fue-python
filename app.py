"""fue diagnostic viewer — Streamlit app.

Run with:
    streamlit run app.py
"""

import io
import os
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

import fue
from fue.plots import plot_residuals_ts, plot_acf_pacf, plot_histogram


# ── Page config ────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="fue diagnostics",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("fue — residual diagnostics")


# ── Sidebar: file selection ────────────────────────────────────────────────────

st.sidebar.header("Model")

mode = st.sidebar.radio("Source", ["Upload .inp file", "Browse test cases"])

model = None
ts = None

if mode == "Upload .inp file":
    uploaded = st.sidebar.file_uploader("Select .inp file", type=["inp", "pre"])
    if uploaded:
        with tempfile.NamedTemporaryFile(suffix=".inp", delete=False) as tmp:
            tmp.write(uploaded.read())
            tmp_path = tmp.name
        try:
            ts, model = fue.load(tmp_path)
        except Exception as e:
            st.sidebar.error(f"Load error: {e}")
        finally:
            os.unlink(tmp_path)

else:
    test_root = os.path.join(os.path.dirname(__file__), "tests", "real_cases")
    inp_files = []
    for root, _, files in os.walk(test_root):
        for f in files:
            if f.endswith(".inp") and not f.startswith("forecast"):
                rel = os.path.relpath(os.path.join(root, f), test_root)
                inp_files.append(rel)
    inp_files.sort()

    selected = st.sidebar.selectbox("Select model", inp_files)
    if selected:
        path = os.path.join(test_root, selected)
        try:
            ts, model = fue.load(path)
        except Exception as e:
            st.sidebar.error(f"Load error: {e}")


# ── Sidebar: plot options ──────────────────────────────────────────────────────

st.sidebar.header("Options")
lags_input = st.sidebar.number_input("ACF/PACF lags (0 = auto)", min_value=0, max_value=100, value=0)


# ── Main: fit + plot ───────────────────────────────────────────────────────────

if model is None:
    st.info("Select or upload a .inp model file in the sidebar.")
    st.stop()

# Model info
freq_label = {1: "Annual", 4: "Quarterly", 12: "Monthly"}.get(ts.freq, f"freq={ts.freq}")
st.sidebar.markdown(
    f"**{ts.name}** · {freq_label} · {ts.nobs} obs  \n"
    f"Start: {ts.start[0]}/{ts.start[1]}"
)

# Fit
with st.spinner("Fitting model…"):
    try:
        model.fit()
    except RuntimeError as e:
        st.error(f"Estimation failed: {e}")
        st.stop()

r = model._result.residuals
npar = model._result.npar
freq = ts.freq

lags = int(lags_input) if lags_input > 0 else None
if lags is None:
    from fue.plots import _default_lags
    lags = _default_lags(len(r), freq)

# Fit summary strip
col1, col2, col3, col4, col5 = st.columns(5)
col1.metric("nobs", ts.nobs)
col2.metric("npar", npar)
col3.metric("σ²", f"{model._result.sigma2:.4f}")
col4.metric("log L", f"{model._result.loglik:.2f}")
col5.metric("AIC", f"{model._result.aic:.2f}")

st.divider()

# ── Row 1: residuals (wide) + ACF/PACF (narrow) ───────────────────────────────

col_ser, col_acf = st.columns([1.65, 1.0])

with col_ser:
    fig_ser, ax_ser = plt.subplots(figsize=(10, 4))
    plot_residuals_ts(r, model=model, title=ts.name, ax=ax_ser)
    fig_ser.tight_layout()
    st.pyplot(fig_ser)
    plt.close(fig_ser)

with col_acf:
    from fue.plots import _default_lags, _round_cmax
    import numpy as np
    from fue.diagnostics import acf as _acf, pacf as _pacf, ljung_box as _lb

    r_arr = np.asarray(r, dtype=float)
    rc = _acf(r_arr, lags=lags)
    pc = _pacf(r_arr, lags=lags)
    band = 2.0 / np.sqrt(len(r_arr))
    raw_max = max(float(np.abs(rc).max()), float(np.abs(pc).max()))
    cmax = _round_cmax(max(raw_max, band * 1.2))

    fig_acf, (ax_acf, ax_pacf) = plt.subplots(
        2, 1, figsize=(5, 6), layout="constrained"
    )
    plot_acf_pacf(r, npar=npar, freq=freq, lags=lags,
                  ax_acf=ax_acf, ax_pacf=ax_pacf)
    st.pyplot(fig_acf)
    plt.close(fig_acf)

# ── Row 2: histogram (centred, not too wide) ───────────────────────────────────

st.divider()
_, col_hist, _ = st.columns([1, 1.2, 1])
with col_hist:
    fig_hist, ax_hist = plt.subplots(figsize=(5, 5))
    plot_histogram(r, title=ts.name, ax=ax_hist)
    fig_hist.tight_layout()
    st.pyplot(fig_hist)
    plt.close(fig_hist)
