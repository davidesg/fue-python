"""Estimate RIPC.2 with pure Python and show diagnostic plots."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fue
from fue.plots import plot_model_diagnostics
from fue.cast_us import estimate_py
from fue.model import FitResult

# ── Load model ─────────────────────────────────────────────────────────────────

inp = os.path.join(
    os.path.dirname(__file__), "..",
    "tests/real_cases/PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.2.inp",
)
ts, model = fue.load(inp)

print(f"Series  : {ts.name}  ({ts.nobs} obs, freq={ts.freq})")
print(f"Model   : AR{[len(f) for f in model.ar]}, d={model.d}, D={model.D}, "
      f"{len(model.interventions)} interventions, boxlam={model.boxlam}")

# ── Estimate (pure Python raxopt) ──────────────────────────────────────────────

print("\nEstimating (pure Python)…")
raw = estimate_py(model)
model._result = FitResult(raw)
if not model._result.converged:
    raise RuntimeError(f"Estimation failed: ifault={model._result.ifault}")
r = model._result

print(f"\nConverged : {r.converged}  (ifault={r.ifault})")
print(f"npar      : {r.npar}")
print(f"sigma²    : {r.sigma2:.6f}")
print(f"log L     : {r.loglik:.4f}")
print(f"AIC       : {r.aic:.4f}")
print()
for i, (p, se) in enumerate(zip(r.params, r.std_errors)):
    print(f"  [{i+1:2d}]  {p:12.6f}  (se {se:.6f})")

# ── Output report ──────────────────────────────────────────────────────────────

here = os.path.dirname(os.path.abspath(__file__))

model.write_out(os.path.join(here, "RIPC.2.out"), inp_name="RIPC.2.inp")
print("\nSaved examples/RIPC.2.out")

# ── Plots ──────────────────────────────────────────────────────────────────────

fig_diag, fig_hist = plot_model_diagnostics(model)
fig_diag.savefig(os.path.join(here, "RIPC.2_diag.png"), dpi=120, bbox_inches="tight")
print("Saved examples/RIPC.2_diag.png")
fig_hist.savefig(os.path.join(here, "RIPC.2_hist.png"), dpi=120, bbox_inches="tight")
print("Saved examples/RIPC.2_hist.png")

