# Why the wheel, and what it costs to do without it

*The question anyone installing this package eventually asks: there are two
engines — what do I lose by running the one that needs no compiler? Measured,
not estimated, over the 23-case benchmark that ships with the tests.*

The short answer: **you lose speed — one to two orders of magnitude — and you
lose nothing else.** The two engines run the same algorithm with the same
optimiser and reach the same optimum to the eighth decimal.

---

## 1. Efficacy first, because it is the one that matters

If the answers differed, speed would be irrelevant. They do not:

| | |
|---|---|
| largest \|Δ log-likelihood\| over 23 real models, Python vs C | **0.0002** |
| in all but one case | **< 0.00005** |

Same likelihood (`elfvarma.py` translates `elfvarma.c` step by step, and both
are checked against the published algorithms — `docs/PROVENANCE.md`), same
optimiser (`raxopt` in both, `docs/PORT.md` §3.1), same stopping rules. The
one case that moves, `IPC-T/Coint/R.4`, is the fit that never converged in the
first place — termcode 2 with a gradient of 1e5 — so a difference in the eighth
decimal there is the flat valley, not the engine.

## 2. Efficiency, over the 23-case benchmark

Measured 13 August 2026, Linux x86-64, CPython 3.12, gcc -O2. Reproduce with

```bash
python -m pytest tests/test_performance.py::test_summary -s
```

| | C (the wheel) | Python + raxopt | Python + L-BFGS-B |
|---|---|---|---|
| slowdown vs C, median | — | **×90** | ×282 |
| slowdown vs C, range | — | ×29 … ×384 | ×29 … ×1178 |

In wall-clock terms, on the models people actually fit:

| case | freq | n | npar | C | Python + raxopt |
|---|---|---|---|---|---|
| PCE/SF/R.1 | 4 | 68 | 1 | 0.1 ms | 4.1 ms |
| IPC-T/R.1 | 4 | 68 | 4 | 0.3 ms | 21.2 ms |
| IPC-T/Coint/R.4 | 4 | 68 | 5 | 1.6 ms | 134 ms |
| RIPC.4 | 12 | 78 | 15 | 8.8 ms | 1.8 s |
| RIPC.1 | 12 | 72 | 14 | 38 ms | 7.0 s |
| RIPC.3 | 12 | 78 | 17 | 14 ms | 2.7 s |

So: a single quarterly model is imperceptible either way. A monthly model with
fifteen parameters is 9 ms against 2 seconds — and that is the difference
between an interactive session and a wait. Multiply by the model ladder of a
real analysis, or by a Monte Carlo, and the wheel stops being a convenience.

**Where the time goes**, profiled on SFNY.2: `cast_us_py` 63%, `flikam_scalar`
19%, the rest below 4% each. It is the casting, not the likelihood, and not
scipy: the Python loops that rebuild the model structure at every function
evaluation are what the C does in compiled form.

## 3. The other study: raxopt against L-BFGS-B

Both are available in the Python engine (`estimate_py` uses `raxopt`;
`optimizer="lbfgsb"` switches). The comparison was worth making because the
intuition — a modern library optimiser must beat a 1995 translation — is wrong
on both counts here:

* **Speed.** raxopt is faster in 21 of the 23 cases, by a median factor of
  **2.4×** and up to 7.3×. L-BFGS-B evaluates the objective more times, and
  every evaluation carries the Python casting overhead of §2.
* **The optimum.** Δ log-likelihood is ≈ 0 for raxopt against the C on every
  case. L-BFGS-B is not uniformly worse — in June 2026 it escaped a local
  optimum on `IPC-T/Coint/R.4` where raxopt stayed (251.68 against 211.21) —
  which is exactly why it is kept as an option rather than removed.

That asymmetry is the honest summary: **raxopt for estimation, L-BFGS-B as a
second opinion when a fit looks trapped.** `bugs/BUG-0005` is the open defect
about spurious optima on multimodal seasonal-AR likelihoods, and it is the
situation where trying the other optimiser is the right reflex.

## 4. What the wheel actually contains

`pip install fue` fetches a `manylinux` wheel with the C engine compiled in and
**GSL linked statically**, so there is nothing else to install. Building from
source needs a C compiler and GSL headers; if neither is available:

```bash
FUE_SKIP_C=1 pip install fue      # skip the extension: pure Python, permanently
```

Everything except the estimator is pure Python in both cases — the `.inp`
parser, the reports, the forecasts, the diagnostics, the plots. The extension is
the estimator and nothing else.

## 5. What was considered and not done

Three hybrid architectures were costed before settling on "ship the C, keep the
translation":

| option | projected speed | dependency | verdict |
|---|---|---|---|
| Python optimiser + C objective function (`fue_objcfunc`) | ≈ C | still GSL | the interesting one, if the optimiser ever needs extending |
| Python optimiser + Python casting + C `flikam` only | ×2.5 over pure Python | still GSL | rejected: `flikam` is 19% of the time, the casting is 63% |
| vectorised pure Python (`np.convolve`, `lfilter`) | ×20-50 vs C | none | not done; it is a rewrite of the loops, and the algorithms are published |

The first remains the candidate if the day comes when the optimiser has to grow
constraints or alternative criteria — the C would keep the inner loop and the
search would move to Python without the ×90.

---

*The full measurement log, including every migration stage since June 2026, is
in `PERFORMANCE.md` at the repository root; this document is the part a user
needs in order to decide.*
