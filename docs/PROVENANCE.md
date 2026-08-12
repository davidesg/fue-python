# Provenance — which algorithm, from which paper, verified how

*This document exists so that a reader who does not trust the software can check
it. Every numerical routine in `fue` is traced to a published algorithm, and every
claim below is a command you can run.*

**Last audited: 12 August 2026.**

---

## 1. The chain

`fue` is not a reimplementation of published algorithms. It is a port of code
written by the authors of those algorithms, and the chain has four links:

```
  published algorithm          DRVUS              FUE               fue (Python)
  ───────────────────          ─────              ───               ────────────
  Melard (1984) AS 197    José Alberto        Treadway renames   port, 2026
  Mauricio (1995) JASA    Mauricio, 2000      DRVUS to FUE       (this package)
  Mauricio (1997) AS 311  SRC/drvus           ("free")
  Ansley (1979)
  Dennis & Schnabel (1983)
```

**DRVUS 1.01** (`SRC/drvus`) is Mauricio's own C, with earlier versions in
`SRC/drvus-source` (1.0 … 1.2.03). Its `readme` maps each module to its paper;
that mapping is the first row of every table below and it was written by the
authors, not reconstructed here.

---

## 2. The C core is Mauricio's, verbatim

The C embedded in this package (`csrc/internal/`) is the DRVUS source with a
licence header and two mechanical edits. Verify it yourself:

```bash
for f in elfvarma usmelard drvmlest nlatools; do
  n=$(diff <(tr -d '\r' < ~/Dropbox/SRC/drvus/src/$f.c) \
           <(tr -d '\r' < csrc/internal/$f.c) | awk '/^[<>]/ {n++} END {print n+0}')
  echo "$f.c: $n differing lines"
done
```

Measured 12 August 2026 (the DRVUS files are CRLF and Latin-1; these are LF and
UTF-8, hence the `tr`):

| file | differing lines | what they are |
|---|---|---|
| `elfvarma.c` | 23 | GPL header (14), `José` in UTF-8, the `#include`, **one functional change**, one closing comment |
| `usmelard.c` | 21 | GPL header, encoding, `#include`. **No functional change** |
| `drvmlest.c` | 22 | idem. **No functional change** |
| `nlatools.c` | 1129 | a separate cleanup — this module *was* rewritten (GSL, 764 lines against 1355) |

The single functional change in the whole likelihood core, `elfvarma.c:513`:

```c
-    eigenqr( a, n, wr, wi );
+    if ( n>1 ) {gsl_eigenqr( a, n, wr, wi );}
```

The Numerical Recipes eigenvalue routine replaced by GSL's, guarded for `n>1`.
**That is the entire delta.** The code implementing AS 311 and AS 197 is the code
Mauricio published.

> ✅ This is an invariant, not a claim: `tests/test_c_core_matches_drvus.py`
> fails if any line differs in a way not declared there. See §6.

---

## 3. The table

Status legend: **✅ analytic** = checked against a closed-form value derived
independently; **✅ oracle** = checked against TASTE, an independent
implementation (§3.6); **✅ cross** = C engine against pure-Python engine;
**⚠ binary** = checked against the reference binary only, so it verifies the port
and not the algorithm; **🔴 none** = no test addresses it.

### 3.1 The likelihood

| routine | what it computes | source | verified by | status |
|---|---|---|---|---|
| `elfvarma.elf_scalar` (py)<br>`elfvarma.c: elf()` (C) | exact Gaussian log-likelihood of a VARMA(p,q), specialised to m=1 | **Mauricio (1997), Algorithm AS 311**, *Appl. Statist.* 46, 157-171; method of **Mauricio (1995)**, *JASA* 90, 282-291; innovations form of **Ansley (1979)** | `test_elfvarma.py::test_elf_scalar_ar1` against the closed-form exact AR(1) Gaussian likelihood (tol 1e-4); `::test_elf_scalar_ma1_iid` against the white-noise reduction at θ=0 (tol 1e-6) | **✅ analytic** |
| `elfvarma.flikam_scalar` (py)<br>`usmelard.c: flikam()` (C) | fast approximate likelihood, Kalman recursions with switch to quick recursions | **Melard (1984), Algorithm AS 197**, *Appl. Statist.* 33, 104-114 | `test_elfvarma.py::test_flikam_scalar_ar1` — only that it is finite and within **2.0** of the exact value, because Melard's normalisation differs | **⚠ weak** |
| `elfvarma.chekma_scalar` | invertibility check on the MA polynomial | companion eigenvalues | `test_elfvarma.py` — five boundary cases, including θ=1.00004 vs 1.00005 | **✅ analytic** |

### 3.2 Estimation

| routine | what it does | source | verified by | status |
|---|---|---|---|---|
| `cast_us.cast_us_py` | builds the stationary series `w` from the ARMAX spec: Box-Cox, differencing, deterministic transfer functions | fue's own cast, mirroring the C chain `populate_globals → cast_us` | `test_cast_us.py` — `w` against the C engine's residual path; and §3.6 `cpi_deterministics` / `synd_estimate` against TASTE, which is what checks that interventions enter at the same level | **⚠ binary** + **✅ oracle** |
| `cast_us.estimate_py` | full exact-ML estimation loop | **Mauricio (1995) §3**, incl. the objective scaled so `F(x₀)=1`, eq. (3.5) | `test_cast_us.py`, `test_reliability{,2,3,4}.py` — loglik, σ², parameters, AIC/BIC, standard errors, residuals against the C engine and against `fue-1.13.1` reference output | **⚠ binary** + **✅ cross** |
| `qnewtopt.raxopt` (py)<br>`qnewtopt.c: raxopt()` (C) | factorised BFGS quasi-Newton with Cholesky factor of the Hessian maintained explicitly | **Dennis & Schnabel (1983)**, ch. 9: A9.4.1 (BFGS), A9.4.2 (rank-1 QR update), A6.3.1 (cubic backtracking line search), §7.2.1 (stopping) | `test_qnewtopt.py` — minimises Rosenbrock and quadratics; `cdgrad` against analytic gradients; `qrupdate`/`bfgsfac` against hand-computed 3×3 cases | **✅ analytic** |
| `nlatools.c` | numerical linear algebra, dynamic allocation | Numerical Recipes, rewritten around GSL | via the routines that use it | **⚠ indirect** |

### 3.3 Forecasting

| routine | what it does | source | verified by | status |
|---|---|---|---|---|
| `forecast.varphi` | combines the stationary AR polynomial with the non-stationary operator | mirrors `usfo.c` / `fuf.c` | *(no test by name)* | **🔴 none** |
| `forecast.forecast` | recursive level forecasts and ψ-weight variances | idem | `test_forecast.py::test_forecast_sfny2_vs_fuf` against `fuf` output; plus shape, convergence to the mean, and increasing σ. §3.6 `wti_forecast` against TASTE at six horizons, forecasts **and** standard errors, tol 1e-03 | **⚠ binary** + **✅ oracle** |
| Box-Cox back-transformation | asymmetric intervals on the original scale | **Box & Cox (1964)** | `test_forecast.py::test_forecast_level_positive_boxcox0` | **⚠ partial** |

### 3.4 Diagnostics

| routine | source | verified by | status |
|---|---|---|---|
| `diagnostics.acf`, `pacf` | standard definitions | `test_api.py` | **⚠ smoke** |
| `diagnostics.jarque_bera`, `ljung_box` | Jarque & Bera (1980); Ljung & Box (1978) | `test_api.py` | **⚠ smoke** |

### 3.6 Against an independent implementation — TASTE

The strongest check available, and the one that does not run inside this package.
**TASTE** (Pascal, 1987-2001) solves the same problem by another route; the
oracle battery in `SRC/atws/Taste/oracle` drives it and compares.

```sh
cd ~/Dropbox/SRC/atws/Taste/oracle && ./battery.py --datos data
```

Measured 12 August 2026 — **20 of 20 cases pass**, of which four are `fue`'s:

| case | what it exercises | values | tolerance |
|---|---|---|---|
| `cpi_deterministics` | 11 deterministics (5 cos, 5 sin, alternator) on the Spanish CPI — that interventions enter at the same level of the series in both programs | μ, φ, three ω | 5e-03 |
| `synd_estimate` | the simplest case in which an intervention enters: one step, d=0 — isolates whether the step is built in the same place | ω, φ | 1e-02 |
| `wti_forecast` | univariate forecast at 6 horizons, AR(1) on ∇log — the forecast recursion and the band | 3 forecasts, 3 s.e. | 1e-03 |
| `fue_es_cpi_airline` | **added 12 Aug 2026** — the first case with SEASONAL differencing. Validates the transformation, not the estimate: see below | θ, Θ | 1e-03 |

### TASTE is a different ESTIMATOR, and that limits what the oracle can say

**TASTE estimates by nonlinear least squares with backforecasting** — Box &
Jenkins (1976) §7.1.4 — and `fue` by exact unconditional maximum likelihood. On
a series whose seasonal MA sits near the non-invertibility boundary the two
criteria pull apart. Measured 12 August 2026 on the airline model
(0,1,1)(0,1,1)₁₂ over the Spanish CPI:

| | θ | Θ |
|---|---|---|
| `fue`, exact ML | −0.421156 | **+0.814706** |
| TASTE, BJ NLS | −0.426710 | **+0.906060** |

**Neither is wrong.** Each sits at its own optimum: on fue's exact likelihood
TASTE's point is 1.494 worse; on the backforecasting sum of squares fue's point
is 14.418 against TASTE's 14.186.

The explanation was established by implementing the Box-Jenkins criterion on the
series `fue` produces after ∇∇₁₂ — it optimises at θ=−0.426706, Θ=+0.906150,
**TASTE to 4e-06 and 9e-05**. Two candidates were ruled out first: it is not
convergence (TASTE gives the same at 60 and 300 iterations) and it is not the
plain conditional sum of squares, which optimises at Θ=+0.7848, *below* fue and
in the opposite direction. `tests/test_taste_nls_criterion.py`.

**What this means for the oracle.** On models where the two criteria agree
closely — short samples, parameters away from the boundary — the comparison
validates the estimate. On models near the boundary, which is exactly where this
suite does its most interesting work, it validates **the transformation, not the
estimate**: that both programs are looking at the same series through the same
operator. Reading such a case as a defect was the first thing that happened, and
it is recorded here so it is not read that way again.

**How independent is it, exactly?** TASTE was *written by José Alberto Mauricio*
too, directed by Treadway and Serrano (`Taste/PROCEDENCIA.md`). So this is an
independent **implementation**, not an independent **author**: it shares no code
with the DRVUS/FUE line and was written in a different language and era, but a
shared author can carry a shared misunderstanding. It rules out implementation
bugs, not conceptual ones. Stated plainly because overselling it would defeat the
purpose of this document.

**And the per-frequency half cannot be validated here at all.** TASTE's model
specification carries `{λ, d, ds, sp}` and nothing else (`USMOD.PAS`): regular and
seasonal differences, no per-frequency factors. There is no way to express a
hybrid MEG in it, so the half of the thesis that makes `fue` unique is, by
construction, the half no independent implementation can check.

What *can* be chained to the oracle is the reduction: since

    ∇₁₂ = (1−B)(1+B)·∏_{f=1}^{5}(1−2cos ω_f·B + B²)

a model with **every** `ifadf` factor active is identically the SARIMA with D=1 —
verified in `fue` to **1.1e-13** in log-likelihood and 1.7e-13 in φ. So if the
per-frequency machinery ever stopped reproducing the seasonal difference, the
`fue_es_cpi_airline` case would move. That does not validate a *mixed* MEG, and
nothing external can.

**Coverage remains thin for `fue`:** four cases against sixteen for `drtran`.
Still unexercised on the univariate side: the Box-Cox family, impulse and ramp
interventions, and seasonal ARMA operators beyond the airline.

### 3.5 Not implemented, and deliberately

| routine | source | status |
|---|---|---|
| `multshea.c` | **Shea (1989), Algorithm AS 242**, *Appl. Statist.* 38, 161-184 | Declared *«no implementado»* in the DRVUS 1.0 `README.txt` (2002) and still out of scope in `drvarma` today. **Twenty-six years of the same decision, taken twice.** |

---

## 4. What the model itself rests on

Not algorithms, but the specification the software implements:

| element | source |
|---|---|
| ARIMA(p,d,q)(P,D,Q)ₛ, identification and diagnosis | **Box & Jenkins (1970, 1976); Box, Jenkins & Reinsel** |
| Box-Cox transformation | **Box & Cox (1964)**, *JRSS B* 26, 211-252 |
| Hybrid / frequency-by-frequency seasonality (MEG) | **Abraham & Box (1978)**; **Gallego (1995, 1996)** |
| Boundary LR critical values for the DCD/MEG witness | **Davis & Dunsmuir (1996)**; **Davis, Chen & Dunsmuir (1995)**; and, for the per-frequency case, the SF_MEG paper of this project |
| Unconditional-ML unit-root test | **Shin & Fuller (2001)** |
| Practice and the BJT procedure | **Treadway (1994)**; **Treadway, Guerrero & Mauricio (2009)** |

---

## 5. The rescaling convention, and where it comes from

`refactor = 100` runs through the whole suite. Its origin is not a modelling
choice and was, until now, written nowhere legible: it is advice from **Arthur B.
Treadway, May 2001**, in `drvus-source/1.0/Drvus/ABTreadway-Drvus.doc`:

> *Es deseable que la norma del gradiente sea cero hasta toda la precisión que
> ofrece el programa. Cuando esto no ocurre […] es muchas veces útil escalar los
> datos para que tengan más precisión de entrada. Esto se hace multiplicando
> todos los valores de la variable lnY […] por, p.e., 100 antes de introducirlos
> como entrada al programa. Por supuesto, también se multiplicará cada parámetro
> de intervención por el mismo factor, y la salida presentará una sigma
> multiplicado por el mismo factor.*

So the ×100 is **a numerical-precision remedy aimed at the gradient norm**, and
the compensating rules — every intervention parameter scales, and σ comes out
scaled — are part of the same advice. That is why the two defects it has produced
(`bugs/BUG-0001`, μ collapsing under rescaling; `bugs/BUG-0007`, μ read 100× off
scale) were both about *scale in the report*, not about estimation. See
`art-python/docs/RESCALING_ARCHITECTURE.md`.

The same document states the `.inp` parser's contract, in 2001:

> *El programa no interpreta los comentarios escritos entre asteriscos, solamente
> lee los números que espera encontrar en la posición correcta.*

A **positional** parser that ignores comments and tolerates no deviation. That is
the specification of the format, and it is why a stray blank line breaks a file.

---

## 6. What is missing

Stated plainly, because a provenance document that only lists what works is an
advertisement.

- [ ] **No test executes the published test cases of AS 311 or AS 197.** Both
      papers ship data and expected results. Today `test_estimation.py` uses
      *"hard-coded reference values obtained from the reference binary"*, which
      verifies **port against binary**, not **binary against paper**. Given the C
      *is* Mauricio's code this is close to a formality — but it is the formality
      that turns "trust" into "check".
- [x] ~~The `diff` of §2 is not a test.~~ **Done, 12 Aug 2026** —
      `tests/test_c_core_matches_drvus.py`, 9 tests. It does **not** count
      differing lines (a count would pass if twenty-three benign lines were
      replaced by twenty-three malicious ones): every differing line must match a
      **declared** exception, and anything else fails with the line printed. The
      one functional change is written out in full as its own exception, so
      altering it further breaks the test instead of sliding through.
      Mutation-checked in both directions: changing `LOG2PI` by one digit fails
      the classification test, and reverting the GSL substitution fails the test
      that guards it.
- [ ] **`flikam_scalar` is only loosely checked** (within 2.0 of the exact value).
      It runs inside the BFGS inner loop, so it drives every iteration.
- [ ] **`forecast.varphi` has no test of its own** — though `wti_forecast` (§3.6)
      exercises it against TASTE indirectly.
- [ ] **The oracle covers four `fue` cases out of twenty.** Still missing on the
      univariate side: Box-Cox, impulse and ramp interventions, seasonal ARMA
      beyond the airline. The per-frequency factors cannot be added — TASTE
      cannot express them (§3.6).
- [ ] **Diagnostics are smoke-tested only.** `acf`, `pacf`, `jarque_bera` and
      `ljung_box` have closed forms and reference implementations; checking them
      against `scipy`/`statsmodels` is cheap.

---

## 7. Known limits of the estimator

Documenting the limits is part of allowing verification; hiding them prevents it.

- **The likelihood profile jumps at the invertibility boundary.** The production
      estimator shows an erratic upward jump exactly at `r=1`, which inflates the
      apparent pile-up and distorts the tail. The exact banded-Cholesky likelihood
      (`art-python/research/sf_meg/dcd_mc.py`) is continuous there. Any decision
      taken *at* the boundary should use the latter. Source: SF_MEG, appendix
      *"Note on the exact likelihood near the boundary"*.
- **A spurious optimum that depends on the build.** Same source, same data, same
      starting values: the Windows wheel settles on `μ̂=−0.144` (AIC −2511) and the
      Linux wheel on `μ̂=0.0021` (AIC −2613), both reporting `converged=True`,
      `ifault=0`. `bugs/BUG-0005`.
- **`converged=True` is not a certificate.** `ifault` reports model adequacy, not
      convergence; the optimiser's own termination code is a different thing. See
      `drvarma`'s note on `termcode`.
- **The optimiser is published work and is not modified.** `raxopt`/`qnewtopt`
      keep Mauricio's stopping criteria and announcements. A defect there is
      answered with a documented diagnosis, not a patch — see `drvarma/TODO.md`
      §*ESTUDIO — el criterio de parada del optimizador*.

---

## 8. How to re-audit this document

```bash
# 1 — the C core is still Mauricio's
for f in elfvarma usmelard drvmlest; do
  echo "$f: $(diff <(tr -d '\r' < ~/Dropbox/SRC/drvus/src/$f.c) \
                   <(tr -d '\r' < csrc/internal/$f.c) \
              | awk '/^[<>]/ {n++} END {print n+0}') lines"
done

# 2 — the citations in the code
grep -rhoE "Mauricio \(19[0-9]{2}|Melard \(1984|Ansley \(19[0-9]{2}|AS 311|AS 197" src/fue/*.py | sort | uniq -c
# measured 12-Aug-2026: Mauricio (1995) ×8, Mauricio (1997) ×6, AS 311 ×4,
#                       AS 197 ×3, Ansley (1979) ×1, Ansley (1982) ×2

# 3 — the tests behind each row
python -m pytest tests/test_elfvarma.py tests/test_qnewtopt.py tests/test_cast_us.py -v
```

If any of the three disagrees with the tables above, this document is stale and
the discrepancy is the finding.
