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
the DRVUS source archive (versions 1.0 … 1.2.03), held privately. Its `readme`
maps each module to its paper;
that mapping is the first row of every table below and it was written by the
authors, not reconstructed here.

---

## 2. The C core is Mauricio's, verbatim

The C embedded in this package (`csrc/internal/`) is the DRVUS source with a
licence header and two mechanical edits. Verify it yourself:

The archives these commands read are not distributed with the package —
Mauricio's DRVUS sources and the TASTE oracle are private — so they are named by
environment variable rather than by anyone's home directory:

```bash
export DRVUS_SRC=/path/to/drvus/src           # Mauricio's C, as published to us
export DRVUS_SOURCE=/path/to/drvus-source     # the versioned archive, 1.0 … 1.2.03
export TASTE_ORACLE=/path/to/Taste/oracle     # the oracle harness

for f in elfvarma usmelard drvmlest nlatools; do
  n=$(diff <(tr -d '\r' < "$DRVUS_SRC"/$f.c) \
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

### 2.1 The optimiser, and the one edit it carries

`qnewtopt.c` — `raxopt`, Mauricio (1995) JASA 90, 282-291 — sat outside that
invariant until 12 August 2026, which is how these things go: it was then the
one file BUG-0012 required editing. Its ancestor is `fue-1.13.1/src/qnewtopt.c`
rather than DRVUS, because the fue line already wrote its progress to `outputv`
instead of stdout.

Comparing **code only** (comments stripped), the delta is three things:

| change | why |
|---|---|
| `José` in UTF-8 | encoding |
| `printf("%4d F: …")` → `if (outputv) fprintf(outputv, …)` | the binding has no stdout to write to |
| `qn_last_termcode`, `qn_last_nit`, `qn_last_gnorm` | **BUG-0012**: record the verdict `report()` computes and then writes to `/dev/null` |

The third is the only addition, and it is a **recording**, not a decision: no
criterion, no announcement, no numerical behaviour. Two tests hold the line —
`test_the_optimizer_is_still_mauricios` (declared exceptions on code lines) and
`test_the_stopping_criteria_themselves_are_untouched`, which compares `umstop0()`
and `umstop()` — the two routines that decide when to stop — character for
character. Changing a stopping criterion is a **study**, not a bug fix; see
`TODO.md`.

The recording pays for itself immediately: on the 28 real-case `.out` files that
the C produced during actual work, the port now reproduces **termcode, iteration
count and gradient norm exactly** — including `Coint/R.4.out`, whose gradient
norm is 116330.0394. Reproducing an iteration count on twenty-eight real fits is
a much stronger statement than reproducing a likelihood
(`tests/test_optimizer_termcode.py`).

⚠ `drvarma` carries the same copy of this file and records none of it yet.

### 2.2 The archived `.out` files are 80-bit runs — and are reproducible

The reference runs preserved from the 2000s were produced by **32-bit x86**
binaries, whose FPU held intermediates in **80-bit** registers. Today's x86-64
builds use SSE2 and work in 64. On well-conditioned problems this is invisible:
8 of the 9 Box-Jenkins cases agree to 1e-10, and 28 of 28 real cases agree down
to the iteration count. On `a1` — a valley so flat that the AR sits at 0.99998 —
it decides the outcome, and `bugs/BUG-0012` is the whole measurement.

That case was run to the end, because the alternative was to leave the port
under suspicion. In summary: DRVUS's own source, rebuilt today at 64 bits,
reproduces **fue's** answer and not its own 2001 one; rebuilt with
`-m32 -O0` it reproduces the 2001 trace iteration by iteration; the preserved
2006 32-bit binary, run today, reproduces the archive exactly; and the preserved
**2001 Borland `.exe` of DRVUS 1.0**, under `wine`, agrees bit for bit with a
rebuild of the 1.0 source at 80 bits. Bisecting the versions with one compiler
and one set of flags puts the only functional change between 1.0 and 1.01 in two
lines of `drvus.c` — `steptol` tightened from `1.0e-5` to `macheps^(2/3)` — and
the estimation path itself (`qnewtop.c`, `usmelard.c`, `elfvarma.c`,
`drvmlest.c`) identical, **0 differing lines of code**. `statsmodels`, which owes
nothing to any of it, confirms the 2001 optimum is the real one, and fue reaches
it once mu is seeded sanely.

What matters for provenance is that the archive is **reproducible**, so a
reference can always be re-derived rather than trusted:

```bash
cp -r "$DRVUS_SOURCE"/1.2.01/drvus/src /tmp/drvus && cd /tmp/drvus
sed -i 's/\bround\s*(/drvus_round(/g' nlatools.c diagnose.c drvus.h
gcc -O0 -mfpmath=387 -o drvus drvus.c drvmlest.c elfvarma.c usmelard.c \
    qnewtopt.c nlatools.c diagnose.c -lm
./drvus a1 eml chk        # reproduces the 2001 a1.out iteration by iteration
```

`-mfpmath=387` restores the 80-bit intermediates and `-O0` keeps the optimiser
from spilling them to 64-bit memory — both are needed; with `-O2 -mfpmath=387`
the result returns to the SSE2 one. The only edit to the source is renaming
`round()`, which did not clash with `math.h` in 2001 and does now.

**Read this as a limit on what "agrees with the reference binary" can mean.**
`⚠ binary` in the table below is weaker than it looks whenever the likelihood is
ill-conditioned: it certifies agreement with one build on one architecture. The
checks that survive the change of arithmetic are the analytic ones, the oracle,
and the independent implementations (`statsmodels` on Series A, TASTE).

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
| `elfvarma.elf_scalar` (py)<br>`elfvarma.c: elf()` (C) | exact Gaussian log-likelihood of a VARMA(p,q), specialised to m=1 | **Mauricio (1997), Algorithm AS 311**, *Appl. Statist.* 46, 157-171; method of **Mauricio (1995)**, *JASA* 90, 282-291; innovations form of **Ansley (1979)** | `test_as311_published_identities.py` — **the paper's own equations**: (2) the exact log-likelihood rebuilt from `f1`/`f2` for several σ², (3) and (4) against Melard's published FORTRAN (S to **1e-14**, |ΛᵀΛ|^(1/n) to **4.4e-16**), σ̂²=S/n, and the equivalence of minimising S·|ΛᵀΛ|^(1/n) along a grid. Plus the ten steps of WP 9316 traced to the `[1]`…`[9]` blocks of `elfvarma.c` and the `(a)`…`(k)` markers of the port. `test_elfvarma.py::test_elf_scalar_ar1` against the closed-form exact AR(1) Gaussian likelihood (tol 1e-4); `::test_elf_scalar_ma1_iid` against the white-noise reduction at θ=0 (tol 1e-6)| **✅ published identities** |
| `elfvarma.flikam_scalar` (py)<br>`usmelard.c: flikam()` (C) | fast approximate likelihood, Kalman recursions with switch to quick recursions | **Melard (1984), Algorithm AS 197**, *Appl. Statist.* 33, 104-114 | `test_as197_published_fortran.py` — **against Melard's own FORTRAN**, transcribed from pages 110-114 and compiled: nine Box-Jenkins specifications, agreement to **5e-08** in log-likelihood. Plus `test_elfvarma.py::test_flikam_scalar_ar1` for the Python path | **✅ published source** |
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
cd "$TASTE_ORACLE" && ./battery.py --datos data      # the oracle harness (private)
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
| Boundary LR critical values for the DCD/MEG witness | **Davis & Dunsmuir (1996)**; **Davis, Chen & Dunsmuir (1995)**; and, for the per-frequency case, Guerrero (2026) |
| Unconditional-ML unit-root test | **Shin & Fuller (2001)** |
| Practice and the BJT procedure | **Treadway (1994)**; **Treadway, Guerrero & Mauricio (2009)** |

---

## 5. The rescaling convention, and where it comes from

`refactor = 100` runs through the whole suite. Its origin is not a modelling
choice and was, until now, written nowhere legible: it is advice from **Arthur B.
Treadway, May 2001** (Treadway, A. B. (2001), *DRVUS: manual de usuario*, unpublished user manual):

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
`RESCALING_ARCHITECTURE.md` in the `art-tseries` repository.

The same document states the `.inp` parser's contract, in 2001:

> *El programa no interpreta los comentarios escritos entre asteriscos, solamente
> lee los números que espera encontrar en la posición correcta.*

A **positional** parser that ignores comments and tolerates no deviation. That is
the specification of the format, and it is why a stray blank line breaks a file.

---

## 6. What is missing

Stated plainly, because a provenance document that only lists what works is an
advertisement.

- [x] ~~No test executes a published benchmark.~~ **Done, 12 Aug 2026** —
      `tests/test_published_benchmark_series_a.py`. The AS 311 and AS 197 papers
      were not to hand, but a better fixture was: **Box & Jenkins Series A**, the
      canonical benchmark of the field, which travels inside Mauricio's own DRVUS
      1.2.01 (`src/Box_y_Jenkings/`) with his `.out` files from ~2001. On the
      IMA(0,1,1), n=197:

      | | θ̂ | σ | logL |
      |---|---|---|---|
      | Box & Jenkins (book) | 0.70 | — | — |
      | DRVUS (Mauricio's C) | 0.699384 | 0.317382 | −53.508690 |
      | `fue` (2026) | 0.699384 | 0.317382 | −53.508690 |
      | **statsmodels** | **0.699384** | **0.317382** | **−53.508686** |

      `statsmodels` is what closes the link: different authors, a different
      algorithm, no shared code with DRVUS or FUE. **fue↔statsmodels agree to
      4.5e-09 on θ**, fue↔DRVUS to 2.7e-07, and both match the textbook 0.70.

      Series B–E are the obvious extension and are blocked today by
      `bugs/BUG-0011` — the DRVUS-era files do not load.
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
- [x] ~~**`flikam_scalar` is only loosely checked** (within 2.0 of the exact
      value).~~ **Done, 13 Aug 2026** — `tests/test_as197_published_fortran.py`.
      Melard's listing is printed in full in the paper (pp. 110-114) and is
      transcribed verbatim in `tests/fortran/as197.f`; compiled with
      `-fdefault-real-8` (the paper's own Precision note) it agrees with the
      engine to **5e-08** in log-likelihood on the nine Box-Jenkins
      specifications, and its `IFAULT = -m` switching point is asserted case by
      case. **This is the first check in the suite against a source that is
      neither ours nor Mauricio's**: everything else — the C, the Python port,
      the `.out` archive — descends from one implementation.

      Two things the exercise produced beyond the check itself. The `TOLER`
      contract of AS 197 — *"it should be negative if the exact likelihood is
      desired"* — is the reason `fue` signs `xitol` as it does
      (`fue_api.c:951-956`), until now justified only by `fue.c:1087`; there is
      a test for both branches. And the first transcription put label `170` one
      line early, which left every MA model right to 1e-8 and every AR model
      wrong by 2 to 7 in log-likelihood — the numbers caught what reading the
      listing had not.
- [ ] **`forecast.varphi` has no test of its own** — though `wti_forecast` (§3.6)
      exercises it against TASTE indirectly.
- [x] ~~**AS 311 is verified only against implementations that descend from
      it.**~~ **Done, 13 Aug 2026** — `tests/test_as311_published_identities.py`.
      AS 311 published no listing, so the check is by its published identities:
      equation (2) rebuilt from the engine's own `S` and `|ΛᵀΛ|`, and — the
      part that is genuinely external — equations (3) and (4) against **Melard's
      FORTRAN**, which shares neither author nor derivation. Six Box-Jenkins
      specifications: S to 1e-14, |ΛᵀΛ|^(1/n) to 4.4e-16.

      ⚠ Both must be evaluated with **negative `xitol`**. With the default
      positive value the Ξ sequence is truncated and S moves by 1e-5 relative —
      correct behaviour, different quantity, and it is exactly the `TOLER`
      contract of AS 197.

      What remains out of reach is the paper's own numerical example (WP 9316,
      Tables 4-5): the estimates are printed to two decimals, the series is not
      printed at all.

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
      (the simulation code accompanying Guerrero 2026) is continuous there. Any decision
      taken *at* the boundary should use the latter. Source: Guerrero (2026), appendix
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
  echo "$f: $(diff <(tr -d '\r' < "$DRVUS_SRC"/$f.c) \
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

---

## 9. References

The published algorithms, in the form a reader can obtain them:

* **Mauricio, J. A. (1997)**, "Algorithm AS 311: the exact likelihood function of
  a vector autoregressive moving average model", *Applied Statistics* **46**(1),
  157-171.
* **Mauricio, J. A. (1995)**, "Exact maximum likelihood estimation of stationary
  vector ARMA models", *Journal of the American Statistical Association* **90**,
  282-291.
* **Mauricio, J. A. (2002)**, "An algorithm for the exact likelihood of a
  stationary vector autoregressive-moving average model", *Journal of Time
  Series Analysis* **23**(4), 473-486.
* **Melard, G. (1984)**, "Algorithm AS 197: a fast algorithm for the exact
  likelihood of autoregressive-moving average models", *Applied Statistics*
  **33**(1), 104-114.
* **Ansley, C. F. (1979)**, "An algorithm for the exact likelihood of a mixed
  autoregressive-moving average process", *Biometrika* **66**, 59-65.
* **Dennis, J. E. and Schnabel, R. B. (1983)**, *Numerical methods for
  unconstrained optimization and nonlinear equations*, Prentice-Hall.
* **Box, G. E. P. and Jenkins, G. M. (1976)**, *Time series analysis:
  forecasting and control*, revised edition, Holden-Day.
* **Shin, D. W. and Fuller, W. A. (1998)**, "Unit root tests based on
  unconditional maximum likelihood estimation for the autoregressive moving
  average", *Journal of Time Series Analysis* **19**(5), 591-599.
* **Davis, R. A. and Dunsmuir, W. T. M. (1996)**, "Maximum likelihood estimation
  for MA(1) processes with a root on or near the unit circle", *Econometric
  Theory* **12**, 1-29.

Unpublished, and cited as such because they are not obtainable from a library:

* **Guerrero, D. E. (2026)**, *Hybrid seasonal models: critical values for
  testing deterministic versus stochastic seasonality frequency by frequency*,
  unpublished manuscript, Universidad Complutense de Madrid. — the source of the
  per-frequency DCD/MEG critical values and of the seasonal AR_f table.
* **Treadway, A. B. (2011)**, *FUE: manual de usuario*, unpublished user manual.
* **Treadway, A. B. (2001)**, *DRVUS: manual de usuario*, unpublished user
  manual. — the source of the positional `.inp` grammar and of the `refactor`
  advice.
* **Mauricio, J. A. (1993)**, *Exact maximum likelihood estimation of stationary
  vector ARMA models*, working paper 9316, Instituto Complutense de Análisis
  Económico. — the ten numbered steps the C implements, with equations
  (2.15)-(2.22).
