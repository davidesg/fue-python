# The cast

*A study of `cast_us()`: what it is, what it guarantees, what it costs, and where
it ends — from the free parameter vector to the driver that hands the result to
the likelihood.*

---

## 0. Attribution

The cast is part of the code signed by **José Alberto Mauricio** in DRVUS; the
credit is his. `csrc/fue_api.c:371` records the extraction point — «cast_us:
translate Tm+Ts+DataMat into Tvarma for the engine. Extracted from `fue.c:3645`».
The methodology around it —the frequency-by-frequency reading of seasonality, the
MEG's ordering, the practice the specification language serves— is Treadway's,
and it is plausible that the design of what the cast must *express* is his too.
**Nothing in the record establishes that**, and this document does not assert it.
See `PROVENANCE.md`.

---

## 1. What the cast is

**The cast is a compiler.** It is not a parameter unpacker, and reading it as one
is what makes it look like plumbing.

It takes a *declarative specification* — operators in factored form, a free/fixed
flag per coefficient, factors pinned to a frequency, the annual difference as a
list of factors, transfer functions ω(B)/δ(B) on inputs — and lowers it to the
*flat canonical form* an ARMA likelihood can evaluate: two expanded polynomials
φ(B), θ(B), a mean, and a stationary series `w`.

That separation between **specification language** and **evaluation kernel** is
what the rest of the suite is built on:

* the transfer function exists because the cast can subtract a filtered input
  before the operator is applied;
* the MEG can state a hypothesis *at one frequency* because the cast can pin an
  AR/MA factor to that frequency;
* `ifadf` — the annual difference as a list of factors — is expressible because
  the non-stationary operator is a list, not the integer pair `(d, D)`;
* `drtran` can embed a fue model inside a VARMA because the cast delivers exactly
  what a VARMA row needs.

The practical consequence is that **fue's model space is larger than the standard
`(p,d,q)(P,D,Q)s` API's** — not because of numerical power, but because there is
somewhere to *say* those models. The cast is that somewhere.

## 2. Where the cast sits: it is an argument, not a step

The architectural keystone is in `csrc/internal/fue.h:167`:

```c
void est( void (*cast)( real *, struct Tvarma *, int *, int, int ),
          int npar, real *par, real *dev, real **cov, int maxits, int nrits,
          real grtol, real sptol, real xitol, int chkma,
          real **a, real *sigma2, real *logelf, int *ifault );
```

**The cast is a function pointer passed to the estimation driver.** `est` does not
know what a fue model is; it knows that something will turn a parameter vector
into a `Tvarma`. The optimiser proposes vectors, the cast translates, `elf`
evaluates, and none of the three knows the other two.

This is why `drtran` can substitute its own casts (`cast_diagonal`,
`cast_embedded`) and reuse everything else, and it is the single most important
property of the design: **the cast is replaceable by construction.** Any proposal
to give a consumer its own cast is working *with* the architecture, not against
it.

## 3. Lifecycle: `firstx` and `lastx`

```c
static void cast_us(real *x, struct Tvarma *armax,
                    int *ifaultx, int firstx, int lastx)
```

The last two arguments are **not** a partial-recomputation hint — a natural and
wrong reading. They are the session boundary:

* `firstx` → block [3] allocates `mu`, `phi`, `theta`, `qq`, `w`, `a`
  (`fue_api.c:434`);
* `lastx == 1` → block [6] frees them all (`fue_api.c:513`);
* every call in between reuses the buffers and rewrites them in place.

So the C cast **already has an open/eval/close shape**, expressed as two flags
rather than as a handle. Anyone considering exposing `cast_us` through the CFFI
should know this: the session API is not invented from nothing, it is `firstx`/
`lastx` promoted to the surface. What is *not* solved by those flags is the state
that lives outside `armax` — see §9.

## 4. The parameter vector: the suite's real interface contract

`x` carries **only the free parameters**, in the order `count_npar_build_par()`
lays them down. From `cast_us_py`'s own docstring, and matching `cast_us()`
step [1]:

| # | group |
|---|---|
| 1 | ω free, per intervention |
| 2 | δ free, per intervention |
| 3 | regular AR free coefficients, factor by factor |
| 4 | seasonal AR free coefficients |
| 5 | regular MA free coefficients (with the MA(1) \|θ\|>1 flip applied) |
| 6 | seasonal MA free coefficients |
| 7 | fixed-frequency AR: the free `coef2`, one per factor |
| 8 | fixed-frequency MA: the free `coef2` (with the \|θ₂\|>1 flip) |
| 9 | μ, if `estimate_mu` |

**This order is a public contract even though it is written nowhere public.**
`drtran` slices the joint vector by `sc.npar` and hands each slice to the cast; if
the order changed, nothing would raise — every consumer would silently compute a
different model. It belongs in this document precisely because the code cannot
enforce it.

Two properties that follow and that consumers rely on:

* **`p` and `q` depend only on the ORDERS, never on the values.** The expanded
  polynomial lengths are constant throughout an optimisation, which is what lets
  `drtran` build the VARMA structure once.
* **A coefficient pinned at zero is still a coefficient**, occupying its slot in
  the factor while contributing nothing to `x`. The distinction between «no
  factor» and «a factor fixed at zero» is real, and it is the difference between
  the path that crashes and the path that does not (BUG-0013).

## 5. The translation, block by block

Following the pure-Python mirror `fue.cast_us.cast_us_py`, which is a literal
mirror of the C:

**[1] Unpack.** `x` → the current values of ω, δ, the AR/MA factors (regular,
seasonal, fixed-frequency) and μ, respecting the free/fixed flags: a fixed
coefficient keeps its declared value and never reads from `x`.

**[2] MA(1) invertibility flip.** A one-coefficient MA factor with |θ₁| > 1 is
replaced by 1/θ₁. This is a **design decision with consequences**: the
optimiser is left *unconstrained* — it may wander outside the invertible region —
and the cast maps the point back to the observationally equivalent invertible one.
The likelihood surface therefore has a fold, which is the mechanism behind
`BUG-0005` (spurious optimum). The alternative — constraining the optimiser —
would have coupled `est` to the model class it is deliberately ignorant of.

**[3] Fixed-frequency AR → a regular 2-lag factor.** The free parameter is
`c2 < 0`, and

    c1 = 2·cos(2π·f/s)·√(−c2)

so the factor `(1 − c1 B − c2 B²)` has its roots at **exactly** the frequency
`f/s`, whatever the modulus. *This is the mechanism that lets a hypothesis be
stated at one frequency*: the frequency is not estimated, it is pinned, and only
the modulus is free. `c2 ≥ 0` is not representable and returns `ifault = 1`.

**[4] Fixed-frequency MA.** The same, plus the |θ₂| > 1 flip (`c2 → 1/c2`) for the
same invertibility reason as [2].

**[5] Expansion.** The factored operators — regular factors, fixed-frequency
factors appended as ordinary 2-lag factors, and the seasonal ones lifted to lag
`s` — are convolved into the flat `phi` and `theta`. This is `_unscramble`, and it
is where «factored» stops and «canonical» begins.

**[6] The deterministic path.** For each intervention, the impulse response
ν(B) = ω(B)/δ(B) is built (`calcnu`) and convolved with the input, and the result
is **subtracted** from the Box-Cox'd data:

```python
lags = len(ω) - 1 if δ is empty else 40      # pure FIR needs no truncation
nu   = calcnu_py(ω_j, δ_j, lags)
for t in 1..nobs:
    z[t] -= Σ_k nu[k] · ind[t-k]
```

Note the truncation rule: a **pure FIR** transfer (δ empty) is exact at `len(ω)−1`
lags; a rational one is truncated at **40**. That 40 is a modelling constant with
no error bound attached to it in the code, and it deserves its own study.

**[7] The non-stationary operator.** `w` is produced by the `rnsop` recursion,

    w[t] = z[j] − Σ_i rnsop[i]·z[j−1−i],     j = ornsop + 1 + t

where `rnsop` is the expanded coefficient list of ∇^d ∇_s^D and the individual
`ifadf` factors. **The operator is a list of coefficients, never the pair
`(d, D)`** — which is exactly why `∇∇₄` written as `d=2, ifadf=[0,1,1]` and
written any other way yield the same polynomial. Consumers that rebuild it from
`d + D·s` are wrong by construction; that is BUG-9 and BUG-10 in `drtran`.

The sample shrinks here: `n_eff = nobs − ornsop`. Every «common sample» question
downstream is a question about `ornsop`, not about `d` and `D`.

## 6. What is constant across an optimisation, and what is not

This table is the one an optimiser-facing consumer needs, and it did not exist
before this document. «Constant» means: identical in every evaluation of a single
`fit()`.

| piece | constant? | why it matters |
|---|---|---|
| `p`, `q` | **yes** | structure can be built once |
| the intervention design columns `ind_data[j]` | **yes** | they are data |
| `rnsop`, `ornsop`, `n_eff` | **yes** | the differencing does not depend on `x` |
| `data0` (Box-Cox'd series) | **yes** | λ is decided before estimation |
| ν(B) for a **pure FIR** intervention | **no**, but it *is* ω itself | the deterministic path is then **linear in ω**: ξ = X·ω with X constant |
| ν(B) for a **rational** intervention | no | δ changes the whole response |
| `phi`, `theta` | no | the point of the optimisation |
| `w` | no | it carries the deterministic subtraction |

The row that matters most in practice: for pure-FIR deterministics — which is
what the full harmonic set of a deterministic-seasonality model is — **block [6]
recomputes, on every evaluation, a matrix-vector product whose matrix never
changes**. Measured on IPC_ES (11 interventions, all pure FIR, one ω each, 288
observations): the block costs ~1.9 ms per call in the Python mirror, against 81 µs
for the exact VARMA likelihood it feeds. See `drtran-python/docs/STUDY_efficiency_vs_c.md`.

This is *not* a defect of the cast. It is a consequence of the C's per-call
contract (§3): the cast is called with a fresh `x` and must produce the full
answer. A consumer that holds state across evaluations can do better, and that is
an argument about **contracts**, not about correctness.

## 7. Guarantees and failure modes

What the cast guarantees to whatever evaluates the likelihood:

* `phi` has exactly `p` entries and `theta` exactly `q`, both in plain expanded
  form, both consistent with the declared orders;
* `w` has `n_eff = nobs − ornsop` entries, with the deterministic path already
  subtracted and the non-stationary operator already applied;
* μ is the mean **of the differenced series**, not an intercept of the level —
  a distinction that has caused confusion outside this file and is worth
  restating: since `ifadf` is differencing, μ scales with `A_f(1)`;
* every returned MA factor is invertible, by flips [2] and [4].

What it does not guarantee, and the failure modes:

* `ifault = 1` — a fixed-frequency coefficient with `coef2 ≥ 0`, i.e. a request
  that has no representation as a real 2-lag factor at that frequency. It is a
  *specification* failure surfaced at translation time;
* `n_eff ≤ 0` returns an empty `w` with `ifault = 0` — the caller must check;
* **no bounds are enforced on the AR side.** Stationarity is not the cast's
  business; `drtran` checks it separately (`ar_is_stationary`), which is itself a
  sign that the guarantee is missing at this layer;
* with `p = q = 0` the C writes out of range (`BUG-0013`). The empty case was not
  defended. The Python mirror handles it, which is why the current fix diverts
  rather than repairs.

## 8. The driver: from the cast to the likelihood

```
    est( cast_us, npar, par, ... )        fue_api.c:958
      └─ optimiser (raxopt / qnewtopt) proposes x
           └─ cast_us(x, &armax, &ifault, first, last)     ← this document
                └─ elf( m, n, p, q, mu, phi, theta, qq, w, ... )
```

`elf` receives `w` **as it is** and computes the *exact unconditional* likelihood,
initialisation included: nothing is truncated outside it, which is why the
embedded cast is preferred for forecasting in `drtran`. The likelihood is
**concentrated** — `elf` returns the pieces from which σ² comes out separately, so
a consumer that reads `qq` as if it were Σ gets the right shape and the wrong
magnitude.

The optimiser side is `raxopt`/`qnewtopt`, and it is **published work that is not
to be modified**: its stopping tests are the relative-gradient test (termcode 1)
and the step test (termcode 2, the stuck detector), with `cmacheps()` adapting the
tolerances to the arithmetic. Anything to be said about `typx ≡ 1` being
hardcoded belongs in a study, not in a patch.

The boundary of this document is exactly here: **the cast ends where `elf`
begins.** What `elf` does with `w` is the subject of AS 197/AS 311 and of
`PROVENANCE.md`.

## 9. Global state, and why it constrains everything downstream

```c
struct Tusmodel Tm;
struct Tseries  Ts;
double        **DataMat;
```

`cast_us()` reads the model, the series and the data from **module-level
globals** (`fue_api.c:37`), matching `fue.c`'s layout — deliberately, so the
extraction stays comparable with the original. `fue_estimate()` populates them on
entry and tears them down on exit, so the public API is one call in, one call out
and the globals are never observable.

The consequences are not academic:

* the cast is **not reentrant**, so it cannot be called concurrently — which
  matters the moment a live MCP server is in the picture;
* **two casts cannot be alive at once**, which is exactly what a two-series
  transfer model needs;
* therefore exposing `cast_us` through the CFFI requires *de-globalising* it —
  and that breaks the deliberate correspondence with `fue.c:3645` that makes the
  extraction auditable.

That trade — auditability against reentrancy — is the central open design
question about the cast, and it should be decided explicitly rather than as a
side effect of an optimisation.

## 10. The Python mirror

`fue.cast_us.cast_us_py` is a **literal** mirror: same seven blocks, same loop
structure, same 1-based indexing where the C has it. That is the right choice and
it should stay that way, for a reason that has nothing to do with style: the
mirror reproducing the C binary's log-likelihood to 1e-9 is what validates the
port. It is an oracle, and an oracle that has been rewritten for speed is no
longer independent evidence.

The cost of being a faithful mirror is the ~90× median slowdown documented in
`PERFORMANCE.md`, concentrated in blocks [6] and [7]. **A consumer that needs
speed should get it somewhere other than by rewriting the mirror.**

## 11. Open questions

1. **The 40-lag truncation** of a rational transfer's impulse response (block [6])
   carries no error bound in the code. What is the worst case as δ approaches the
   unit circle?
2. **De-globalise or not** (§9): auditability against reentrancy.
3. **`p = q = 0`** — the out-of-range write is diverted, not fixed (`BUG-0013`).
4. **The MA flip's fold** and its role in `BUG-0005` (spurious optimum): is the
   fold reachable by the optimiser under the current stopping tests?
5. **Whether the packing order should be machine-checkable** rather than
   documented — a consumer asserting `npar` against a declared layout would turn
   a silent miscomputation into a loud failure.
6. **Design attribution** (§0): whether the specification language the cast
   implements was designed jointly. Answerable from sources outside the code.
