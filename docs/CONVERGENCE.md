# Convergence — how to read a fit that stopped

*What `raxopt` reports when it stops, why the two stopping tests are not
interchangeable, and what the tolerances are actually made of.*

Until August 2026 this document could not have been written, because the engine
did not return the answer: `raxopt` announces its verdict through `outputv`,
which the Python binding sends to `/dev/null`. `FitResult.converged` meant
`ifault == 0` — *"the engine did not crash"* — and a fit that stopped 6.86 in
log-likelihood below the optimum came back as good. That is `bugs/BUG-0012`, and
the reading below is what came out of it.

---

## 1. The two tests are not the same question

From `csrc/internal/qnewtopt.c:230-247`, unchanged from Mauricio's source:

```c
max1 = |g_i| · (|x_i| + 1) / (|f| + 1)      ≤ gradtol   →  termcode 1
max2 = |Δx_i| / (|x_i| + 1)                 ≤ steptol   →  termcode 2
```

* **`max1` asks "am I at a minimum?"** — the scaled gradient. This is the
  convergence test.
* **`max2` asks "can I still move?"** — the relative step. This is a **stuck
  detector**, not a convergence test.

Both return a termination code, and that is the whole trap: they are not two
ways of saying the same thing. A run that ends on `max2` has stopped *because
the iterates froze*, wherever they froze.

| termcode | what raxopt means | is it a maximum? |
|---|---|---|
| 1 | the scaled gradient fell below `gradtol` | **yes** |
| 2 | the iterates stopped moving | not necessarily |
| 3 | the last line search found no lower point | no |
| 4 | the iteration limit was reached | no |
| 5 | five consecutive steps of maximum length | no — likely unbounded |

`fue` reports all of this: `FitResult.termcode`, `.niter`, `.gnorm`, and
`.termination` in words. `converged` is now `ifault == 0 and termcode in (0, 1)`,
and anything else raises a `RuntimeWarning` naming the reason, the iteration
count and ‖g‖. Engine faults are still exceptions; a fit that exists but is not a
maximum is not — it is something you must be able to decide about.

```python
m.fit()
r = m._result
r.termcode        # 1
r.termination     # 'criterio del gradiente satisfecho'
r.niter, r.gnorm  # 7, 3.9e-07
```

The `.out` report carries the same block the C wrote, wording included.

---

## 2. What the tolerances are made of

`fue` uses Mauricio's values from DRVUS 1.01 onwards (`fue_api.c:70-71`):

```c
gradtol = pow( DBL_EPSILON, 1.1 / 3.0 );    /* 1.82e-06 */
steptol = pow( DBL_EPSILON, 2.0 / 3.0 );    /* 3.67e-11 */
```

They are **statements about the arithmetic, not about the problem**, and the
three quantities that fix the scale are:

| | value | what it is |
|---|---|---|
| `macheps^(2/3)` | 3.7e-11 | the noise floor of the gradient: `cdgrad` uses **central differences** with step `macheps^(1/3)·max(\|x\|,1)`, so truncation and roundoff meet here |
| `√macheps` | 1.5e-8 | the accuracy with which a minimiser can be *located* at all: near a minimum `f` is quadratic, so an error ε in `f` is √ε in `x` |
| `macheps` | 2.2e-16 | below this a relative step does not exist — `x + Δx == x` |

Read the two tolerances against that scale and the design is plain:

* **`gradtol` = 1.82e-06 sits ~5·10⁴ above the gradient's noise floor.** It asks
  for a gradient that a finite-difference gradient can actually deliver.
* **`steptol` = 3.67e-11 sits three orders *below* the resolution limit**
  `√macheps`. It cannot fire while the search is still making meaningful
  progress — so when it does fire, the search really is stuck.

Dennis & Schnabel (1983) §7.2 is the source of both forms; the `1.1/3` exponent
is Mauricio's, slightly stricter than the textbook `1/3` (1.82e-06 against
6.06e-06).

### Why not fixed constants

DRVUS 1.0 used fixed values, and they were the wrong way round:

```c
gradtol = 1.0e-7;      /* stricter than the formula */
steptol = 1.0e-5;      /* six orders looser */
```

With `steptol = 1e-5` the stuck detector fires while the iterates are still
moving, and it becomes the de facto convergence test. Measured on Box-Jenkins
Series A, ARMA(1,1) on levels: DRVUS 1.0 stops at iteration 25, announces
**"PARAMETER STOPPING CRITERIUM SATISFIED"**, and leaves φ = 0.9999784856 with
‖g‖ = 2e-4 — a declaration of convergence on the AR boundary, 6.86 in
log-likelihood below the optimum. Changing those two lines, and nothing else,
turns 1.0 into 1.01: 64 iterations and −50.7450915148.

And there is a second reason, which is decisive and measurable. `cmacheps()` is
evaluated at run time, so the formula **adapts to the arithmetic it is running
on**. The same source, compiled two ways on the same machine:

| build | macheps | gradtol | steptol |
|---|---|---|---|
| 64-bit (SSE2) | 2.220e-16 | 1.82e-06 | 3.67e-11 |
| 32-bit x87 (80-bit registers) | **1.084e-19** | **1.11e-07** | **2.27e-13** |

On the wider arithmetic it tightens the gradient test 16× and the step test 161×,
by itself. A constant cannot do that: it silently means something different on
every build — and note that DRVUS 1.0's `1.0e-7` is almost exactly the 80-bit
value `1.11e-07`, which is what those constants were: hand-calibrated on the
hardware of the day. The formula is the correct generalisation of a number
Mauricio had already tuned.

> This is not an abstraction. `bugs/BUG-0012` is a case where the same source,
> at 64 and at 80 bits, reaches two different answers — and
> `docs/PROVENANCE.md` §2.2 explains what that costs a reference binary as
> evidence.

---

## 3. What to do when a fit stops on the step criterion

A `termcode 2` says the iterates froze. Three causes, in order of how often they
turn out to be the one:

1. **A stale starting value.** The commonest by far, and the cheapest to test.
   `a1.inp` of the Box-Jenkins bank seeds μ = 2.5 on a series whose mean is
   17.06; from that seed `fue` stops on the boundary, and from μ = 17 it reaches
   the published optimum by the gradient criterion in seven iterations. Seed μ
   from the mean of the **differenced** variable — see
   `docs/FILE_CONTRACT.md` §2.4 — which is what `art` does.

2. **Scale.** Treadway's advice from May 2001 was about exactly this symptom —
   *"es deseable que la norma del gradiente sea cero hasta toda la precisión que
   ofrece el programa; cuando esto no ocurre […] es muchas veces útil escalar los
   datos"* — and it is where `refactor` comes from. Multiply the series by 100
   and remember that every intervention parameter and σ scale with it
   (`docs/FILE_CONTRACT.md` §2.7). It is not always a cure: on Series A it makes
   things worse, because the problem there is the seed, not the scale.

3. **The model.** Near-cancelling AR and MA factors, an operator sitting on the
   invertibility boundary, or two parameters chasing the same direction — the
   parameter correlation matrix in the `.out` will show ±1.000 when that is what
   is happening. Simplify before re-estimating.

What **not** to do is raise the iteration limit. On Series A the search stops at
iteration 23 of 500, and with `steptol = 0` — the test disabled — a 64-bit build
runs all 500 iterations and ends at the same point. The step test is reporting
the situation, not causing it.

---

## 4. The known weakness, stated

Look again at the two tests:

```c
max1 = |g_i| · (|x_i| + 1) / (|f| + 1)
max2 = |Δx_i| / (|x_i| + 1)
```

The `+ 1.0` is `typx = typf = 1` — Dennis & Schnabel's *typical size* vectors,
hard-coded to one. When the parameters live on one scale that is harmless. When
they do not, it is not: on Series A, μ ≈ 17 and φ ≈ 0.9 enter the same test with
weights 18 and 1.9, so the gradient criterion is dominated by the mean — which
is precisely the ill-conditioned direction of that problem.

**This is a limitation of the tolerances as used, not of the formulas**, and it
is the one place where a real improvement is available. It is also not a bug fix:
`raxopt` is Mauricio (1995), *JASA* 90, 282-291, refereed and published, and
supplying a scaling vector changes the algorithm's behaviour on every model ever
fitted with it. It is carried in `TODO.md` as a **study**, with the burden of
proof that a study carries.

---

## 5. Where the numbers in this document come from

Every figure above is measured, and reproducible:

```bash
# the two builds and their tolerances, and the reference runs they produce
tools/reproduce_drvus_reference.sh 1.2.01 a1

# the port's verdict against 28 preserved C runs
python -m pytest tests/test_optimizer_termcode.py -v

# Series A: where it stops, and where it goes from a sane seed
python -m pytest tests/test_box_jenkins_series.py -v
```

See also `bugs/BUG-0012` for the full measurement, `docs/PROVENANCE.md` §2.1 for
what exactly was added to `qnewtopt.c` (three globals that record; no criterion,
no announcement, no numerical behaviour) and §2.2 for what all this means for
"agrees with the reference binary" as a form of evidence.
