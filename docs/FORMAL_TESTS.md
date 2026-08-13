# Formal hypothesis tests

*Translation and update of §"Contrastes Formales de Hipótesis" of Treadway
(2011). Full references are in `docs/PROVENANCE.md` §9; "Guerrero (2026)" below
is the unpublished manuscript on hybrid seasonal models, which is where the
per-frequency critical values come from. The tests are the same; what is new is that the critical values
the manual left as an empty table now exist, and that this document says which
ones are tabulated, which are simulated and which are interpolated.*

---

## 1. The rule that comes before every test

Treadway states it first, and it is not a preamble:

> A formal test should be applied **only on a model that is efficiently
> estimated, statistically adequate and parsimoniously parameterised**, and
> should be used to simplify the parameterisation only when the message of the
> test is very clear — not when the null is accepted or rejected marginally at
> one conventional level.
>
> From an over-parameterised model one can obtain any answer; parameterise
> enough and one obtains the answer one wants.

The operational form of that rule: **a formal test on an inadequate model is not
a weak test — it is not a test at all.** ADF, KPSS and the rest are *initial
specification* tools; the tests below belong at the end of the process, on a
model that has passed diagnosis.

`fue` computes none of these automatically. Each is a **likelihood ratio between
two runs of the program**: fit the model free, fit it again with the parameter
fixed at the boundary, and compare `logelf` in the two `.out` files. That is why
every section below shows the `.inp` fragment that imposes the restriction — a
coefficient line whose flag is `0` instead of `1`.

Not covered here: Jarque-Bera and Ljung-Box, whose statistics `fue` already
prints in the `.out`.

---

## 2. Non-stationarity of an AR(1) factor — Shin-Fuller

**H₀:** φ₁ = 1  **H₁:** φ₁ < 1

```
             ⎧ 2[ l(φ̂₁) − l(φ̄₁) ]     if φ̂₁ ≤ φ̄₁
  SF   =     ⎨
             ⎩ 0                       if φ̂₁ > φ̄₁

  φ̄₁ = 1 − 4/n
```

with n the number of observations of the potentially under-differenced model
under H₁, and φ₁ the parameter of the **largest** AR(1) factor in the model.
`l(φ̂₁)` is `logelf` from the unrestricted fit; `l(φ̄₁)` requires re-estimating
with φ₁ fixed at φ̄₁ — for n = 100, at 0.96:

```
** Number and orders of regular AR operators:
1 1
**
.96 0
```

Shin & Fuller (1998) show the asymptotics do not differ between an ARMA(p+1,q)
with a strictly stationary and invertible ARMA(p,q) factor and a plain AR(1),
and give evidence that the test dominates the alternatives in size and power.

**Critical values** — Shin & Fuller (1998), Table II. Reject H₀ when the
statistic exceeds them:

| n | 10% | 5% | 1% |
|---|---|---|---|
| 25 | 1.02 | 1.68 | 3.33 |
| 50 | 1.06 | 1.75 | 3.41 |
| 100 | 1.07 | 1.75 | 3.41 |
| 250 | 1.07 | 1.76 | 3.44 |
| 500 | 1.08 | 1.77 | 3.46 |

⚠ **The manual printed this table empty.** These are the published values, and
`art`'s `formal_tests.shin_fuller` interpolates them linearly between the five
sample sizes — an approximation for intermediate n, not a tabulation.

⚠ **The statistic is Φ̂₁ᵤ = L_free − L_constrained, not 2·ΔL.** The manual writes
the 2[·] form; `art` implements the paper's Φ̂₁ᵤ, and the critical values above
belong to Φ̂₁ᵤ. Using one convention with the other's table doubles or halves the
statistic.

---

## 3. Non-invertibility of an MA — DCD

For MA(1) and MA_f operators estimated inside the invertibility region, the
generalised likelihood ratio of Davis & Dunsmuir and Davis, Chen & Dunsmuir.

**Regular MA(1), positive parameter.** H₀: θ₁ = 1, H₁: θ₁ < 1

```
  DCD = 2[ l(θ̂₁) − l(θ₁ = 1) ]
```

```
** Number and orders of regular MA operators:
1 1
**
1 0
```

**Frequency-fixed MA_f.** H₀: λ_f = −1, H₁: λ_f > −1

```
  DCD = 2[ l(λ̂_f) − l(λ_f = −1) ]
```

The same fragment applies, with the parameter fixed at −1.

### Critical values, and which regime a frequency belongs to

This is where the update matters most, because the law is **not** the same at
every frequency. It is governed by the **order of the factor**, not by which
frequency it sits at — Monte Carlo evidence in Guerrero (2026); implemented
in `art.formal_tests`:

| regime | which frequencies | d.o.f. | pile-up | 10% | 5% | 1% |
|---|---|---|---|---|---|---|
| **real root** (s=1) | regular MA(1), f = 0, Nyquist f = s/2 | 1 | 0.6575 | 1.00 | 1.94 | 4.41 |
| **complex pair** (s=2) | interior f = 1 … s/2−1 | 2 | 0.616 | 1.11 | 2.04 | 4.52 |

The real-root values are ≈ n-invariant. The complex-pair values carry a mild
finite-sample dependence, tabulated by simulation and **interpolated between
these sample sizes**:

| n | 10% | 5% | 1% |
|---|---|---|---|
| 120 | 1.12 | 2.06 | 4.64 |
| 240 | 1.13 | 2.07 | 4.52 |
| 480 | 1.10 | 2.04 | 4.53 |
| 960 | 1.11 | 2.03 | 4.52 |
| → ∞ | 1.11 | 2.04 | 4.52 |

These supersede the interpolated values of the thesis (1.07 / 2.02 / 4.52).

### The published simulations, quoted

The tables above are what `art` implements; Guerrero (2026, Table 1)
publishes the whole simulation, with Monte Carlo standard errors, and the two
are close enough that the interpolation costs little — which is worth knowing
before anyone re-runs it:

| regime | n | pile-up | 10% | 5% | 1% |
|---|---|---|---|---|---|
| real root (f = 0, 6) | 120 | 0.651 (.002) | 1.00 (.02) | 1.95 (.03) | 4.42 (.09) |
| | 240 | 0.645 (.002) | 1.03 (.02) | 1.97 (.03) | 4.47 (.07) |
| | 480 | 0.651 (.003) | 1.01 (.02) | 1.91 (.03) | 4.38 (.09) |
| | 960 | 0.650 (.003) | 1.02 (.02) | 1.96 (.05) | 4.52 (.11) |
| complex pair (f = 1…5) | 120 | 0.615 (.002) | 1.12 (.01) | 2.06 (.02) | 4.64 (.04) |
| | 240 | 0.615 (.002) | 1.13 (.01) | 2.07 (.02) | 4.52 (.05) |
| | 480 | 0.617 (.002) | 1.10 (.01) | 2.04 (.02) | 4.53 (.06) |
| | 960 | 0.617 (.002) | 1.11 (.02) | 2.03 (.03) | 4.52 (.07) |
| Davis MA(1) | ∞ | 0.6575 | 1.00 | 1.94 | 4.41 |
| thesis interpolation | | | 1.07 | 2.02 | 4.52 |

And the realistic case (Guerrero 2026, Table 2) — the witness together with a mean and the
nine surviving harmonics, k = 10 regressors — which is the one to use on a real
model:

| n | pile-up | 10% | 5% | 1% |
|---|---|---|---|---|
| 120 | 0.559 (.002) | 1.63 (.02) | 2.87 (.04) | 5.81 (.10) |
| 240 | 0.586 (.002) | 1.34 (.02) | 2.39 (.03) | 5.16 (.11) |
| 480 | 0.601 (.002) | 1.23 (.02) | 2.20 (.03) | 4.67 (.08) |
| 960 | 0.611 (.003) | 1.13 (.02) | 2.10 (.04) | 4.70 (.10) |
| bare, n→∞ | 0.616 | 1.11 | 2.04 | 4.52 |

The gap between the two tables is the whole point: at n = 120 the 5% value moves
from 2.06 to **2.87** once the model carries what a real model carries.

### The seasonal AR_f statistic

The MEG sweep of §5 needs the other half — the analogue of Shin-Fuller at a
seasonal frequency — and Guerrero (2026, Table 3) tabulates it, for the statistic
Φ̂_f with fixed truncation ρ_m = 1 − c/n, 2×10⁴ replications:

| factor | deterministic | s | c | 10% | 5% | 1% |
|---|---|---|---|---|---|---|
| interior f = 1…5 | harmonic (2-dim) | 2 | **3** | 1.34 | 2.12 | 3.90 |
| | | | | [1.32, 1.37] | [2.09, 2.15] | [3.85, 3.94] |
| real (f = 0, Nyquist) | constant / alternator | 1 | **4** | 1.06 | 1.75 | 3.47 |
| | | | | [1.05, 1.06] | [1.70, 1.79] | [3.46, 3.48] |
| Shin-Fuller Table II, AR(1), c=4 | | | | 1.07 | 1.75 | 3.41 |

Two things to take from it. The real-frequency row reproduces Shin-Fuller's
Table II — as it must, since it is the same law — which is a check on the
simulation rather than a new result. And **the interior frequencies use c = 3,
not 4**: the truncation constant of the two-dimensional case is not the one of
the AR(1), and using ρ_m = 1 − 4/n there is testing a different point.

⚠ `art` does not implement this statistic yet: `formal_tests` carries the DCD
regimes and Shin-Fuller, and the MEG sweep decides on the MA_f witness. The
table is quoted here because the values exist and the sweep is where they
belong.

⚠ **Two production caveats, both measured, both easy to get wrong.**

1. **Compute the boundary likelihood exactly.** The restricted fit must profile
   over a fixed grid rather than let the free MA optimiser walk to the boundary:
   the optimiser is biased at the second-order non-invertibility boundary and
   produces a spurious pile-up of ≈0.82 against the correct ≈0.62.

2. **In a realistic model — one carrying a mean and deterministic harmonics —
   the correct finite-sample critical values are HIGHER**: at n = 120, 1.63 /
   2.87 / 5.81 at 10/5/1%, an effect that vanishes as n grows. Using the bare
   values there over-rejects.

---

## 4. Fixed frequency for an AR(2) with imaginary roots

Worth testing when the estimated frequency of an AR(2) factor is close to one of
f = 1 … s/2−1 relative to its standard error.

**H₀:** f = 1 (say)  **H₁:** f ≠ 1

```
  LR = 2[ l(H₁) − l(H₀) ]     ~ χ²₁ under H₀
```

`l(H₁)` from the free AR(2) with imaginary roots, `l(H₀)` from the same model
with the factor constrained to the fixed frequency (an AR_f operator). Reject at
1−α when the statistic exceeds the χ²₁ quantile.

This one is an ordinary LR with a standard distribution, because the restriction
is interior — unlike §2 and §3, which sit on a boundary and therefore have
non-standard laws.

---

## 5. Stochastic seasonality — the MEG sweep

Starting from a model with **fully deterministic** seasonality, evaluate
frequency by frequency whether the seasonality is stochastic instead.

Two ordering rules, and they matter:

* **Evaluate stochastic seasonality BEFORE simplifying the deterministic
  seasonal terms**, so the deterministic representation is as flexible as
  possible while you test.
* **Simplify the intervention terms and the ARMA structure FIRST**, so the model
  is as parsimonious as possible when you start.

For each frequency f = 1 … s/2, estimate one model that adds **both**:

* the homogeneously non-stationary AR_f operator at that frequency, and
* an MA_f operator at the same frequency as an **over-differencing witness**
  (seed λ_f at, say, −0.9).

Applying the non-stationary AR_f annihilates the deterministic seasonal terms at
that frequency, which is what makes the comparison clean. For monthly data that
is six models, one per frequency, each with its DCD test on the MA_f.

**Decision:** if the MA_f comes out literally non-invertible — or the DCD says
λ̂_f does not differ significantly from the non-invertible value — reject
stochastic seasonality at that frequency.

⚠ **Do not apply this rigidly**, and Treadway is explicit:

> Even when the DCD tests say an MA_f is invertible, with and without the
> corresponding AR_f overfit, it is sometimes worth integrating at that
> frequency anyway, to see whether the representation loses anything noticeable.
> Faced with ambiguity — the data not discriminating between stochastic and
> deterministic — the analyst should choose the representation that suits a
> pre-established criterion.
>
> **No program and no algorithm can replace the analyst's judgement in building
> models from data.**

There is a companion to this in the confirmatory direction: at f = 0 the
Shin-Fuller test (§2) and the DCD test (§3) have **opposite nulls** and bracket
the quasi-cancellation band, so their disagreement is itself the diagnostic.
See Guerrero (2026), the comparison table.

---

## 6. Simplification of intervention terms

When an intervention carries more than one parameter, test the hypothesis of
zero **long-run gain** with a Student t.

```
  g ≡ ω_s(1) / δ_r(1)
```

the long-run effect on the output of a permanent unit increase in the input.
Imposing g = 0 on two consecutive steps, for example, is the same as specifying
that intervention as a single impulse in the level.

The point is not only saving parameters: it is that the simplified form is the
one that says what actually happened. An intervention you can interpret against
the extra-sample information about the incident is worth more than one that fits
marginally better.

---

## 7. Where these are implemented

`fue` computes none of them automatically — it provides `logelf` and the
parameter table, and the tests are ratios between runs.

| test | implemented in | notes |
|---|---|---|
| Shin-Fuller (§2) | `art.formal_tests.shin_fuller` | Φ̂₁ᵤ convention, Table II interpolated |
| DCD, DCD_f (§3) | `art.formal_tests.dcd`, `dcd_overdiff_regular` | regime chosen by factor order |
| MEG sweep (§5) | `art.formal_tests.meg` | sweep or explicit frequency list |
| fixed frequency (§4) | — | LR by hand; χ²₁ |
| gain t-test (§6) | `drtran` reports the gain and its s.e. | |

The division is deliberate and is the architecture of the suite: **`fue` is the
engine, `art` is the criterion.** A likelihood ratio needs two fits and a
decision rule; the engine supplies the fits and stays out of the decision.
