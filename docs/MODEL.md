# The model, formally

*The class of models `fue` estimates, written without ambiguity. It is the
translation and update of §"Nociones de Series Temporales" of Treadway, A. B. (2011), *FUE: manual de usuario*, unpublished user manual,
with the notation aligned to Guerrero, D. E. (2026), *Hybrid seasonal models: critical values for testing deterministic versus stochastic seasonality frequency by frequency*, unpublished manuscript, Universidad Complutense de Madrid. Where this document and the
parser disagree, see `docs/FILE_CONTRACT.md`, which is generated from the
reading order of the code.*

---

## 1. Transformation

Let `Yₜ` be the original series, `t = 1 … N`, and `B` the backshift operator,
`Bxₜ ≡ xₜ₋₁`. The one-parameter Box-Cox family:

```
        ⎧ (Yₜ^λ − 1)/λ     λ ≠ 0
  zₜ =  ⎨
        ⎩ ln Yₜ            λ = 0
```

`zₜ` is the **transformed variable**. In economics λ=0 is the common case and
λ=1 (the identity, up to a constant) appears occasionally.

⚠ **λ is specified, not estimated.** To estimate it, run `fue` repeatedly on a
grid of λ and compare likelihoods in the Box & Cox (1964) sense — the
likelihood values are comparable only after the Jacobian correction, which is
why the program does not do it silently.

There is a second scale factor, `refactor`, which is **not** part of the model:
it is a numerical remedy, and `docs/FILE_CONTRACT.md` §2.7 explains what it
scales along with the series.

## 2. Two components, additive in the transformed level

```
  zₜ ≡ ξₜ + Nₜ
```

with `ξₜ` purely **deterministic** and `Nₜ` **stochastic**. `fue` specifies and
estimates the parameters of both. The additivity is in the *transformed* level:
with λ=0 it is multiplicative in the original units, which is the usual reason
for taking logs in the first place.

### 2.1 The deterministic component

```
  ξₜ ≡ Σᵢ ξᵢₜ ,   i = 1 … n_ξ
```

including two components (`cos f`, `sin f`) for each seasonal frequency
`f = 1 … s/2 − 1` of deterministic seasonality, with s = 4 or 12.

Each `ξᵢₜ` is a **linear transfer function** (s,b,r) with b=0 acting on a
non-parametric deterministic input:

```
  ω_s(B)/δ_r(B) ,   ω_s(B) = ω₀ − ω₁B − … − ω_sB^s
                    δ_r(B) = 1 − δ₁B − … − δ_rB^r
```

with the stability condition `|δ_r(B)| = 0 ⇒ |B| > 1`. In most cases the
transfer function is a single constant ω₀ (s = r = 0).

**This is the second thing `fue` has that the alternatives do not.** With
`r > 0` the response to an input is *dynamic*: the denominator makes the effect
build up or decay geometrically instead of switching on and staying put. A
static regression coefficient — which is what an `exog` column gives you
elsewhere — is the special case δ(B) = 1, and the difference is not cosmetic:
the long-run gain

```
  g = ω(1)/δ(1)
```

is a parameter you can test (`docs/FORMAL_TESTS.md` §6), and it is what says
whether an incident moved the level permanently or was absorbed.

⚠ **The `s` of ω_s(B) is not the `s` of the seasonal period.** The manual warns
about this and it is worth repeating: the same symbol carries two meanings and
only the context separates them.

The inputs, with `t*` a date:

| keyword | ξₜ |
|---|---|
| `impulse` | 1 at t = t*, 0 elsewhere |
| `step` | 0 for t < t*, 1 for t ≥ t* |
| `ramp` | 0 for t < t*, 1 + (t − t*) for t ≥ t* |
| `compimp` | +1 at t*, −1 at t*+1, 0 elsewhere |
| `trend` | the time index |
| `easter` | a specification of the Easter effect (monthly data) |
| `cos f`, `sin f` | the harmonic pair at frequency f |
| `alter` | the alternator, (−1)ᵗ |
| *custom* | a non-standard input the user supplies as an extra data column |

`compimp` is a **different regressor** from `impulse` — folding them together
silently estimates another model, which is `bugs/BUG-0006`.

### 2.2 The stochastic component: ARIMA(p,d,q)(P,D,Q)ₛ

The starting point for `Nₜ`. Let

```
  wₜ ≡ ∇^d ∇ₛ^D Nₜ ,   ∇ ≡ 1 − B ,   ∇ₛ ≡ 1 − B^s
```

be stationary, with d = 0,1,2 and D = 0,1 in practice. `wₜ` need not have zero
mean; **μ is the mean of wₜ**, so that `wₜ − μ` does. Then

```
  φ_p(B) Φ_P(B^s) (wₜ − μ) = θ_q(B) Θ_Q(B^s) aₜ
```

| operator | |
|---|---|
| `φ_p(B) = 1 − φ₁B − … − φ_pB^p` | AR(p) |
| `Φ_P(B^s) = 1 − Φ₁B^s − … − Φ_PB^{sP}` | AR(P)ₛ |
| `θ_q(B) = 1 − θ₁B − … − θ_qB^q` | MA(q) |
| `Θ_Q(B^s) = 1 − Θ₁B^s − … − Θ_QB^{sQ}` | MA(Q)ₛ |
| `aₜ ~ iid N(0, σ²)` | Gaussian white noise |

with no common factors between the AR and MA structures, and

```
  φ_p(B) Φ_P(B^s) = 0 ⇒ |B| > 1      stationarity
  θ_q(B) Θ_Q(B^s) = 0 ⇒ |B| > 1      invertibility
```

⚠ **The likelihood is also defined on the non-invertibility boundary**, and
`fue` can return MA estimates with |B| = 1 as legitimate maximum-likelihood
values. Treadway's manual says the user must know how to handle that
competently; this documentation adds the machinery for it —
`docs/CONVERGENCE.md` for how a fit reports where it stopped, and
`docs/FORMAL_TESTS.md` §3 for the DCD test that decides whether the boundary is
where the parameter belongs.

**μ scales with the whole differencing operator.** Because `ifadf` (§2.4) is
differencing too, μ is the mean of the *fully* differenced variable. Printing
the factors outside the μ parenthesis once made a correct model look
inconsistent — `art`'s BUG-0012.

### 2.3 Generalised seasonality (MEG)

Define the **annual moving sum**

```
  S(B) ≡ Σ_{k=0}^{s−1} B^k
```

If `S(B)Nₜ` is what has a stationary and invertible representation, `Nₜ` has
*stochastic seasonal non-stationarity at every frequency f = 1 … s/2*.

Treating that by applying the annual difference is a blunt instrument, and the
manual is precise about why:

```
  ∇ₛ ≡ (1 − B) · S(B)
```

so ∇ₛ applies, **on top of the annual moving sum, a regular difference that has
nothing to do with seasonality of any kind** — and if `wₜ` has an invertible
representation, `∇wₜ` does not. Hence the name: ∇ₛ is the **annual difference**,
not the seasonal difference.

The annual moving sum factors into s/2 irreducible factors, one per frequency
`f = 1 … s/2` in cycles per year:

```
  S(B) = (1 + B) · ∏_{f=1}^{s/2−1} (1 − 2cos(2πf/s)·B + B²)
```

The factor at f = s/2, `(1 + B)`, is undamped alternation; the rest are
second-order with imaginary roots — undamped oscillation at their frequency.

**`fue` lets you apply these one frequency at a time.** That is the field
`ifadf` in the input file, and it is what no other program in this family
offers: if `Nₜ` (s=12) has stochastic seasonality at f = 1 and 5 only, then

```
  (1 − √3·B + B²)(1 + √3·B + B²) Nₜ
```

is what has a stationary and invertible representation — and the other
frequencies stay deterministic.

⚠ **Count the degree, not the flags.** An interior factor costs two
observations, the Nyquist one costs one. Computing the loss as `d + D·s` and
ignoring `ifadf` is a real defect that has occurred twice (`drtran` BUG-5 and
BUG-9).

#### Why MEG rather than (P,1,Q)ₛ

The MEG generalises ARIMA(p,d,q)(P,1,Q)ₛ. The seasonal operators factor as

```
  Φ₁(B^s) = (1 − Φ₁^{1/s}B) · ( Σ_{k=0}^{s−1} Φ₁^{k/s} B^k )
  Θ₁(B^s) = (1 − Θ₁^{1/s}B) · ( Σ_{k=0}^{s−1} Θ₁^{k/s} B^k )
```

where the first factor on the right is a **regular** AR(1)/MA(1) — which has
nothing to do with seasonality — and the second is a weighted annual moving sum.

So `(1 − Θ₁B^s)` carries **one** parameter, and that single parameter appears as
`Θ₁^{1/s}` in a regular MA(1) *and* s−1 more times, at different powers, across
the seasonal factors. The MEG relaxes exactly that: it drops the restriction
that all the roots of the MA(1)ₛ share a modulus, while the MA_f operators keep
the frequencies of the irreducible factors. The MA(1)ₛ is equivalent to a
regular MA(1) times s/2 irreducible factors; at f = s/2 that factor is an MA(1)
with positive parameter, the rest are second-order with imaginary roots. The
same reading applies to the AR(1)ₛ.

The advantage is one restriction fewer — a restriction that the (P,1,Q)ₛ form
imposes arbitrarily on every seasonal frequency at once.

### 2.4 The individual factors of the annual difference

In the input file this is a list of `s/2 + 1` flags, indices `0 … s/2`:

| index | factor | what it is | observations lost |
|---|---|---|---|
| 0 | (1 − B) | the regular difference, f = 0 | 1 |
| 1 … s/2−1 | (1 − 2cos ω_f B + B²) | oscillation at ω_f = 2πf/s | 2 |
| s/2 | (1 + B) | Nyquist alternation | 1 |

With every flag set the operator **is** ∇ₛ — verified against `D=1` to 1.1e-13.
For annual data the list is empty.

**Compare polynomials, not (d,D) tuples.** Two specifications with different
(d, D, ifadf) can be the same operator, and the only reliable comparison is the
product.

---

## 2.5 The cast, which is the part with no equivalent

Everything above is a *specification*, and the piece that turns it into
something estimable is `cast_us` — the routine that maps a free parameter vector
into the full structure and back. It is worth naming because it is what has no
counterpart elsewhere:

* **operators enter factored, not expanded.** `2 1 1` is `(1−φ₁B)(1−φ₂B)`, two
  first-order factors, and not one AR(2). The factorisation is the hypothesis;
  expanding it would estimate a different model;
* **any coefficient can be fixed** at a value — that is the `0` flag beside it —
  which is what makes every likelihood-ratio test in `docs/FORMAL_TESTS.md`
  possible: the restricted fit is the same model with one flag flipped;
* **operators can be constrained to a frequency** (`AR_f`, `MA_f`), so a
  second-order factor is estimated with its root *pinned* to ω_f and only its
  modulus free;
* **the annual difference is a list of factors**, not a power (§2.4);
* **each deterministic input carries its own rational transfer function**, with
  its own orders and its own free/fixed pattern.

All of that is one parameter vector handed to the optimiser, and one exact
likelihood evaluated on the result. The usual alternative — a state-space
package with an `exog` matrix — can express none of the five: it estimates
expanded polynomials with free coefficients and static regressors.

This is the reason the two applications in `docs/README.md` could be written the
way they were: a convergence *transition path* is a restricted transfer function
whose restrictions are the economics, and the cast is what lets those
restrictions be imposed and tested rather than assumed.

## 3. What fue estimates, and how

Given the specification above, `fue` maximises the **exact unconditional
Gaussian likelihood** of the stationary variable — not a conditional sum of
squares, and not an approximation:

```
  l(Φ,Θ,μ,σ²,Q|w) = −½{ n·log(2πσ²) + log|ΛᵀΛ| + S(Φ,Θ,μ,Q|w)/σ² }
```

which is equation (2) of Mauricio (1997), Algorithm AS 311. `docs/PROVENANCE.md`
§3 says which routine computes which piece and what verifies it; the short
version is that the likelihood is Mauricio's published algorithm, checked
against Melard's published FORTRAN (AS 197) to 1e-14 on the quadratic form.

Two consequences worth stating for anyone comparing `fue` with other software:

* **Exact ML is not the Box-Jenkins criterion.** Box & Jenkins estimate by
  nonlinear least squares with backforecasting, and on a series whose MA sits
  near the non-invertibility boundary the two criteria pull apart — measured, on
  the Spanish CPI airline model: Θ = 0.8147 (exact ML) against 0.9061 (BJ NLS).
  Neither is wrong; each sits at its own optimum. See
  `tests/test_taste_nls_criterion.py`.

* **ARMAX by exact ML, with seasonality frequency by frequency, is the thesis of
  this program.** `statsmodels.SARIMAX` gives whole-operator seasonality or
  none. That is the reason `fue` exists and the reason `ifadf` is a list rather
  than a flag.
