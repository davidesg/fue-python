# fue

**ARMAX by exact unconditional maximum likelihood, with rational transfer
functions on the inputs and seasonality resolved one frequency at a time.**

Three things, and no other program in this space has all three:

1. **The exact unconditional likelihood** — not a conditional sum of squares and
   not an approximation — of a model that carries deterministic inputs.

2. **A rational transfer function on each input**, `ω_s(B)/δ_r(B)`, so an
   intervention has a *dynamic response*: an effect that builds up, decays or
   settles at a long-run gain `g = ω(1)/δ(1)`. `statsmodels`' `exog` gives a
   static coefficient — the effect is the same in the period of the event and
   ten years later. That is a different model, not a simplified one.

3. **Seasonality frequency by frequency.** `SARIMAX` gives you the whole
   seasonal operator or none of it; `fue` lets a series be stochastic at f=1 and
   deterministic at f=2, which is what real series usually are.

### What that versatility is for

Two published applications, both of which need the transfer function to be a
*model of a process* rather than a dummy variable:

* **García-Hiernaux & Guerrero (2021)**, "Price convergence: representation and
  testing", *Economic Modelling* **104**, 105641. A general notion of price
  convergence — steady-state and catching-up, weak and strong — in which the
  **transition phase is an exogenous deterministic input passed through a
  transfer function**: `ω` is how far the process travels, `δ` its speed, and the
  starting date is the input's date. The shape and starting point are then
  *estimated from the data* rather than assumed, and the definitions imply
  parameter restrictions that are tested. Applied to the aggregate price levels
  of Germany, France and Italy after the euro, and to transatlantic 19th-century
  wheat prices.

* **García-Hiernaux, González-Pérez & Guerrero (2023)**, "Eurozone prices: a tale
  of convergence and divergence", *Economic Modelling* **126**, 106418. The same
  framework on relative prices across the EMU, 2001-2020, identifying the date,
  shape and velocity of each convergent or divergent process: convergence for
  over 80% of the EMU between 2001 and 2011, and price-level divergence from
  2012 onwards.

Neither study is a regression with a break dummy. The object being estimated is
a **transition path** — a dynamic response with a shape, a speed and a long-run
gain — inside a model whose stochastic part is estimated jointly by exact ML.
That is the class of question `fue`'s casting opens.

The engine is José Alberto Mauricio's C, embedded and unmodified: Algorithm
AS 311 (1997) for the likelihood, AS 197 (Melard 1984) for the scalar case, and
his own quasi-Newton optimiser from *JASA* 90, 282-291. See
[PROVENANCE.md](PROVENANCE.md) for who wrote what and what checks it.

---

## What it fits

```
  zₜ = ξₜ + Nₜ                                     (Box-Cox transformed level)

  ξₜ = Σ  ω_s(B)/δ_r(B) · xᵢₜ                      inputs with a DYNAMIC response
                                                   (step, impulse, ramp, Easter,
                                                    harmonics, your own column)

  φ_p(B) Φ_P(B^s) (wₜ − μ) = θ_q(B) Θ_Q(B^s) aₜ ,  wₜ = ∇^d ∇ₛ^D Nₜ
```

with ∇ₛ replaceable, factor by factor, by the frequencies you actually want to
difference. [MODEL.md](MODEL.md) writes it out without ambiguity.

## Ten lines that run

```python
import fue
from fue.datasets import ripc

ts = ripc()                       # Spanish CPI, monthly, 2002-2007
m = fue.Model(ts, d=1, boxlam=0.0, refactor=100.0,
              ar=[[0.3]], ar_free=[[True]],
              mu=0.1, estimate_mu=True)
m.fit()

print(m.ar[0][0], m.loglik)                    # 0.430784  -48.995195
print(m._result.termination, m._result.gnorm)  # gradient criterion, 3.9e-06
```

Five graded examples are in [`examples/`](../examples), from this to a mixed MEG;
they run in CI, so they still work.

## Why the numbers can be trusted

Because that is checkable rather than asserted, and the checks are of four
different kinds:

| | what it proves |
|---|---|
| **the paper's own FORTRAN** | Melard's AS 197 listing, transcribed and compiled: agreement to **5e-08** in log-likelihood on nine Box-Jenkins specifications |
| **the paper's own identities** | AS 311's equations (2)-(4) on the engine's outputs; its quadratic form matches Melard's to **1e-14** |
| **an implementation with no shared ancestry** | `statsmodels` agrees on all nine to ~**1e-7** in the parameters |
| **the original program's preserved runs** | 28 archived `.out` files reproduced exactly — termination code, iteration count and gradient norm |

And the limits are documented too, because a document that lists only what works
is an advertisement: [PROVENANCE.md §6](PROVENANCE.md), and
[`bugs/`](../bugs), which is public.

## The documentation

| | |
|---|---|
| [GETTING_STARTED.md](GETTING_STARTED.md) | install, first model, reading the `.out` |
| [MODEL.md](MODEL.md) | the class of models, formally |
| [FILE_CONTRACT.md](FILE_CONTRACT.md) | `.inp` / `.out` / `.pre`, field by field |
| [FORMAL_TESTS.md](FORMAL_TESTS.md) | non-stationarity, non-invertibility, MEG; with the critical values |
| [CONVERGENCE.md](CONVERGENCE.md) | what the optimiser reports and what to do about it |
| [PROVENANCE.md](PROVENANCE.md) | which algorithm, from which paper, verified how |
| [PORT.md](PORT.md) | how the port was done, and what it found |
| [API.md](API.md) | every public symbol, generated from the docstrings |
| [MIGRATION.md](MIGRATION.md) | for users of the FUE in C: what changes and what does not |

`fue` is the **engine**. The Box-Jenkins-Treadway orchestration — identification,
diagnosis, the decision rules — is [`art`](https://pypi.org/project/art-tseries/);
transfer functions are `drtran`, and multivariate is `drvarma`. The separation is
deliberate: *the engine supplies the fits and stays out of the decisions.*

## Licence and lineage

GPL. The C is Mauricio's, the file format and the user manual are Arthur B.
Treadway's, and both authorised the licence. The port is 2026.
