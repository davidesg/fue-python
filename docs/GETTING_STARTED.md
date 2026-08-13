# Getting started

*Install, fit one model, read the output. Treadway's §"Cómo aprender a usar FUE"
made two points that have survived the port and are worth keeping at the top:*

> You need a minimum knowledge of univariate time-series analysis. You will
> learn more through practice.
>
> FUE is not very communicative about the causes of failure, but it is itself
> very reliable.

The second is less true than it was — `docs/CONVERGENCE.md` exists because the
engine now says why it stopped — but the first is unchanged, and no amount of
documentation replaces it.

---

## 1. Install

```bash
pip install fue
```

That fetches a wheel with the C engine compiled in and GSL linked statically:
nothing else to install. If no wheel matches your platform, `pip` builds from
source and you need a C compiler and GSL headers.

### Which engine am I running?

This matters more than it looks — see §5 — so check rather than assume:

```python
try:
    import fue._fue_engine          # noqa: F401
    print("C engine: Mauricio's, compiled")
except ImportError:
    print("pure-Python engine: the translation, ~10× slower")
```

Both compute the same likelihood and both optimise with `raxopt`; the Python one
is the fallback when the extension is not available. `docs/PORT.md` §1 explains
why the package carries both.

## 2. The first model

```python
import fue
from fue.datasets import ripc

ts = ripc()
print(ts)          # TimeSeries('RIPC', nobs=72, freq=12, start=(2002, 1))
```

A price index is I(1) — the level wanders, the monthly change does not — and it
is modelled in logs, ×100 so the differences read as percentages:

```python
m = fue.Model(ts,
              d=1,                      # one regular difference
              boxlam=0.0,               # 0 = log, 1 = identity
              refactor=100.0,           # the ×100
              ar=[[0.3]], ar_free=[[True]],
              mu=0.1, estimate_mu=True) # mu is the mean of the DIFFERENCED series
m.fit()
```

`fit()` writes the estimates back into the model, so after it `m.ar[0][0]` is the
estimate and not the seed. Read the result:

```python
r = m._result
r.params, r.std_errors      # in the order they were specified, mu last
m.loglik, r.aic, r.bic
r.termination, r.niter, r.gnorm
r.converged
```

⚠ **Always read `converged` and `termination`, not only the parameters.**
`converged` means the optimiser reached a maximum by the gradient criterion —
not merely that nothing crashed. A fit that stopped because the iterates froze
returns normally, with `converged = False` and a warning naming the reason. On
Box-Jenkins Series A that difference was 6.86 in log-likelihood
(`bugs/BUG-0012`).

## 3. Three seeds worth getting right

Most fits that go wrong go wrong here, not in the engine:

1. **μ.** It is the mean of the *fully differenced* variable — seed it from the
   data (`np.diff(y, d).mean()`), not from zero and not from the level. A stale
   μ is what sends Series A to the boundary.
2. **The differencing.** `d`, `D` and `ifadf` are a *specification*, not a
   search. `art` decides them; `fue` estimates what you specify.
3. **The intervention dates.** `at` is a **0-based observation index**, and an
   intervention is extra-sample information about something that happened — the
   program does not look for it.

## 4. Reading the `.out`

```python
fue.write_out(m, "model.out")
```

The layout is the C's, so a `.out` from either program can be diffed against the
other. Four blocks:

```
Estimation method      : exact maximum likelihood      ← what was maximised
Check for invertibility: constrained search

**** GRADIENT STOPPING CRITERIUM SATISFIED …           ← the verdict; read it
**** CONVERGENCE OBTAINED AFTER 22 ITERATIONS [GRADIENT NORM = 0.0000]

Coefficients for regular AR factor 1:                  ← estimate (s.e.) [slot]
      0.402839  (0.097146) [ 1]
Mean parameter (mu):
      0.154472

     logelf: -53.5086902793                            ← what the tests compare
```

The parameter **correlation matrix** further down is worth a look every time:
a ±1.000 off the diagonal says two parameters are chasing the same direction,
and that is a specification problem no optimiser can fix.

## 5. What to do when it does not work

| symptom | first thing to check |
|---|---|
| `RuntimeWarning: … criterio de paso` | the seeds, especially μ — `docs/CONVERGENCE.md` §3 |
| `ValueError: Unexpected end of .inp file` | the file is positional; a missing section shifts everything — `docs/FILE_CONTRACT.md` §1 |
| an MA parameter pinned at 1 or −1 | possibly correct: the boundary is a legitimate ML estimate. Test it — `docs/FORMAL_TESTS.md` §3 |
| numbers differ from another machine | if the likelihood is ill-conditioned they can, and it is not a bug — `docs/PROVENANCE.md` §2.2 |
| results differ from Box-Jenkins' book | different estimator: the book is nonlinear least squares with backforecasting, `fue` is exact ML. Two decimals is the agreement to expect |

Treadway's advice on asking for help still applies, and is still the right
protocol: **bring the `.inp` and the `.out`**, and the question you formulated.

## 6. Where to go next

* [`examples/`](../examples) — five graded examples, from the minimal flow to a
  mixed MEG. Examples 3-5 are simulated with a fixed seed so you can see what
  the program recovers against a known truth.
* [MODEL.md](MODEL.md) — what the class of models actually is.
* [FILE_CONTRACT.md](FILE_CONTRACT.md) — when you need to read or write a file
  by hand, and the reasons not to.
* [API.md](API.md) — every public symbol with its signature and docstring.
* [`art`](https://pypi.org/project/art-tseries/) — if what you want is not to
  specify the model yourself. `fue` is the engine; `art` is the criterion.
