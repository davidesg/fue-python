# The file contract — `.inp`, `.out`, `.pre`

*The authoritative grammar, taken from the parser (`src/fue/inp.py`), which
reproduces the reading order of `fue-1.13.1/src/fue.c` §3.0–3.7. Where this
document and the parser disagree, the parser is right and this document is a
defect.*

A user-facing companion, about which program reads which file and why, is in
[`atsw-suite/docs/FILE_FORMATS.md`](https://davidesg.github.io/art-python/FILE_FORMATS/).
This one is the reference: field by field, in order.

---

## 1. The rule that governs everything

Stated by Arthur B. Treadway in 2001, in the DRVUS manual, and still exactly
true:

> *El programa no interpreta los comentarios escritos entre asteriscos,
> solamente lee los números que espera encontrar en la posición correcta. Esto
> quiere decir que si se cometiera el error de introducir una línea en blanco,
> por ejemplo, al comienzo del fichero de entrada, el programa bloquearía al
> ordenador y el proceso de estimación no se llevaría a cabo.*

**The parser is positional.** Lines beginning with `**` are separators: their
text is ignored, their *presence* is not. Everything else is a data line, read in
a fixed order. Reordering sections, dropping a separator or inserting a blank
line does not produce an error message — it produces a different model, or a
crash.

Two consequences worth stating plainly:

* **The comments are for you, not for the program.** You may rewrite them freely;
  you may not move them.
* **There is no schema validation.** A `.inp` that parses is not thereby a
  correct `.inp`. This is the strongest argument for never writing one by hand
  when a program can write it — see §5.

---

## 2. `.inp` — the specification

The reading order below is the parser's. `**` marks a separator line.

### 2.0 Header

The first lines are skipped until a separator whose text contains `frequency`.
That is the only place where the parser looks at what a comment *says*, and it is
how the fixed preamble is tolerated.

```
************************************************
* Input file for program FUE                   *
* DOCTYPE ATSW-interface SYSTEM                *
************************************************
```

### 2.1 Frequency, sample and name

```
** Frequency of time series: either 1(A), 4(Q) or 12(M):
 12
** Number of observations and starting date of time series:
 216  1 2002 IPC_ES
```

| field | meaning |
|---|---|
| frequency | `1` annual, `4` quarterly, `12` monthly. The literal `number` means unnumbered data and is read as frequency 1 |
| line 2 | `nobs`, start period, start year, name |

⚠ **Annual files put the year twice.** For `freq = 1` the parser reads
`nobs, <ignored>, begyear, name` and forces `begtime = 1`. Writing the year in
the second field is the historical form; `1` is the current one and both load.
This is where BUG-0018 lived.

### 2.1b Forecast horizon — `fuf` files only

```
** Forecast horizon and estimated innovation variance
 24  0.0625
```

Present only in `fuf` inputs, inserted between the sample line and the
deterministics. `fue.load()` detects it by key and `fue.load_fuf()` requires it.

### 2.2 Deterministic variables

```
** Number of deterministic variables (including seasonal components):
11
**
cos 1
sin 1
...
alter
**
0 0 0 0 0 0 0 0 0 0 0
**
0.000000  1
...
**
0 0 0 0 0 0 0 0 0 0 0
```

In order: the **count**; then one line per variable with its **type and date**;
then the **orders of ω(B)** on a single line; then one `**`-separated block of
**ω coefficients** per variable; then the **orders of δ(B)**; then, if any order
is non-zero, the **δ coefficients** in the same shape.

Types the parser knows:

| keyword | argument | what it is |
|---|---|---|
| `impulse` | date | a single spike |
| `compimp` | date | the **compensated** impulse: +1 then −1. **A different regressor from `impulse`** — folding the two together silently estimates another model (fue BUG-0006) |
| `step` | date | permanent level shift |
| `ramp` | date | linear from the date on |
| `cos` / `sin` | harmonic index | seasonal pair at frequency *f* |
| `alter` | — | the Nyquist alternator (−1)ᵗ |
| `easter`, `trend` | — | calendar and deterministic trend |
| *anything else* | — | a **custom** regressor, whose data arrives as an extra column in the series block (§2.7) |

`pulse` is accepted as a synonym of `impulse` because this package once spelled
it that way; **the C does not know that keyword**. Do not write it in files meant
to travel.

Every coefficient is a pair: `value  flag`, where the flag is `1` to estimate and
`0` to hold fixed. The same convention runs through every section below.

### 2.3 The ARMA structure

Six sections, always present, always in this order, each `count [orders…]`
followed by one `**` block of coefficients per factor:

1. regular AR
2. annual/seasonal AR
3. regular MA
4. annual/seasonal MA
5. fixed-frequency regular AR(2) — `count [frequencies…]`
6. fixed-frequency regular MA(2)

```
** Number and orders of regular AR operators:
 1 1
**
0.402839  1
```

A zero count is written as a bare `0` and consumes no further block. **Some
structure must be specified even if null** — the sections cannot be omitted.

Between them the parser tolerates four optional sections that `fue-1.13.1`
carries commented out (seasonal AR / MA with fixed frequency, and the annual
f-fixed pair). They are detected by key and skipped.

**Factorised operators.** `2 1 1` means *two* first-order factors, i.e.
`(1−φ₁B)(1−φ₂B)`, not one AR(2). `1 2` means one second-order factor. The
distinction is the whole point of the notation and it changes the model.

### 2.4 The mean

```
** Mean parameter (mu):
0.154472  1
```

μ is the mean of the **fully differenced** variable. That matters: `ifadf` is
differencing too, so μ scales with the gain of the whole operator at B=1 — see
`art`'s BUG-0012, where printing the factors outside the μ parenthesis made a
correct model look inconsistent.


⚠ **The seed matters, and a stale one can cost you the optimum.** μ is a
starting value like any other, and `a1.inp` of the Box-Jenkins bank ships
μ₀ = 2.5 against a series whose mean is 17.06. From that seed fue stops on the
AR boundary 6.86 in log-likelihood below the published estimate; from μ₀ = 17
—or 0, or anything in 6…20— it lands on it in seven iterations. That is
`bugs/BUG-0012`, and the lesson is general: **seed μ from the mean of the
differenced variable**, which is what `art` does and why the case never appears
along its path.

### 2.5 Box-Cox and differencing

```
** Box-Cox lambda, regular differences and complete annual differences:
0.00 1 0
```

λ, `d`, `D`. λ=0 is the log; λ=1 the identity.

### 2.6 The individual factors of the annual difference

```
** Individual factors of the annual difference (from freq 0.0):
 0 0 0 0 0 0 0
```

`freq//2 + 1` flags, indices `0 … s/2`, each turning on one factor of the
factorised annual difference:

```
∇₁₂ = (1−B)(1+B)·∏_{f=1}^{5}(1 − 2cos ω_f·B + B²),   ω_f = 2πf/12
```

Index 0 is `(1−B)`, index `s/2` the Nyquist `(1+B)`, the rest the
complex-conjugate pairs. **This is the frequency-by-frequency seasonality**, and
it is the field no other program in this family has. With every flag set the
operator *is* ∇₁₂ — verified against `D=1` to 1.1e-13.

For annual data the list is empty.

⚠ **Count the degree, not the flags.** An interior factor costs **two**
observations, the Nyquist one. Computing the loss as `d + D·s` and ignoring
`ifadf` is a real defect that has occurred twice: `drtran` BUG-5 and BUG-9.

### 2.7 Bands, rescaling and the data

```
** ACF/PACF bands (0 Automatic) and reescaling factor:
 0 100.00
** Time series (stochastic and non-standard deterministic variables):
70.6600000000
70.9400000000
...
```

The second field of the bands line is `refactor`. **A zero is read as 1.0.**

On `refactor`, and it is not decoration: it is advice from Treadway, May 2001,
about the *gradient norm* —

> *Es deseable que la norma del gradiente sea cero hasta toda la precisión que
> ofrece el programa. Cuando esto no ocurre […] es muchas veces útil escalar los
> datos […] multiplicando todos los valores de la variable lnY […] por, p.e.,
> 100 […] también se multiplicará cada parámetro de intervención por el mismo
> factor, y la salida presentará una sigma multiplicado por el mismo factor.*

— so it is a numerical-precision remedy, and the compensating rules (every
intervention parameter scales; σ comes out scaled) are part of it. See
`docs/PROVENANCE.md` §5 and `RESCALING_ARCHITECTURE.md` in the `art-tseries`
repository.

Then `nobs` data lines. **Column 0 is the series**; any further columns are the
data of the custom deterministic variables of §2.2, in the order those variables
were declared.

### 2.8 Encoding

The parser reads the file with the platform default. Files written by the
original C on a Latin-1 system fail with `UnicodeDecodeError` — `bugs/BUG-0010`,
open. Until it is fixed, convert:

```python
open(dst, "w", encoding="utf-8").write(open(src, encoding="latin-1").read())
```

---

## 3. `.out` — the estimation record

Text report of a fit: the parameter table with standard errors and t-ratios, σ̂,
the log-likelihood, AIC and BIC, residual diagnostics and the ACF/PACF listing.
The Python port reproduces the C's layout, which is why `.out` files from either
can be diffed against each other.

**It is a record, not an input.** Nothing reads it back.

---

## 4. `.pre` — an optimum in re-runnable form

Same grammar as `.inp`. The difference is not syntactic but a **claim**:

> A `.pre` says *these values are an optimum*. An `.inp` says *these values are a
> starting point*.

`write_pre()` takes a fitted model and writes the estimates as the new initial
values, preserving the free/fixed flags and the original data. So a `.pre` is
simultaneously the record of a fit and the seed of the next one, which is what
makes the ladder work:

```
.inp  →  fue  →  .out  (the record)
                 .pre  (the optimum, re-runnable)
                   ↓
                 .inp for the next rung
```

**The invariant, and it is testable:** run `fue` on a `.pre` and the numbers do
not move. If they do, the file was not an optimum — either it was hand-edited,
or it came from a different specification. `drtran`'s `load_pre` measures exactly
this and reports the optimality gap.

**Only the program that estimated may assert an optimum.** That is why
`drtran.write_inp` writes `.inp` and not `.pre`: it emits a specification, and a
specification is not a claim about optimality.

---

## 5. Do not write these by hand

The format has no validation: it is positional, and a file that parses is not
thereby correct. Every program in the suite can write one —

```python
import fue
from fue.report import write_pre          # not re-exported at package level

ts, m = fue.load("model.inp")             # .inp or .pre, same reader
m.fit()
fue.write_out(m, "model.out")             # the record
write_pre(m, "model.pre")                 # the optimum, re-runnable
```

`write_out` is on the package; `write_pre` is not, and lives in `fue.report`.
The asymmetry has no reason behind it — it is worth fixing when the API
reference of the documentation plan is written.

— and `art`'s `create_inp` / `confirm_and_estimate` build them from a
specification. The one case that genuinely needs a hand-written file is a custom
deterministic regressor (§2.2), and even there the safe route is to write the
file with a program and add the column.

---

## 6. Historical note

The format is Mauricio's, from DRVUS (2000); the six-part description of the
input file is Treadway's, from the 2001 manual
(Treadway, A. B. (2001), *DRVUS: manual de usuario*, unpublished user manual). The extended sections of
`fue-1.13.1` (`compimp`, `easter`, `trend`, the fixed frequency factors,
`ifadf`) came later and are additive.

⚠ **But no DRVUS-era file loads today**, and this section claimed the opposite
until the claim was tried. Measured on three of them — `S.inp` (the example of
Treadway's manual, an ARMAX with two consecutive ramps), `M2.inp`, and the
Box-Jenkins series that ship with DRVUS 1.2.01 — every one fails with

```
ValueError: Unexpected end of .inp file
```

because the bands/`refactor` section of §2.7 **did not exist then** and the
parser reads it unconditionally, consuming the first observation in its place.
Insert those two lines and the file loads and estimates, so what is incompatible
is one section rather than the format. `bugs/BUG-0011`, together with
`bugs/BUG-0010` (the Latin-1 encoding), which is the other barrier on the same
files. Between the two, the historical archive — the series Mauricio and Treadway
worked with, and their reference `.out`s — is unreachable without hand-editing.
