# Porting fue to Python — the record of the process

*What was done, what changed and what did not. Written because it exists nowhere
else and because a port that cannot be audited is a rewrite with a familiar
name. The model is `drtran/docs/PORTE.md`.*

---

## 1. The shape of the thing: two engines, one API

`fue` is not a translation of the C, and it is not a wrapper around it. It is
**both at once**, and that is the central decision of the port:

```
  fue.Model(...).fit()
        │
        ├── csrc/            Mauricio's C, embedded and compiled into a
        │                    cffi extension (_fue_engine).  3.742 lines.
        │                    Used whenever it is available.
        │
        └── src/fue/         a pure-Python re-implementation of the same
                             engine — elfvarma.py, cast_us.py, qnewtopt.py:
                             1.895 lines.  Used when the extension is not.
```

The selection is one `try/except ImportError` in `src/fue/_engine.py:34`. There
is no flag and no configuration: if the extension imported, it runs; if it did
not, the Python engine runs and the answers stay comparable.

**Why carry both.** A wrapper cannot be read — the algorithm stays behind a
`.so` and nobody can check that the port is faithful, which is the whole purpose
of this package. A pure rewrite loses the reference: there is nothing left to
compare against. Carrying both makes the comparison *routine* — it runs in the
test suite as `✅ cross` in `docs/PROVENANCE.md` §3 — and it costs one
`try/except`.

| | C engine | Python engine |
|---|---|---|
| likelihood | `csrc/internal/elfvarma.c` (AS 311) + `usmelard.c` (AS 197) | `src/fue/elfvarma.py` |
| casting | `csrc/fue_api.c: cast_us()` | `src/fue/cast_us.py` |
| optimiser | `csrc/internal/qnewtopt.c` (raxopt) | `src/fue/qnewtopt.py` — **raxopt too**, translated; L-BFGS-B optional |
| speed | reference | ~10× slower |

## 2. What is reused and what is ported

**Reused verbatim** — `csrc/internal/`: `elfvarma.c`, `usmelard.c`,
`drvmlest.c`, `qnewtopt.c`. These are Mauricio's files with a GPL header, an
encoding change and the `#include` rename. That claim is an **invariant, not a
statement**: `tests/test_c_core_matches_drvus.py` classifies every differing line
against a declared exception and fails on anything else. The single functional
change in the whole likelihood core is written out in full:

```c
-    eigenqr( a, n, wr, wi );
+    if ( n>1 ) {gsl_eigenqr( a, n, wr, wi );}
```

**Rewritten** — `nlatools.c`, around GSL: 764 lines against Mauricio's 1355.
It is the one module the provenance test makes no claim about, deliberately and
in writing.

**Newly written** — `csrc/fue_api.c` (1000 lines): the transport struct, the
`.inp`-to-model casting and the result packing. This is the layer that did not
exist in the original, because the original was a program and this is a library.

**Ported to Python** — the 1.895 lines above, translated from the C with the
AS 311 step markers `(a)`…`(k)` carried across so the two can be read side by
side (`tests/test_as311_published_identities.py` asserts the markers survive).

## 3. Porting decisions that are not translation

Things that were decided differently, on purpose:

1. **Both engines run raxopt, and that was the right call.** `qnewtopt.py` is a
   line-by-line translation of Mauricio's optimiser — Dennis & Schnabel (1983)
   A9.4.1, with `cdgrad`, `lnsrch`, `qrupdate` and the two stopping tests — and
   `estimate_py` calls it by default (`cast_us.py:586`). scipy's L-BFGS-B is
   available as `optimizer="lbfgsb"` and is **not** the default: it is there to
   escape a poor local optimum on purpose, not to estimate behind your back.

   Translating a published optimiser doubles the surface on which two engines
   can silently differ, which is the argument against; what settles it is that
   the alternative is worse. With a different optimiser, any disagreement
   between the engines is ambiguous — likelihood or search? — and the whole
   point of carrying two engines is that the comparison means something. It
   does: on Series A both stop at the same boundary point, with the same
   termination code and the same iteration count.

2. **`fit()` writes the estimate back into the model attributes.** In the C the
   model record and the result are separate; here, after a fit, `m.ar[0][0]` is
   the estimate, not the seed. Anything else guarantees the bug that
   `bugs/BUG-0004` records: forecasting from the seeds.

3. **The `.inp` reader is the single source of truth for the format.** Nothing
   else parses it — not the oracle harness, not `drtran`, not the examples. A
   second reader is a second grammar that drifts.

4. **Errors are exceptions, and non-convergence is not an error.** `ifault != 0`
   raises; a fit that stopped without reaching a maximum returns, with
   `converged = False` and a `RuntimeWarning` naming the reason. See
   `docs/CONVERGENCE.md`; the distinction was `bugs/BUG-0012`.

5. **`refactor` is carried, not silently applied.** It is a numerical remedy
   from Treadway (2001), and everything it scales — every intervention
   parameter, σ — scales with it. `docs/FILE_CONTRACT.md` §2.7 and
   `RESCALING_ARCHITECTURE.md` in the `art-tseries` repository.

## 4. What the port found in the C

A port that finds nothing has not looked. Twelve defect reports, of which these
are the ones that were **in the original program**, not in the port:

| | what | where it hurt |
|---|---|---|
| `BUG-0006` | `compimp` read as a plain `impulse`; `easter` and `trend` unsupported by the binding | a `.pre` estimated a **different model in silence** |
| `BUG-0007` | fue C's `.pre` writer omits `easter`, `trend` and custom variables | **fue C cannot re-read its own `.pre`**, and segfaults on it |
| `BUG-0008` | the reporting plots segfault on a degenerate (zero-variance) series | the program **dies instead of reporting** |
| `BUG-0009` | the embedded GSL had no error handler | a failed eigensolve **aborted the Python interpreter** |
| `BUG-0010` | the `.inp` reader assumed one encoding | **no DRVUS-era file loaded at all** |
| `BUG-0011` | the bands/`refactor` section did not exist in the DRVUS format, and the format carries no version | the reader consumed the first observation in its place |
| `BUG-0012` | the engine never returned the optimiser's verdict | `converged=True` on a fit 6.86 below the optimum |

The `nlatools.c` `tensor()` crash (fixed 15 June 2026, recorded in `TODO.md`) is
the sharpest of them, because it shows what a port is *for*: the standalone
binary never hit it, since the fixed-frequency sections of `cast_us()` in
`fue.c` are commented out and `q` therefore never exceeded the regular MA order.
The hybrid does expand them — and the latent out-of-bounds write in
`tensor(-q+1, …)` became a heap corruption. **The defect was always there; only
the port reached it.**

⚠ And the counterexample, which matters just as much: `BUG-0012` looked like a
defect of the port for a month and was not one. Mauricio's own C, rebuilt today,
fails identically. `docs/PROVENANCE.md` §2.2.

## 5. What the port verified about the original

The traffic is not one-way. Because the port could be run against things the
binary never was:

* **28 preserved `.out` files** reproduce exactly — termcode, iteration count and
  gradient norm, including one whose norm is 116330.0394.
* **Melard's published FORTRAN** agrees with the engine to 5e-08 in
  log-likelihood over the nine Box-Jenkins specifications
  (`tests/test_as197_published_fortran.py`).
* **The published identities of AS 311** hold on the engine's own outputs, and
  its quadratic form matches Melard's `SUMSQ` to 1e-14.
* **`statsmodels`** — no shared ancestry — agrees on all nine to ~1e-7 in the
  parameters.

## 6. Wheels, and why they are built the way they are

`cibuildwheel` (`.github/workflows/wheels.yml`, `[tool.cibuildwheel]` in
`pyproject.toml`), skipping `win32`, `manylinux_i686` and PyPy. GSL is linked
statically into the wheel so that `pip install fue` does not require a system
GSL — verified in a clean virtualenv on 0.1.8.

The wheels run `tests/test_smoke.py` and **not** the golden battery, and that is
deliberate: `tests/test_real_cases.py` pins values that are platform-dependent
where the likelihood is ill-conditioned (`Coint/R.4` is the example, and
`bugs/BUG-0012` is the explanation). A per-wheel gate on those numbers would fail
for reasons that have nothing to do with the wheel.

## 7. How to reproduce the comparison

```bash
# the two engines, on the same model
python -m pytest tests/test_cast_us.py -v

# the C core against Mauricio's DRVUS sources
python -m pytest tests/test_c_core_matches_drvus.py -v

# the algorithms against their publications
python -m pytest tests/test_as197_published_fortran.py \
                tests/test_as311_published_identities.py -v

# and the original binary, rebuilt
tools/reproduce_drvus_reference.sh 1.2.01 a1
```

## 8. What is missing

* **`✅ cross` is still a weaker claim than `✅ oracle`**, and for a reason that
  survives the two engines running the same optimiser: the Python engine is a
  *translation*, so a shared misreading of the C would agree with itself. That
  is what the checks against Melard's FORTRAN and against `statsmodels` are for.
* **`bugs/BUG-0005` is open**: the optimiser reaches a spurious optimum on
  multimodal seasonal-AR likelihoods, and which basin it lands in is
  platform-dependent. Now that the termcode is propagated, the diagnosis is at
  least visible; the fix is not written.
* **`nlatools.c` has no verbatim claim and no independent check.** It was
  rewritten, so the provenance test excludes it; nothing else covers it.
