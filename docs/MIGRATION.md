# Migrating from the FUE in C

*For people who already use Treadway's FUE. What is the same, what is called
differently, and what genuinely differs in the output.*

The short version: **your `.inp` files work unchanged, the `.out` is the same
layout, and the command line takes the same arguments.** The differences are at
the edges, and they are listed here rather than left to be discovered.

---

## 1. The command line

```bash
fue      model  [eml|aml] [chk|nochk] [-f [N]]      # the C program
fue-py   model  [eml|aml] [chk|nochk] [-f [N]]      # this package
```

Same arguments, same defaults (`eml`, `chk`, horizon 24), same convention of
omitting the `.inp` extension. Likewise `fuf` → `fuf-py`.

⚠ **Why the `-py` suffix, and why it will not be dropped.** `fue` and `fuf` are
the C programs and they live in `/usr/local/bin`; `pip install --user` puts its
scripts in `~/.local/bin`, which comes **first** in the PATH. Installing this
package as `fue` would silently leave you running the port when you believe you
are running the original — on the same file, in the same directory, with the
same output name. The two coexist on purpose, and being able to compare them is
the point of the whole exercise (`docs/PORT.md`).

| C | port | note |
|---|---|---|
| `eml` / `aml` | same | exact / approximate ML |
| `chk` / `nochk` | same | check MA invertibility |
| `-f [N]` | same | write the forecast `.inp` for `fuf` |
| `geom` | **not implemented** | the geometric transform of `fue.c` has no counterpart yet |

## 2. Files in, files out

**In.** The `.inp` is the same format, read by the same grammar — the parser
reproduces the reading order of `fue.c` §3.0-3.7 (`docs/FILE_CONTRACT.md`).
Since August 2026 the DRVUS-era files load too, without editing: the missing
bands/`refactor` section is tolerated (`bugs/BUG-0011`) and a Latin-1 file is
read as Latin-1 (`bugs/BUG-0010`). The nine Box-Jenkins specifications that ship
with DRVUS 1.2.01 are in the test suite precisely so this keeps being true.

**Out.**

| file | C | port |
|---|---|---|
| `.out` | yes | yes — same layout, diffable line by line |
| `.pre` | **disabled** in `fue-1.13.1` (`preputf` is commented out) | yes |
| forecast `.inp` with `-f` | yes | yes |
| `.tex`, `_dist.tex`, `_res.tex`, `.dvi` | yes | **not implemented** — the port warns and skips |
| plots | gnuplot | matplotlib, from the Python API |

The `.pre` is worth a paragraph. In the C it is commented out, and there is a
reason recorded as `bugs/BUG-0007`: the C's `.pre` writer omitted `easter`,
`trend` and non-standard variables, so **fue C could not re-read its own
`.pre`** and segfaulted on it. The port writes a complete one, and the invariant
that makes it useful is testable: run `fue` on a `.pre` and the numbers must not
move (`docs/FILE_CONTRACT.md` §4).

## 3. What differs in the numbers

Nothing, to the precision that matters — and that claim is measured rather than
asserted:

* on the **28 preserved `.out` files** of real work, the port reproduces the
  estimates, the termination code, the iteration count and the gradient norm;
* on the **nine Box-Jenkins specifications**, log-likelihoods agree with the
  archived DRVUS runs to 1e-10 or better in eight of them.

The ninth is `a1`, and it is the one thing worth knowing before comparing an old
`.out` with a new run: **archived output from the 2000s was produced by 32-bit
binaries carrying 80-bit intermediates.** On a well-conditioned likelihood this
is invisible; on `a1` it decides the outcome, and Mauricio's own C rebuilt today
gives what the port gives. `docs/PROVENANCE.md` §2.2, and
`tools/reproduce_drvus_reference.sh` if you want to rebuild the old arithmetic
and see it for yourself.

## 4. What the port adds

* **The optimiser's verdict reaches you.** `converged` now means the gradient
  criterion was satisfied, not merely that nothing crashed; a fit that stopped
  because the iterates froze says so (`docs/CONVERGENCE.md`).
* **A library, not only a program.** `fue.Model(...).fit()` — the `.inp` is one
  way in, not the only one.
* **It does not die on you.** The C aborts the process on a degenerate series
  and on a failed eigensolve (`bugs/BUG-0008`, `bugs/BUG-0009`); the port raises.
* **A pure-Python engine** when no compiled extension is available, running the
  same algorithm and the same optimiser.

## 5. What the port does not have

Said plainly, because a migration guide that omits this wastes your afternoon:

* the LaTeX reports (`.tex`, `_dist.tex`, `_res.tex`) and the `.dvi`;
* the `geom` transform;
* the gnuplot graphics as such — there are matplotlib equivalents through the
  Python API, but they are not the same files.

## 6. A five-minute check before you switch

Treadway's advice for learning FUE was to run it on the shipped examples and
compare with the shipped output; the same protocol works for the migration:

```bash
fue     model            # the original
fue-py  model            # this package
diff model.out model.out.new
```

If the two `.out` differ in more than the last digits, that is a finding and it
should be reported with both files attached — which is exactly the protocol the
manual asks for, and the one that produced most of `bugs/`.
