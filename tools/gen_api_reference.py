#!/usr/bin/env python3
"""Generate `docs/API.md` from the docstrings of fue's public API.

The reference is GENERATED, never edited: a hand-written API document is a
second source of truth that drifts from the first one, and the drift is silent.
`tests/test_docs_match_the_code.py` regenerates it and fails if the file on disk
differs, so the docstring is the only place to fix anything.

    python tools/gen_api_reference.py            # write docs/API.md
    python tools/gen_api_reference.py --check    # exit 1 if it is out of date
"""
import inspect
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_OUT = os.path.join(_ROOT, "docs", "API.md")

# Document THIS working copy, not whatever is installed.
#
# The generator imported `fue` from wherever Python found it, which is
# site-packages: a docstring fixed in `src/` did not reach the reference, and
# `--check` agreed, because both sides read the same installed package. That is
# how `datasets.ripc` published a description of the data that was wrong in the
# opposite direction — and nothing complained.
sys.path.insert(0, os.path.join(_ROOT, "src"))

#: The public surface, grouped as a reader meets it rather than alphabetically.
_GROUPS = [
    ("Building a model", [
        "TimeSeries", "Model", "Intervention", "FixedFreqFactor",
    ]),
    ("Reading and writing files", [
        "load", "load_fuf", "write_out", "write_fuf", "write_fuf_out",
        "write_forecast_report",
    ]),
    ("Diagnostics", [
        "acf", "pacf", "ljung_box", "jarque_bera",
    ]),
    ("Results", [
        "ForecastResult",
    ]),
]

#: Datasets. Not attributes of the package — they live in `fue.datasets` — but
#: every example starts with one, and their docstrings say what the numbers are,
#: which is exactly what a reader needs and what `ripc` got wrong until
#: 2026-08-13.
_DATASETS = ["ripc", "sfny"]

#: Modules re-exported at package level. Listed, not expanded: their contents
#: are reached through the objects above.
_MODULES = ["datasets", "diagnostics", "forecast", "inp", "intervention",
            "model", "report", "report_forecast", "series"]


def _signature(obj):
    try:
        if inspect.isclass(obj):
            return f"{obj.__name__}{inspect.signature(obj.__init__)}".replace(
                "(self, ", "(").replace("(self)", "()")
        sig = f"{obj.__name__}{inspect.signature(obj)}"
        # methods are printed as Class.method(...), without the receiver
        return sig.replace("(self, ", "(").replace("(self)", "()")
    except (TypeError, ValueError):          # builtins without introspection
        return obj.__name__


def _doc(obj):
    return inspect.getdoc(obj) or "*(no docstring)*"


def _methods(cls):
    """Public methods of a class, in definition order."""
    out = []
    for name, fn in inspect.getmembers(cls, inspect.isfunction):
        if name.startswith("_"):
            continue
        if fn.__qualname__.split(".")[0] != cls.__name__:
            continue                          # inherited, documented elsewhere
        out.append((name, fn))
    return sorted(out, key=lambda kv: kv[1].__code__.co_firstlineno)


def _version():
    """The version of the REPOSITORY, not of whatever is installed.

    `fue.__version__` comes from site-packages, so generating the reference in
    a working copy whose version has been bumped produced a document declaring
    the previous release — and the --check guard agreed with it, because it
    compared against the same installed package. The published API page of
    0.1.10 said "fue 0.1.9" until this was fixed.
    """
    import re

    p = os.path.join(_ROOT, "pyproject.toml")
    m = re.search(r'^version\s*=\s*"([^"]+)"', open(p, encoding="utf-8").read(),
                  re.M)
    return m.group(1) if m else ""


def build():
    import fue
    import fue.datasets

    lines = [
        "# API reference",
        "",
        "*Generated from the docstrings by `tools/gen_api_reference.py`. "
        "Do not edit: fix the docstring and regenerate. "
        "`tests/test_docs_match_the_code.py` checks that this file is current.*",
        "",
        f"`fue` {_version()}".rstrip(),
        "",
        "---",
        "",
    ]

    for titulo, nombres in _GROUPS:
        lines += [f"## {titulo}", ""]
        for nombre in nombres:
            obj = getattr(fue, nombre, None)
            if obj is None:
                continue
            lines += [f"### `{_signature(obj)}`", "", _doc(obj), ""]
            if inspect.isclass(obj):
                for mname, fn in _methods(obj):
                    lines += [f"#### `{obj.__name__}.{_signature(fn)}`", "",
                              _doc(fn), ""]
        lines.append("")

    lines += ["## Datasets", "",
              "Shipped with the package; `from fue.datasets import ripc`.", ""]
    for nombre in _DATASETS:
        fn = getattr(fue.datasets, nombre, None)
        if fn is None:
            continue
        lines += [f"### `fue.datasets.{_signature(fn)}`", "", _doc(fn), ""]
    lines.append("")

    lines += ["## Modules", "",
              "Re-exported at package level; their contents are reached through "
              "the objects above.", ""]
    for nombre in _MODULES:
        mod = getattr(fue, nombre, None)
        if mod is None:
            continue
        primera = (_doc(mod).split("\n")[0]) if mod.__doc__ else ""
        lines.append(f"* **`fue.{nombre}`** — {primera}")
    lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def main():
    texto = build()
    if "--check" in sys.argv:
        actual = open(_OUT, encoding="utf-8").read() if os.path.exists(_OUT) else ""
        if actual != texto:
            print("docs/API.md is out of date — run tools/gen_api_reference.py")
            return 1
        print("docs/API.md is current")
        return 0
    with open(_OUT, "w", encoding="utf-8") as fh:
        fh.write(texto)
    print(f"{_OUT}: {len(texto.splitlines())} lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
