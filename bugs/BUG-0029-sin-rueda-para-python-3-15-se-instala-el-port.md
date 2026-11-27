---
id: BUG-0029
title: Binary wheels stop at cp314; on CPython 3.15 pip and uv install the py3-none-any wheel and every estimate runs on the port
status: fixed
severity: high
component: packaging
found_in: 0.1.22
fixed_in: unreleased (after 0.1.22)
reported: 2026-11-27
reporter: David / Claude — the P08 run of the course (UCM), with drvarma BUG-0031
tags: [packaging, wheels, engine, python-3.15]
references:
  - pyproject.toml [tool.cibuildwheel] build / skip
  - .github/workflows/wheels.yml
  - bugs/BUG-0029-repro/repro.py
  - drvarma BUG-0031 (the same defect), BUG-0024 (the fallback warning)
---

## Summary

`[tool.cibuildwheel]` builds `cp310-* … cp314-*` and a `py3-none-any` wheel is
published next to them. On CPython 3.15, final since 2026-10-09, pip and uv
take the pure wheel and fue estimates on its Python port. The course installs
the suite with `uv tool install atsw==…`, which takes the newest Python.

## Impact

In the P08 the art scripts of the groups whose drvarma ran on the port show
fue on «motor Python» too: slower, and on a flat likelihood possibly another
optimum (BUG-0005). The warning (BUG-0024) goes to stderr, which an MCP
client does not show.

## Reproduction

`python3 bugs/BUG-0029-repro/repro.py` reads pyproject.toml and setup.py
offline and exits 1 while a final CPython is outside the tested list or the
wheel is not abi3.

## Root cause

The build list is a fixed set of minor versions.

## Fix

The extension already compiles against the limited API (cffi's setuptools
integration sets `py_limited_api`); the wheel is now tagged `cp310-abi3`
(setup.py), one binary per platform for every CPython from 3.10 on.
cibuildwheel 4.3.5 builds it on cp310 and smoke-tests it on cp310–cp315, with
abi3audit. The pure wheel stays (the port is fue on every platform without a
binary). Free-threaded builds and 32-bit Windows still get the port.

## Validation

The repro exits 0; a local build gives `fue-…-cp310-abi3-linux_x86_64.whl`
with 0 abi3audit violations.
