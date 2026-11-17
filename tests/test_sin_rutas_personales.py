"""Nothing that ships names a personal folder.

The sdist carries src/, csrc/, tests/, bugs/*.md and README.md
(MANIFEST.in). Until 0.1.22 five of its files named folders of the
author's machine -- one of them a private analysis. A test that needs a
local file reads its path from an environment variable (ATSW_GUI,
FUE_C_BIN, FUE_SPAIN_FUF) or finds it next to the repository, and skips
without it (docs/PUBLISHING.md: no personal path in the sdist)."""
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHIPS = ("src", "csrc", "tests", "bugs", "README.md")
#: written split, so that this file does not find itself
PATTERNS = ("/" + "home/", "/" + "Users/", "C:" + "\\\\Users")
TEXT = (".py", ".md", ".c", ".h", ".txt", ".inp", ".pre", ".csv", ".json", ".toml", ".cfg")


def _files():
    for top in SHIPS:
        p = os.path.join(ROOT, top)
        if os.path.isfile(p):
            yield p
        for d, _, fs in os.walk(p):
            for f in fs:
                if f.endswith(TEXT):
                    yield os.path.join(d, f)


def test_nothing_that_ships_names_a_personal_folder():
    if not os.path.isdir(os.path.join(ROOT, "src")):
        pytest.skip("not a source tree (the suite runs against an installed fue)")
    hits = []
    for f in _files():
        try:
            text = open(f, encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        for n, line in enumerate(text.splitlines(), 1):
            if any(pat in line for pat in PATTERNS):
                hits.append(f"{os.path.relpath(f, ROOT)}:{n}: {line.strip()[:100]}")
    assert not hits, "personal paths in what ships:\n" + "\n".join(hits)
