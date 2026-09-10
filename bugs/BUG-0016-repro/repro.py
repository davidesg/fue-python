"""BUG-0016 — `fue.__version__` estaba escrito a mano y se quedó atrás.

Determinista y sin motor: compara lo que declara el paquete con lo que dice el
`pyproject.toml`, y comprueba además que no vuelva a ser un literal.

    python bugs/BUG-0016-repro/repro.py
"""
import os
import re
import sys

RAIZ = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(RAIZ, "src"))

import fue  # noqa: E402

pyproject = open(os.path.join(RAIZ, "pyproject.toml"), encoding="utf-8").read()
esperada = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, re.M).group(1)

fallos = 0
print(f"pyproject.toml    : {esperada}")
print(f"fue.__version__   : {fue.__version__}")
if fue.__version__ != esperada:
    print("FALLO  la version esta escrita dos veces y las copias divergen")
    fallos += 1
else:
    print("OK     coinciden")

src = open(os.path.join(RAIZ, "src", "fue", "__init__.py"), encoding="utf-8").read()
if re.search(r'^__version__\s*=\s*"[\d.]', src, re.M):
    print("FALLO  __version__ es un literal: volvera a quedarse atras")
    fallos += 1
else:
    print("OK     __version__ se deriva, no se escribe")

sys.exit(1 if fallos else 0)
