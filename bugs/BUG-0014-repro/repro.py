"""BUG-0014 — la previsión de Python daba efecto NULO a tres deterministas.

`forecast._build_xi` reimplementaba la construcción del indicador —indexada por
`type_code` en vez de por nombre— con ramas para 8 de los 11 tipos y sin `else`
final. `compimp` (8), `easter` (9) y `trend` (10) caían fuera y su indicador
salía idénticamente nulo, sin aviso.

    python bugs/BUG-0014-repro/repro.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

import numpy as np  # noqa: E402

import fue  # noqa: E402
from fue.forecast import _build_xi  # noqa: E402

ts = fue.TimeSeries([100.0] * 160, freq=12, start=(2005, 1), name="X")
fallos = 0

for tipo in ("step", "impulse", "compimp", "easter", "trend", "alter"):
    at = 0 if tipo in ("easter", "trend", "alter") else 100
    m = fue.Model(ts, d=1, ar=[[0.0]], ar_free=[[False]], mu=0.0,
                  estimate_mu=False, refactor=100.0,
                  interventions=[fue.Intervention(tipo, at=at, omega=[1.0],
                                                  omega_free=[False])])
    xi = _build_xi(m, 160, 12, 12, [[1.0]], [[]])
    total = float(np.abs(xi).sum())
    if total == 0.0:
        print(f"FALLO  {tipo:9} regresor IDÉNTICAMENTE NULO en la previsión")
        fallos += 1
    else:
        print(f"OK     {tipo:9} Σ|xi| = {total:.3f}")

# Y un tipo que no existe no puede caer en silencio en NINGUNO de los dos
# generadores: devolver ceros es prever un modelo distinto del pedido.
from fue.cast_us import _build_indicator  # noqa: E402


class _Inventado:
    type = "no_existe"
    at = 0
    data = None
    harmonic = 1


try:
    _build_indicator(_Inventado(), 60, 12, 1, 2005)
    print("FALLO  un tipo desconocido devuelve ceros en silencio")
    fallos += 1
except ValueError:
    print("OK     un tipo desconocido se rechaza en vez de valer cero")

sys.exit(1 if fallos else 0)
