"""BUG-0006 — deterministic types that this port did not build like fue C.

Three of fue C's nine deterministic regressors were wrong or missing here, and
the worst of them failed **in silence**: `compimp` was read as a plain `impulse`,
so the model estimated cleanly — just not the model that was asked for.

The arbiter is fue C itself. Each test states the indicator fue C builds
(`fue.c:295-435`) and checks this package builds the same one. Where the fue C
binary is available the check is made end to end: it fits, rewrites the `.pre`
with its estimates, and this package is evaluated **at that same point** — if the
regressor matches, the two likelihoods match.

The fourth finding is about *vocabulary*, not arithmetic: fue C does not reject a
keyword it does not know, it takes it for a non-standard variable whose data
comes as extra columns. So any word the two programs do not share is a silent
misreading waiting to happen, and `.pre`/`.inp` files are meant to be shared
between the two interpreters.
"""

import math
import os
import subprocess

import numpy as np
import pytest

from fue.cast_us import _build_indicator, easter_date, obs_to_date
from fue.intervention import Intervention
from fue.report import _itv_name_line

FUEC = "/home/david/Dropbox/SRC/atws/fue/fue-1.13.1/bin/fue"


def _ind(type, nobs=12, freq=12, begtime=1, begyear=2002, at=0, **kw):
    return _build_indicator(Intervention(type, at=at, **kw),
                            nobs, freq, begtime, begyear)


# ── compimp: el fallo silencioso ─────────────────────────────────────────────
def test_compimp_es_un_impulso_COMPENSADO_no_un_impulso():
    """+1 en la fecha y −1 en la SIGUIENTE (fue.c:317).

    Es la diferencia que importa: sobre una serie diferenciada un impulso
    desplaza el nivel y un impulso compensado no. Leer uno por el otro estima
    otro modelo sin decirlo.
    """
    ind = _ind("compimp", at=4)
    assert ind[5] == 1.0 and ind[6] == -1.0
    assert ind.sum() == 0.0, "un impulso compensado suma cero; un pulso, uno"

    pulso = _ind("impulse", at=4)
    assert pulso[5] == 1.0 and pulso[6] == 0.0
    assert not np.array_equal(ind, pulso)


def test_compimp_al_final_de_la_muestra_no_se_sale():
    ind = _ind("compimp", nobs=6, at=5)          # última observación
    assert ind[6] == 1.0 and len(ind) == 7


def test_pulse_sigue_aceptandose_pero_se_normaliza_a_impulse():
    """`impulse` es el nombre de la escuela y el del formato; `pulse` era el de
    este paquete. Se aceptan los dos y se guarda UNO, para que no exista un
    fichero con una palabra que el otro intérprete no conoce."""
    itv = Intervention("pulse", at=3)
    assert itv.type == "impulse"
    assert itv.type_code == Intervention.TYPES["impulse"] == 0
    assert _itv_name_line(itv, 2002, 1, 12).startswith("impulse ")


# ── easter y trend: los que faltaban ─────────────────────────────────────────
def test_easter_es_el_algoritmo_de_fue_c():
    """`Easter` (nlatools.c:693), portado literalmente: el indicador tiene que
    ser el que construye fue C, no uno más correcto."""
    # comprobado contra el propio Easter() de fue C, compilado, 2002-2020
    assert easter_date(2002) == (31, 3)
    assert easter_date(2005) == (27, 3)
    assert easter_date(2008) == (23, 3)
    assert easter_date(2011) == (24, 4)
    assert easter_date(2018) == (1, 4)
    for y in range(1990, 2031):
        d, m = easter_date(y)
        assert m in (3, 4) and 1 <= d <= 31


def test_easter_reparte_con_marzo_cuando_cae_a_primeros_de_abril():
    """Si el domingo de Pascua cae en los tres primeros días de abril, el peso se
    reparte 0.5/0.5 con marzo: la semana que mueve la actividad cayó en marzo."""
    # 2013: Pascua el 31 de marzo -> todo el peso en marzo
    ind = _ind("easter", nobs=24, begyear=2013)
    assert easter_date(2013) == (31, 3)
    assert ind[3] == 1.0 and ind[4] == 0.0

    # 2002: Pascua el 31 de marzo tambien; 2018 cae el 1 de abril -> mitades
    assert easter_date(2018) == (1, 4)
    ind = _ind("easter", nobs=24, begyear=2018)
    assert ind[4] == 0.5 and ind[3] == 0.5


def test_easter_solo_es_mensual():
    """fue C exige freq == 12 (fue.c:374); con otra frecuencia no hay variable."""
    assert _ind("easter", nobs=12, freq=4, begyear=2002).sum() == 0.0


def test_trend_es_1_2_3():
    ind = _ind("trend", nobs=5)
    assert list(ind[1:]) == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_obs_to_date_es_el_del_c():
    assert obs_to_date(2002, 1, 1, 12) == (2002, 1)
    assert obs_to_date(2002, 1, 12, 12) == (2002, 12)
    assert obs_to_date(2002, 1, 13, 12) == (2003, 1)
    assert obs_to_date(2002, 6, 8, 12) == (2003, 1)      # empieza en junio


# ── el vocabulario compartido ────────────────────────────────────────────────
def test_no_se_escribe_una_palabra_que_fue_c_no_conoce():
    """fue C no da error ante una palabra desconocida: la toma por variable NO
    ESTÁNDAR y estima otra cosa. Así que escribir `seasonal` —que no existe en
    fue C, porque la estacionalidad determinista va con armónicos— produciría un
    fichero que miente. Mejor negarse."""
    itv = Intervention("seasonal", at=2)
    with pytest.raises(ValueError, match="no .pre/.inp representation"):
        _itv_name_line(itv, 2002, 1, 12)


@pytest.mark.parametrize("tipo", ["impulse", "compimp", "step", "ramp"])
def test_los_fechados_se_escriben_con_su_palabra_y_su_fecha(tipo):
    linea = _itv_name_line(Intervention(tipo, at=5), 2002, 1, 12)
    assert linea == f"{tipo} 6 2002"


@pytest.mark.parametrize("tipo", ["easter", "trend", "alter"])
def test_los_sin_fecha_se_escriben_solos(tipo):
    assert _itv_name_line(Intervention(tipo), 2002, 1, 12) == tipo


# ── el árbitro: fue C ────────────────────────────────────────────────────────
PLANTILLA = """************************************************
*        Input file for program DRVUS          *
*   BUG-0006: {tipo}
************************************************

** Frequency of time series: either 1(A), 4(Q) or 12(M):
 12
** Number of observations and starting date of time series:
 {nobs}  1 2002 TEST
** Number of deterministic variables (including seasonal components):
1
**
{linea}
**
0
**
1.000000  1
**
0
**Number and orders of regular AR operators:
1 1
**
0.4000 1
** Number and orders of annual AR operators:
0
** Number and orders of regular MA operators:
0
** Number and orders of anual MA operators:
0
** Number and frequencies of regular AR(2) operators with fixed frequency:
0
** Number and frequencies of regular MA(2) operators with fixed frequency:
0
** Mean parameter (mu):
0
** Box-Cox lambda, regular differences and complete annual differences:
1.00 1 0
** Individual factors of the annual difference (from freq 0.0):
 0 0 0
** ACF/PACF bands (0 Automatic) and reescaling factor:
 0.00 1.00
** Time series (stochastic and non-standard deterministic variables):
{datos}
"""


def _ll_en_semillas(path):
    """Verosimilitud exacta concentrada en las semillas del `.pre` (sin ajustar)."""
    import fue
    from drvarma.estimate_py import _elf_f1f2
    from fue.cast_us import _build_initial_x, build_est_spec, cast_us_py

    _ts, mod = fue.load(path)
    p, q, phi, theta, mu, w, _ifa = cast_us_py(_build_initial_x(mod),
                                               build_est_spec(mod))
    w = np.asarray(w, float).reshape(-1, 1)
    n = len(w)
    ph = np.asarray(phi, float).reshape(-1, 1, 1) if p else np.zeros((0, 1, 1))
    th = np.asarray(theta, float).reshape(-1, 1, 1) if q else np.zeros((0, 1, 1))
    f1, f2, _ = _elf_f1f2(w, np.array([float(mu)]), ph, th, np.ones((1, 1)), -1e-3)
    return (-0.5 * n * (math.log(2 * math.pi) - math.log(n) + 1.0)
            - 0.5 * n * (math.log(f1) + math.log(f2)))


@pytest.mark.skipif(not os.path.exists(FUEC), reason="falta el binario de fue C")
@pytest.mark.parametrize("tipo,linea", [
    ("impulse", "impulse 6 2005"),
    ("compimp", "compimp 6 2005"),
    ("step",    "step 6 2005"),
    ("ramp",    "ramp 6 2005"),
    ("easter",  "easter"),
    ("trend",   "trend"),
    ("cos",     "cos 1"),
    ("sin",     "sin 1"),
    ("alter",   "alter"),
])
def test_los_nueve_deterministas_homologan_con_fue_c(tipo, linea, tmp_path):
    """fue C ajusta y reescribe el `.pre`; este paquete se evalúa EN ESE PUNTO.

    Comparar en un punto DADO es lo que aísla el regresor: si los dos construyen
    el mismo indicador, las dos verosimilitudes coinciden, sin que intervenga el
    optimizador.

    De paso guarda el `.pre` de fue C, que hasta el arreglo de BUG-0007 perdía la
    palabra de `easter` y `trend` y no reproducía su propio modelo. Si estas dos
    fallan y las otras siete pasan, el binario de fue C es anterior al arreglo.
    """
    rng = np.random.default_rng(20260730)
    datos = "\n".join(f"{v:.10f} " for v in 100.0 + np.cumsum(rng.normal(0, 1.0, 216)))
    base = tmp_path / f"T_{tipo}"
    base.with_suffix(".inp").write_text(
        PLANTILLA.format(tipo=tipo, nobs=216, linea=linea, datos=datos))

    r = subprocess.run([FUEC, f"T_{tipo}", "eml"], cwd=str(tmp_path),
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-400:]

    ll_c = None
    for ln in base.with_suffix(".out").read_text().splitlines():
        if ln.lower().startswith("logelf"):
            ll_c = float(ln.split(":")[1])
            break
    assert ll_c is not None, "fue C no reportó logelf"

    pre = base.with_suffix(".pre").read_text()
    assert f"\n{tipo}" in pre, (
        f"el .pre que escribió fue C no contiene la palabra {tipo!r}: "
        "binario anterior al arreglo de BUG-0007")
    assert _ll_en_semillas(str(base.with_suffix(".pre"))) == pytest.approx(
        ll_c, abs=1e-4)
