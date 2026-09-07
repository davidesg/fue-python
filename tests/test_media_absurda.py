"""Un óptimo con una media absurda no puede reportarse como éxito (BUG-0005).

El optimizador es una búsqueda LOCAL sobre una superficie que puede ser
multimodal, y no comprobaba si el óptimo al que llegaba tenía sentido. En el
caso del informe la rueda de Windows se quedaba en una cuenca espuria con
μ̂ = −0,144 sobre una serie cuya media diferenciada es +0,0022, y lo reportaba
como `converged=True, ifault=0` sin un solo aviso.

La no-reproducibilidad entre plataformas no se puede arreglar desde aquí —es
aritmética de coma flotante de dos compiladores— pero **que el disparate pase
callando, sí**.
"""
import os
import warnings

import numpy as np
import pytest

import fue
from fue.cast_us import build_est_spec, cast_us_py
from fue.model import UMBRAL_MEDIA_ABSURDA, _avisa_si_la_media_es_absurda


def _serie(seed=7, n=200, deriva=0.02):
    rng = np.random.default_rng(seed)
    y = 100.0 + np.cumsum(rng.standard_normal(n) * 0.3 + deriva)
    return fue.TimeSeries(y.tolist(), freq=12, start=(2005, 1), name="S")


def _modelo(ts):
    return fue.Model(ts, d=1, ar=[[0.0]], ar_free=[[True]], mu=0.0,
                     estimate_mu=True, refactor=1.0)


def _avisos(fn):
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter("always")
        fn()
    return [str(w.message) for w in W if "NO es plausible" in str(w.message)]


# ── el caso sano no molesta ───────────────────────────────────────────

def test_un_ajuste_normal_no_avisa():
    m = _modelo(_serie())
    assert _avisos(m.fit) == []


def test_el_umbral_esta_muy_por_encima_de_lo_observado():
    """Medido sobre 1.523 ajustes reales con μ identificada: máximo 1,38,
    percentil 99,9 = 0,73, ninguno por encima de 2."""
    assert UMBRAL_MEDIA_ABSURDA >= 3.0


# ── el caso espurio sí ────────────────────────────────────────────────

def test_una_media_absurda_avisa():
    ts = _serie()
    m = _modelo(ts)
    m.fit()
    _p, _q, _phi, _th, _mu, w, _f = cast_us_py(np.asarray(m._result.params, float),
                                               build_est_spec(m))
    lejos = w.mean() + 50.0 * w.std(ddof=1)
    par = np.asarray(m._result.params, float).copy()
    par[-1] = lejos                       # μ va la última
    m._result.params = par
    avisos = _avisos(lambda: _avisa_si_la_media_es_absurda(m))
    assert avisos, "un μ a 50 sd pasa en silencio"
    assert "búsqueda LOCAL" in avisos[0]
    assert "BUG-0005" in avisos[0]


def test_el_aviso_da_las_dos_cifras():
    """Sin el μ̂ y la media contra la que se compara, el aviso no se puede
    juzgar: sería una alarma sin evidencia."""
    ts = _serie()
    m = _modelo(ts)
    m.fit()
    par = np.asarray(m._result.params, float).copy()
    par[-1] = 99.0
    m._result.params = par
    a = _avisos(lambda: _avisa_si_la_media_es_absurda(m))[0]
    assert "99" in a and "desviaciones típicas" in a


# ── lo que NO se juzga, y es la mitad importante ──────────────────────

def test_una_mu_SIN_identificar_no_se_juzga():
    """El término de deriva es μ·φ(1). Con un AR de raíz unitaria escrito con
    d=0, φ(1)=0: μ no entra en la verosimilitud y puede valer cualquier cosa.

    Medido: 47 de los 1.570 modelos del ecosistema están en ese caso —los VIX
    con λ extrema y AR(1) con φ=1— y avisar de ellos serían 47 falsos
    positivos."""
    rng = np.random.default_rng(3)
    y = 100.0 + np.cumsum(rng.standard_normal(300) * 0.3)
    ts = fue.TimeSeries(y.tolist(), freq=1, start=(1990, 1), name="RW")
    m = fue.Model(ts, d=0, ar=[[1.0]], ar_free=[[False]], mu=-50.0,
                  estimate_mu=True, refactor=1.0)
    m._result = type("R", (), {"params": np.array([-50.0])})()
    assert _avisos(lambda: _avisa_si_la_media_es_absurda(m)) == []


def test_sin_mu_estimada_no_hay_nada_que_comprobar():
    ts = _serie()
    m = fue.Model(ts, d=1, ar=[[0.0]], ar_free=[[True]], mu=0.0,
                  estimate_mu=False, refactor=1.0)
    assert _avisos(m.fit) == []


# ── la comprobación no puede hacer daño ───────────────────────────────

def test_la_comprobacion_nunca_tumba_una_estimacion():
    """Una guarda de sanidad que revienta es peor que no tenerla."""
    class _Roto:
        estimate_mu = True
        _result = type("R", (), {"params": np.array([1.0, 2.0])})()
    _avisa_si_la_media_es_absurda(_Roto())      # no debe levantar


def test_no_reimplementa_la_serie_del_motor():
    """`w` la da `cast_us_py`. Reconstruirla aquí sería repetir BUG-0014, donde
    un segundo generador del mismo objeto se quedó atrás."""
    import inspect
    src = inspect.getsource(_avisa_si_la_media_es_absurda)
    assert "cast_us_py" in src
    assert "np.diff" not in src and "boxcox" not in src.lower()


# ── el caso del informe, medido (BUG-0005, revisión 2026-09-07) ───────

_US_CPI = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "..", "..", "..", "ART", "art-python", "bugs",
                       "BUG-0006-repro", "US_CPI.pre")


@pytest.mark.skipif(not os.path.exists(_US_CPI), reason="sin el .pre del repro")
def test_el_caso_del_informe_es_robusto_a_la_semilla():
    """La causa raíz del informe —«superficie multimodal cuya cuenca decide el
    compilador»— no la sostiene la evidencia. Lo que la sostiene es la línea que
    el propio informe dejó al margen: una semilla con el SIGNO CAMBIADO metía la
    búsqueda en la región equivocada.

    La convención de `fue` es la de Box y Jenkins para todo operador
    —ω(B) = ω₀ − ω₁B − …— y sembrar con el signo al revés no es empezar «algo
    peor»: es empezar en otro sitio. Arreglado eso (art/BUG-0006), el caso no se
    mueve desde ningún arranque.
    """
    import fue as _fue
    logliks = set()
    for a, b in [(-0.109, -0.0935), (0.109, 0.0935), (0.5, 0.3), (-0.5, -0.3),
                 (0.9, 0.0), (0.0, 0.0)]:
        ld = _fue.load(_US_CPI)
        m = ld[1] if isinstance(ld, tuple) else ld
        m.ar_s = [[a, b]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit()
        logliks.add(round(m._result.loglik, 2))
    assert len(logliks) == 1, f"varias cuencas: {logliks}"


@pytest.mark.skipif(not os.path.exists(_US_CPI), reason="sin el .pre del repro")
def test_el_punto_espurio_no_es_un_optimo():
    """Segunda prueba de que no es multimodalidad: a lo largo de μ el objetivo
    decrece monótonamente desde el óptimo y su derivada en −0,144 vale +19,6.
    Un punto con gradiente no nulo NO es una cuenca."""
    import fue as _fue
    from fue.elfvarma import elf_scalar
    ld = _fue.load(_US_CPI)
    ts, m = (ld if isinstance(ld, tuple) else (ld.series, ld))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit()
    spec = build_est_spec(m)
    par = np.asarray(m._result.params, float)

    def obj(x):
        p, q, phi, th, mu, w, f = cast_us_py(np.asarray(x, float), spec)
        r = elf_scalar(len(w), p, q, phi, th, w, 1.0, mu)
        return float(r[0] if isinstance(r, (tuple, list)) else r)

    base = obj(par)
    previo = base
    for mu in (-0.005, -0.05, -0.144, -0.30):
        x = par.copy(); x[-1] = mu
        v = obj(x)
        assert v < previo, "el objetivo no decrece: habría otra cima"
        previo = v
    h = 1e-6
    x1, x2 = par.copy(), par.copy()
    x1[-1], x2[-1] = -0.144 + h, -0.144 - h
    g = (obj(x1) - obj(x2)) / (2 * h)
    assert abs(g) > 1.0, f"gradiente {g}: sería un punto estacionario"
