"""Barrido del ecosistema: qué .pre/.inp tocan los deterministas afectados.

Tres poblaciones distintas, que conviene no mezclar:

  CORRUPTO   el bloque de deterministas tiene MENOS líneas de tipo que las
             declaradas -> BUG-0007 se comió una palabra al escribir el .pre.
             El fichero ya no reproduce su modelo, y no se puede saber desde él
             cuál era la variable perdida.
  compimp    el fichero está BIEN; lo que estaba mal era leerlo con fue Python
             anterior a 0.1.9. Hay que rehacer lo estimado con él.
  easter/    en riesgo: fue Python anterior a 0.1.9 ni los cargaba, y fue C
  trend      anterior al arreglo los perdía al reescribir el .pre.
  pulse      palabra que fue Python aceptaba y fue C no conoce: si ese fichero
             pasa por fue C, se lee como variable no estándar, en silencio.
"""
import os
import sys
from collections import defaultdict

CONOCIDAS = {"impulse", "compimp", "step", "ramp", "easter", "trend",
             "cos", "sin", "alter"}


def bloque_deterministas(path):
    """(n_declarado, [lineas de tipo]) o None si el fichero no tiene el bloque."""
    try:
        with open(path, encoding="latin-1") as f:
            lineas = f.read().splitlines()
    except OSError:
        return None

    i = None
    for k, ln in enumerate(lineas):
        if "number of deterministic" in ln.lower():
            i = k
            break
    if i is None:
        return None

    n = None
    j = i + 1
    while j < len(lineas):
        s = lineas[j].strip()
        if s and not s.startswith("*"):
            try:
                n = int(s.split()[0])
            except ValueError:
                return None
            break
        j += 1
    if n is None:
        return None
    if n == 0:
        return 0, []

    while j < len(lineas) and not lineas[j].strip().startswith("**"):
        j += 1
    j += 1                                    # tras el '**' de apertura
    tipos = []
    while j < len(lineas) and not lineas[j].strip().startswith("**"):
        tipos.append(lineas[j].strip())
        j += 1
    return n, tipos


def main(raiz):
    corruptos, usa = [], defaultdict(list)
    total = con_bloque = 0

    for dirpath, _dirs, files in os.walk(raiz):
        for fn in files:
            if not fn.endswith((".pre", ".inp")):
                continue
            total += 1
            p = os.path.join(dirpath, fn)
            r = bloque_deterministas(p)
            if r is None:
                continue
            con_bloque += 1
            n, tipos = r
            vivos = [t for t in tipos if t]
            if len(vivos) < n:
                corruptos.append((p, n, len(vivos)))
            for t in vivos:
                kw = t.split()[0].lower()
                if kw in ("compimp", "easter", "trend", "pulse"):
                    usa[kw].append(p)

    rel = lambda p: os.path.relpath(p, raiz)                      # noqa: E731
    print(f"ficheros .pre/.inp: {total}   con bloque de deterministas: {con_bloque}\n")

    print(f"== CORRUPTOS (faltan lineas de tipo): {len(corruptos)}")
    for p, n, v in sorted(corruptos)[:60]:
        print(f"   {rel(p)}   declara {n}, hay {v}")
    if len(corruptos) > 60:
        print(f"   ... y {len(corruptos)-60} mas")

    for kw in ("compimp", "easter", "trend", "pulse"):
        ps = usa[kw]
        print(f"\n== usan '{kw}': {len(ps)}")
        for p in sorted(ps)[:40]:
            print(f"   {rel(p)}")
        if len(ps) > 40:
            print(f"   ... y {len(ps)-40} mas")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/home/david/Dropbox")
