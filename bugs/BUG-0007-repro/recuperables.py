"""De los .pre corruptos: cuáles se pueden recuperar y de dónde."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from barrido import bloque_deterministas  # noqa: E402

RAIZ = "/home/david/Dropbox"


def main():
    corruptos = []
    for dirpath, _d, files in os.walk(RAIZ):
        for fn in files:
            if not fn.endswith((".pre", ".inp")):
                continue
            p = os.path.join(dirpath, fn)
            r = bloque_deterministas(p)
            if r is None:
                continue
            n, tipos = r
            vivos = [t for t in tipos if t]
            if len(vivos) < n:
                corruptos.append((p, n, len(vivos)))

    por_proy, recup, perdidos, inp_rotos = {}, [], [], []
    for p, n, v in corruptos:
        rel = os.path.relpath(p, RAIZ)
        proy = "/".join(rel.split("/")[:2])
        por_proy[proy] = por_proy.get(proy, 0) + 1
        if p.endswith(".inp"):
            inp_rotos.append(rel)
            continue
        # ¿hay un .inp hermano con el bloque intacto?
        cand = p[:-4] + ".inp"
        ok = False
        if os.path.exists(cand):
            r2 = bloque_deterministas(cand)
            if r2 and len([t for t in r2[1] if t]) >= r2[0] > 0:
                ok = True
        (recup if ok else perdidos).append(rel)

    print(f"corruptos: {len(corruptos)}   "
          f"({len(recup)} con .inp hermano intacto, {len(perdidos)} sin él, "
          f"{len(inp_rotos)} son .inp ya rotos)\n")
    print("por proyecto:")
    for proy, k in sorted(por_proy.items(), key=lambda kv: -kv[1]):
        print(f"  {k:4d}  {proy}")

    print(f"\nSIN .inp hermano intacto ({len(perdidos)}):")
    for r in sorted(perdidos)[:50]:
        print(f"  {r}")
    if len(perdidos) > 50:
        print(f"  ... y {len(perdidos)-50} más")

    if inp_rotos:
        print(f"\n.inp ya rotos ({len(inp_rotos)}) — son fuente, no derivados:")
        for r in sorted(inp_rotos)[:20]:
            print(f"  {r}")


if __name__ == "__main__":
    main()
