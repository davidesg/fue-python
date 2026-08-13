#!/usr/bin/env bash
#
# Rebuild DRVUS from Mauricio's source and reproduce the archived reference runs.
#
# Written while resolving bugs/BUG-0012, where fue disagreed with the `.out` that
# DRVUS left on Box-Jenkins Series A and the port was, reasonably, under
# suspicion. It was not the port: the same source rebuilt today at 64 bits gives
# fue's answer, and only a 32-bit unoptimised build — 80-bit x87 intermediates —
# gives the archived one.
#
# Two builds, and the difference between them IS the point:
#
#   64-bit (any -O)    → a1: PARAMETER stop, 23 iterations, logelf -57.6038617504
#   -m32 -O0           → a1: GRADIENT stop,  64 iterations, logelf -50.7450915148
#
# The second reproduces the 2001/2006 archive iteration by iteration. -O1 is
# already enough to lose it: the optimiser spills the 80-bit registers to 64-bit
# memory.
#
# Usage:  tools/reproduce_drvus_reference.sh [version] [case]
#         tools/reproduce_drvus_reference.sh 1.2.01 a1
#
set -euo pipefail

VERSION="${1:-1.2.01}"
CASE="${2:-a1}"
SRC="${DRVUS_SOURCE:-$HOME/Dropbox/SRC/drvus-source}/$VERSION/drvus/src"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

[ -d "$SRC" ] || { echo "DRVUS $VERSION not found at $SRC"; exit 1; }

cp "$SRC"/*.c "$SRC"/*.h "$WORK/" 2>/dev/null

# The only edits to Mauricio's source, both for things the platform changed:
#   round()   collided with C99's math.h in 2001 it did not
#   \x1a      the DOS end-of-file mark, in the oldest trees
#   strcmpi   Borland, absent from glibc
python3 - "$WORK" <<'PY'
import glob, os, re, sys
work = sys.argv[1]
for p in glob.glob(os.path.join(work, "*.c")) + glob.glob(os.path.join(work, "*.h")):
    raw = open(p, "rb").read().replace(b"\x1a", b"")
    txt = re.sub(r"\bround\s*\(", "drvus_round(", raw.decode("latin-1"))
    open(p, "wb").write(txt.encode("latin-1"))
PY

cat > "$WORK/shim.c" <<'EOF'
/* itoa() was Borland's; glibc has no such function. Nothing here computes. */
#include <stdio.h>
char *itoa( int value, char *str, int base )
{
    if ( base == 16 ) sprintf( str, "%x", value ); else sprintf( str, "%d", value );
    return str;
}
EOF

cd "$WORK"
SOURCES=$(ls *.c | tr '\n' ' ')

for flags in "-O2" "-m32 -O0"; do
    out="drvus_$(echo "$flags" | tr -d ' -')"
    if ! gcc $flags -Dstrcmpi=strcasecmp -o "$out" $SOURCES -lm 2>/dev/null; then
        echo "  $flags: does not build here (32-bit multilib missing?)"
        continue
    fi
    cp "$SRC/Box_y_Jenkings"/*/"$CASE.inp" . 2>/dev/null || \
        { echo "  case $CASE not found"; exit 1; }
    rm -f "$CASE.out"
    "./$out" "$CASE" eml chk >/dev/null 2>&1 || true
    python3 - "$CASE.out" "$flags" <<'PY'
import re, sys
t = open(sys.argv[1], encoding="latin-1", errors="replace").read()
m = re.search(r"AFTER (\d+) ITERATIONS \[GRADIENT NORM = *([\d.]+)\]", t)
l = re.search(r"logelf:\s*([-\d.]+)", t)
c = ("GRADIENT" if "GRADIENT STOPPING" in t else
     "STEP" if "PARAMETER STOP" in t else "?")
print(f"  {sys.argv[2]:10}  {c:9} iter={m.group(1) if m else '?':>4}  "
      f"logelf={l.group(1) if l else '?'}")
PY
done

echo
echo "The archived reference for comparison:"
grep -h "logelf:" "$SRC/Box_y_Jenkings"/*/"$CASE.out" 2>/dev/null | head -1
