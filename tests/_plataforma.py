"""Where the suite can pin a trajectory, and where it cannot.

Some tests compare against what the C engine did on Linux: the bytes of a
.out, an iteration count, a termcode. Those are facts about the optimiser's
TRAJECTORY, and the trajectory moves with the libm's last bits. Measured
(fue 0.1.19, the preflight and the podman sandbox, Nov 2026):

  * the C engine on Linux x86_64 with glibc (AlmaLinux 8's 2.28 through
    Fedora's): identical, every test passes;
  * the C engine with musl (Alpine), on macOS arm64 and on Windows: RIPC.1
    takes 118 iterations instead of 122 and stops 2e-4 away;
  * the pure-Python engine (no wheel with the extension): another optimiser,
    and on a flat likelihood possibly another optimum (BUG-0005).

Off the first case those tests are skipped, not loosened: a byte golden with
a tolerance is not a byte golden. What does not depend on the trajectory
(the estimates to their tolerance, the report's structure) runs everywhere.
"""
import platform
import sys

import pytest

import fue

C_ENGINE = fue.engine_backend() == "c"

EXACT = (C_ENGINE and sys.platform.startswith("linux")
         and platform.machine() in ("x86_64", "AMD64")
         and platform.libc_ver()[0] == "glibc")

exact_trajectory = pytest.mark.skipif(
    not EXACT, reason="pins the C engine's trajectory on Linux x86_64 with "
                      "glibc; another libm or engine takes another path")

c_engine = pytest.mark.skipif(
    not C_ENGINE, reason="pinned at the C engine's precision; the pure-Python "
                         "engine may stop elsewhere on a flat likelihood (BUG-0005)")
