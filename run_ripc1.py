#!/usr/bin/env python3
"""
RIPC.1 con un factor AR(2) adicional.
"""
import fue
from fue import Model

INP = "tests/real_cases/PRICES/IPC/Mensual/sample_1.2002_12.2007/RIPC.2.inp"

ts, m = fue.load(INP)


m.fit()
print(m.write_out())
