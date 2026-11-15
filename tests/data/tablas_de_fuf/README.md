# The forecast table of fuf C, as the port must write it (BUG-0027)

Each case is a fuf input (`<id>.inp`) and the forecast table that atsw-gui's
fuf writes for it (`<id>.tabla`): from `FORECAST REPORT:` to the closing
border, cut from the `.out`. fuf fixed its table on 15-nov-2026 (its
BUG-0001, BUG-0002 and BUG-0003, `engines/fuf/BUGS.md`, atsw-gui 9636562).

- The fourteen cases of `engines/fuf/tests/golden`, from the inputs of
  `engines/fuf/tests/corpus`. Six of them (`CPI_USA_model.fuf`,
  `CPI_USA_model_AR.fuf`, `IPC_DE_model.fuf`, `IPC_DE_model_AR.fuf`,
  `IPC_ES_m01_ar1.fuf`, `IPC_FR_model.fuf`) carry a rescaling factor of 1
  with parameters estimated at 100: the forecasts explode, and the port and
  the C part in the last digit (at most 4e-5 relative) and in the sign of a
  zero. fuf keeps them as its byte baseline, not as models to believe.
- `CPI_USA_model_AR_refactor100` and `IPC_ES_m01_ar1_refactor100`: the
  inputs of fuf's `bugs/BUG-0001-repro`, written by fue with the right
  factor; one in levels (lambda = 1), one in logarithms.
- `MINK_m01_anual`: Jenkins and Alavi's mink (atsw-gui
  `docs/ejemplos/jenkins_alavi/MINK/work/MINK_m01.pre`), an annual series;
  the input written by `fue.report.write_fuf` with 10 periods.

The tables were written by fuf built from atsw-gui 9636562 on Linux x86_64
(gcc, glibc), which passes its own battery (`tests/run_tests.sh`: 15 runs,
0 failures).
