# Test datasets

Real LC-MS feature tables used by the test suite (fixtures in
[`../conftest.py`](../conftest.py)). All are MassCube exports, positive mode, from
runs of the 15-standard panel (`mix15`). Each folder holds the two runs a
calibration needs: `samples/` (the table to calibrate) and `standards/` (the
standards-mixture run from the same column).

| folder | instrument | contents | used by |
|---|---|---|---|
| `Orbitrap/` | Orbitrap | `m/z` + `RT` only — the minimal input (22,314 / 8,768 features) | most end-to-end tests: calibration, CLI, errors, I/O, reference, viz |
| `QTOF/` | QTOF | `m/z` + `RT` only (22,559 / 7,220 features) | `test_calibrate.py` — extrapolation flags |
| `QTOF_full/` | QTOF | the same QTOF runs with isotope / in-source-fragment flags, MS2 and per-injection intensities, plus `samples/single_files/` — one peak list per injection (6 injections) | `test_calibrate.py`, `test_io.py`, `test_viz.py` — the per-sample tier (`single_files` → `RI_spread`) |

The second axis of every calibration — the reference column — is not here: it
is bundled with the package in `src/rt_anchor/reference_data/col35/`.

Copies elsewhere in the repository: `Orbitrap/` is byte-identical to
[`examples/mix15/`](../../../examples/mix15/) and to the package's demo data
(`src/rt_anchor/example/`); the `QTOF_full/` aligned tables are the desktop app's
**Load example** dataset (`RT Anchor Desktop/example/`).
