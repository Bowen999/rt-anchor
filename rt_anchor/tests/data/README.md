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

`expected_series_repo_datasets.json` holds the reference implementation's
stage-1b (homologous-series term) numbers for the three datasets above —
series counts, gate statistics, and per-row raw corrections — and is asserted
against in [`../test_series.py`](../test_series.py).

`expected_v121_orbitrap_golden.json` holds the numeric goldens for the two
v1.2.1-equivalent configurations (series term off, anchors on / off) on the
Orbitrap dataset: the sampled result columns, value counts of the text
columns, the anchors table's numbers, the pair counts, and the `curve` /
`anchors` / `irt` model blocks. Also asserted in `test_series.py` (to 1e-8;
generated from this tree with the series term off, the configuration the
reviewer verified reproduces v1.2.1 exactly).

Copies elsewhere in the repository: `Orbitrap/` is byte-identical to
[`examples/mix15/`](../../../examples/mix15/), and the `QTOF_full/` aligned
tables are the desktop app's **Load example** dataset (`RT Anchor Desktop/example/`).
