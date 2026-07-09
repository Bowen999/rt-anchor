# RT Anchor

Standard-panel **retention-index (iRT) calibration** for LC-MS lipidomics — a
Python package (`rt-anchor`) and a macOS desktop app.

## Purpose

Convert the retention time of every feature in an LC-MS table into a
dimensionless, portable **retention index (RI/iRT)** anchored on a spiked
standard panel, so retention is comparable across injections, batches, and
instruments. A monotone warp is fitted from the standards and applied to all
features; the original table is never altered — calibration only appends columns.

## Download

- **Python package:** `pip install rt-anchor` (`pip install "rt-anchor[report,structures]"` for the HTML/PDF report + structure hover)
- **macOS desktop app:** download `RT Anchor.dmg` from the Releases page, or build it from [`RT Anchor Desktop/`](RT%20Anchor%20Desktop/).

## Input

Required: a **sample feature table** (MS-DIAL / MZmine / MassCube / LipidScreener,
auto-detected; needs an **m/z** and a **retention-time** column), a **standards
run** of the panel, and the **polarity** (`positive`/`negative`).

Optional: per-injection files (per-sample tier + repeatability QC), a custom
manifest / reference baseline, and matching parameters (`mz_tol_ppm`,
`rt_window_min`, `min_anchors`, `extrapolate`).

## Output

- **`*_calibrated.csv`** — original table + `RI`, `RI_uncertainty`, `RI_reliability`,
  `RI_spread`, `n_contributing`, `is_extrapolated`, `calibration_scope`, `warp_source`.
- **`*_model.json`**, **`*_anchors.csv`**, **`*_log.txt`** — the fitted model, anchors, and log.
- **`*_report.html`** + **`*_report.pdf`** — interactive + static report (on by default).

## Repository layout

- [`rt_anchor/`](rt_anchor/) — the Python package (PyPI: **rt-anchor**).
- [`RT Anchor Desktop/`](RT%20Anchor%20Desktop/) — the macOS desktop app (pywebview + PyInstaller).

## License

MIT
