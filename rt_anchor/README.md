# rt-anchor

Standard-panel **retention-index (iRT) calibration** for LC-MS lipidomics.

## Purpose

Convert the retention time of every feature in an LC-MS table into a
dimensionless, portable **retention index (RI/iRT)** anchored on a spiked
standard panel — so retention is comparable across injections, batches, and
instruments. A monotone warp is fitted from the standards and applied to all
features; nothing in the original table is dropped, calibration only appends
columns.

## Download

```bash
pip install rt-anchor
pip install "rt-anchor[report,structures]"   # + HTML/PDF report and structure hover
```

## Input

Required:

| Input | What |
|-------|------|
| **Sample feature table** | The table to calibrate. MS-DIAL / MZmine / MassCube / LipidScreener — auto-detected. Must have an **m/z** and a **retention-time** column. |
| **Standards run** | A run of the standard panel (same formats) — builds the native anchor template. **Required.** |
| **Polarity** | `positive` or `negative`. |

Optional: per-injection files (per-sample tier + repeatability QC), a custom
standard manifest / reference baseline, and the matching parameters
(`mz_tol_ppm`, `rt_window_min`, `min_anchors`, `extrapolate`).

```python
from rt_anchor import calibrate, write_results
res = calibrate("samples.txt", polarity="positive", standards_table="standards.txt")
write_results(res, "out/run1")
```
```bash
rt-anchor calibrate --samples samples.txt --standards standards.txt --polarity positive --out out/run1
```

## Output

- **`*_calibrated.csv`** — the original table + appended columns: `RI` (iRT-style
  index), `RI_uncertainty`, `RI_reliability` (high/medium/low/none), `RI_spread`,
  `n_contributing`, `is_extrapolated`, `calibration_scope`, `warp_source`.
- **`*_model.json`** — the fitted warp, anchors, and QC metrics.
- **`*_anchors.csv`**, **`*_log.txt`** — anchors used and a run log.
- **`*_report.html`** + **`*_report.pdf`** — an interactive + static report
  (detection, calibration warp, repeatability). On by default; `--no-report` to skip.

## License

MIT
