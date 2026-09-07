# rt-anchor

**Cross-column retention-time calibration** for LC-MS lipidomics — plus a
portable retention index (iRT).

## Purpose

Retention time is a property of the column and the gradient, not of the
molecule. The same lipid elutes at 4.5 min on a 15-minute method and at 17 min
on a 90-minute one, so RT cannot be compared across methods and RT-based
annotation does not transfer between labs.

rt-anchor puts every feature of your run onto a **single shared time axis**: the
retention time it *would have had* on a reference column. From a run of your
standards mixture and the same mixture on the reference column, it learns a
monotone map `f: RT_yours → RT_reference`, then applies it to every feature of
your biological sample. Two outputs come out of it:

* **`Cal_RT_min`** — your feature's RT on the reference column's minute scale.
* **`iRT`** — the same position as a dimensionless 1–100 index, pinned to the
  standards detected on the reference column.

Nothing in your table is dropped or reordered; calibration only appends columns.

### How it works

**Stage 1 — the anchor-free curve.** Every feature of your standards run is
matched to the reference standards run by accurate m/z alone (reciprocal best
match, 0.008 Da), and a robust monotone curve is fitted through the pairs
(LOESS → isotonic → PCHIP, with MAD outlier trimming). m/z-matched features from
the two *sample* runs are merged in, because the biological sample covers the
very early and very late elution regions where the mixture is sparse.

No standard identities are used at this stage. That is deliberate: the
calibration still works when your mixture is unknown, or when nothing in the
panel is detectable at all.

**Stage 2 — sample anchors, gated.** Seventeen high-abundance endogenous plasma
lipids (LPC/PC/SM/CE/TG) are located in both sample runs, validated by adduct
consistency, within-class elution order, and a curve-assisted isomer re-pick.
Their residuals against the Stage-1 curve are class-systematic on long
gradients, so the correction is decomposed into a class offset plus a
class-detrended piecewise-linear term, both shrunk by leave-one-anchor-out
cross-validation.

The correction is then **gated**: it is applied only if it reduces
leave-one-anchor-out error by at least 20%. Otherwise the pure Stage-1 curve is
used. On a non-plasma matrix, too few anchors validate and Stage 2 simply does
not engage — the run succeeds on the curve alone, which is the intended
behaviour, not a failure.

**iRT.** The chosen panel's standards are located in the *reference* standards
run; the earliest and latest set the 1 and 100 ends of the scale. Because the
scale is spaced linearly in RT, it is an **affine** function of `Cal_RT_min` —
the interior landmarks certify that the panel was found, they do not bend the
scale. iRT therefore depends only on the reference column and the panel, never
on whether your own run contains the standards.

## Install

```bash
pip install rt-anchor
pip install "rt-anchor[report,structures]"   # + HTML/PDF report and structure hover
```

## Input

| Input | What | Required |
|-------|------|----------|
| **Sample feature table** | The biological sample to calibrate. MS-DIAL / MZmine / MassCube / LipidScreener — auto-detected. Must carry a **precursor m/z** and a **retention-time** column. | yes |
| **Standards run** | The standards mixture, **from the same column and gradient as the sample**. Same formats. | yes |
| **Panel** | `mix15` (15-standard Caley mix), `mix21` (Mix 4.4, default), or `none`. | yes |
| **Polarity** | `positive` or `negative`. | yes |
| **Reference pair** | The reference sample + standards runs. Defaults to the bundled **35-minute column** dataset. | no |

The same-column requirement is the one assumption the method rests on: the
standards run is what tells us how *your* time axis maps onto the reference
one, and the sample inherits that map.

`panel="none"` means "my mixture is not one of the known panels". Stage 1 needs
no identities, so calibration proceeds normally; the iRT ruler falls back to
`irt_landmark_panel` on the reference run, and the model records that it did.

Optional: per-injection files (per-sample tier + `RI_spread` repeatability QC), a
custom manifest, and the matching parameters (`match_mz_tol_da`, `mz_tol_ppm`,
`curve_frac`, `min_anchors`, `extrapolate`).

```python
from rt_anchor import calibrate, write_results
res = calibrate("samples.csv", polarity="positive",
                standards_table="standards.csv", panel="mix21")
write_results(res, "out/run1")
```
```bash
rt-anchor calibrate --samples samples.csv --standards standards.csv \
                    --polarity positive --panel mix21 --out out/run1
rt-anchor references          # list the bundled reference datasets
rt-anchor describe table.csv  # detect format + summarise, no calibration
```

### Using your own reference column

```bash
rt-anchor calibrate ... --reference-sample ref_serum.csv \
                        --reference-standards ref_stds.csv
```

Both must be given together. Outputs then live on *that* column's axis and are
no longer comparable with default-reference runs — the model JSON records which
reference was used, so a result can always be traced back to its axis.

## Output

**`*_calibrated.csv`** — your original table with these columns appended:

| Column | Meaning |
|---|---|
| `Cal_RT_min` | RT on the reference column's minute scale |
| `Cal_RT_uncertainty_min` | 1σ, from the local scatter of the matched pairs about the curve, combined with the anchor LOO error |
| `iRT` | dimensionless 1–100 retention index |
| `iRT_uncertainty` | 1σ of `iRT` |
| `iRT_reliability` | `high` / `medium` / `low` / `none` |
| `is_extrapolated` | RT outside the span the matched pairs cover |
| `calibration_scope` | `project` or `sample` |
| `warp_source` | `curve` or `curve+anchors` |
| `RI_spread`, `n_contributing` | per-injection dispersion (per-sample tier only) |

Companion files: **`*_model.json`** (curve, anchors, gate decision and its
reason, iRT definition, reference provenance, detection QC),
**`*_anchors.csv`** (the Stage-2 anchors and their residuals),
**`*_pairs.csv`** (every matched pair and whether it survived trimming),
**`*_landmarks.csv`** (the panel standards located on the reference run),
**`*_log.txt`**, and **`*_report.html`** + **`*_report.pdf`**.

## Accuracy

Against the five-column validation set (15 / 22 / 25 / 45 / 90-minute methods,
human serum, QTOF), calibrating onto the 35-minute reference column:

| source method | median abs. error before | Stage 1 only | + Stage 2 |
|---|---|---|---|
| 15 min | 9.44 min | 0.10 min | **0.09 min** |
| 22 min | 1.72 min | 0.11 min | **0.11 min** |
| 45 min | 17.89 min | 0.17 min | **0.15 min** |
| 90 min | 39.77 min | 0.19 min | **0.19 min** — gate off |

Measured on held-out serum features that were never used to fit the curve. The
90-minute column is where the gate earns its keep: its anchor residuals are
class-incoherent, the correction would not have generalised, and the gate
declines it, leaving the Stage-1 result untouched.
Calibrating the reference column against itself returns the input RT to within
1e-11 min.

Two honest caveats. `Cal_RT_min` can show hairline non-monotonicity (observed
worst case 0.033 min, an order of magnitude below the method's own accuracy)
where the Stage-2 correction is active, because curve + piecewise-linear
correction is not monotone by construction. And Stage 2's anchor panel is
human plasma/serum: other matrices get the Stage-1 curve only.

## Desktop client

A native-window client (pywebview) wraps the engine with a sectioned UI —
Overview · Detection · Profile · Warp · Anchors · Repeatability · Table ·
Export. See `RT Anchor Desktop/` for the macOS and Windows builds.

```bash
pip install "rt-anchor[app,report]"
rt-anchor-gui          # or:  python -m rt_anchor.gui
```

## License

MIT
