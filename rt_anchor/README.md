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
match within 15 ppm), and a robust monotone curve is fitted through the pairs
(LOESS → isotonic → PCHIP, with MAD outlier trimming). m/z-matched features from
the two *sample* runs are merged in, because the biological sample covers the
very early and very late elution regions where the mixture is sparse.

No standard identities are used at this stage. That is deliberate: the
calibration still works when your mixture is unknown, or when nothing in the
panel is detectable at all.

**Stage 1b — the homologous-series term.** Two methods can differ in how
strongly they retain a double bond relative to a CH2 group, and the error the
curve leaves behind is then systematic within a *homologous series* (same
class, same unsaturation, different chain length): two lipids that co-elute on
your column can sit more than a minute apart on the reference one, and no
single curve can place both. Members of a series differ by CH2, so they are
recognised from *m/z* alone through the Kendrick mass — a shared position on
the 14-Da mass-defect circle, **not an identification**. Each feature that
belongs to a validated series is corrected with the curve residuals of the
*other* matched pairs in its series, interpolated along the series' elution
order. Three safeguards keep it honest:

* **Carbon-order validation on both columns** — a cluster of same-phase pairs
  becomes a series only if it elutes in carbon-number order on your column
  *and* on the reference column; a cluster whose order disagrees is ignored.
* **Leave-own-family-out** — a feature is never corrected with pairs that
  co-elute with it, so a pair's own residual (or a mis-assigned co-eluter's)
  cannot leak into its own correction.
* **A 20% gate** — the term is applied only if predicting every matched pair
  this way cuts the pairs' mean squared error by at least 20% below the curve
  alone, on at least 20 covered pairs; otherwise it leaves the stage-1 curve
  untouched.

A feature outside every validated series is not corrected by this term, and on
a matrix with few homologous series the gate will simply not engage; the
lattice term, next, fills in some of what it leaves.

**Stage 1c — the lattice term (since 1.2.3).** The series term needs at least
four members in one CH2 ladder (same class, same number of double bonds), and
some classes never offer that: serum cholesteryl esters have at most three per
double-bond count, so the series term never touches them and their class-wide
offset against the curve stays in the calibrated RT. Yet the class as a whole
has a dozen members, spread over *several* ladders that differ by H2 — one
double bond. The lattice term uses that. CH2 ladders whose Kendrick phases are
one H2 step apart (two, across a missing row) are joined into a *lattice
family*, and every member gets integer coordinates (carbons, H2 count) read off
its mass. A feature whose m/z sits on a family's lattice point is corrected
with the weighted mean curve residual of its nearest family members: at most
four, at least two, no further than three lattice steps away (one double bond
counts like two carbons). Like the series term it works from *m/z* alone — a
lattice family is **not an identification**. Four safeguards keep it honest:

* **Order validation on both columns** — a family is kept only if its members
  elute in the order reversed-phase chromatography demands (more carbons and/or
  fewer double bonds elutes later) on your column *and* on the reference column;
  members that break the order are dropped, and a family needs at least six
  members over at least two double-bond counts.
* **Leave-own-family-out, with a reach limit** — a feature is never corrected
  with pairs that co-elute with it, and it borrows only from members within
  three lattice steps.
* **A fallback behind the series term, never a replacement** — a feature the
  series term corrects keeps that correction; the lattice term fills in only
  where the series term is silent.
* **A 20% gate, judged behind the series term** — the term is applied only if
  predicting every matched pair it covers this way, among the pairs the series
  term leaves alone, cuts their mean squared error by at least 20% below the
  curve alone, on at least 20 such pairs; otherwise it leaves the stage-1 curve
  untouched.

Lattice families share a weak spot: a species from another class that shares a
family's formula series lands on the same lattice and is corrected as if it
belonged. LPC and ether PC share one, and so does the sodium adduct of a lipid
two carbons shorter and three double bonds more saturated (2.4 mDa from the
protonated lattice point); see Accuracy below. And on a matrix or method where
the families carry no consistent selectivity difference, the gate does not
engage.

**Stage 2 — sample anchors, gated (opt-in since 1.2.2).** Seventeen
high-abundance endogenous plasma lipids (LPC/PC/SM/CE/TG) are located in both
sample runs, validated by adduct consistency, within-class elution order, and
a curve-assisted isomer re-pick. Their residuals against the Stage-1 result
(curve, plus the series and lattice terms where they engaged) are
class-systematic on long gradients, so the correction is decomposed into a
class offset plus a class-detrended piecewise-linear term, both shrunk by
leave-one-anchor-out cross-validation.

The correction is then **gated**: it is applied only if it reduces
leave-one-anchor-out error by at least 20%. Otherwise the pure Stage-1 curve is
used. Since 1.2.2 the whole stage is **off by default** — enable it with
`use_sample_anchors=True` or the CLI's `--sample-anchors`. On a non-plasma
matrix, too few anchors validate and Stage 2 simply does not engage — the run
succeeds on the curve alone, which is the intended behaviour, not a failure.

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
custom manifest, and the matching parameters (`match_mz_tol_ppm` for feature
matching and `mz_tol_ppm` for panel standards — both 15 ppm by default —
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
rt-anchor calibrate ... --no-series-term   # skip stage 1b (on by default)
rt-anchor calibrate ... --no-lattice-term  # skip stage 1c (on by default since
                                           # 1.2.3); alone, it reproduces v1.2.2's
                                           # default output
rt-anchor calibrate ... --sample-anchors   # opt in to stage 2 (plasma/serum only;
                                           # off by default since 1.2.2)
rt-anchor calibrate ... --no-series-term --no-lattice-term --sample-anchors
                                           # reproduce v1.2.1's default output
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
| `warp_source` | per row: `curve`, plus `+series` or `+lattice` where that term actually corrected the row, plus `+anchors` on every row when stage 2 engaged — e.g. `curve+lattice`, `curve+series+anchors`; a row carries both `+series` and `+lattice` only in the per-sample tier |
| `RI_spread`, `n_contributing` | per-injection dispersion (per-sample tier only) |
| `series_correction_min` | minutes the series term added to the curve's prediction (0 = no correction) |
| `series_n_members` | number of series members behind the correction (0 = not corrected) |
| `lattice_correction_min` | minutes the lattice term added to the curve's prediction (0 = no correction) |
| `lattice_n_members` | number of lattice-family members in reach behind the correction (0 = not corrected) |

Companion files: **`*_model.json`** (curve, the series and lattice terms' fit
and gate decisions, stage-2 anchors, iRT definition, reference provenance,
detection QC), **`*_anchors.csv`** (the Stage-2 anchors and their residuals),
**`*_pairs.csv`** (every matched pair, whether it survived trimming, and
whether it sits in a validated homologous series or in a lattice family),
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

These figures were measured with the validation configuration, a flat 0.008 Da
matching window (`match_mz_tol_ppm=0, match_mz_tol_da=0.008` reproduces it).
The default is now 15 ppm, and the table has not been re-measured at that
width. The spec records that `max(15 ppm, 0.008 Da)` — the same window above
m/z 533 — turns column 90's gate on (28%) while lowering its held-out error, so
that row is the one most likely to move.

The table was measured before 1.2.2, with stage 2 enabled and without the
series term; `--no-series-term --no-lattice-term --sample-anchors` reproduces
that configuration.

The series term was measured separately, on a six-method human-serum
validation set (five source methods calibrated onto the reference column at
the default 15 ppm matching window; 134 lipids located from exact mass
independently of the calibration, 576 lipid-by-method points): it lowered the
median |Cal_RT − reference RT| from 0.106 to 0.064 min and the 90th percentile
from 0.341 to 0.271 min, and the share of points within 0.2 min rose from 74%
to 82%. The gate engaged on four methods and declined on the one whose
selectivity already matched the reference. Validated series covered about half
of the points; a feature outside any validated series keeps the stage-1 curve.
Every prediction is leave-own-family-out by construction, so these are not
in-sample figures.

The lattice term (1.2.3) was measured on the same six-method human-serum
validation set that was used for the series term (five source methods
calibrated onto the reference column at the default 15 ppm matching window; 134
lipids located from exact mass independently of the calibration, 576
lipid-by-method points). Adding it to 1.2.2 lowered the median |Cal_RT −
reference RT| from 0.064 to 0.059 min and the 90th percentile from 0.271 to
0.218 min; the share of points within 0.2 min rose from 82% to 87%. The lattice
gate engaged on the same four methods as the series gate (leave-own-out MSE
46–51% below the curve alone, on the pairs the series term leaves alone) and
declined on the one whose selectivity already matched the reference. On those
four methods the series term corrected 301 of the points and the lattice term
139 more. The gain sits where the series term could not reach: cholesteryl
esters (55 points, none of them touched by the series term) went from a median
of 0.175 to 0.066 min, from 64% to 87% within 0.2 min, and their class-wide
offset (median signed error) from −0.147 to −0.004 min. Point by point, 55
improved by more than 0.1 min and 8 got worse by more than 0.1 min.

One class got worse: lysophosphatidylcholines (52 points) went from a median of
0.055 to 0.074 min and from 83% to 79% within 0.2 min; five of the eight points
that worsened by more than 0.1 min are LPC. LPC shares its lattice with ether
PC (the same formula series), and the sodium adducts of lighter LPCs land on
LPC lattice points, so its neighbours on the lattice are not all LPC. Lattice
families are recognised from m/z alone and are not identifications, and where
the families carry no consistent selectivity difference the gate does not
engage. Every prediction is leave-own-family-out by construction, so these are
not in-sample figures; but the term's design choices (the neighbour rule and
its reach) were made on this same data set, and it has not been tested on
another matrix or instrument.

Two honest caveats. `Cal_RT_min` can show hairline non-monotonicity (observed
worst case 0.033 min, an order of magnitude below the method's own accuracy)
where the Stage-2 correction is active, because curve + piecewise-linear
correction is not monotone by construction. And Stage 2's anchor panel is
human plasma/serum: other matrices get the Stage-1 curve only.

## Desktop app

The graphical client is the RT Anchor desktop app — Overview · Detection ·
Profile · Curve · Table · Export over this engine. Installers for macOS and
Windows are on the [releases page](https://github.com/Bowen999/rt-anchor/releases/latest);
the source is `RT Anchor Desktop/` in the
[repository](https://github.com/Bowen999/rt-anchor), and runs from a checkout:

```bash
git clone https://github.com/Bowen999/rt-anchor && cd rt-anchor
pip install -e "./rt_anchor[app,report]"     # this engine + pywebview
python "RT Anchor Desktop/app.py"
```

## License

MIT
