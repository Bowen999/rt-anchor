# RT Anchor v2 — cross-column RT calibration: port specification

Status: **implementation contract**. Supersedes the method described in
`RI_CALIBRATION_SPEC.md` (v1, standard-panel warp). Written 2026-09-06.

The v1 engine warped a sample's RT onto a dimensionless iRT scale using the
spiked standard panel found *in that same sample*. v2 replaces that method
wholesale with the cross-column calibration validated in
`v2/5 chroms/results/rt_calibration/` (reference implementation, reproduced
2026-09-06 — see §9 for the accepted numbers).

---

## 1. What the user gives us, what they get back

### Inputs (required)

| Input | Meaning | Minimum content |
|---|---|---|
| `sample_table` | biological sample run, the table to calibrate | precursor **m/z** + **RT** |
| `standards_table` | standards-mixture run **from the same column / same method** | precursor **m/z** + **RT** |
| `panel` | `"mix15"` \| `"mix21"` \| `"none"` | see §6 |
| `polarity` | `positive` \| `negative` | |

`sample_table` and `standards_table` **must come from the same LC column and
gradient**. That is the single assumption the whole method rests on: the
standards run tells us how *this column's* time axis maps onto the reference
column's, and the sample run inherits that map.

### Inputs (advanced / defaults)

| Input | Default |
|---|---|
| `reference_sample` | bundled **35-min column** serum run (`reference_data/col35/reference_sample.csv`) |
| `reference_standards` | bundled **35-min column** standards run (`reference_data/col35/reference_standards.csv`) |
| everything in `CalibrationConfig` | §7 |

The bundled reference pair is the *Carly 35-min* dataset (internally column
`25`), the run against which the v2 method was developed and validated. It
defines both the `Cal_RT` time axis and the iRT ruler. Users may substitute
their own reference pair; the outputs then live on *that* column's axis and
are no longer comparable with default-reference runs (the model JSON records
which reference was used — §5).

### Output

The sample's **original table, unmodified, with columns appended** (the engine
never drops or reorders an input column):

| Column | Meaning |
|---|---|
| `Cal_RT_min` | the feature's RT expressed on the reference column's time axis (minutes) |
| `Cal_RT_uncertainty_min` | 1-sigma local uncertainty of `Cal_RT_min` |
| `iRT` | dimensionless 1–100 retention index (§4) |
| `iRT_uncertainty` | 1-sigma uncertainty of `iRT` |
| `iRT_reliability` | `high` \| `medium` \| `low` \| `none` |
| `is_extrapolated` | RT outside the span covered by matched pairs |
| `calibration_scope` | `project` \| `sample` |
| `warp_source` | `curve` \| `curve+anchors` |
| `RI_spread` | per-injection dispersion of `iRT` (per-sample tier only, else NaN) |
| `n_contributing` | number of injections contributing (per-sample tier only, else 1) |

Name collision rule is unchanged: if the input table already has a column of
that name, the appended one is prefixed `rtanchor_`.

Companion outputs (`write_results`) — §5.

---

## 2. The algorithm

Two stages, ported from `v2/5 chroms/results/rt_calibration/rt_calibrate.py`
and `plasma_lipids.py`. **Read those two files before implementing.** The port
must reproduce their numbers (§9), so port the arithmetic literally; the only
permitted changes are the generalisations in §3.

### Stage 1 — anchor-free monotone curve `f: RT_source → RT_reference`

1. **Reciprocal m/z matching** of the *standards* runs — every feature of the
   user's standards run against every feature of the reference standards run,
   window `match_mz_tol_da` (§7). For each feature the in-window candidate with
   the highest **S/N** wins — `FeatureTable.sn()`, i.e. the format's S/N column
   where one exists and summed abundance where it does not — in both
   directions; only mutual best matches are kept, giving a 1-to-1 pairing. No
   standard identities are used, which is why the method still works when the
   mixture is unknown or undetectable.

   > Ranking by S/N rather than by summed abundance is load-bearing, not
   > cosmetic: abundance scales differ between runs and injections, S/N is a
   > within-run quality measure. Ranking anonymous match candidates by summed
   > abundance was measured to break parity outright — column 90's gate flips
   > from OFF to engaged, and only 11 of 40 check lipids stay within 0.01 min.
2. **Sample pairs merged in** — the same reciprocal matching between the user's
   *sample* run and the reference *sample* run. Serum features densely cover
   the very early and very late elution regions where the mixture is sparse,
   and this is what makes the curve trustworthy at the edges. Pairs whose
   source m/z falls within `mz_tol` of any **stage-2 anchor candidate** m/z are
   excluded, so anchor residuals cannot be contaminated by the anchors' own
   feature pairs.
3. **Robust monotone fit** — LOESS (`statsmodels.lowess`, `frac=curve_frac`,
   `it=3`) → `sklearn.isotonic.IsotonicRegression(increasing=True)` →
   `scipy.interpolate.PchipInterpolator` through the isotonic knots, with up to
   `curve_iter` rounds of MAD-based outlier removal (`|resid − median| >
   curve_mad_k · 1.4826 · MAD` dropped; removal only, never resurrection).
   Below `curve_min_points` pairs it degrades to a straight line fit.
4. **Ends** — linear extension at the terminal slope beyond the knot range;
   every RT outside `[x0, x1]` is flagged `is_extrapolated`.

`curve_frac = 0.1` is not a guess: it was chosen by 5-fold CV on the
panel-masked training pairs of all five columns. Do not change it.

### Stage 2 — class-aware sample-anchor refinement (optional, gated)

Anchors are **endogenous high-abundance plasma lipids**, not mixture
standards — the v2 hit-rate table showed most mixture standards are simply
absent from serum. The 17-lipid panel, its exact masses computed from sum
composition, and its within-class RT-order rules are in `plasma_lipids.py`.

1. **Validate** each candidate in the user's sample run *and* the reference
   sample run: m/z window → adduct-consistency filter → Fill%/intensity pick →
   within-class RT-order rules → isomer-candidate census.
2. **Curve-assisted isomer refinement** — with the Stage-1 curve in hand,
   re-pick each source-run candidate as the in-window, adduct-consistent
   feature whose *curve-calibrated* RT is closest to the reference-run RT
   (ties within `tie_tol_min` broken by Fill% then intensity). This is what
   caught the LPC 18:1 mis-pick on column 15. The curve is built from
   anonymous m/z-matched features only, so no lipid identity leaks in.
3. **Residual decomposition** — anchor residuals on long gradients are
   class-systematic (SM/CE sit high, PC/TG low) and the classes interleave
   along the RT axis, so a single piecewise-linear correction both misses the
   sign flips and gets shrunk by cross-class contamination:

   ```
   resid_i               = m[class_i] + g(x_i) + eps
   correction(x, cls)    = lam_c · m[cls] + lam_g · g(x)
   Cal_RT                = f(RT) + correction(RT, class)
   ```

   `m_c` = median residual within class `c`; `g` = piecewise-linear through
   the class-detrended residuals; `(lam_g, lam_c)` chosen on a 2-D grid by
   leave-one-anchor-out MSE. Features of unknown class get the global median
   offset. Without class labels, fall back to the plain shrunk
   `AnchorRefiner` with LOO-chosen `lam`.
4. **The gate** — engage the correction only when its LOO MSE is at least
   `anchor_gate_min_mse_reduction` (20%) below the zero-correction MSE.
   Otherwise fall back to the pure Stage-1 curve. This is the "do no harm"
   rule that correctly switches column 90 off (14% gain), where the anchor
   residuals are class-incoherent single-lipid deviations that do not
   generalise.
5. Skipped entirely when fewer than `min_anchors` (3) anchors validate, or
   when `use_sample_anchors=False`.

**Matrix caveat.** The 17-lipid panel is human plasma/serum. For any other
matrix Stage 2 will find too few anchors and the engine silently uses the
Stage-1 curve — which is by design (the archaeal-lipid case): log it, report
it, do not fail.

---

## 3. Required generalisations over the reference implementation

The v2 scripts read MS-DIAL exports directly and hard-code MS-DIAL column
names. The engine supports MS-DIAL / MZmine / MassCube / LipidScreener, so the
port must go through `FeatureTable`. Add these accessors to
`io/schema.py:FeatureTable` and use them everywhere instead of literal column
names:

| accessor | MS-DIAL source | fallback |
|---|---|---|
| `intensity()` | sum of `sample_cols`, else `S/N average` | `peak_height`/`peak_area`/`maxo`/`height`/`area`/`intensity`; else NaN |
| `sn()` | `S/N average` | `intensity()` |
| `adducts()` | `Adduct type` | `adduct` / `main_adduct`; else `None` (skip the adduct filter) |
| `fill()` | `Fill %` | `detection_rate`; else `None` (skip the Fill% preference) |
| `feature_id()` | `Alignment ID` | `feature_id` / `row ID`; else the row index |

Reuse `identify._intensity` for `intensity()` rather than writing a third copy.

Two behaviours must be preserved exactly where the column exists, and skipped
(not faked) where it does not: the adduct-consistency filter and the
`Fill % >= 1.0` preference.

**Loader fix (small, do it):** `io/adapters.load_msdial` currently treats the
MS-DIAL `Average` / `Stdev` trailing columns as sample columns, because their
class-header label is the literal string `NA`. Exclude columns whose class
label is empty *or* case-insensitively `na`/`nan`. This inflates `intensity()`
today; the fix must be made before parity testing (§9) so the baseline is
compared against the fixed behaviour.

---

## 4. The iRT scale

```
landmarks = RTs at which the chosen panel's standards are detected
            in the REFERENCE standards run
iRT(rt_ref) = 1 + 99 · (rt_ref − landmark_first) / (landmark_last − landmark_first)
```

Ported from `rt_calibrate.IRTMapper`. Note honestly, in the docstring and in
the report: because the reference implementation spaces iRT *linearly in RT*
across the landmarks, the landmarks are collinear and the mapping is
**affine** — the interior landmarks do not bend the scale, they only certify
that the panel was found. Only the earliest and latest detected landmark set
the endpoints. (This is the deliberate v2 definition and the user has chosen
it; §10 records the alternative.)

Consequences to respect:

* iRT is a property of the **reference column plus the landmark panel** alone.
  It does not require the user's own column to contain the panel. Compute it
  as `iRT(Cal_RT_min)`.
* Detect landmarks with the existing `identify.build_native_template` against
  the reference standards run (it already does m/z window → intensity pick →
  longest-increasing-subsequence monotonicity enforcement).
* If fewer than 2 landmarks are found, `iRT` is NaN for every feature,
  `iRT_reliability = "none"`, and the log/model say why. `Cal_RT_min` is still
  produced — the two outputs are independent.

---

## 5. Companion outputs

`write_results(result, prefix)` writes, unchanged in spirit from v1:

| file | content |
|---|---|
| `<p>_calibrated.csv` | the output table of §1 |
| `<p>_model.json` | §5.1 |
| `<p>_anchors.csv` | the Stage-2 anchors actually used: `label`, `lipid_class`, `rt_src`, `rt_ref`, `residual_min`, `n_isomer_candidates`, `isomer_rts`, `pick_refined`, `dropped_by_sanity_filter` |
| `<p>_pairs.csv` | **new** — the Stage-1 matched pairs: `mz_src`, `rt_src`, `rt_ref`, `source` (`standards`\|`sample`), `kept` (survived MAD trimming) |
| `<p>_landmarks.csv` | **new** — panel standards located in the reference standards run: `name`, `class`, `mz`, `rt_ref_run_min`, `iRT` |
| `<p>_log.txt` | the run log |
| `<p>_report.html` / `.pdf` | the visual report |

### 5.1 `model.json` — required keys

```jsonc
{
  "method": "cross-column-v2",
  "calibration_scope": "project",
  "reference": { "key": "col35", "label": "...", "sample": "...", "standards": "...",
                 "is_default": true, "provenance": { ... } },
  "panel": { "key": "mix21", "label": "...", "n_standards": 21 },
  "curve": { "n_pairs": 0, "n_pairs_kept": 0, "n_pairs_standards": 0, "n_pairs_sample": 0,
             "frac": 0.1, "rt_span_src_min": [0, 0], "rt_span_ref_min": [0, 0],
             "residual_min": { "median_abs": 0, "p90_abs": 0 } },
  "anchors": { "engaged": true, "n_validated": 0, "n_used": 0, "lam_g": 0, "lam_c": 0,
               "gate_mse_reduction": 0, "gate_threshold": 0.2,
               "loo_residual_min": { "median": 0, "p90": 0, "max": 0 },
               "table": [ ... ] },
  "irt": { "n_landmarks": 0, "landmark_rt_span_min": [0, 0], "irt_range": [1, 100],
           "definition": "affine on the reference-column RT axis", "panel_used": "mix21" },
  "detection_qc": { "user_standards_run": { "n_panel": 21, "n_detected": 0, "native_rt": {} } },
  "n_features": 0, "n_features_extrapolated": 0,
  "confidence_counts": { "high": 0, "medium": 0, "low": 0, "none": 0 },
  "config": { ... }
}
```

`detection_qc` is what keeps the app's **Detection** section meaningful: it is
where the chosen panel's standards were found in the *user's own* standards
run. It is reporting only — it no longer drives the calibration.

### 5.2 Uncertainty (new, replaces the v1 LOO-warp sigma)

```
sigma_curve(rt) = local robust scatter of the matched-pair residuals about the
                  fitted curve, from the nearest `sigma_window_pairs` (default 50)
                  pairs in RT, as 1.4826 · MAD  [minutes]
sigma_anchor    = RMS of the Stage-2 leave-one-anchor-out residuals, or 0 when
                  the correction is gated off                       [minutes]
Cal_RT_uncertainty_min = sqrt(sigma_curve^2 + sigma_anchor^2)
iRT_uncertainty        = Cal_RT_uncertainty_min · d(iRT)/d(rt_ref)
```

`iRT_reliability`: `high` if `iRT_uncertainty < conf_high_irt`, `low` if
`> conf_low_irt` **or** `is_extrapolated`, else `medium`; `none` when `iRT` is
NaN. Thresholds keep their v1 values and meaning.

---

## 6. Panel choice

Panel definitions move out of the desktop app (`rtad/mixtures.py`) **into the
engine** (`rt_anchor/mixtures.py`) so the two can never drift. `rtad.mixtures`
becomes a thin re-export for backwards compatibility.

| key | UI label | effect |
|---|---|---|
| `mix15` | **15** — Caley lipid RT-calibration mix | landmark panel + detection QC |
| `mix21` | **21** — Mix 4.4 extended panel | landmark panel + detection QC |
| `none` | **Others** — unknown / no panel | no identity claimed for the user's mixture |

`mix21lpc` is retained as a non-default key so `v2/5 chroms/iRT/run_iRT.py`
and existing results stay reproducible, but it is **not** offered as a card.

**`none` behaviour.** Stage 1 needs no identities, so calibration proceeds
normally. For the iRT ruler the engine falls back to `irt_landmark_panel`
(config, default `mix21`) detected on the *reference* standards run, and
records `irt.panel_used = "mix21 (fallback: user panel = none)"` plus a log
line. Detection QC is skipped. If a user-supplied reference standards run
yields <2 landmarks, iRT is NaN per §4.

---

## 7. `CalibrationConfig` — v2 fields

**Two m/z tolerances, two jobs.** They must not be conflated:

* `match_mz_tol_da = 0.008` — **absolute**, for *anonymous* work: reciprocal
  feature matching (§2.1), the panel-mask exclusion, and the plasma-lipid m/z
  windows. This is the value the v2 method was validated with, so it is what
  ships by default.
* `mz_tol_ppm = 15`, `mz_tol_min_da = 0.0` — for *targeted* panel-standard
  identification (landmarks, detection QC) via `identify`. This is the v1
  meaning, unchanged.

Do **not** apply the ppm window to anonymous matching. Measured consequence of
doing so: at `max(15 ppm, 0.008 Da)` the window widens above m/z 533, column
90's gate flips from OFF (14%) to engaged (28%), and 7 of 40 check lipids leave
the 0.01 min band. Widening is not obviously *worse* there — column 90's
held-out error improves — but the shipping default must reproduce the validated
method, and a decision that flips on the tolerance is exactly the decision to
pin rather than to drift.

Keep: `mz_tol_ppm`, `mz_tol_min_da`, `rt_window_min`, `rt_window_seed_min`,
`min_anchors`, `conf_high_irt`, `conf_low_irt`, `min_gaussian_similarity`,
`asymmetry_range`, `rt_unit`, `verbose`, and the `qtof()` / `orbitrap()`
presets. Also keep `tie_epsilon_min` and `max_drop_anchors`: they are read by
`identify._drop_nonmonotone`, which survives as part of the detection-QC and
landmark path.

Remove (v1 warp only): `scale_rt_lo_min`, `scale_rt_hi_min`,
`irt_from_rt()`, `slope_irt_per_min()`, `coverage_gap_min`,
`min_anchors_hard`, `sigma_gap_factor`, `stds_fallback`,
`max_extrapolation_min`. `from_dict` must raise `ConfigError` naming the v2
replacement when it sees a removed key, so old config files fail loudly
instead of silently.

Removing `irt_from_rt` breaks `panel.build_panel`, which stamps a v1 affine
`irt` onto every target. Fix it there: `Panel.targets.irt` becomes NaN and the
column is kept only for shape compatibility — the v2 iRT scale is derived from
the landmarks, never from scale constants. `irt.landmark_panel` already
bypasses `build_panel` for exactly this reason; make `build_panel` agree with
it rather than the reverse.

Add:

| field | default | note |
|---|---|---|
| `match_mz_tol_da` | `0.008` | absolute window for anonymous matching — see above |
| `curve_frac` | `0.1` | LOESS fraction — CV-chosen, do not change |
| `curve_iter` | `3` | MAD outlier rounds |
| `curve_mad_k` | `3.0` | |
| `curve_min_points` | `20` | below this the curve degrades to a line |
| `use_sample_pairs` | `True` | merge sample-run pairs into the curve |
| `use_sample_anchors` | `True` | Stage 2 on/off |
| `class_aware` | `True` | class-offset decomposition |
| `anchor_gate_min_mse_reduction` | `0.2` | the gate |
| `anchor_tie_tol_min` | `0.15` | isomer re-pick tie window |
| `isomer_sn_frac` | `0.10` | isomer-census threshold |
| `sigma_window_pairs` | `50` | local scatter window |
| `irt_landmark_panel` | `"mix21"` | fallback landmarks when `panel="none"` |
| `extrapolate` | `True` | v2 always assigns a value and flags it |
| `extrapolate_mode` | `"linear"` | `"linear"` = v2 terminal-slope extension; `"clamp"` = hold the edge value |

`extrapolate` defaulting to `True` is a **deliberate change from v1** (which
defaulted to NaN beyond span): the v2 curve extends linearly by construction
and every out-of-span row carries `is_extrapolated=True` and `low`
reliability. Say so in the changelog.

---

## 8. Module layout

```
rt_anchor/
  crosscolumn.py   NEW  match_features_by_mz · MonotoneCurve · fit_robust_curve ·
                        AnchorRefiner · ClassAwareRefiner · loo_mse_* ·
                        choose_lambda_loo* · ColumnCalibrator · build_calibrator
  plasma_lipids.py NEW  plasma_lipid_candidates · RT_ORDER_RULES ·
                        validate_candidates · summarize_picks · refine_picks_with_curves
  irt.py           NEW  IRTMapper · detect_landmarks
  reference.py     NEW  REFERENCE_SETS · resolve_reference · load_reference
  mixtures.py      NEW  moved from rtad/mixtures.py (+ the `none` key)
  reference_data/col35/{reference_sample.csv,reference_standards.csv,provenance.json}
                   DONE (already created; RT/mz byte-identical to
                        v2/5 chroms/data/25/, two spectrum columns dropped)

  pipeline.py      REWRITTEN  calibrate() orchestrates the above
  config.py        REWRITTEN  §7
  io/schema.py     EXTENDED   §3 accessors + RESULT_COLUMNS
  io/adapters.py   FIXED      §3 loader fix
  identify.py      KEPT       build_native_template (landmarks + detection QC);
                              identify_anchors kept for detection QC only
  calibrate.py     DELETED    MonotoneWarp / fit_warp / apply_warp (v1 warp)
  panel.py         KEPT       manifest/adduct machinery; DEFAULT_MANIFEST stays
  helpers.py       UPDATED    §5 writers
  cli.py           UPDATED    §11
  viz/             UPDATED    §12
  __init__.py      UPDATED    exports
```

`calibrate()` signature:

```python
def calibrate(sample_table: str,
              polarity: str,
              standards_table: str,
              panel: str = "mix21",
              reference_sample: str | None = None,     # None -> bundled col35
              reference_standards: str | None = None,  # None -> bundled col35
              reference_key: str = "col35",
              single_files: list[str] | None = None,
              config: CalibrationConfig | None = None,
              manifest: pd.DataFrame | None = None,    # overrides `panel`
              source_format: str | None = None,
              rt_unit: str | None = None) -> CalibrationResult
```

Guard rails to implement as clean `RtAnchorError` subclasses:

* `standards_table` missing → `ConfigError` (unchanged wording intent).
* fewer than `curve_min_points` matched pairs and no sample pairs →
  `CalibrationError` naming m/z tolerance and polarity as the likely causes.
* user's standards run and the reference standards run share <5 matched pairs →
  `CalibrationError`: almost always a polarity mismatch or a wrong reference.
* `reference_sample`/`reference_standards` supplied one without the other →
  `ConfigError`.

**Per-injection tier** (`single_files`) is retained: each injection is
m/z-matched against the reference *sample* run, gets its own curve, and the
per-feature `iRT` is the median with `RI_spread` = 1.4826·MAD across
injections; `calibration_scope="sample"`. Features no injection covers fall
back to the project-level curve.

---

## 9. Parity acceptance (non-negotiable)

A verified reference run of the v2 implementation is at
`<scratchpad>/v2_baseline/` (regenerate with the command in
`v2/5 chroms/results/rt_calibration/README.md`). The port must reproduce it:

| check | tolerance |
|---|---|
| `Cal_RT` for the 10 check lipids × 4 source columns vs `manual_check_table_with_anchors.csv` | **± 0.01 min** |
| `(lam_g, lam_c)` and gate decision per column | **exact** (15: 0.0/0.8 engaged 48% · 22: 0.8/0.7 engaged 30% · 45: 1.0/0.7 engaged 59% · 90: gated OFF at 14%) |
| number of validated plasma lipids | 17/17 |
| held-out serum median abs error per column | within **± 0.02 min** of 0.102 / 0.113 / 0.175 / 0.193 (no-anchor) |
| iRT landmark span on the reference run | 1.46–21.61 min with the 15-standard list |

Differences are expected in two places and must be *explained, not
tolerated silently*: the §3 loader fix changes intensity-based candidate
picking, and production uses **all** sample pairs where the v2 evaluation used
an 80% training split. Run the parity check both ways
(`sample_pair_frac=0.8, seed=0` reproduces the baseline exactly) and report
both numbers.

---

## 10. Deliberate decisions worth revisiting

Recorded here so they are not silently re-litigated. Discussed with the user
at the end of the run.

1. **iRT is affine, not piecewise-linear** (§4). Interior landmarks are
   decorative. A rank-spaced definition (landmark *k* of *n* → iRT
   `1 + 99k/(n−1)`) would make the interior landmarks load-bearing and the
   scale robust to gradient shape, at the cost of breaking comparability with
   every number produced so far.
2. **Panel choice moves the iRT scale**, because mix15 and mix21 have
   different earliest/latest detected landmarks on the reference run. iRT
   values are therefore comparable only within a panel choice. The alternative
   — always ruling with a fixed panel — is one config line
   (`irt_landmark_panel`) away.
3. **The 17-lipid Stage-2 panel is plasma-specific.** Non-plasma matrices get
   the Stage-1 curve only. Correct, but it means two users with the same data
   and different matrices get corrections of different quality with no error.
4. **`extrapolate` now defaults to `True`** (§7).
5. **The per-injection tier's semantics changed**: injections are now matched
   against the reference sample run rather than self-anchored on the panel.

---

## 11. CLI

```
rt-anchor calibrate --samples S --standards T --polarity positive \
                    --panel {mix15,mix21,none} \
                    [--reference-sample R --reference-standards RT] \
                    [--single-files ...] [--no-sample-anchors] [--no-sample-pairs] \
                    [--curve-frac F] [--mz-tol-ppm P] [--mz-tol-da D] \
                    [--min-anchors N] [--no-extrapolate] \
                    [--extrapolate-mode {linear,clamp}] \
                    --out PREFIX
rt-anchor describe <table>            # unchanged
rt-anchor references                  # NEW: list bundled reference datasets
```

Removed flags (`--manifest` keeps working as a `panel` override;
`--reference` for the old baseline CSV, `--no-stds-fallback` are gone) must
exit 2 with a message naming the replacement.

---

## 12. Visualisation / report / desktop app

The report and the in-app charts must reflect the new method. Sections:

| section | v1 content | v2 content |
|---|---|---|
| Overview | KPI scorecard + radar | same shape; KPIs become n pairs, curve residual median/P90, anchor gate state, n landmarks, % extrapolated |
| Detection | panel standards found in the sample | panel standards found in the **user's standards run** (`detection_qc`), explicitly labelled "QC only — does not drive the calibration" |
| Profile | RT/iRT feature profile | unchanged, driven by `iRT` |
| Warp | the panel warp curve | **the Stage-1 curve**: matched pairs scatter (standards vs sample coloured differently), the fitted curve, the anchor-refined curve when engaged, the Stage-2 anchors as rings, and the residual panel below — i.e. the two panels of `rt_calibration_diagnostics.png` |
| Anchors | *(new)* | the Stage-2 anchor table with residuals, class offsets, `(lam_g, lam_c)`, gate state and reason |
| Repeatability | per-sample only | unchanged |
| Export | CSV / report | + the new `_pairs.csv` / `_landmarks.csv` |

Desktop app (`RT Anchor Desktop`, **macOS source only** — do not touch
`app_win.py` or `build/windows/`):

* Input screen: panel cards become **15 / 21 / Others**, from
  `rt_anchor.mixtures`.
* Advanced: two new file pickers, **Reference sample** and **Reference
  standards**, both showing "bundled 35-min column" until overridden, plus a
  "Reset to bundled reference" control. Add `use_sample_anchors`,
  `use_sample_pairs`, `curve_frac`. Remove the `stds_fallback` control (§7).
* `build/macos/RTAnchor.spec` must bundle `rt_anchor/reference_data/**` as
  data files (PyInstaller will not pick up package data automatically), and
  `statsmodels` + `sklearn` must be present in the frozen bundle — verify with
  `app.py --selftest`, extending `_selftest()` to run a real 2-run calibration
  on the bundled reference against itself.
* `rtad/api.py`: `_lpc_seed_window` and the mix21lpc special-case are dead
  under the new method — remove them; `run_calibration` passes `panel`,
  `reference_sample`, `reference_standards`.

---

## 13. Definition of done

1. `pytest rt_anchor/tests` green, with the v1 warp tests replaced by v2
   equivalents (curve monotonicity, reciprocal matching, gate behaviour,
   class-aware LOO, iRT affinity, extrapolation flags, reference resolution,
   `none` panel, non-plasma matrix falls back to Stage 1 cleanly).
2. §9 parity table reproduced and reported.
3. `rt-anchor calibrate` runs end-to-end on all four source columns
   (`v2/5 chroms/data/{15,22,45,90}`) against the bundled reference, writing
   every §5 output.
4. Desktop app launches from source on macOS, runs a calibration, and exports.
5. `rt_anchor/README.md`, `RT Anchor Desktop/README.md` and this spec agree
   with the code.
