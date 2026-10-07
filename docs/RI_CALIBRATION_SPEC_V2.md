# RT Anchor v2 — cross-column RT calibration: port specification

Status: **implementation contract**. Supersedes the method described in
`RI_CALIBRATION_SPEC.md` (v1, standard-panel warp). Written 2026-09-06.
The v1 spec has been removed from the tree; it is in git history
(`git show 8269da9:RI_CALIBRATION_SPEC.md`).

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
| `warp_source` | **per row** since 1.2.2: `curve`, then `+series` where the series term actually corrected the row (§14), then — since 1.2.3 — `+lattice` where the lattice term did (§15), then `+anchors` on every row when stage 2 engaged. The values: `curve` \| `curve+series` \| `curve+lattice` \| `curve+series+lattice`, each also with `+anchors`; a row carries both `+series` and `+lattice` only in the per-sample tier (§15.9) |
| `RI_spread` | per-injection dispersion of `iRT` (per-sample tier only, else NaN) |
| `n_contributing` | number of injections contributing (per-sample tier only, else 1) |
| `series_correction_min` | minutes the homologous-series term added to the curve's prediction; 0 where none was applied (§14) |
| `series_n_members` | number of series members behind the correction; 0 = not corrected (§14) |
| `lattice_correction_min` | since 1.2.3: minutes the lattice term added to the curve's prediction; 0 where none was applied (§15) |
| `lattice_n_members` | since 1.2.3: number of lattice-family members in reach behind the correction; 0 = not corrected (§15) |

Name collision rule is unchanged: if the input table already has a column of
that name, the appended one is prefixed `rtanchor_`.

Companion outputs (`write_results`) — §5.

---

## 2. The algorithm

Two stages, ported from `v2/5 chroms/results/rt_calibration/rt_calibrate.py`
and `plasma_lipids.py`. **Read those two files before implementing.** The port
must reproduce their numbers (§9), so port the arithmetic literally; the only
permitted changes are the generalisations in §3.

v1.2.2 adds a stage between the two — the homologous-series term (stage 1b,
§14) — and makes stage 2 opt-in (`use_sample_anchors=False` by default, §7).
v1.2.3 adds one more right after it — the lattice term (stage 1c, §15), the
fallback behind the series term: it corrects, from lattice neighbours, the
features the series term leaves alone.

### Stage 1 — anchor-free monotone curve `f: RT_source → RT_reference`

1. **Reciprocal m/z matching** of the *standards* runs — every feature of the
   user's standards run against every feature of the reference standards run,
   window `match_mz_tol_ppm` (§7). For each feature the in-window candidate with
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

### Stage 2 — class-aware sample-anchor refinement (optional, gated; opt-in since 1.2.2)

Anchors are **endogenous high-abundance plasma lipids**, not mixture
standards — the v2 hit-rate table showed most mixture standards are simply
absent from serum. The 17-lipid panel, its exact masses computed from sum
composition, and its within-class RT-order rules are in `plasma_lipids.py`.
Since 1.2.2 the whole stage is opt-in: `use_sample_anchors` defaults to
`False` (§7). When it runs, the anchor residuals are taken against the
stage-1 result *after* the series and lattice terms (curve + series + lattice,
§14.7, §15.8) — the anchor m/z mask on the sample pairs (stage 1, step 2) is
unchanged and applies whether or not stage 2 runs.

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
   when `use_sample_anchors=False` — the default since 1.2.2.

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
| `<p>_anchors.csv` | the Stage-2 anchors actually used: `label`, `lipid_class`, `rt_src`, `rt_ref`, `residual_min`, `loo_residual_min`, `n_isomer_candidates`, `isomer_rts`, `pick_refined`, `dropped_by_sanity_filter` |
| `<p>_pairs.csv` | **new** — the Stage-1 matched pairs: `mz_src`, `rt_src`, `rt_ref`, `source` (`standards`\|`sample`), `kept` (survived MAD trimming), `in_series` (in a validated homologous series — §14), `in_lattice` (in an accepted lattice family — §15) |
| `<p>_landmarks.csv` | **new** — panel standards located in the reference standards run: `name`, `class`, `mz`, `rt_ref_run_min`, `iRT` |
| `<p>_log.txt` | the run log |
| `<p>_report.html` / `.pdf` | the visual report |

Since 1.2.2 `residual_min` (and the gate's own arithmetic) is measured against
the stage-1 result *after* the series term — curve + series when the term
engaged, the bare curve otherwise (§14.7); since 1.2.3 after the lattice term
as well — curve + series + lattice, each as far as it engaged (§15.8). The log
gains a `series term: ...` line right after the stage-1 line (§14.6) and, since
1.2.3, a `lattice term: ...` line right after that (§15.7).

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
  "series": { "enabled": true, "engaged": false, "gate_mse_reduction": null,
              "gate_threshold": 0.2, "gate_reason": "...", "n_series": 0,
              "n_pairs_covered": 0, "n_features_corrected": 0,
              "kmd_tol": 0.008, "min_members": 4, "end_reach_ch2": 2,
              "exclude_rt_min": 0.03 },
  "lattice": { "enabled": true, "engaged": false, "gate_mse_reduction": null,
               "gate_threshold": 0.2, "gate_reason": "...", "n_families": 0,
               "n_pairs_in_families": 0, "n_pairs_covered": 0,
               "n_features_corrected": 0, "kmd_tol": 0.008, "min_members": 6,
               "max_distance": 3.0, "min_neighbours": 2, "max_neighbours": 4,
               "exclude_rt_min": 0.03 },
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

**Two m/z tolerances, two jobs — one default.** They are tuned separately, and
both default to **15 ppm with no absolute floor** (the window is
`max(ppm, floor)`):

* `match_mz_tol_ppm = 15`, `match_mz_tol_da = 0.0` — for *anonymous* work:
  reciprocal feature matching (§2.1), the panel-mask exclusion, and the
  plasma-lipid m/z windows.
* `mz_tol_ppm = 15`, `mz_tol_min_da = 0.0` — for *targeted* panel-standard
  identification (landmarks, detection QC) via `identify`. This is the v1
  meaning, unchanged.

**Changed 2026-10-04 — the anonymous default is 15 ppm.** The v2 method was
validated with a flat 0.008 Da window (`match_mz_tol_ppm = 0`,
`match_mz_tol_da = 0.008`), and the first port shipped `max(15 ppm, 0.008 Da)`;
both remain one setting away. The decision is one 15 ppm threshold for every
match. For the record, what was measured before it: at `max(15 ppm, 0.008 Da)`
the window widens above m/z 533, column 90's gate flips from OFF (14%) to
engaged (28%), and 7 of 40 check lipids leave the 0.01 min band — widening is
not obviously *worse* there, column 90's held-out error improves. The §9
numbers were measured with the flat 0.008 Da window and have not been
re-measured at 15 ppm.

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
| `match_mz_tol_ppm` | `15.0` | anonymous matching window — see above |
| `match_mz_tol_da` | `0.0` | optional absolute floor under it; `0.008` with ppm `0` is the validation window |
| `curve_frac` | `0.1` | LOESS fraction — CV-chosen, do not change |
| `curve_iter` | `3` | MAD outlier rounds |
| `curve_mad_k` | `3.0` | |
| `curve_min_points` | `20` | below this the curve degrades to a line |
| `use_sample_pairs` | `True` | merge sample-run pairs into the curve |
| `use_series_term` | `True` | stage 1b on/off (§14) |
| `series_kmd_tol` | `0.008` | tolerance on the Kendrick phase |
| `series_min_members` | `4` | smallest cluster / chain that counts as a series |
| `series_end_reach_ch2` | `2` | how many CH2 past a series end a query may sit |
| `series_exclude_rt_frac` | `0.0015` | co-elution exclusion, fraction of the curve's knot span |
| `series_exclude_rt_floor_min` | `0.03` | ... with this absolute floor, minutes |
| `series_gate_min_mse_reduction` | `0.2` | the series term's "do no harm" gate |
| `series_min_covered_pairs` | `20` | below this the series gate declines unread |
| `use_lattice_term` | `True` | stage 1c on/off (§15) |
| `lattice_kmd_tol` | `0.008` | tolerance: cluster gap, H2 link, lattice-point match |
| `lattice_min_members` | `6` | smallest lattice family |
| `lattice_max_distance` | `3.0` | how far (lattice steps) a correction may be borrowed from |
| `lattice_min_neighbours` | `2` | members needed within that distance (at least 1) |
| `lattice_max_neighbours` | `4` | nearest members averaged (at least 1) |
| `lattice_gate_min_mse_reduction` | `0.2` | the lattice term's "do no harm" gate |
| `lattice_min_covered_pairs` | `20` | below this the lattice gate declines unread |
| `use_sample_anchors` | `False` | Stage 2 on/off — opt-in since 1.2.2 (was `True`) |
| `class_aware` | `True` | class-offset decomposition |
| `anchor_gate_min_mse_reduction` | `0.2` | the gate |
| `anchor_tie_tol_min` | `0.15` | isomer re-pick tie window |
| `isomer_sn_frac` | `0.10` | isomer-census threshold |
| `sigma_window_pairs` | `50` | local scatter window |
| `irt_landmark_panel` | `"mix21"` | fallback landmarks when `panel="none"` |
| `extrapolate` | `True` | v2 always assigns a value and flags it |
| `extrapolate_mode` | `"linear"` | `"linear"` = v2 terminal-slope extension; `"clamp"` = hold the edge value |

The lattice term has no co-elution fields of its own: it uses the series term's
`series_exclude_rt_frac` / `series_exclude_rt_floor_min` — one notion of
"co-eluting", one setting (§15.6).

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
  series.py        NEW  stage 1b, v1.2.2 (§14): SeriesTerm · kendrick_phase ·
                        KENDRICK_FACTOR · CH2_MASS · PERIOD
  lattice.py       NEW  stage 1c, v1.2.3 (§15): LatticeTerm · H2_MASS · H2_PHASE
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
back to the project-level prediction (§14.8, §15.9).

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

Parity is defined for the validation window: run it with
`match_mz_tol_ppm=0, match_mz_tol_da=0.008`, not the 15 ppm default (§7).

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
6. **The default m/z window is 15 ppm** (§7), one threshold for anonymous
   matching and panel identification alike — not the flat 0.008 Da window the
   method was validated with.
7. **Stage 2 is opt-in since 1.2.2** (`use_sample_anchors=False` by default).
   The plasma-lipid anchor panel is only meaningful on human plasma/serum, and
   the series term (§14) is the safer default: it is m/z-only, makes no
   identity claim, and is gated the same way. `--sample-anchors` restores
   stage 2.
8. **Validated series cover only part of the features.** A feature outside
   every validated series keeps the stage-1 curve, and on a matrix with few
   homologous series the gate does not engage at all — on the validation set
   about half the points sat in a validated series, and the rest were not
   corrected. (Under 1.2.2; since 1.2.3 the lattice term, item 9, can reach
   some of what the series term leaves.)
9. **The lattice term is a fallback behind the series term** (§15.5), not a
   peer. A feature the series term corrects keeps that correction untouched,
   and the lattice gate is judged only on the pairs the series term leaves
   alone — so `use_lattice_term=False` gives exactly the 1.2.2 output, and
   with stage 2 off every feature the lattice term does not correct keeps its
   1.2.2 value. The precedence is fixed: nothing compares, per feature, which
   of the two would have predicted better, so a feature the series term
   corrects keeps the series correction even where lattice neighbours might
   have predicted it better.
10. **Shared lattices and sodium adducts are corrected as members** (§15.10).
    The lattice is a property of a formula series, not of a lipid class, and
    the term is m/z-only: LPC and ether PC share one lattice, and the sodium
    adduct of a lipid two carbons shorter and three double bonds more saturated
    sits 2.4 mDa from the protonated lattice point, inside `lattice_kmd_tol`.
    Such a species is corrected as if it belonged to the family. It is the
    known weak spot — LPC is the one class that got worse on the validation set
    (§15.11) — and m/z alone cannot tell these species apart.
11. **The lattice distance and the reach are fixed, not tuned per data set.**
    `d = |Δa| / 2 + |Δb|` (one double bond counts like two carbons) is a
    constant of the method, not a setting, and `lattice_max_distance`,
    `lattice_min_neighbours` and `lattice_max_neighbours` are defaults, not
    values fitted to a run. The neighbour rule and its reach were chosen on the
    same six-method validation set that measures the term's evidence (§15.11),
    so that evidence is not an out-of-sample test of those choices; the term
    has not been tested on another matrix or instrument.

---

## 11. CLI

```
rt-anchor calibrate --samples S --standards T --polarity positive \
                    --panel {mix15,mix21,none} \
                    [--reference-sample R --reference-standards RT] \
                    [--single-files ...] [--sample-anchors] [--no-sample-anchors] \
                    [--no-series-term] [--no-lattice-term] [--no-sample-pairs] \
                    [--curve-frac F] [--match-mz-tol-ppm P] [--mz-tol-da D] \
                    [--mz-tol-ppm P] \
                    [--min-anchors N] [--no-extrapolate] \
                    [--extrapolate-mode {linear,clamp}] \
                    --out PREFIX
rt-anchor describe <table>            # unchanged
rt-anchor references                  # NEW: list bundled reference datasets
```

Since 1.2.2 stage 2 is opt-in: `--sample-anchors` enables it, and
`--no-sample-anchors` — now the default — is kept for explicitness and wins if
both are given. The stage-1b series term (§14) is on by default;
`--no-series-term` switches it off. Since 1.2.3 the stage-1c lattice term
(§15) is on by default too; `--no-lattice-term` switches it off, and alone
reproduces v1.2.2's default output. `--no-series-term --no-lattice-term
--sample-anchors` reproduces v1.2.1's default output.

Removed flags (`--manifest` keeps working as a `panel` override;
`--reference` for the old baseline CSV, `--no-stds-fallback` are gone) must
exit 2 with a message naming the replacement.

---

## 12. Visualisation / report / desktop app

The report and the in-app charts must reflect the new method. Sections:

| section | v1 content | v2 content |
|---|---|---|
| Overview | KPI scorecard + radar | same shape; KPIs become n pairs, curve residual median/P90, the series and lattice terms' gate states (one tile, two stacked rows; the anchor gate rides on that tile's note when stage 2 is requested), n landmarks, % extrapolated |
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

Since 1.2.2 (series term) and 1.2.3 (lattice term) the report and the app also
show the stage-1 terms, from the `series` and `lattice` blocks of `model.json` —
they read the model, they never recompute (§14.6, §15.7):

* Report: the fourth KPI tile is **Series + lattice** — two stacked rows,
  `series` and `lattice`, each the term's headline (`on · {p}%`,
  `off · {p}%` for a shortfall against the gate, `off` for a measured loss, for
  too few covered pairs or when disabled, `not needed` for self-calibration),
  with the note `{k_series} + {k_lattice} of {n} features corrected` and, when
  stage 2 was requested, the anchor gate's state appended. The radar axis
  **Series + lattice** is the larger of the two gates' leave-own-out MSE
  reductions among the terms that are engaged, clipped to [0, 1], and 0 when
  neither is engaged — zero is not "bad": where the selectivity already matches
  the reference there is nothing to correct, and the gates declining is the
  method working. The method-facts table has a "Series term" and a "Lattice
  term" row (the gate reason, or "not fitted"), and the curve figure of the PDF
  prints both decisions, one line each, under its title.
* Desktop app: the Advanced section has a **Series term** and a **Lattice term**
  checkbox (both on by default; they send `use_series_term` and
  `use_lattice_term`); the Curve section lists both gate reasons; the Table
  preview adds `series corr (min)`, `series n`, `lattice corr (min)` and
  `lattice n`; `run_info.json` records `results.series` and `results.lattice`
  (`enabled`, `engaged`, `gate_mse_reduction`, `gate_threshold`, `gate_reason`,
  the term's count of series or families, `n_pairs_covered`,
  `n_features_corrected`).

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

---

## 14. v1.2.2 — the homologous-series term

Stage 1b, fitted after the stage-1 curve and before the (now opt-in) stage 2.
Implemented in `series.py`; on by default (`use_series_term=True`, CLI
`--no-series-term` switches it off).

**Why.** Two methods can differ in how strongly they retain a double bond
relative to a CH2 group, so two lipids that co-elute on the source column can
sit more than a minute apart on the reference column, and no single curve can
place both. The error the curve leaves behind is systematic within a
**homologous series** (same class, same number of double bonds, different
chain length): members share it, smoothly along RT. Members of a homologous
series differ by CH2, so they can be recognised from m/z alone — no
identities — through the Kendrick mass.

### 14.1 Constants and the Kendrick phase

```
KENDRICK_FACTOR = 14.0 / 14.01565006
CH2_MASS        = 14.01565006
PERIOD          = 14.0
phase(mz)       = (mz * KENDRICK_FACTOR) mod PERIOD        # in [0, 14)
circ(p, q)      = min(|p - q|, PERIOD - |p - q|)           # circular distance
```

Rescaling the mass axis so CH2 is exactly 14.0 makes every member of one
series share a mass defect — one position on the 14-Da circle (the *phase*).

### 14.2 Fit

Input: the stage-1 pairs that survived MAD trimming — source m/z `mz`, source
RT `xa`, reference RT `xb`, curve residual `r = xb − f(xa)` — plus the curve's
source-RT knot span `x0, x1`.

1. Pairs with any non-finite value are dropped; the rest keep their given
   order, indexed `0..n−1`.
2. `excl = max(series_exclude_rt_floor_min, series_exclude_rt_frac · (x1 − x0))`
   is the co-elution exclusion radius; `reach = series_end_reach_ch2 · CH2_MASS
   + 1.0` is how far past a series end a query may hang (the +1 Da is slack for
   the mass-defect spread within a cluster).
3. **Phase clusters.** Pair indices are sorted by `(phase, index)`; a new
   cluster starts whenever the phase gap exceeds `series_kmd_tol`. If there
   are at least two clusters and `phase[first of first cluster] + PERIOD −
   phase[last of last cluster] ≤ series_kmd_tol`, the last cluster is
   prepended to the first — a cluster straddling the 0/14 seam is rejoined.
4. **Validated series.** Each cluster with at least `series_min_members`
   pairs is ordered by `(mz, index)` → `o[0..m−1]`, and its longest chain in
   which RT rises with carbon number on **both** columns is found by dynamic
   programming:

   ```
   best[i] = 1, prev[i] = -1 for all i
   for i in 0..m-1:
       for j in 0..i-1:
           if mz[o[i]] - mz[o[j]] > 7.0 and xa[o[i]] > xa[o[j]] and xb[o[i]] > xb[o[j]] \
                  and best[j] + 1 > best[i]:
               best[i] = best[j] + 1; prev[i] = j
   k     = first index attaining max(best)        # numpy argmax semantics
   chain = backtrack from k through prev, reversed   # ascending m/z
   ```

   The `> 7.0` step is "gained at least one carbon" (a CH2 is 14.01565 Da; the
   half-step admits the mass-defect spread while excluding non-homologues).
   Tie-breaking: the strict `best[j] + 1 > best[i]` keeps the earliest
   improving predecessor, and `argmax` the earliest of the longest chains. A
   chain of at least `series_min_members` pairs is a validated series; its
   **centre** is `(p0 + median(d)) mod PERIOD`, where `p0` is the phase of the
   chain's first member and `d` the members' phase offsets from `p0`, each
   wrapped into `(−7, 7]` (`d > 7 → d − 14`; `d ≤ −7 → d + 14`), so a series
   straddling the seam centres correctly. One cluster yields at most one
   series; pairs outside the chain are simply not used; series keep the order
   they were found in.

### 14.3 Prediction of one query `(mq, xq)`

The **raw** correction — `raw_correction()`, available whether or not the
gate engages — returns `(correction_min, n_members)` per query; "no
correction" is `(0.0, 0)`:

1. Non-finite `mq`/`xq`, or no validated series → `(0.0, 0)`.
2. `pq = phase(mq)`; take the series whose centre has the smallest
   `circ(centre, pq)` (first on ties). If that distance is
   `> series_kmd_tol` → `(0.0, 0)`.
3. `M` = that series' members with `|xa − xq| > excl` — the
   **leave-own-family-out** rule: a feature is never corrected with pairs that
   co-elute with it. If `len(M) < series_min_members − 1` → `(0.0, 0)`.
4. `lo` = members of `M` with `mz < mq − 7.0`; `hi` = members with
   `mz > mq + 7.0`:
   - **Both non-empty** (bracketed interpolation): `a` = last of `lo`,
     `b` = first of `hi`. If not `xa[a] < xq < xa[b]`, the query does not sit
     in the series' elution order → `(0.0, 0)`. Otherwise
     `t = (xq − xa[a]) / (xa[b] − xa[a])`, correction
     `= r[a] + t · (r[b] − r[a])`, `n_members = len(M)`.
   - **Only `lo`** (past the upper end): `a` = last of `lo`; if
     `mq − mz[a] ≤ reach` and `xq > xa[a]` → `(r[a], len(M))` — the end
     member's residual, carried at most `series_end_reach_ch2` CH2 out.
   - **Only `hi`** (below the first member): symmetric — `b` = first of `hi`;
     if `mz[b] − mq ≤ reach` and `xq < xa[b]` → `(r[b], len(M))`.
   - Otherwise `(0.0, 0)`.

The applied `correction()` is the raw correction when the gate engaged and
all zeros (and `n_members` all zero) otherwise. Queries are grouped by their
nearest series centre, so correcting a whole feature table is vectorised, not
a Python loop over features against series.

### 14.4 The gate

Every fitted pair is predicted as a query at its own `(mz, xa)` — the
co-elution exclusion turns that into a leave-own-family-out prediction —
giving `pred[i]`, `n[i]`; `cov = n > 0`.

Two counts describe the coverage: the pairs in a validated series (`in_series`
in `*_pairs.csv`) and the covered pairs (`n_pairs_covered`), those of them
that get a leave-own-out prediction, which the gate is read on; the covered
count is the smaller one, because a pair whose co-eluting neighbours are
excluded can be left with too few members, or not sit in the series' elution
order (§14.3, rules 3 and 4).

- `n_pairs_covered = cov.sum()`; if `< series_min_covered_pairs` → not
  engaged, `gate_mse_reduction = NaN`, reason
  `"only {n} matched pairs sit in a validated homologous series (need {min}) — stage-1 curve only"`.
- `mse0 = mean(r[cov]²)`; `mse1 = mean((r[cov] − pred[cov])²)`.
- `mse0 < 1e-12` → not engaged, `gate_mse_reduction = 0.0`, reason
  `"no correction needed: the stage-1 curve already reproduces the series pairs"`
  (the self-calibration case; there is genuinely nothing to correct, so no
  percentage is reported).
- Otherwise `reduction = 1 − mse1/mse0`, engaged iff
  `reduction ≥ series_gate_min_mse_reduction`. With `p = round(100 ·
  reduction)` and `t = round(100 · threshold)`:
  - engaged: `"series term engaged: leave-own-out MSE {p}% below curve-only (threshold {t}%)"`;
  - gated off, `p ≥ 0`: `"series term gated off: leave-own-out MSE only {p}% below curve-only (threshold {t}%)"`;
  - gated off, `p < 0` (a measured *loss*): `"series term gated off: leave-own-out MSE {−p}% above curve-only (threshold: {t}% below)"`.

### 14.5 Parameters and defaults

| field | default | note |
|---|---|---|
| `use_series_term` | `True` | stage 1b on/off |
| `series_kmd_tol` | `0.008` | tolerance on the Kendrick phase |
| `series_min_members` | `4` | smallest cluster / chain that counts as a series |
| `series_end_reach_ch2` | `2` | how many CH2 past a series end a query may sit |
| `series_exclude_rt_frac` | `0.0015` | co-elution exclusion, fraction of the knot span |
| `series_exclude_rt_floor_min` | `0.03` | ... with this absolute floor, minutes |
| `series_gate_min_mse_reduction` | `0.2` | the gate |
| `series_min_covered_pairs` | `20` | below this the gate declines unread |

### 14.6 New outputs

* `*_calibrated.csv` gains, as its last two columns, `series_correction_min`
  (minutes added to the curve's prediction; 0 when no correction was applied)
  and `series_n_members` (number of series members behind the correction;
  0 = not corrected). `warp_source` is per row: `curve`, `curve+series` (only
  where a correction was actually applied), and with stage 2 engaged
  `curve+anchors` / `curve+series+anchors` (`+anchors` on every row).
* `*_pairs.csv` gains `in_series` — the pair sits in a validated series
  (False for trimmed pairs and non-members; all False when the term is
  disabled).
* `*_model.json` gains a `series` block, always present: `enabled`,
  `engaged`, `gate_mse_reduction`, `gate_threshold`, `gate_reason`,
  `n_series`, `n_pairs_covered`, `n_features_corrected`, `kmd_tol`,
  `min_members`, `end_reach_ch2`, `exclude_rt_min`. When the term is disabled:
  `enabled=False`, `engaged=False`, `gate_reason="disabled
  (use_series_term=False)"`, counts 0, reduction null.
* The log gains, right after the stage-1 line,
  `"series term: {gate_reason}; {n_series} homologous series ({n_in_series} pairs, {n_pairs_covered} with a leave-own-out prediction); {k}/{n} features corrected"`,
  where `n_in_series` is the number of pairs in a validated series
  (`int(member_mask.sum())`, §14.4), or
  `"series term: disabled (use_series_term=False) — stage-1 curve only"`.
* The report: the fourth KPI tile is "Series term" (with the anchor gate
  state appended to its note when stage 2 was requested), the radar axis
  "Anchors" is now "Series term", the method-facts table gains a "Series
  term" row, and the curve figure annotates the series-term decision.
  (Since 1.2.3 the tile and the axis read "Series + lattice", §15.7.)
* API: `rt_anchor` exports `SeriesTerm`, `kendrick_phase`,
  `KENDRICK_FACTOR`, `CH2_MASS`, `PERIOD`. `ColumnCalibrator.predict(rt,
  classes=None, mz=None)` without `mz` returns the curve (plus anchors) only.

### 14.7 Interaction with stage 2, and the `use_sample_anchors` default change

Stage 2 is **opt-in** since 1.2.2: `use_sample_anchors` defaults to `False`;
CLI `--sample-anchors` turns it on; `--no-sample-anchors` is still accepted
and wins. When stage 2 runs, the anchor residuals — and with them the sanity
filter, the `(lam_g, lam_c)` search, the LOO residuals and the gate — are
measured against `stage1_predict`: curve **plus the series term** when the
term engaged and the anchor frame carries the anchors' source m/z (`mz_src`,
which the pipeline's anchor table always does). The series term's work is
never "corrected" a second time. The anchor m/z mask on the sample pairs
(stage 1, step 2) is unchanged, and applies whether or not stage 2 runs.
(Since 1.2.3 `stage1_predict` adds the applied lattice correction too, §15.8.)

Reproducing 1.2.1: `use_series_term=False, use_sample_anchors=True` (CLI
`--no-series-term --sample-anchors`) gives v1.2.1's default output;
`use_series_term=False` alone gives v1.2.1's `--no-sample-anchors` output.
Since 1.2.3 both recipes also need `use_lattice_term=False` (CLI
`--no-lattice-term`, §15).

### 14.8 Per-sample tier

When `single_files` are given, each injection fits its own series term from
its own kept pairs, with the same config and its own gate; the injection's
prediction is its curve plus its series correction. Rows no injection covers
fall back to the project-level prediction. The reported
`series_correction_min` is the median over injections of the applied
corrections and `series_n_members` the maximum member count. (Since 1.2.3 each
injection also fits its own lattice term behind its series term, §15.9.)

### 14.9 What it does not do

* The series are recognised from m/z alone (the Kendrick phase) and are
  **not identifications**; nothing is claimed about a feature's class or
  identity.
* Coverage is partial by construction: a feature outside every validated
  series keeps the stage-1 curve, and on a matrix with few homologous series
  the gate will simply not engage. On the shipped Mix 15 example (Orbitrap),
  27 validated series holding 182 of the 433 kept pairs were found; the gate,
  read on the 148 of them that get a leave-own-out prediction, **declined** —
  the leave-own-out error would have been 35% *higher* than the curve alone —
  so the term is not applied; that is the designed behaviour on data where the
  series carry no consistent selectivity difference, not a failure. (Since
  1.2.3 a feature the series term leaves alone may be corrected by the lattice
  term behind it, §15.)
* The term depends on m/z *and* RT, so it is not a curve: the curve figure
  keeps drawing the stage-1 curve, and `predict(rt)` without `mz` returns the
  curve (plus anchors) only, exactly as before.
* It never corrects a feature with its co-eluting family members
  (leave-own-family-out), and it extrapolates at most
  `series_end_reach_ch2` CH2 past a series end.

### 14.10 Measured evidence

On a six-method human-serum validation set (five source methods calibrated
onto the reference column at the default 15 ppm matching window; 134 lipids
located from exact mass independently
of the calibration, 576 lipid-by-method points), the series term lowered the
median |Cal_RT − reference RT| from 0.106 to 0.064 min and the 90th
percentile from 0.341 to 0.271 min; the share of points within 0.2 min rose
from 74% to 82%. The gate engaged on four methods and declined on the one
whose selectivity already matched the reference. Validated series covered
about half of the points; a feature outside any validated series keeps the
stage-1 curve. Every prediction is leave-own-family-out by construction, so
these are not in-sample figures.

---

## 15. v1.2.3 — the lattice term

Stage 1c, fitted right after the series term (stage 1b) and before the
(opt-in) stage 2. Implemented in `lattice.py`; on by default
(`use_lattice_term=True`, CLI `--no-lattice-term` switches it off). It is the
two-dimensional sibling of the series term (§14): the series term interpolates
along one CH2 ladder, the lattice term borrows across ladders that differ by
one double bond.

**Why.** The series term needs at least `series_min_members` (4) members in
one CH2 ladder — same class, same number of double bonds — and some lipid
classes never offer that: serum cholesteryl esters have at most three members
per double-bond count, so the series term never touches them and their
class-wide offset against the stage-1 curve stays in the calibrated RT. Yet the
class as a whole has a dozen members, spread over *several* ladders that differ
by H2 (one double bond). The lattice term uses that: CH2 ladders whose Kendrick
phases are one H2 step apart (two, across a missing row — rule 4 of §15.2) are
joined into a **lattice family**, and every member gets integer coordinates
`(a, b)` — carbons and H2 count, relative to the family's lightest member —
read straight off its mass,
`mz = m0 + a · CH2_MASS + b · H2_MASS`. A feature's correction is then a
distance-weighted mean of the curve residuals of the family's members near its
own lattice point. Like the series term it is identity-free: it works from m/z
alone.

**A fallback behind the series term, never a replacement.** A feature the
series term corrects keeps that correction and is not touched by the lattice
term; the lattice term fills in where the series term is silent, and its gate
is judged only on the pairs the series term leaves alone (§15.5). The series
term interpolates along one ladder, which is the better evidence wherever a
ladder is long enough to trust; the lattice term borrows across ladders where
it is not.

### 15.1 Constants

```
H2_MASS  = CH2_MASS - 12.0                    # 2.01565006: one H2, one double bond
H2_PHASE = H2_MASS * KENDRICK_FACTOR          # one H2 step on the Kendrick circle, ~2.0134
circ(d)  = min(d mod PERIOD, PERIOD - (d mod PERIOD))     # circular distance of a phase difference
wrap(d)  = d - 14 if d > 7;  d + 14 if d <= -7;  else d   # into (-7, 7]
d        = |Δa| / 2 + |Δb|                    # lattice distance between two lattice points
```

`KENDRICK_FACTOR`, `CH2_MASS`, `PERIOD` and `phase(mz)` are §14.1's; `circ(d)`
is §14.1's `circ(p, q)` written for the difference `d = p − q`. `b` counts H2,
so one more `b` is one double bond fewer, and two species one double bond apart
differ by `H2_MASS` in m/z and by `H2_PHASE` on the Kendrick circle. The lattice
distance counts one double bond like two carbons; it is a constant of the
method, not a setting. `lattice_kmd_tol` is applied as it stands on both axes —
to the Kendrick phase (cluster gap, H2 link) and to the m/z error in Da (member
error, lattice-point match); the implementation does not rescale between them.

### 15.2 Fit

Input: the stage-1 pairs that survived MAD trimming — source m/z `mz`, source
RT `xa`, reference RT `xb`, curve residual `r = xb − f(xa)` — plus the curve's
source-RT knot span `x0, x1` and `series_covered`, a boolean per pair that is
True where the **engaged** series term gives that pair a leave-own-out
correction (`series.correction(mz, xa)[1] > 0`; all False when the series term
is disabled or gated off). Pairs in `series_covered` may sit in a family and
donate; only the gate (§15.4) skips them.

1. Pairs with any non-finite `mz`, `xa`, `xb` or `r` are dropped (with the same
   entries of `series_covered`) and are never members; the rest keep their
   given order, indexed `0..n−1`.
2. `excl = max(series_exclude_rt_floor_min, series_exclude_rt_frac · (x1 − x0))`
   is the co-elution exclusion radius — the series term's (§14.2, rule 2), one
   setting for both terms.
3. **Clusters and centres.** Exactly the series term's clustering (§14.2, rule
   3) with `lattice_kmd_tol`: pair indices are sorted by `(phase, index)`, a new
   cluster starts whenever the phase gap exceeds `lattice_kmd_tol`, and if there
   are at least two clusters and `phase[first of first cluster] + PERIOD −
   phase[last of last cluster] ≤ lattice_kmd_tol`, the last cluster is prepended
   to the first and removed. Unlike the series term, **every** cluster is kept,
   whatever its size: a ladder too short to be a series can still be one row of
   a lattice family. Clusters are numbered `0..K−1` in the resulting order (a
   seam-rejoined cluster first). The **centre** of a cluster is
   `(p0 + median(wrap(phase[members] − p0))) mod PERIOD` over all its members,
   with `p0` the phase of its first member in its stored order — sorted by
   `(phase, index)`; for a seam-rejoined cluster the former last cluster's
   members come first. (The series term's centre is taken over a validated
   chain; there is none yet here.)
4. **H2 links.** `step(u, v)` is the first `s` in `(+1, +2, −1, −2)` with
   `circ(centre[v] − centre[u] − s · H2_PHASE) ≤ lattice_kmd_tol`, else none:
   `±1` links neighbouring ladders of one class, `±2` bridges a ladder whose
   middle row is missing from the data.
5. **Components and the H2 index, breadth-first.** `seen` = all False. For `s0`
   in `0..K−1` ascending, if not seen: mark it seen, `b_cluster[s0] = 0`,
   `queue = [s0]` (FIFO). While the queue is not empty: pop `u` from the front;
   for `v` in `0..K−1` ascending, if `v` is not seen and `step(u, v)` is some
   `s`: mark `v` seen, `b_cluster[v] = b_cluster[u] + s`, append `v` to the
   queue. The clusters visited from one `s0` are a **component**, and a
   cluster's H2 index is the one the first path to it assigned. A component with
   a single cluster is skipped — a lone ladder is the series term's business;
   every other component goes through rules 6–8, in the order the components
   are found.
6. **Members and coordinates.** All pairs of the component's clusters, ordered
   by `(mz, index)` ascending → `o[0..m−1]`. Each has `b` = its cluster's
   `b_cluster`; the first (lightest) member's `b` is subtracted from all. Then

   ```
   root = mz[o[0]]
   a    = rint((mz[o] - root - b * H2_MASS) / CH2_MASS)       # numpy rint, as int
   err  = mz[o] - root - a * CH2_MASS - b * H2_MASS
   keep the members with |err - median(err)| <= lattice_kmd_tol    # order preserved
   ```

   The coordinates are relative to that lightest pair, which can itself be
   dropped here or in rule 7: a family need not contain `(0, 0)`, and `b` can
   be negative.
7. **Order pruning.** On a reversed-phase column more carbons and/or fewer
   double bonds elutes later, and `b` counts H2. For the kept members `i, j`:

   ```
   above(i, j)    = a[i] >= a[j] and b[i] >= b[j] and (a[i] > a[j] or b[i] > b[j])
   conflict(i, j) = above(i, j) and not (xa[i] > xa[j] and xb[i] > xb[j])   # made symmetric
   alive = all members
   loop:
       v[i] = number of alive members in conflict with alive member i
       if max(v) == 0: stop
       cand = the alive members with v == max(v), in member order
       med  = median of r over the alive members
       remove the member of cand with the largest |r - med|        # the first on ties
   ```

   Two members are in conflict when one sits above the other on the lattice but
   does not elute later than it on **both** columns (the comparisons are strict:
   an equal RT is a conflict). Members on the same lattice point never conflict
   with each other — isomers are allowed.
8. **Acceptance.** The family is kept iff at least `lattice_min_members` members
   remain **and** they span at least two distinct `b`.
   `m0 = median(mz − a · CH2_MASS − b · H2_MASS)` over the remaining members. A
   family stores its members (still in ascending `(mz, index)` order) with
   their `a`, `b` and `m0`; families keep the order they were found in.

Tie-breaking: `step` takes the first matching `s` in the order given; the
breadth-first search visits clusters in ascending index; members are in
`(mz, index)` order; the pruning removes the first of the equally far-off
candidates. A pair belongs to at most one family; the pairs of the accepted
families are the **members** (`member_mask`, `in_lattice` in `*_pairs.csv`), and
pairs outside every family are never donors.

### 15.3 Prediction of one query `(mq, xq)`

The **raw** correction — `raw_correction()`, available whether or not the gate
engages — returns `(correction_min, n_members)` per query; "no correction" is
`(0.0, 0)`. `D = lattice_max_distance`; `⌊·⌋` is the floor.

1. Non-finite `mq`/`xq`, or no accepted family → `(0.0, 0)`.
2. **Locate the lattice point.** For each family `k` in order, and each integer
   `bq` from `b_min − ⌊D⌋` to `b_max + ⌊D⌋` ascending (`b_min`, `b_max` over the
   family's members):

   ```
   aq = rint((mq - m0 - bq * H2_MASS) / CH2_MASS)
   skip if aq < a_min - floor(2*D) or aq > a_max + floor(2*D)     # a_min, a_max over the members
   e  = |mq - m0 - aq * CH2_MASS - bq * H2_MASS|
   a candidate if e <= lattice_kmd_tol
   ```

   The query's lattice point is the candidate with the **smallest** `e` over all
   families and rows; on exactly equal `e` the first one found (family order,
   then ascending `bq`). No candidate → `(0.0, 0)`.
3. **Donors.** In that family, `d = |a − aq| / 2 + |b − bq|` for each member.
   The donors are the members with `|xa − xq| > excl` — the
   **leave-own-family-out** rule: a feature is never corrected with pairs that
   co-elute with it — **and** `d ≤ D`, in member order. Fewer than
   `lattice_min_neighbours` donors → `(0.0, 0)`.
4. **Weights.** The `lattice_max_neighbours` donors with the smallest `d` — a
   stable sort, so ties keep member order — are averaged with
   `w = 1 / (d + 0.5)²`: `correction = Σ w · r / Σ w`. `n_members` is the number
   of donors — all of them in reach, not only the ones averaged.

The applied `correction()` is the raw correction when the gate engaged and all
zeros (and `n_members` all zero) otherwise. The location step is vectorised over
the queries per family and H2 row, and the donor step over each family's
queries, so correcting a whole feature table never loops in Python over
features.

### 15.4 The gate

Every fitted pair is predicted as a query at its own `(mz, xa)` — the
co-elution exclusion turns that into a leave-own-family-out prediction —
giving `pred[i]`, `n[i]`; `cov = (n > 0) and not series_covered`: the pairs the
lattice term would actually correct.

Two counts describe the coverage, and neither bounds the other: the pairs in an
accepted family (`in_lattice` in `*_pairs.csv`, `n_pairs_in_families` in the
model) and the covered pairs (`n_pairs_covered`), which the gate is read on. A
member can be uncovered — too few donors are left once its co-eluting neighbours
are excluded, or the engaged series term already corrects it — and a pair that
is not a member (pruned, or outside the member error band) can still sit on a
family's lattice point with enough donors, and so be covered.

- `n_pairs_covered = cov.sum()`; if `< lattice_min_covered_pairs` → not
  engaged, `gate_mse_reduction = NaN`, reason
  `"only {n} matched pairs the series term leaves alone sit in a lattice family (need {min}) — no lattice correction"`.
- `mse0 = mean(r[cov]²)`; `mse1 = mean((r[cov] − pred[cov])²)`.
- `mse0 < 1e-12` → not engaged, `gate_mse_reduction = 0.0`, reason
  `"no correction needed: the stage-1 curve already reproduces the lattice pairs"`
  (the self-calibration case; there is genuinely nothing to correct, so no
  percentage is reported).
- Otherwise `reduction = 1 − mse1/mse0`, engaged iff
  `reduction ≥ lattice_gate_min_mse_reduction`. With `p = round(100 ·
  reduction)` and `t = round(100 · threshold)`:
  - engaged: `"lattice term engaged: leave-own-out MSE {p}% below curve-only (threshold {t}%)"`;
  - gated off, `p ≥ 0`: `"lattice term gated off: leave-own-out MSE only {p}% below curve-only (threshold {t}%)"`;
  - gated off, `p < 0` (a measured *loss*): `"lattice term gated off: leave-own-out MSE {−p}% above curve-only (threshold: {t}% below)"`.

The percentages are rounded for the reason string only; the decision is made on
the unrounded reduction, so a reduction just under the threshold can read "only
20% below curve-only (threshold 20%)".

### 15.5 Composition with the series term

The calibrator composes the two terms, not the lattice object. For a feature
`(rt, mz)`:

```
s_corr, s_n = series.correction(mz, rt)     # zeros when the series term is absent or gated off
l_corr, l_n = lattice.correction(mz, rt)    # zeros when the lattice term is absent or gated off
where s_n > 0:  l_corr = 0.0, l_n = 0       # the series term has precedence
Cal_RT (anchor-free part) = curve(rt) + s_corr + l_corr
```

`ColumnCalibrator.stage1_terms(rt, mz)` returns the **applied**
`(series_corr, series_n, lattice_corr, lattice_n)`: zeros for a term that is
absent or gated off, and the lattice values zeroed wherever `series_n > 0`.
`stage1_predict(rt, mz=None)` is the curve plus both applied corrections when
`mz` is given, and the bare curve without it, exactly as before the terms
existed; `predict(rt, classes=None, mz=None)` adds the stage-2 refiner on top.
`lattice_correction_min` and `lattice_n_members` in the output are the applied
values, so in the per-project tier a row never carries both a series and a
lattice correction.

The same rule sets `series_covered` at fit time. After the series term is fitted
on the kept pairs, the lattice term is fitted on the same kept pairs and curve
residuals with `series_covered = series.correction(mz_kept, rt_src_kept)[1] > 0`
— all False when there is no series term (disabled) or it was gated off. If the
series term corrects every pair the lattice term would cover,
`n_pairs_covered = 0` and the lattice gate cannot engage.

### 15.6 Parameters and defaults

| field | default | note |
|---|---|---|
| `use_lattice_term` | `True` | stage 1c on/off |
| `lattice_kmd_tol` | `0.008` | tolerance: cluster gap, H2 link, member error, lattice-point match |
| `lattice_min_members` | `6` | smallest lattice family (after pruning) |
| `lattice_max_distance` | `3.0` | how far (lattice steps, §15.1) a correction may be borrowed from |
| `lattice_min_neighbours` | `2` | donors needed within that distance; clamped to at least 1 |
| `lattice_max_neighbours` | `4` | nearest donors averaged; clamped to at least 1 |
| `lattice_gate_min_mse_reduction` | `0.2` | the gate |
| `lattice_min_covered_pairs` | `20` | below this the gate declines unread |

A parameter is an explicit override, else the config value, else the default in
`crosscolumn.DEFAULTS` — the series term's resolution. The co-elution exclusion
is the series term's `series_exclude_rt_frac` and `series_exclude_rt_floor_min`
(§14.5): one notion of "co-eluting", one setting. The lattice distance (§15.1)
is not a setting. `lattice_min_neighbours` and `lattice_max_neighbours` are
clamped to at least 1 where the term stores them: a correction needs at least
one donor, and a weighted mean over none would be 0/0; the defaults are
unaffected, and the model block of a fitted term reports the values it ran
with. `from_dict` keeps rejecting unknown keys.

### 15.7 New outputs

* `*_calibrated.csv` gains, as its last two columns (after the series ones),
  `lattice_correction_min` (minutes added to the curve's prediction; 0 when no
  correction was applied) and `lattice_n_members` (the number of donors behind
  the correction — family members in reach that do not co-elute with the
  feature; 0 = not corrected). `warp_source` is per row: `curve`, then
  `+series` where `series_n_members > 0`, then `+lattice` where
  `lattice_n_members > 0`, then `+anchors` on every row when stage 2 engaged —
  so `curve+lattice`, and with stage 2 engaged `curve+lattice+anchors`. In the
  per-project tier a row has at most one of `+series` and `+lattice` (§15.5); in
  the per-sample tier, where both columns aggregate over injections, a row may
  read `curve+series+lattice` (§15.9). At run level
  `ColumnCalibrator.warp_source` is `curve` + `+series` + `+lattice` +
  `+anchors`, each part only when that stage is engaged, in that order.
* `*_pairs.csv` gains `in_lattice`, after `in_series` — the pair is a member of
  an accepted lattice family (False for trimmed pairs and non-members; all False
  when the term is disabled).
* `*_model.json` gains a `lattice` block, always present, right after `series`:
  `enabled`, `engaged`, `gate_mse_reduction`, `gate_threshold`, `gate_reason`,
  `n_families`, `n_pairs_in_families`, `n_pairs_covered`, `n_features_corrected`,
  `kmd_tol`, `min_members`, `max_distance`, `min_neighbours`, `max_neighbours`,
  `exclude_rt_min`. `n_pairs_in_families` is the number of member pairs
  (`member_mask.sum()`), `n_pairs_covered` the pairs the gate was read on
  (§15.4), `n_features_corrected` the number of rows of the result table with
  `lattice_n_members > 0`. When the term is disabled: `enabled=False`,
  `engaged=False`, `gate_reason="disabled (use_lattice_term=False)"`, counts 0,
  reduction null, the parameters from the config. The `config` block records
  the new fields like any other.
* The log gains, right after the series-term line,
  `"lattice term: {gate_reason}; {n_families} lattice families ({n_in_families} pairs, {n_pairs_covered} judged by the gate); {k}/{n} features corrected"`,
  where `n_in_families` is `int(member_mask.sum())` and `k` the number of the
  sample table's features that receive a lattice correction, or
  `"lattice term: disabled (use_lattice_term=False)"`.
* The report and the app (§12). Report: the fourth KPI tile and the radar axis
  are "Series + lattice", the method-facts table gains a "Lattice term" row
  right after "Series term", and the curve figure of the PDF prints the
  lattice-term decision on the line after the series-term one. App: the radar
  axis is "Series + lattice" too (the app does not show the tile), the Curve
  section lists the lattice-term reason after the series-term one, and the app
  gains the Advanced checkbox "Lattice term" (on by default), the table-preview
  columns `lattice corr (min)` and `lattice n`, and the `results.lattice` block
  of `run_info.json`. Both read the model's `lattice` block; neither recomputes
  it.
* API: `rt_anchor` exports `LatticeTerm` and `H2_MASS`. `LatticeTerm.fit(mz,
  rt_src, rt_ref, residual, x0, x1, config=None, series_covered=None,
  **overrides)` builds the term; it exposes `engaged`, `gate_mse_reduction`,
  `gate_threshold`, `gate_reason`, `n_families`, `family_sizes`,
  `family_offsets` (`m0` per family), `n_pairs_covered`, `exclude_rt_min`,
  `member_mask` (aligned with the arrays passed to `fit`), `raw_correction`,
  `correction` and `to_model()` (the fit-level part of the model block; the
  pipeline adds `n_pairs_in_families` and `n_features_corrected`).
  `ColumnCalibrator` gains `lattice`, `lattice_used` and `stage1_terms(rt, mz)`;
  `CalibrationResult` gains `lattice_used`.
* Unchanged: `Cal_RT_uncertainty_min` (§5.2) and `is_extrapolated` do not
  depend on the lattice term.

### 15.8 Interaction with stage 2

Stage 2 stays opt-in (§14.7). When it runs, the anchor residuals — and with
them the sanity filter, the `(lam_g, lam_c)` search, the LOO residuals and the
gate — are measured against `stage1_predict`: the curve **plus the applied
series and lattice corrections** when at least one of the two terms engaged and
the anchor frame carries the anchors' source m/z (`mz_src`, which the
pipeline's anchor table always does); otherwise the bare curve, exactly as
before 1.2.2. The terms' work is never "corrected" a second time, and
`residual_min` in `*_anchors.csv` is on the same basis (§5). The anchor m/z mask
on the sample pairs (stage 1, step 2) is unchanged.

Reproducing 1.2.2: `use_lattice_term=False` (CLI `--no-lattice-term`) gives
v1.2.2's output; the test suite pins this against numeric snapshots of 1.2.2
runs. With `use_series_term=False` as well, and `use_sample_anchors=True`, it is
v1.2.1's default output (§14.7).

### 15.9 Per-sample tier

When `single_files` are given, each injection fits its own series term (when
enabled) from its own kept pairs and then its own lattice term (when enabled)
from the same kept pairs and residuals, with the same config, its own
`series_covered` (the pairs its own series term corrects; all False when that
term is disabled or gated off) and its own gate. The injection's prediction is
its curve plus its applied series and lattice corrections, under the same
precedence: the lattice correction is zeroed wherever that injection's series
term corrected the feature. The reported `lattice_correction_min` is the median
over injections of the applied lattice corrections and `lattice_n_members` the
maximum member count — built exactly like the series columns — so a row may
show both terms (`curve+series+lattice`), because each pair of columns
aggregates over injections on its own. Rows no injection covers fall back to the
project-level values (`stage1_terms` of the project calibrator).

In this tier `model["lattice"]` is the project-level term's fit, while its
`n_features_corrected` counts the rows of the final table with
`lattice_n_members > 0`; the log line's `{k}` counts the project-level
calibrator's corrections.

### 15.10 What it does not do

* The families are recognised from m/z alone — the Kendrick phase and H2 steps —
  and are **not identifications**; nothing is claimed about a feature's class or
  identity.
* A lattice belongs to a formula series, not to a lipid class. A species from
  another class that shares the family's formula series lands on the same
  lattice and is corrected as if it belonged: LPC and ether PC do, and so does
  the sodium adduct of a lipid two carbons shorter and three double bonds more
  saturated (2.4 mDa from the protonated lattice point, inside
  `lattice_kmd_tol`). This is the known weak spot (see §15.11, LPC).
* Coverage is partial by construction. A feature is not corrected by the term if
  its m/z is not on a family's lattice point, if fewer than
  `lattice_min_neighbours` family members within `lattice_max_distance` lattice
  steps are left once its co-eluting neighbours are excluded, or if the series
  term corrects it; and on a matrix or method where the families carry no
  consistent selectivity difference the gate does not engage. On the shipped
  Mix 15 example (Orbitrap) both gates decline — the designed behaviour on data
  where neither the series nor the lattice families carry a consistent
  selectivity difference, not a failure (the figures are in `examples/README.md`).
* The term depends on m/z *and* RT, so it is not a curve: the curve figure keeps
  drawing the stage-1 curve, and `predict(rt)` without `mz` returns the curve
  (plus anchors) only, exactly as before.
* It never corrects a feature with its co-eluting family members
  (leave-own-family-out), borrows only from members within
  `lattice_max_distance` lattice steps, and never replaces a series correction.
* Its distance metric, tolerance and reach are fixed, not tuned per data set
  (§10, item 11).

### 15.11 Measured evidence

On the same six-method human-serum validation set that was used for the series
term (five source methods calibrated onto the reference column at the default
15 ppm matching window; 134 lipids located from exact mass independently of the
calibration, 576 lipid-by-method points), adding the lattice term to 1.2.2
lowered the median |Cal_RT − reference RT| from 0.064 to 0.059 min and the 90th
percentile from 0.271 to 0.218 min; the share of points within 0.2 min rose from
82% to 87%. The lattice gate engaged on the same four methods as the series gate
(leave-own-out MSE 46–51% below the curve alone, on the pairs the series term
leaves alone) and declined on the one whose selectivity already matched the
reference. On those four methods the series term corrected 301 of the points and
the lattice term 139 more. The gain sits where the series term could not reach:
cholesteryl esters (55 points, none of them touched by the series term) went
from a median of 0.175 to 0.066 min, from 64% to 87% within 0.2 min, and their
class-wide offset (median signed error) from −0.147 to −0.004 min. Point by
point, 55 improved by more than 0.1 min and 8 got worse by more than 0.1 min.

Per class (median |Cal_RT − reference RT| in min, and the share within 0.2 min,
1.2.2 → 1.2.3):

| class | points | median | within 0.2 min |
|---|---|---|---|
| CE | 55 | 0.175 → 0.066 | 64% → 87% |
| DG | 16 | 0.093 → 0.090 | 81% → 94% |
| LPC | 52 | 0.055 → 0.074 | 83% → 79% |
| PC | 110 | 0.098 → 0.080 | 76% → 81% |
| ether PC | 57 | 0.066 → 0.052 | 70% → 79% |
| SM | 97 | 0.047 → 0.047 | 92% → 93% |
| TG | 189 | 0.058 → 0.055 | 89% → 92% |

One class got worse: lysophosphatidylcholines (52 points) went from a median of
0.055 to 0.074 min and from 83% to 79% within 0.2 min; five of the eight points
that worsened by more than 0.1 min are LPC. LPC shares its lattice with ether PC
(the same formula series), and the sodium adducts of lighter LPCs land on LPC
lattice points, so its neighbours on the lattice are not all LPC. Every
prediction is leave-own-family-out by construction, so these are not in-sample
figures; but the term's design choices (the neighbour rule and its reach) were
made on this same data set, and it has not been tested on another matrix or
instrument.
