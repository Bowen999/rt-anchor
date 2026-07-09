# Retention‑Index (RI) Calibration — Specification v3

Status: **implementation‑ready draft**. Revised through **two adversarial review rounds** (round‑1: 4 lenses
on the concept; round‑2: verify‑fixes + stress‑test the iRT decision + consistency + a numerical
worked‑example on real data). All config is now locked with defaults; no open design decisions remain.

Scope: convert raw LC‑MS retention time (RT) into a retention index (RI) using a spiked lipid‑standard
panel, **per sample**. Built from first principles + the real `mzML/` data; no reuse of prior code.

Changelog v2→v3: iRT scale redefined from fixed method constants (was anchor‑pinned [0,100] — non‑extensible,
pos/neg‑incomparable, anchored on the worst landmarks); native anchor template promoted to a **required**
input (empirically load‑bearing); §10 split into panel‑CV vs per‑sample deployment error; Step‑0 filter and
anchor weights rewritten in real `single_files` columns; warp model + numeric defaults locked; naming/units
unified (`iRT` throughout).
Changelog v3→v3.1: **`single_files` made optional.** Default calibration is now **per‑project** from the
aligned consensus RT (already MassCube‑aligned across the batch); **per‑sample self‑anchoring auto‑engages**
when per‑sample RT (`single_files`/pkl) is present — adding per‑injection drift correction + the inter‑sample
`RI_spread` QC. Minimal required input is now just the aligned table (`m/z`+`RT`) + a standard source + polarity.

---

## 0. Locked decisions

| # | Decision | Choice |
|---|---|---|
| 1 | Purpose | Portable RI **for annotation/library matching**, scoped to *within‑method* (§1) |
| 2 | RI scale | **Dimensionless iRT** (not minutes), defined by fixed method constants (§3) |
| 3 | Granularity | **Per‑project by default** (from the aligned consensus RT); **auto‑upgrades to per‑sample self‑anchoring** when per‑sample RT (`single_files`/pkl) is present (§2) |
| 4 | Reference axis | **Defined on the Orbitrap reference run** (honest; no platform‑neutral claim) |

`iRT(RT) = 100 · (RT − RT_lo) / (RT_hi − RT_lo)`, with `RT_lo, RT_hi` **frozen method constants**
(default 1.0 / 18.4 min — documented, arbitrary scale labels, *not* pinned to any anchor). One map, **shared
across polarities**, decoupled from the panel roster. The affine is **pure relabeling that does no
methodological work** — its only purpose is to avoid the "minutes‑but‑not‑minutes" footgun; *all* drift
removal and portability come from the warp (§7). 1 min ≈ 5.75 iRT on the default scale.

---

## 1. Principle & scope

RT carries identification info orthogonal to m/z but drifts (column age, temperature, mobile phase) and
differs strongly and **non‑linearly** across instruments (measured Orbitrap→QTOF divergence: 0 → −1.4 →
−3.5 min). A retention index removes this by anchoring RT to reference compounds run under the same
conditions and interpolating.

**Portability boundary.** Warping raw RT onto a frozen index is portable across runs sharing the **same
chromatography** (column, gradient, solvents, flow, temperature, dead time); it is **robust to
instrument/detector, column age, and batch** (the Orbitrap↔QTOF evidence here is cross‑*detector* on the
same nominal LC method). It is **not** portable to a different gradient — the fundamental RI limitation
(Kováts indices transfer only within a stationary phase). **The reference method's gradient table must
travel with the library**; RI must not be compared to raw‑RT libraries from other methods without
recalibration.

The 15 standards are **RT landmarks only** (class‑agnostic). They are *not* a clean homologous series
(PC 6:0…14:0 are saturated diacyl PCs with uneven carbon steps; the upper five vary carbon number **and**
double‑bond count, ordered by ECN = C − 2·DB), so the Kováts/"homologue ladder" analogy is dropped.

---

## 2. Data model & I/O plumbing

MassCube emits, per project (method × polarity):
- `aligned_feature_table.txt` — one row per aligned feature; **one consensus `RT`** (MassCube **RT‑corrected**
  space; see `project_files/rt_correction_models.pkl`); per‑sample columns are **intensities only**. Also
  carries `isotope_state`, `detection_rate`, `detection_rate_gap_filled` (aligned‑stage only).
- `single_files/<sample>.txt` — per‑sample peak list with a **real per‑sample `RT`** (raw, pre‑alignment),
  `RT_start/RT_end`, and per‑peak `is_isotope`, `is_in_source_fragment`, `charge`, `Gaussian_similarity`,
  `asymmetry_factor`, `noise_score`, `peak_shape`. **No** `isotope_state`/`detection_rate`/gap‑fill columns —
  filters (§5/§6) must use the per‑file columns that exist. Single‑file rows are genuine detections by
  construction (gap‑filling is an alignment‑stage concept, not applicable here).
- `project_files/aligned_features.pkl` — MassCube objects linking each aligned feature to its contributing
  per‑file peaks (requires `import masscube`).
- `project_files/sample_table_with_time.csv` — injection timestamps (injection‑order drift QC).

**Two calibration tiers (auto‑selected by what is provided).**

**Default — per‑project (aligned table only).** The consensus `RT` is already MassCube RT‑corrected/aligned
across the batch, so it is a valid single axis. Identify anchors in the aligned table by m/z (→ consensus RT),
fit **one** warp `consensus_RT → iRT` (§7) **in corrected consensus space**, and apply it to every feature's
consensus RT → **one `RI` per feature**. No per‑sample join or aggregation. This removes the cross‑batch /
cross‑instrument offset (the main goal — the 3.5 min Orbitrap↔QTOF divergence is cross‑*project* and is
absorbed here), but **not** residual within‑batch per‑injection drift beyond what MassCube alignment already
did, and yields no inter‑sample `RI_spread` QC.

**Upgrade — per‑sample self‑anchoring (auto‑engaged when `single_files/<sample>.txt` OR `aligned_features.pkl`
is present).** Per‑sample RT is then available, so:
1. Read each anchor's and feature's **per‑sample RT** (from `single_files`, raw space; or per‑file RTs in the pkl).
2. Fit + apply the warp **entirely in raw per‑file RT space** (never mixed with the corrected consensus RT —
   MassCube's RT correction is upstream and must not be double‑counted).
3. **Join** per‑sample peaks → aligned features: `aligned_features.pkl` per‑file indices primary; m/z+RT
   fallback *deliberately crosses RT spaces*, so use a **≥0.2 min** window (median offset ~0.02, p90 ~0.10 min)
   or map raw RT through `rt_correction_models.pkl` first. Flag fallback‑only joins.
4. **Aggregate** per‑sample iRTs → each feature: inverse‑variance / anchor‑support‑weighted (median fallback);
   emit `RI_spread` (MAD) + `n_contributing`. Never warp the consensus RT through one sample's warp.

Both tiers set `calibration_scope ∈ {project, sample}` in the output.

Assumption: standards are spiked into every sample (verified: **13/15** in all 6 Orbitrap samples). The 2
not consistently present in Orbitrap samples are **PC 6:0_6:0 and PC 11:0_11:0** (their effective low‑end
anchor is PC 9:0 ~1.8 min; earlier‑eluting features are extrapolated). RT unit = minutes — **assert at load**.

---

## 3. iRT scale & reference axis (locked)

**Scale.** `iRT[anchor] = iRT(RT_ref[anchor])` using the §0 fixed‑constant map. Dimensionless, method‑defined,
**shared across polarities** (so a feature at a given RT gets the same iRT in pos and neg). Void features may be
<0, super‑late >100 — legal and explicitly flagged.

**Reference RTs (`RT_ref`).** From the **Orbitrap reference panel run** (xlsx "Exp RT Orbitrap"; "Exp RT
Agilent" is empty). No platform‑neutral claim: RI is the reference LC method **as realized on the Orbitrap
panel**; every run — incl. Orbitrap — warps onto it. (A two‑platform mean was available — QTOF‑native anchor
RTs were extracted — but is not used, by decision #4.)

**Uncertainty & provenance.** Build `RT_ref` + per‑anchor **σ** from the **replicate** panel injections (≥2)
— do not treat a single injection as exact. Reference RTs are quoted to 0.01 min and pos/neg agree to ~0.03 min
→ **accuracy floor ≈ 0.01–0.03 min (≈ 0.06–0.17 iRT)**; never report finer. **Version + checksum** the frozen
`(anchor → RT_ref, σ, iRT; RT_lo, RT_hi)` table with provenance — it is a **required** artifact (§9), not optional.

**Local recalibration (reusing the library on a new same‑method run).** Detect the spiked anchors, fit the
σ‑weighted monotone warp `measured_RT → iRT` (§7), convert. Assumes same‑gradient chromatography **and**
reference‑anchor stability (QC each run). For the strongly non‑linear cross‑detector case, a 2‑point line is
insufficient — the full monotone warp is required.

---

## 4. Standard panel (the 15 landmarks)

| Class | Standards | Pos adduct | Neg adduct | RT_ref range (min) | Notes |
|---|---|---|---|---|---|
| PC ×6 (synthetic) | 6:0_6:0 … 14:0_14:0 | [M+H]+ | [M+HCOO]− | 1.06 – 5.63 | clean; **scale‑defining pair drawn from here** |
| PC ×5 (native‑overlapping) | 16:0_18:3 … 18:0_18:0 | [M+H]+ | [M+HCOO]− | 6.63 – 11.07 | endogenous → isomer‑bias; down‑weight |
| DG ×2 | 16:0_18:0, 18:0_18:0 | [M+NH4]+ | [M+HCOO]− | 12.65 – 14.99 | endogenous + TG‑ISF risk |
| CE ×2 | 20:5, 18:1 | [M+NH4]+ | **none (neg)** | 17.37 – 18.39 | broad tail; not scale‑defining |
| t0 marker | *unretained* | — | — | ~void | **recommend adding** (now safe — scale is roster‑decoupled) |

Verified detection: pos **15/15** both instruments; neg **13** (CE genuinely not neg‑ionizable — the only
apparent CE‑neg hit is the DG 18:0_18:0 [M+HCOO]− collision at ~10 ppm, QTOF‑resolution only). Target m/z from
elemental formulas with **correct cation/anion masses** (e.g. [M+H]+ = M + 1.007276 − mₑ), not the xlsx
"Added Mass" (neutral atom masses → +1.2 ppm bias). Scale endpoints are **not** pinned to panel extremes (§0/§3).

---

## 5. Standard identification — panel run (native template; **required output**)

**Step 0 — filter (per‑file columns only):** keep `is_isotope == False` (M+0), `is_in_source_fragment == False`,
`charge == 1`; enforce a peak‑quality floor (`Gaussian_similarity ≥ 0.7`, `asymmetry_factor ∈ [0.5, 2]` — tunable).

For each target (standard × expected adduct at the run's polarity):
1. **m/z match** within ppm (Orbitrap 8; QTOF 15 — tightened from 30 because 2‑Da‑spaced PCs' ¹³C₂
   isotopologue overlaps the next PC's monoisotopic target at ~11.8 ppm).
2. **Disambiguate** by peak quality + intensity + **polarity‑keyed MS2**.
3. Aggregate the **replicate** injections → per‑anchor RT = median, σ = MAD.

Output = the **method‑native anchor template** (RT + σ per anchor, per method × polarity). This is a **hard
prerequisite** (§9): the worked example showed that without it Orbitrap mis‑picks anchors (isobaric
interference at 10–12 min out‑intensities the true low‑RT anchor) and QTOF recovers only 7/15 (a
reference‑RT‑centered window misses the 3.5‑min‑shifted DG/CE).

**Polarity‑keyed MS2 confirmation:** PC pos [M+H]+ → **184.073** (phosphocholine); PC neg [M+HCOO]− →
**[M−CH3]−** + sn‑1/sn‑2 **RCOO⁻**; DG pos [M+NH4]+ → **[M+H−RCOOH]+** acyl losses; DG neg [M+HCOO]− → RCOO⁻;
CE pos [M+NH4]+ → **369.352** (cholestene); CE neg → n/a.

---

## 6. Standard identification — per sample

For each `single_files/<sample>.txt`, each anchor: apply the Step‑0 per‑file filter, m/z match, then search in
an RT window **centred on the §5 native template RT** (default Orbitrap ±0.3, QTOF ±0.5 min) — *not* the
reference RT (essential on compressing methods). Deterministic tie‑break: (1) MS2‑confirmed, (2) not
isotope/ISF, (3) closest RT to the native template, (4) highest intensity, (5) lowest feature_ID.

- **Down‑weight** endogenous anchors (§4) and void‑region anchors (PC 6:0/9:0 near t0); require MS2/peak‑shape
  agreement with the panel template before accepting an endogenous anchor.
- Record per anchor: found?, sample RT, intensity, quality, MS2‑confirmed?, weight; log zero‑candidate anchors;
  report the **effective first *retained* anchor** per sample.

---

## 7. Calibration model (per‑sample robust monotone warp) — locked

Anchor set: `(RT_sample_i, iRT_i, weight_i)`, `weight_i ∝ Gaussian_similarity / asymmetry_factor` (per‑file)
`× 1/σ_ref_i²` (replicate reference σ). No `detection_rate` term (it is a cross‑sample statistic, unavailable
per‑file).

1. **Robust outlier/monotonicity** (replaces greedy scan): keep the **longest strictly‑increasing subsequence**
   (minimal deletions; tie epsilon 0.005 min), cap deletions at 3, re‑check coverage after. Report every
   dropped/down‑weighted anchor with a reason.
2. **Fit** `f: RT_sample → iRT`, **primary = weighted shape‑constrained monotone smoother** (monotone I‑spline /
   `scam` GAM, or weighted isotonic regression with IRLS Huber reweighting) — does **not** interpolate noisy
   anchors exactly. **Fallback = PCHIP with `extrapolate=False`** when anchors < 8 or the smoother is
   unavailable. Enforce a **minimum positive derivative** (RT→iRT invertible; guard duplicate/near‑equal iRT).
   *Implementation trap:* `scipy PchipInterpolator` defaults to `extrapolate=True` — must set `False` or mask.
3. **Partial pooling** across a batch's samples: fit a **shared warp shape with small per‑sample offsets**
   (hierarchical shape + per‑injection shift) rather than N fully‑independent nonlinear warps — keeps the drift
   correction while stopping low‑anchor samples from manufacturing spread at the ends.
4. **No unbounded extrapolation:** beyond `[anchor_min, anchor_max]` → `iRT = NaN`, `reason = beyond_anchor_support`.
   If a value is required: robust terminal slope over last k ≥ 3 anchors, capped at `max_extrapolation = 1 min`.
5. **Coverage gate** (not a bare count): require the feature's RT within the anchor span with local inter‑anchor
   gap ≤ **3 min**; min anchors ≥ **6** with coverage, else piecewise‑linear low‑confidence; < **3** → uncalibrated.

**Fallback for anchor‑poor samples:** 2–4 anchors → self‑anchor weighted piecewise‑linear (low conf); < 2 → mark
uncalibrated. If borrowing is unavoidable, borrow the **median warp of sibling biological samples** (shared
matrix + injection sequence), tagged `warp_source = batch`, wider σ_RI — never the standards run.

---

## 8. RI output + numeric uncertainty

Per aligned feature (after §2 aggregation):
- `RI` — the dimensionless **iRT** (§3); `calibration_scope ∈ {project, sample}` (which tier ran, §2).
  `RI_spread`, `n_contributing` — cross‑sample QC, **populated only in per‑sample mode** (NaN under per‑project).
- `sigma_RI` — numeric predicted iRT s.e. = `sqrt(σ_local² + σ_scale²)`, where
  `σ_local(RT) = SLOPE · max( LOO_resid_interp(RT), k_gap · nearest_anchor_gap(RT) )` (RT‑space residual → iRT
  via `SLOPE = 100/(RT_hi−RT_lo)`; `k_gap` default 0.1), plus **`σ_scale`** = the frozen‑map covariance from the
  replicated defining anchors (an approximately global term, reported separately in provenance). A monotone GP,
  if used, gives `σ_local` natively as posterior sd.
- `RI_confidence` ∈ {high, medium, low} from `sigma_RI` thresholds **< 0.6 / 0.6–3 / > 3 iRT** (≈ <0.1 / 0.5 min),
  **not** hardcoded RT regions.
- Provenance flags: `is_extrapolated`, `warp_source ∈ {self, batch}`, `join_source ∈ {pkl, mz_rt_fallback}`.

---

## 9. Inputs & outputs

### A. Must be provided by the user (no default possible — experiment‑specific)
- **Sample data to calibrate** — `aligned_feature_table.txt`; **required values per feature: `m/z` and `RT`**
  (consensus RT) — nothing else strictly required. This alone runs the default **per‑project** calibration (§2).
- **Standard source** — the standard‑panel run's `aligned_feature_table.txt` (**required: each standard's
  `m/z` + `RT`**); or the standards as spiked into the samples' own aligned table. Two manifest modes:
  - **Mode A (file only):** identities resolved against the **built‑in default manifest** (§B).
  - **Mode B (file + manifest):** a user manifest overrides the default (custom/extended mixes).
- **Polarity** of each project.

### B. Necessary but with built‑in defaults (override optional)
- **Standard manifest** — default = the 15‑lipid panel. Required per row: `name`, `class`, `formula` (→
  formula‑derived m/z), `adduct_pos`, `adduct_neg` (or `ND`), endogenous/void flags.
- **Reference baseline** (the frozen "ruler") — default = `RT_ref` + σ per anchor and scale constants
  `RT_lo/RT_hi`, built from this method's Orbitrap panel run, versioned + checksummed (§3). **Tied to the
  reference LC method** — rebuild for a different gradient.
- **Native anchor template** — built automatically from the supplied standard‑panel run; the reference‑baseline
  RTs seed its search windows (falls back to the reference RTs if no panel run is given — safe only for small
  within‑method drift, not large cross‑detector shifts).
- **Numeric config** — ppm, RT window, quality floors, coverage/extrapolation caps, tie epsilon, confidence
  thresholds, warp model (§12 locked defaults; all override‑able).

### C. Optional but helpful
- **Per‑sample RT** — `single_files/<sample>.txt` **or** `aligned_features.pkl`: presence **auto‑upgrades to
  per‑sample self‑anchoring** (per‑injection drift correction + inter‑sample `RI_spread` QC; §2). Required
  per‑peak values when `single_files` is used: `m/z` + `RT`; quality/isotope/MS2 columns used if present.
- MS2 spectra / class‑diagnostic fragment rules (anchor confirmation); replicate panel injections (per‑anchor σ,
  reference stability); `sample_table_with_time.csv` (injection‑order drift QC); a t0/void marker (early‑RT anchoring).

### Outputs
- **Per‑sample calibration models**: warp + anchor table `(anchor, sample RT, σ, weight, iRT, residual,
  in/outlier, MS2)` + diagnostics.
- **Calibrated aligned feature table**: original + `RI, sigma_RI, RI_confidence, RI_spread, n_contributing,
  is_extrapolated, warp_source, join_source`.
- **QC**: anchor coverage per sample; per‑sample RT→iRT curves overlaid; **per‑feature inter‑sample RI spread**
  (the key self‑anchoring QC); LOO + leave‑one‑sample‑out residuals stratified by region; extrapolated fraction
  per sample; cross‑detector agreement of shared features after calibration.

---

## 10. Validation

Two distinct error estimates — do not conflate them:

**(A) Panel‑run leave‑one‑anchor‑out CV** — interpolation error of an already‑correctly‑picked anchor set.
Optimistic lower bound; the Orbitrap row is additionally near‑circular (axis = Orbitrap) and below the 0.01 min
ruler. Residual (min | iRT):

| combo | anchors | median | 90th | max |
|---|---|---|---|---|
| Orbitrap pos | 15 | 0.009 \| 0.05 | 0.029 \| 0.17 | 0.063 \| 0.36 |
| Orbitrap neg | 13 | 0.017 | 0.121 | 0.243 |
| QTOF pos | 15 | 0.065 | 1.178 | 2.332 |
| QTOF neg | 13 | 0.048 | 0.640 | 1.415 |

**(B) Per‑sample deployment error (leave‑one‑sample‑out / per‑injection)** — the honest number for a real
feature; includes interference, per‑injection drift, missing void anchors. Measured on a real Orbitrap sample:
**median ~0.033 min, max ~0.167 min** (4–20× the panel‑CV row). QTOF per‑sample was ~0.12 min median / ~1.3 min
max (in line with the panel QTOF row). **Gate downstream matching on (B), not (A).**

Honest reading: **sub‑0.1 min in the PC‑dense region (1–11 min)**; **1–2 min in the sparse DG/CE tail** on QTOF
→ σ‑gate/censor the tail. Report accuracy no finer than 0.01 min. Also run: synthetic anchor‑mis‑pick sensitivity;
region‑stratified (not pooled) residuals.

---

## 11. Known limitations & recommendations

1. **Sparse hydrophobic tail** (pos 4 anchors 11–18 min; neg ends at DG ~15 min) → real interpolation error,
   σ‑gated. *Panel improvement (now safe — scale is roster‑decoupled):* add 2–4 anchors 12–18 min (DG/TG; a
   neg‑ionizable hydrophobe — long‑chain FFA/PG/PI — for the **neg** tail, since CE can't).
2. **No void/t0 anchor** → early gradient unconstrained; add an unretained marker (lands at iRT≈0, no rescale).
3. **Axis is Orbitrap‑defined** on one panel run — build `RT_ref` + σ from the replicate injections; QC
   reference‑anchor stability each new run.
4. **Endogenous‑anchor isomer bias, TG‑ISF at DG, ¹³C₂ QTOF collisions** → handled by §5 Step‑0 + weighting;
   synthetic short/odd‑chain PCs are the trustworthy core.
5. **Cross‑method (different gradient) use is out of scope** — RI portable within the method only (§1).

---

## 12. Implementation plan

New standalone Python module (pandas/numpy/scipy; masscube for the join; `scam`/GP optional for the smoother):
1. `panel.py` — formula‑derived m/z per polarity/adduct; frozen `RT_ref`+σ, `RT_lo/RT_hi`, endogenous/void flags;
   versioned+checksummed config.
2. `io_masscube.py` — load single_files (raw RT), aligned table, join via `aligned_features.pkl` (m/z+RT fallback
   with ≥0.2 min window); RT‑unit assertion; RT‑space consistency.
3. `identify.py` — §5 native template (replicate‑aggregated, σ) + §6 per‑sample (per‑file filters, tie‑break).
4. `calibrate.py` — §7 weighted robust monotone warp (**PCHIP `extrapolate=False`**) + partial pooling + §8
   σ_RI/confidence/provenance; per‑sample → §2 weighted aggregate.
5. `qc.py` — §10 panel‑CV + leave‑one‑sample‑out, region‑stratified; RI‑spread + coverage figures.
6. `cli.py` — project(s) + config in → calibrated table + models + QC out.

**Locked defaults:** ppm Orbi 8 / QTOF 15; RT window Orbi ±0.3 / QTOF ±0.5 min; join window ≥0.2 min; quality
floor Gaussian ≥0.7, asymmetry ∈[0.5,2]; coverage gap ≤3 min, min anchors ≥6 (piecewise‑linear 2–4, uncal <2);
extrapolation NaN, terminal cap 1 min; confidence <0.6/0.6–3/>3 iRT; `RT_lo/RT_hi` = 1.0/18.4 min. All override‑able.

---

## 13. Prior‑art positioning

- **iRT** (Escher et al., Proteomics 2012) — spiked standards + warp to a **dimensionless** RT; this design
  adopts it. The dimensionless choice is **documentation hygiene only** (avoids the minutes footgun) — all drift
  removal/portability come from the anchor warp, *not* the affine relabeling.
- **Kováts (1958) / Van den Dool–Kratz (1963)** — GC RI; transfer only within a stationary phase → basis for
  the within‑method scoping.
- **ECN / ACN+DB lipid RTI** — structure‑anchored indexing; the field's route to cross‑method transfer and to
  predicting unspiked lipids, but ill‑defined across mixed classes on one axis (a separate task).
- **XCMS peakgroups LOESS (Smith 2006) / obiwarp (Prince & Marcotte 2006)** — full‑profile / many‑landmark
  alignment; a shared‑feature/profile alignment **refinement layer** on top of the anchor warp could out‑perform
  15 anchors for the strongly non‑linear QTOF case — a future extension.
