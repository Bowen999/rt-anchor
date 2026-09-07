"""End-to-end orchestration of the v2 cross-column calibration.

What :func:`calibrate` does, in the order it does it, and why the order matters:

1. **Resolve the reference.** The v2 method calibrates one column *onto
   another*, so a second axis is always required. By default that is the
   bundled Column 25 pair (:mod:`rt_anchor.reference`); a user may supply
   their own, in which case every number comes out on *their* axis and the model
   JSON records it.
2. **Validate the stage-2 anchor candidates** — the 17 endogenous plasma lipids
   — in the user's sample run *and* in the reference sample run. A lipid that is
   not confidently found in both is not an anchor.
3. **Match the sample runs anonymously** by accurate m/z and mask out anything
   sitting on an anchor candidate's m/z. Those pairs extend the curve into the
   early and late elution regions where a standards mixture is sparse; the mask
   is what keeps the anchors from partly fitting themselves.
4. **Fit the stage-1 curve** from the two *standards* runs plus those sample
   pairs. No identities are used, which is why the method survives an unknown
   or undetectable mixture.
5. **Refine the isomer picks with that curve** — each source-run candidate is
   re-picked as the in-window feature whose *curve-calibrated* RT lands closest
   to the reference-run RT. This step needs the curve, which is why it cannot
   happen before step 4.
6. **Re-fit with stage 2 engaged**, gated: the class-aware correction is applied
   only if its leave-one-anchor-out MSE is at least 20% below doing nothing.
7. **Build the iRT ruler** from the chosen panel's landmarks detected on the
   *reference* standards run, and evaluate it at each feature's ``Cal_RT_min``.

Two tiers, as in v1: **per-project** (one curve for the table) and, when
``single_files`` are given, **per-sample** — each injection gets its own curve
against the reference sample run and the per-feature ``iRT`` is the median with
``RI_spread`` the dispersion across injections.

A non-plasma matrix finding too few anchors is a normal outcome, not an error:
stage 2 is skipped, the stage-1 curve is used, and the model and log say so.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import mixtures
from .config import CalibrationConfig
from .crosscolumn import (
    ColumnCalibrator,
    MonotoneCurve,
    build_calibrator,
    exclude_pairs_near_mz,
    fit_robust_curve,
    match_features_by_mz,
)
from .errors import CalibrationError, ConfigError, RtAnchorError
from .identify import build_native_template
from .io.loader import load_feature_table
from .io.schema import RESULT_COLUMNS, FeatureTable
from .irt import IRTMapper, build_irt
from .panel import Panel, _norm_polarity, build_panel, is_default_manifest
from .plasma_lipids import (
    build_anchor_table,
    plasma_lipid_candidates,
    refine_picks_with_curves,
    summarize_picks,
    validate_candidates,
)
from .reference import load_reference

#: Run keys used inside the plasma-lipid validation tables.
SOURCE_KEY = "source"
REF_KEY = "reference"

#: Columns of ``<prefix>_pairs.csv`` (spec §5).
PAIR_COLUMNS = ["mz_src", "rt_src", "rt_ref", "source", "kept"]

#: Columns of ``<prefix>_anchors.csv`` (spec §5).
ANCHOR_COLUMNS = ["label", "lipid_class", "rt_src", "rt_ref", "residual_min",
                  "loo_residual_min", "n_isomer_candidates", "isomer_rts",
                  "pick_refined", "dropped_by_sanity_filter"]


@dataclass
class CalibrationResult:
    """Everything one calibration produced: the table, the model, the evidence."""

    table: pd.DataFrame                 # original columns + the appended §1 columns
    anchors: pd.DataFrame               # stage-2 anchors, used and sanity-dropped
    pairs: pd.DataFrame                 # stage-1 matched pairs (PAIR_COLUMNS)
    landmarks: pd.DataFrame             # iRT landmarks on the reference standards run
    calibrator: Optional[ColumnCalibrator]   # the fitted stage-1 curve (+ stage-2 refiner)
    irt: Optional[IRTMapper]            # None when no iRT scale could be built
    model: Dict                         # serialisable calibration model (spec §5.1)
    log: List[str] = field(default_factory=list)
    # ---- metadata for downstream visualisation / QC ----
    rt_col: str = ""                    # RT column name in `table`
    rt_unit: str = "min"
    mz_col: str = ""
    sample_cols: List[str] = field(default_factory=list)
    polarity: Optional[str] = None
    source_format: str = ""
    columns: Dict[str, str] = field(default_factory=dict)  # result-col -> actual name
    panel: Dict = field(default_factory=dict)   # detection QC on the user's standards run
    reference: Dict = field(default_factory=dict)   # which reference pair was used
    panel_key: str = ""                 # "mix15" | "mix21" | "none" | "custom"
    default_panel: bool = True          # built-in 15-standard panel (structure hover)

    # ---- convenience ----
    @property
    def curve(self) -> Optional[MonotoneCurve]:
        return None if self.calibrator is None else self.calibrator.curve

    @property
    def anchors_used(self) -> bool:
        """True when the stage-2 correction passed the gate and is applied."""
        return bool(self.calibrator is not None and self.calibrator.anchors_used)

    def to_model_json(self) -> Dict:
        return self.model

    def rt_minutes(self) -> pd.Series:
        rt = pd.to_numeric(self.table[self.rt_col], errors="coerce")
        return rt / 60.0 if self.rt_unit == "sec" else rt

    def col(self, result_col: str) -> str:
        """Actual table column name for a result column (handles collision rename)."""
        return self.columns.get(result_col, result_col)

    def values(self, result_col: str) -> pd.Series:
        """Numeric values of one appended result column."""
        return pd.to_numeric(self.table[self.col(result_col)], errors="coerce")


# ---------------------------------------------------------------------------
# The entry point
# ---------------------------------------------------------------------------

def calibrate(sample_table: str,
              polarity: str,
              standards_table: str,
              panel: str = mixtures.DEFAULT_MIXTURE,
              reference_sample: Optional[str] = None,
              reference_standards: Optional[str] = None,
              reference_key: str = "col35",
              single_files: Optional[List[str]] = None,
              config: Optional[CalibrationConfig] = None,
              manifest: Optional[pd.DataFrame] = None,
              source_format: Optional[str] = None,
              rt_unit: Optional[str] = None) -> CalibrationResult:
    """Calibrate a feature table onto the reference column's time axis.

    ``sample_table`` and ``standards_table`` **must come from the same column and
    gradient** — that is the assumption the whole method rests on. ``panel`` is
    one of ``mix15`` / ``mix21`` / ``none`` and selects the iRT landmark panel
    and the detection QC; ``manifest`` overrides it with a custom standard list.
    ``reference_sample`` / ``reference_standards`` replace the bundled Column 25
    reference pair and must be given together.
    """
    if not standards_table:
        raise ConfigError(
            "A standards run is required. Pass standards_table=<the standard-panel "
            "run from the same column and gradient as the sample>; it is one of the "
            "two runs the cross-column curve is built from. If the standards are "
            "spiked into the sample itself, pass that same table as standards_table."
        )
    cfg = config or CalibrationConfig()
    pol = _norm_polarity(polarity)
    log: List[str] = []

    # ---- inputs -----------------------------------------------------------
    ft = load_feature_table(sample_table, source_format=source_format,
                            polarity=pol, rt_unit=rt_unit or cfg.rt_unit)
    log.append(f"loaded sample table: {ft.summary()}")
    std_ft = load_feature_table(standards_table, source_format=source_format,
                                polarity=pol, rt_unit=rt_unit or cfg.rt_unit)
    log.append(f"loaded standards run: {std_ft.n_features()} features, "
               f"RT {_span(std_ft)}")

    # A user-supplied reference pair may be in the user's own export format; the
    # bundled pair is a plain MS-DIAL export and must not be forced into it.
    custom_ref = reference_sample is not None or reference_standards is not None
    ref = load_reference(key=reference_key, sample=reference_sample,
                         standards=reference_standards, polarity=pol,
                         source_format=source_format if custom_ref else None,
                         rt_unit=(rt_unit or cfg.rt_unit) if custom_ref else None)
    log.append(f"reference: {ref.pair.label} [{ref.pair.key}] — sample "
               f"{ref.sample.n_features()} features, standards "
               f"{ref.standards.n_features()} features")

    panel_key, panel_manifest = _resolve_panel(panel, manifest)
    panel_obj = (None if panel_manifest is None
                 else build_panel(pol, cfg, manifest=panel_manifest))
    log.append(f"panel: {panel_key}"
               + (f" — {len(panel_obj)} standards ionisable in {pol} mode"
                  if panel_obj is not None else " — no panel identity claimed"))

    # ---- stages 1 + 2 -----------------------------------------------------
    cal, anchors, stage_log = _fit_calibrator(ft, std_ft, ref, cfg)
    log.extend(stage_log)

    # ---- the iRT ruler ----------------------------------------------------
    mapper, landmarks, panel_used, irt_reason = _build_irt(ref, panel_obj, cfg, pol,
                                                          panel_key=panel_key)
    log.append(f"iRT ({panel_used}): {irt_reason}")

    # ---- evaluate -------------------------------------------------------—-
    rt = ft.rt_minutes().to_numpy(dtype=float)
    if single_files:
        res, n_inj, tier_log = _per_sample(rt, single_files, ref, cal, cfg, mapper,
                                           pol, source_format, rt_unit)
        log.extend(tier_log)
    else:
        res, n_inj = _per_project(rt, cal, cfg, mapper), 0

    out, colmap = _append_results(ft.df, res)

    # ---- detection QC (reporting only — it does not drive the calibration) -
    qc = _detection_qc(std_ft, panel_obj, cfg)
    log.append(f"detection QC on the user's standards run: "
               f"{qc['n_detected']}/{qc['n_panel']} panel standards located "
               f"(reporting only — it does not drive the calibration)")

    model = _model_dict(cal, anchors, mapper, panel_used, irt_reason, ref,
                        panel_key, panel_manifest, qc, res, cfg,
                        scope=res["calibration_scope"].iloc[0] if len(res) else "project",
                        n_injections=n_inj)

    result = CalibrationResult(
        table=out, anchors=anchors, pairs=_pair_frame(cal), landmarks=landmarks,
        calibrator=cal, irt=mapper, model=model, log=log, columns=colmap,
        rt_col=ft.rt_col, rt_unit=ft.rt_unit, mz_col=ft.mz_col,
        sample_cols=list(ft.sample_cols), polarity=pol, source_format=ft.source_format,
        panel=qc, reference=ref.pair.to_dict(), panel_key=panel_key,
        default_panel=is_default_manifest(panel_manifest),
    )
    return result


# ---------------------------------------------------------------------------
# Stage 1 + stage 2
# ---------------------------------------------------------------------------

def _fit_calibrator(ft: FeatureTable, std_ft: FeatureTable, ref, cfg
                    ) -> Tuple[ColumnCalibrator, pd.DataFrame, List[str]]:
    """The two-stage fit, in the order stage 2 depends on."""
    log: List[str] = []
    anon = _anon_tol(cfg)
    cands = plasma_lipid_candidates()
    cand_mz = cands["precursor_mz"].to_numpy()

    # -- sample-run pairs, with the anchor candidates' own m/z masked out ----
    extra_pairs = None
    if cfg.use_sample_pairs:
        sp = match_features_by_mz(ft, ref.sample, cfg, **anon)
        extra_pairs = exclude_pairs_near_mz(sp, cand_mz, cfg, **anon)
        log.append(f"sample-run pairs: {len(sp)} matched, {len(extra_pairs)} kept "
                   f"after masking the stage-2 anchor m/z")
    else:
        log.append("sample-run pairs: disabled (use_sample_pairs=False)")

    # -- stage 1 ------------------------------------------------------------
    cal = build_calibrator(std_ft, ref.standards, cfg, source_label="source",
                           use_anchors=False, extra_pairs=extra_pairs, **anon)
    log.append(f"stage 1: {cal.n_pairs} pairs "
               f"({cal.n_pairs_standards} standards + {cal.n_pairs_sample} sample), "
               f"{cal.n_pairs_kept} kept after MAD trimming; source RT span "
               f"{cal.curve.x0:.2f}-{cal.curve.x1:.2f} -> reference "
               f"{cal.curve.y0:.2f}-{cal.curve.y1:.2f} min")

    if not cfg.use_sample_anchors:
        log.append("stage 2: disabled (use_sample_anchors=False) — stage-1 curve only")
        return cal, _anchor_frame(cal), log

    # -- stage 2: validate, then refine the picks WITH the stage-1 curve ----
    tables = {SOURCE_KEY: ft, REF_KEY: ref.sample}
    picks, _ = validate_candidates(tables, cands, cfg, ref_key=REF_KEY, **anon)
    if picks.empty:
        log.append("stage 2: no plasma-lipid candidates found in the sample run "
                   "— stage-1 curve only (expected for a non-plasma matrix)")
        return cal, _anchor_frame(cal), log

    picks = refine_picks_with_curves(picks, tables, {SOURCE_KEY: cal.curve}, cands,
                                     cfg, ref_key=REF_KEY, **anon)
    n_refined = int(picks["pick_refined"].sum())
    summary = summarize_picks(picks, cands, list(tables), ref_key=REF_KEY)
    validated = summary[summary["validated"]]["lipid"].tolist()
    anchor_tbl = build_anchor_table(picks, SOURCE_KEY, validated, cands, ref_key=REF_KEY)
    log.append(f"stage 2: {len(validated)}/{len(cands)} plasma lipids validated in "
               f"both sample runs, {len(anchor_tbl)} usable as anchors"
               + (f"; {n_refined} pick(s) moved by curve-assisted isomer refinement"
                  if n_refined else "; no pick moved by isomer refinement"))

    if len(anchor_tbl) < cfg.min_anchors:
        # Normal outcome for any non-plasma matrix: the 17-lipid panel is human
        # plasma/serum. Say so and keep the stage-1 curve — never fail.
        cal.n_anchors_validated = len(anchor_tbl)
        cal.anchors = anchor_tbl
        cal.gate_reason = (f"only {len(anchor_tbl)} plasma-lipid anchors validated "
                           f"(min_anchors={cfg.min_anchors}) — stage-1 curve only; "
                           f"expected when the matrix is not human plasma/serum")
        log.append(f"stage 2 skipped: {cal.gate_reason}")
        return cal, _anchor_frame(cal), log

    cal = build_calibrator(std_ft, ref.standards, cfg, source_label="source",
                           use_anchors=True, sample_anchors=anchor_tbl,
                           extra_pairs=extra_pairs, **anon)
    log.append(f"stage 2 gate: {cal.gate_reason}")
    if cal.anchors_used:
        log.append(f"stage 2 engaged: lam_g={cal.lam_g:.1f}, lam_c={cal.lam_c:.1f}, "
                   f"{len(cal.anchors)} anchors, sigma_anchor="
                   f"{cal.sigma_anchor:.3f} min")
    return cal, _anchor_frame(cal), log


def _anon_tol(cfg: CalibrationConfig) -> Dict[str, float]:
    """The m/z window used for every anonymous operation (spec §7).

    ``max(ppm, floor)``: the ppm window (``match_mz_tol_ppm``) with the
    absolute floor (``match_mz_tol_da``) beneath it, so low-m/z pairs keep a
    sensible minimum width. The ppm window belongs to targeted identification
    of named standards and is deliberately not applied here.
    """
    return {"mz_tol_ppm": float(cfg.match_mz_tol_ppm),
            "mz_tol_min_da": float(cfg.match_mz_tol_da)}


# ---------------------------------------------------------------------------
# The iRT ruler
# ---------------------------------------------------------------------------

def _build_irt(ref, panel_obj: Optional[Panel], cfg: CalibrationConfig, pol: str,
               panel_key: str = "") -> Tuple[Optional[IRTMapper], pd.DataFrame, str, str]:
    """Landmarks on the *reference* standards run; ``panel="none"`` falls back.

    Returns ``(mapper, landmarks, panel_used, reason)``. ``panel_used`` names the
    panel that actually ruled the scale — iRT values are only comparable between
    runs that used the same one, so it is recorded rather than left implicit.
    """
    if panel_obj is not None:
        mapper, landmarks, reason = build_irt(ref.standards, panel_obj, cfg, polarity=pol)
        return mapper, landmarks, (panel_key or "user panel"), reason

    # panel="none": stage 1 needed no identities, but the ruler does. Fall back
    # to the configured landmark panel, and be explicit that we did.
    fallback = mixtures.get_manifest(cfg.irt_landmark_panel)
    if fallback is None:
        return None, pd.DataFrame(), f"{cfg.irt_landmark_panel} (fallback)", (
            f"irt_landmark_panel='{cfg.irt_landmark_panel}' has no manifest — "
            f"iRT is NaN; Cal_RT_min is unaffected")
    panel_obj = build_panel(pol, cfg, manifest=fallback)
    mapper, landmarks, reason = build_irt(ref.standards, panel_obj, cfg, polarity=pol)
    return mapper, landmarks, f"{cfg.irt_landmark_panel} (fallback: user panel = none)", reason


# ---------------------------------------------------------------------------
# Per-project / per-sample tiers
# ---------------------------------------------------------------------------

def _per_project(rt: np.ndarray, cal: ColumnCalibrator, cfg: CalibrationConfig,
                 mapper: Optional[IRTMapper]) -> pd.DataFrame:
    cal_rt = cal.predict(rt)
    unc = cal.uncertainty(rt)
    res = _assemble(rt, cal_rt, unc, cal.is_extrapolated(rt), cal.warp_source,
                    "project", cfg, mapper)
    res["RI_spread"] = np.nan
    res["n_contributing"] = 1
    return res


def _per_sample(rt: np.ndarray, single_files: Sequence[str], ref,
                cal: ColumnCalibrator, cfg: CalibrationConfig,
                mapper: Optional[IRTMapper], pol: str,
                source_format: Optional[str], rt_unit: Optional[str]
                ) -> Tuple[pd.DataFrame, int, List[str]]:
    """One curve per injection against the reference *sample* run.

    ``Cal_RT_min`` / ``iRT`` become the median across injections and
    ``RI_spread`` their 1.4826·MAD dispersion — the honest inter-injection QC.
    Uncertainty and the extrapolation flag stay on the project-level calibrator:
    they describe the calibration, and the project curve is the one built from
    the standards runs as well as the samples.
    """
    log: List[str] = []
    anon = _anon_tol(cfg)
    preds: List[np.ndarray] = []
    for p in single_files:
        try:
            sf = load_feature_table(p, source_format=source_format, polarity=pol,
                                    rt_unit=rt_unit or cfg.rt_unit)
            pairs = match_features_by_mz(sf, ref.sample, cfg, **anon)
            if len(pairs) < 2:
                raise CalibrationError(
                    f"only {len(pairs)} m/z-matched pairs against the reference "
                    f"sample run")
            curve, _ = fit_robust_curve(
                pairs["rt_a"].to_numpy(), pairs["rt_b"].to_numpy(),
                frac=cfg.curve_frac, n_iter=cfg.curve_iter, mad_k=cfg.curve_mad_k,
                min_points=cfg.curve_min_points, mode=cfg.extrapolate_mode)
            preds.append(curve.predict(rt))
            log.append(f"  per-injection curve ok: {os.path.basename(p)} "
                       f"({len(pairs)} pairs)")
        except RtAnchorError as e:
            log.append(f"  per-injection curve FAILED for {os.path.basename(p)}: {e} "
                       f"(skipped)")
    if not preds:
        log.append("no per-injection curve succeeded -> falling back to per-project")
        return _per_project(rt, cal, cfg, mapper), 0, log

    import warnings

    mat = np.vstack(preds)                       # (n_injections, n_features)
    n_contrib = np.sum(np.isfinite(mat), axis=0)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        cal_rt = np.nanmedian(mat, axis=0)
    cal_rt[n_contrib == 0] = np.nan

    # features no injection covers -> the project curve, so a row is never blank
    proj = cal.predict(rt)
    fallback = (n_contrib == 0) & np.isfinite(proj)
    if fallback.any():
        cal_rt[fallback] = proj[fallback]
        log.append(f"per-injection fallback: {int(fallback.sum())} features covered "
                   f"by no injection took the project-level curve")

    unc = cal.uncertainty(rt)
    res = _assemble(rt, cal_rt, unc, cal.is_extrapolated(rt), cal.warp_source,
                    "sample", cfg, mapper)
    if mapper is not None:
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            irt_mat = np.vstack([mapper.to_irt(row) for row in mat])
            res["iRT"] = np.where(n_contrib > 0, np.nanmedian(irt_mat, axis=0),
                                  res["iRT"].to_numpy())
            res["RI_spread"] = _nanmad(irt_mat, axis=0)
    else:
        res["RI_spread"] = np.nan
    res["n_contributing"] = np.maximum(n_contrib, 1)
    res.loc[n_contrib == 0, "n_contributing"] = 1
    log.append(f"per-sample tier: {len(preds)} injection curves, median iRT with "
               f"1.4826·MAD spread")
    return res, len(preds), log


def _assemble(rt: np.ndarray, cal_rt: np.ndarray, unc: np.ndarray,
              extrapolated: np.ndarray, warp_source: str, scope: str,
              cfg: CalibrationConfig, mapper: Optional[IRTMapper]) -> pd.DataFrame:
    """The §1 output columns for one set of predictions."""
    rt = np.asarray(rt, dtype=float)
    bad = ~np.isfinite(rt)                       # no RT -> no calibrated RT
    cal_rt = np.asarray(cal_rt, dtype=float).copy()
    unc = np.asarray(unc, dtype=float).copy()
    cal_rt[bad] = np.nan
    unc[bad] = np.nan
    extrapolated = np.asarray(extrapolated, dtype=bool).copy()
    extrapolated[bad] = False                    # unknown, not out-of-span

    if mapper is None:
        irt = np.full(rt.shape, np.nan)
        irt_unc = np.full(rt.shape, np.nan)
    else:
        irt = mapper.to_irt(cal_rt)
        irt_unc = unc * np.abs(mapper.slope(cal_rt))
        irt[~np.isfinite(cal_rt)] = np.nan
        irt_unc[~np.isfinite(cal_rt)] = np.nan

    rel = np.full(rt.shape, "medium", dtype=object)
    with np.errstate(invalid="ignore"):
        rel[np.asarray(irt_unc, dtype=float) < cfg.conf_high_irt] = "high"
        rel[np.asarray(irt_unc, dtype=float) > cfg.conf_low_irt] = "low"
    rel[extrapolated] = "low"
    rel[~np.isfinite(np.asarray(irt, dtype=float))] = "none"

    return pd.DataFrame({
        "Cal_RT_min": cal_rt,
        "Cal_RT_uncertainty_min": unc,
        "iRT": irt,
        "iRT_uncertainty": irt_unc,
        "iRT_reliability": rel,
        "is_extrapolated": extrapolated,
        "calibration_scope": scope,
        "warp_source": warp_source,
        "RI_spread": np.nan,
        "n_contributing": 1,
    })


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve_panel(panel: Optional[str], manifest: Optional[pd.DataFrame]
                   ) -> Tuple[str, Optional[pd.DataFrame]]:
    """``(panel_key, manifest_or_None)``; an explicit manifest overrides the key."""
    if manifest is not None:
        return "custom", manifest
    key = (panel or mixtures.NONE).strip().lower()
    if key not in mixtures.mixture_keys():
        raise ConfigError(
            f"Unknown panel '{panel}'. Offered: {mixtures.offered_keys()} "
            f"(registered: {mixtures.mixture_keys()}). Pass manifest=<DataFrame> "
            f"for a custom standard list."
        )
    return key, mixtures.get_manifest(key)


def _detection_qc(std_ft: FeatureTable, panel_obj: Optional[Panel],
                  cfg: CalibrationConfig) -> Dict:
    """Where the chosen panel's standards sit in the **user's own** standards run.

    Reporting only: under v2 this does not drive the calibration at all. It is
    what keeps the app's Detection section meaningful, and it is the first place
    to look when a run's numbers are surprising.
    """
    srt = std_ft.rt_minutes()
    info: Dict = {
        "rt_range": [float(np.nanmin(srt)), float(np.nanmax(srt))] if len(srt) else [np.nan, np.nan],
        "n_features": int(std_ft.n_features()),
        "n_panel": 0 if panel_obj is None else int(len(panel_obj)),
        "n_detected": 0,
        "native_rt": {},
        "targets": pd.DataFrame() if panel_obj is None else panel_obj.targets.copy(),
        "note": "QC only — does not drive the calibration",
    }
    if panel_obj is None:
        return info
    native = build_native_template(std_ft, panel_obj, cfg)
    info["native_rt"] = {k: float(v) for k, v in native.items()}
    info["n_detected"] = len(native)
    return info


def _pair_frame(cal: Optional[ColumnCalibrator]) -> pd.DataFrame:
    """``<prefix>_pairs.csv``: the stage-1 evidence, one row per matched pair."""
    if cal is None or cal.pairs is None or cal.pairs.empty:
        return pd.DataFrame(columns=PAIR_COLUMNS)
    p = cal.pairs.rename(columns={"mz_a": "mz_src", "rt_a": "rt_src", "rt_b": "rt_ref"})
    for c in PAIR_COLUMNS:
        if c not in p.columns:
            p[c] = np.nan
    return p[PAIR_COLUMNS].reset_index(drop=True)


def _anchor_frame(cal: Optional[ColumnCalibrator]) -> pd.DataFrame:
    """``<prefix>_anchors.csv``: anchors used *and* those the sanity filter dropped."""
    if cal is None:
        return pd.DataFrame(columns=ANCHOR_COLUMNS)
    used = cal.anchor_table()
    dropped = cal.anchors_dropped
    frames = [f for f in (used, dropped) if f is not None and len(f)]
    if not frames:
        return pd.DataFrame(columns=ANCHOR_COLUMNS)
    out = pd.concat(frames, ignore_index=True, sort=False)
    if "dropped_by_sanity_filter" not in out.columns:
        out["dropped_by_sanity_filter"] = False
    out["dropped_by_sanity_filter"] = out["dropped_by_sanity_filter"].fillna(False).astype(bool)
    # recompute for the dropped rows too, so every row's residual is comparable
    out["residual_min"] = (out["rt_ref"].to_numpy(dtype=float)
                           - cal.curve.predict(out["rt_src"].to_numpy(dtype=float)))
    for c in ANCHOR_COLUMNS:
        if c not in out.columns:
            out[c] = np.nan
    rest = [c for c in out.columns if c not in ANCHOR_COLUMNS]
    return out[ANCHOR_COLUMNS + rest].sort_values("rt_src").reset_index(drop=True)


def _append_results(df: pd.DataFrame, res: pd.DataFrame):
    """Append the result columns without dropping or reordering an original column.

    Returns ``(out_df, colmap)`` where colmap maps each logical result column to
    its actual name (``rtanchor_<name>`` on collision).
    """
    out = df.copy()
    res = res.reset_index(drop=True)
    colmap: Dict[str, str] = {}
    for col in RESULT_COLUMNS:
        name = col if col not in out.columns else f"rtanchor_{col}"
        out[name] = res[col].to_numpy()
        colmap[col] = name
    return out, colmap


def _nanmad(mat: np.ndarray, axis: int = 0) -> np.ndarray:
    med = np.nanmedian(mat, axis=axis)
    return np.nanmedian(np.abs(mat - med), axis=axis) * 1.4826


def _span(ft: FeatureTable) -> str:
    rt = ft.rt_minutes()
    if not len(rt):
        return "empty"
    return f"{float(np.nanmin(rt)):.2f}-{float(np.nanmax(rt)):.2f} min"


def _model_dict(cal: ColumnCalibrator, anchors: pd.DataFrame,
                mapper: Optional[IRTMapper], panel_used: str, irt_reason: str,
                ref, panel_key: str, panel_manifest: Optional[pd.DataFrame],
                qc: Dict, res: pd.DataFrame, cfg: CalibrationConfig,
                scope: str, n_injections: int) -> Dict:
    """The model JSON of spec §5.1."""
    blocks = cal.to_model()
    conf = pd.Series(res["iRT_reliability"]).value_counts().to_dict()
    try:
        panel_label = mixtures.get_meta(panel_key)["label"]
    except KeyError:
        panel_label = "custom manifest"
    irt_block = ({"n_landmarks": 0, "landmark_rt_span_min": [None, None],
                  "irt_range": [None, None],
                  "definition": "affine on the reference-column RT axis",
                  "panel_used": panel_used}
                 if mapper is None else mapper.to_model(panel_used))
    irt_block["reason"] = irt_reason

    model = {
        "method": "cross-column-v2",
        "calibration_scope": scope,
        "reference": ref.pair.to_dict(),
        "panel": {
            "key": panel_key,
            "label": panel_label,
            "n_standards": (0 if panel_manifest is None else int(len(panel_manifest))),
        },
        "curve": blocks["curve"],
        "anchors": {**blocks["anchors"], "table": anchors.to_dict(orient="records")},
        "irt": irt_block,
        "detection_qc": {
            "user_standards_run": {
                "n_panel": qc["n_panel"],
                "n_detected": qc["n_detected"],
                "native_rt": qc["native_rt"],
                "rt_range_min": qc["rt_range"],
                "n_features": qc["n_features"],
                "note": qc["note"],
            }
        },
        "n_features": int(len(res)),
        "n_features_extrapolated": int(res["is_extrapolated"].sum()),
        "confidence_counts": {str(k): int(v) for k, v in conf.items()},
        "config": cfg.to_dict(),
    }
    if n_injections:
        model["n_injections_used"] = int(n_injections)
    return model
