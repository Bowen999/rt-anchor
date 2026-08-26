"""End-to-end orchestration: load -> panel -> identify -> warp -> apply -> assemble.

Two tiers (auto-selected):

* **per-project** (default): identify the standards in the table being
  calibrated, fit one warp on their observed RT, apply to every feature.
* **per-sample** (when ``single_files`` are given): fit one warp per injection
  and aggregate (median RI + spread), yielding the inter-sample ``RI_spread`` QC.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .calibrate import MonotoneWarp, apply_warp, fit_warp
from .config import CalibrationConfig
from .errors import RtAnchorError
from .identify import build_native_template, identify_anchors
from .io.loader import load_feature_table
from .io.schema import FeatureTable, RESULT_COLUMNS
from .panel import Panel, build_panel


@dataclass
class CalibrationResult:
    table: pd.DataFrame                 # original columns + appended RI columns
    anchors: pd.DataFrame               # anchors used (per-project) or pooled
    warp: Optional[MonotoneWarp]
    model: Dict                         # serialisable calibration model
    log: List[str] = field(default_factory=list)
    # ---- metadata for downstream visualisation / QC ----
    rt_col: str = ""                    # RT column name in `table`
    rt_unit: str = "min"
    mz_col: str = ""
    sample_cols: List[str] = field(default_factory=list)
    polarity: Optional[str] = None
    source_format: str = ""
    columns: Dict[str, str] = field(default_factory=dict)  # result-col -> actual name in table
    panel: Dict = field(default_factory=dict)              # standards-run summary (RT range, detection, native RT, resolved targets, n_features)
    default_panel: bool = True                             # built-in 15-standard panel (enables structure hover)

    def to_model_json(self) -> Dict:
        return self.model

    def rt_minutes(self):
        rt = pd.to_numeric(self.table[self.rt_col], errors="coerce")
        return rt / 60.0 if self.rt_unit == "sec" else rt

    def col(self, result_col: str) -> str:
        """Actual table column name for a result column (handles collision rename)."""
        return self.columns.get(result_col, result_col)


def calibrate(sample_table: str,
              polarity: str,
              standards_table: str,
              single_files: Optional[List[str]] = None,
              config: Optional[CalibrationConfig] = None,
              manifest: Optional[pd.DataFrame] = None,
              reference: Optional[pd.DataFrame] = None,
              source_format: Optional[str] = None,
              rt_unit: Optional[str] = None) -> CalibrationResult:
    if not standards_table:
        from .errors import ConfigError
        raise ConfigError(
            "A standards run is required. Pass standards_table=<the standard-panel "
            "run>; it builds the native anchor template the calibration depends on "
            "(especially on QTOF, where RT is shifted several minutes from the "
            "reference). If the standards are spiked into the sample itself, pass "
            "that same table as standards_table."
        )
    cfg = config or CalibrationConfig()
    log: List[str] = []

    ft = load_feature_table(sample_table, source_format=source_format,
                            polarity=polarity, rt_unit=rt_unit)
    log.append(f"loaded sample table: {ft.summary()}")
    panel = build_panel(polarity, cfg, manifest=manifest, reference=reference)
    log.append(f"panel: {len(panel)} targets in {polarity} mode")

    # native template from the standard-panel run (required): find where each
    # anchor actually elutes on this method to centre the later search windows
    std_ft = load_feature_table(standards_table, source_format=source_format,
                                polarity=polarity, rt_unit=rt_unit)
    native = build_native_template(std_ft, panel, cfg)
    log.append(f"native template built from standards run: {len(native)}/{len(panel)} anchors")

    if single_files:
        result = _per_sample(ft, single_files, panel, cfg, native, polarity,
                             source_format, rt_unit, log)
    else:
        result = _per_project(ft, panel, cfg, native, log)

    # attach metadata for visualisation / QC
    result.rt_col, result.rt_unit, result.mz_col = ft.rt_col, ft.rt_unit, ft.mz_col
    result.sample_cols = list(ft.sample_cols)
    result.polarity, result.source_format = polarity, ft.source_format
    from .panel import is_default_manifest
    result.default_panel = is_default_manifest(manifest)
    srt = std_ft.rt_minutes()
    result.panel = {
        "rt_range": [float(srt.min()), float(srt.max())],
        "n_detected": len(native), "n_panel": len(panel),
        "native_rt": {k: float(v) for k, v in native.items()},
        "targets": panel.targets.copy(),
        "n_features": int(std_ft.df.shape[0]),
    }
    return result


# ---------------------------------------------------------------- per-project -

def _per_project(ft: FeatureTable, panel: Panel, cfg: CalibrationConfig,
                 native: Optional[Dict[str, float]], log: List[str]) -> CalibrationResult:
    anchors = identify_anchors(ft, panel, cfg, native_template=native)
    n_drop = anchors.attrs.get("n_dropped_nonmonotone", 0)
    log.append(f"identified {len(anchors)} anchors (dropped {n_drop} non-monotone)")
    _guard_anchor_count(anchors, cfg, log)

    warp = fit_warp(anchors["rt_obs_min"].to_numpy(), anchors["irt"].to_numpy(), cfg)
    log.append(f"fitted warp on {warp.rt.size} anchors, span "
               f"{warp.rt_min:.2f}-{warp.rt_max:.2f} min")

    res = apply_warp(ft.rt_minutes().to_numpy(), warp, cfg, scope="project", warp_source="self")
    # per-project has no cross-sample dispersion
    res["RI_spread"] = np.nan
    res["n_contributing"] = 1
    out, colmap = _append_results(ft.df, res)
    model = _model_dict("project", anchors, warp, cfg, res, n_drop)
    return CalibrationResult(table=out, anchors=anchors, warp=warp, model=model,
                             log=log, columns=colmap)


# ---------------------------------------------------------------- per-sample --

def _per_sample(ft: FeatureTable, single_files: List[str], panel: Panel,
                cfg: CalibrationConfig, native: Optional[Dict[str, float]],
                polarity: str, source_format: Optional[str], rt_unit: Optional[str],
                log: List[str]) -> CalibrationResult:
    """Fit one warp per injection; evaluate each on the aligned consensus RT;
    aggregate to median RI + spread. (MassCube consensus RT is already batch-
    aligned, so per-sample RT ~ consensus RT up to the small residual drift the
    warps model.)"""
    consensus_rt = ft.rt_minutes().to_numpy()
    preds = []
    used_anchors = []
    for p in single_files:
        try:
            sf = load_feature_table(p, source_format=source_format or "masscube",
                                    polarity=polarity, rt_unit=rt_unit)
            a = identify_anchors(sf, panel, cfg, native_template=native)
            w = fit_warp(a["rt_obs_min"].to_numpy(), a["irt"].to_numpy(), cfg)
            preds.append(w.predict(consensus_rt, extrapolate=cfg.extrapolate,
                                   max_extrap_min=cfg.max_extrapolation_min))
            a["_sample"] = p
            used_anchors.append(a)
            log.append(f"  per-sample warp ok: {p} ({w.rt.size} anchors)")
        except RtAnchorError as e:
            log.append(f"  per-sample warp FAILED for {p}: {e} (skipped)")
    if not preds:
        log.append("no per-sample warps succeeded -> falling back to per-project")
        return _per_project(ft, panel, cfg, native, log)

    import warnings

    mat = np.vstack(preds)                      # (n_samples, n_features)
    n_contrib = np.sum(np.isfinite(mat), axis=0)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        # features with no contributing sample -> all-NaN column; the NaN is
        # expected and masked below, so silence the RuntimeWarning noise.
        warnings.simplefilter("ignore", RuntimeWarning)
        ri = np.nanmedian(mat, axis=0)
        spread = _nanmad(mat, axis=0)
    ri[n_contrib == 0] = np.nan

    # reuse a representative warp (pooled anchors) for sigma/confidence + span flag
    pooled = pd.concat(used_anchors, ignore_index=True)
    pooled_med = (pooled.groupby("name")
                  .agg(rt_obs_min=("rt_obs_min", "median"), irt=("irt", "first"),
                       rt_ref_min=("rt_ref_min", "first"), mz=("mz", "first"),
                       class_=("class", "first"))
                  .reset_index().sort_values("rt_ref_min"))
    warp = fit_warp(pooled_med["rt_obs_min"].to_numpy(), pooled_med["irt"].to_numpy(), cfg)
    res = apply_warp(consensus_rt, warp, cfg, scope="sample", warp_source="self")
    res["RI"] = ri                              # median across per-sample warps
    res["RI_spread"] = spread
    res["n_contributing"] = n_contrib
    res.loc[n_contrib == 0, "RI_reliability"] = "none"

    out, colmap = _append_results(ft.df, res)
    model = _model_dict("sample", pooled_med.rename(columns={"class_": "class"}),
                        warp, cfg, res, n_dropped=0)
    model["n_injections_used"] = len(preds)
    return CalibrationResult(table=out, anchors=pooled, warp=warp, model=model,
                             log=log, columns=colmap)


# ---------------------------------------------------------------- helpers -----

def _guard_anchor_count(anchors: pd.DataFrame, cfg: CalibrationConfig, log: List[str]) -> None:
    if len(anchors) < cfg.min_anchors_hard:
        from .errors import AnchorIdentificationError
        raise AnchorIdentificationError(
            f"Only {len(anchors)} anchors (< hard minimum {cfg.min_anchors_hard}); "
            f"cannot calibrate. Provide a standards run / relax tolerances / check polarity."
        )
    if len(anchors) < cfg.min_anchors:
        log.append(f"WARNING: only {len(anchors)} anchors (< recommended {cfg.min_anchors}) "
                   f"-> results are low-confidence")


def _append_results(df: pd.DataFrame, res: pd.DataFrame):
    """Append result columns without dropping/altering any original column.

    Returns ``(out_df, colmap)`` where colmap maps each logical result column to
    its actual name (``rtanchor_<name>`` on collision)."""
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


def _model_dict(scope: str, anchors: pd.DataFrame, warp: MonotoneWarp,
                cfg: CalibrationConfig, res: pd.DataFrame, n_dropped: int) -> Dict:
    interior = np.isfinite(warp.loo_resid)
    loo = np.abs(warp.loo_resid[interior]) if interior.any() else np.array([])
    conf = pd.Series(res["RI_reliability"]).value_counts().to_dict()
    return {
        "calibration_scope": scope,
        "n_anchors_used": int(warp.rt.size),
        "n_anchors_dropped_nonmonotone": int(n_dropped),
        "anchor_span_min": [float(warp.rt_min), float(warp.rt_max)],
        "loo_residual_irt": {
            "median": (float(np.median(loo)) if loo.size else None),
            "p90": (float(np.percentile(loo, 90)) if loo.size else None),
            "max": (float(np.max(loo)) if loo.size else None),
        },
        "confidence_counts": {str(k): int(v) for k, v in conf.items()},
        "n_features": int(len(res)),
        "n_features_extrapolated": int(res["is_extrapolated"].sum()),
        "scale_constants": {"RT_lo_min": cfg.scale_rt_lo_min, "RT_hi_min": cfg.scale_rt_hi_min},
        "config": cfg.to_dict(),
        "anchors": anchors.to_dict(orient="records"),
    }
