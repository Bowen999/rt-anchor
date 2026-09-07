"""App-native visualization data (v2 cross-column method).

Extracts the numeric arrays each in-app chart needs from a ``CalibrationResult``
and hands them to the front-end, which draws them as hand-authored SVG
(``web/charts.js``). **Detection** is the package's own Plotly figure; the
**Profile** section is two app-built Plotly mirror pseudo-TICs (before/after
calibration, and calibrated input vs. reference).

Where a number also appears in the exported report, it is taken from the same
engine helper the report uses rather than recomputed here:

* KPI tiles and the radar   -> ``rt_anchor.viz.metrics.kpi_tiles`` / ``radar_axes``
* the cross-column curve    -> ``rt_anchor.viz.performance.compute_curve``

That is deliberate: the app and the PDF a reviewer receives must not be able to
disagree about what the run did. The *drawing* is the app's own (brutalist SVG,
cool-blue palette); only the arithmetic is shared.
"""

from __future__ import annotations

import math
import os
from typing import Dict, List

import numpy as np
import pandas as pd


# --------------------------------------------------------------- json helpers --

def _f(v):
    """A finite float, or None (JSON-safe)."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _arr(a) -> List:
    return [_f(x) for x in np.asarray(a, dtype=float).ravel()]


def _bools(a) -> List[bool]:
    return [bool(x) for x in np.asarray(a).ravel()]


def _s(v) -> str:
    return "" if v is None else str(v)


# --------------------------------------------------------------- shared bits ---

def _intensity(result) -> np.ndarray:
    df = result.table
    if result.sample_cols:
        return (df[result.sample_cols].apply(pd.to_numeric, errors="coerce")
                .fillna(0.0).sum(axis=1).to_numpy())
    cols = {c.lower(): c for c in df.columns}
    for key in ("peak_height", "peak_area", "maxo", "height", "area", "intensity"):
        if key in cols:
            return pd.to_numeric(df[cols[key]], errors="coerce").fillna(0.0).to_numpy()
    return np.ones(len(df), dtype=float)


def _tile(label, value, note) -> Dict:
    """Engine KPI tuple -> the front-end's tile shape.

    ``value`` is either a headline string or a list of ``(caption, value)``
    pairs; the two render differently (bold tile vs. a captioned metric row), so
    the distinction is carried through rather than flattened.
    """
    d = {"label": str(label), "sub": _s(note)}
    if isinstance(value, (list, tuple)):
        d["rows"] = [{"k": str(k), "v": str(v)} for k, v in value]
    else:
        d["value"] = str(value)
    return d


# --------------------------------------------------------------- KPIs ----------

#: Engine KPI tiles that the app deliberately does not show.
_REMOVED_TILES = {"Matched pairs", "Anchor gate", "iRT landmarks"}


def kpis(result) -> List[Dict]:
    """The v2 KPI tiles, straight from the engine's own tile builder.

    A few tiles are dropped in the app: the pair census, the stage-2 gate and
    the iRT landmark count are evidence-level details the results screen no
    longer surfaces.
    """
    from rt_anchor.viz import metrics as _vm
    return [_tile(label, value, note) for label, value, note in _vm.kpi_tiles(result)
            if str(label) not in _REMOVED_TILES]


def radar(result) -> List[Dict]:
    """Normalised quality axes (outward = better), from the engine."""
    from rt_anchor.viz import metrics as _vm
    return [{"axis": a, "value": max(0.0, min(1.0, _f(v) or 0.0))}
            for a, v in _vm.radar_axes(result)]


# --------------------------------------------------------------- curve ---------

def curve(result, ngrid: int = 240) -> Dict:
    """Stage-1 cross-column curve + its residuals (spec §12, the Warp section).

    Mirrors ``rt_calibration_diagnostics.png``: the matched pairs with the two
    sources drawn differently, the MAD-rejected pairs kept visible as
    *excluded*, the anchor-free fit, the anchor-refined fit when the gate
    engaged, the stage-2 anchors as rings, and the residual strip below.
    """
    from rt_anchor.viz import performance as _vp
    d = _vp.compute_curve(result, ngrid=ngrid)
    if d.get("empty"):
        return {"empty": True, "reason": _s(d.get("reason"))}

    rt_src = np.asarray(d["rt_src"], dtype=float)
    rt_ref = np.asarray(d["rt_ref"], dtype=float)
    is_std = np.asarray(d["is_std"], dtype=bool)
    kept = np.asarray(d["kept"], dtype=bool)
    r_curve = np.asarray(d["r_curve"], dtype=float)
    r_ref = d.get("r_refined")
    r_ref = None if r_ref is None else np.asarray(r_ref, dtype=float)

    anc = d.get("anchors")
    anchors = []
    if anc is not None and len(anc):
        for _, r in anc.iterrows():
            anchors.append({"x": _f(r.get("rt_src")), "y": _f(r.get("rt_ref")),
                            "label": _s(r.get("label")), "class": _s(r.get("lipid_class")),
                            "resid": _f(r.get("residual_min")),
                            "loo": _f(r.get("loo_residual_min"))})

    span_src = [_f(np.nanmin(rt_src)), _f(np.nanmax(rt_src))]
    span_ref = [_f(np.nanmin(rt_ref)), _f(np.nanmax(rt_ref))]
    return {
        "empty": False,
        "fit": {"x": _arr(d["grid"]), "y": _arr(d["fit"])},
        "refined": (None if d.get("refined") is None
                    else {"x": _arr(d["grid"]), "y": _arr(d["refined"])}),
        "pairs": {"x": _arr(rt_src), "y": _arr(rt_ref),
                  "std": _bools(is_std), "kept": _bools(kept),
                  "r": _arr(r_curve),
                  "rr": (None if r_ref is None else _arr(r_ref))},
        "anchors": anchors,
        "rt_span": span_src, "ref_span": span_ref,
        "engaged": bool(d["engaged"]),
        "lam_g": _f(d["lam_g"]), "lam_c": _f(d["lam_c"]),
        "class_aware": bool(d["class_aware"]),
        "gate_reduction": _f(d["gate_reduction"]),
        "gate_threshold": _f(d["gate_threshold"]),
        "gate_reason": _s(d["gate_reason"]),
        "legend": {"fit": _vp.fit_legend_label(d) if d["engaged"] else "",
                   "anchors": _vp.anchor_legend_label(d)},
        "counts": {"n_pairs": int(rt_src.size), "n_kept": int(kept.sum()),
                   "n_standards": int(is_std.sum()), "n_sample": int((~is_std).sum())},
        "labels": {"src": _s(d.get("source_label")), "ref": _s(d.get("ref_label"))},
    }


# --------------------------------------------------------------- profile -------

def _pseudo_tic(x, inten):
    """Binned, lightly smoothed pseudo-chromatogram of (x, intensity) rows.

    Feature tables carry no scan-level signal, so a TIC-like profile is
    *reconstructed* by binning feature apex intensities — the same binning the
    engine's own TIC figure uses, so the app and the report agree.
    """
    x = np.asarray(x, dtype=float)
    inten = np.asarray(inten, dtype=float)
    ok = np.isfinite(x) & np.isfinite(inten)
    if not ok.any():
        return None
    lo, hi = float(x[ok].min()), float(x[ok].max())
    from rt_anchor.viz import tic as _vt
    c, h = _vt._reconstruct(x[ok], inten[ok], lo, hi, "clean")
    return c, h


def _reference_tic(result):
    """``(rt_min, intensity)`` of the reference *sample* run.

    The reference run's feature table is not stored on the result, so it is
    re-loaded from the recorded path (bundled or user-supplied). ``None`` when
    it cannot be loaded — the calibrated-comparison figure then shows the input
    alone.
    """
    ref = result.reference or {}
    path = ref.get("sample")
    if not path or not os.path.isfile(path):
        return None
    try:
        from rt_anchor.io.loader import load_feature_table
        ft = load_feature_table(path, polarity=result.polarity)
        rt = ft.rt_minutes().to_numpy(dtype=float)
        inten = pd.to_numeric(ft.intensity(), errors="coerce").fillna(0.0).to_numpy()
    except Exception:
        return None
    return _pseudo_tic(rt, inten)


def _ticks(lo, hi, step):
    start = math.ceil(lo / step) * step
    return [float(v) for v in np.arange(start, hi + 1e-9, step)]


def _tic_figure(top_c, top_h, top_title, top_color,
                bot_c, bot_h, bot_title, bot_color, x_title, height=440):
    """Mirror pseudo-TIC on ONE shared time axis, with a visible Y axis.

    Top trace points up, bottom trace points down, and both are drawn against
    the same X coordinate system (retention time, min) so the two profiles are
    directly comparable. The Y axis is normalised intensity (±1).
    """
    import plotly.graph_objects as go
    from rt_anchor.viz import theme

    top_c = np.asarray(top_c, dtype=float)
    bot_c = np.asarray(bot_c, dtype=float)
    top_h = np.asarray(top_h, dtype=float) / (float(np.max(top_h)) or 1.0)
    bot_h = np.asarray(bot_h, dtype=float) / (float(np.max(bot_h)) or 1.0)

    lo = float(min(top_c[0], bot_c[0]))
    hi = float(max(top_c[-1], bot_c[-1]))
    pad = 0.03 * (hi - lo) if hi > lo else 0.5
    x0, x1 = lo - pad, hi + pad

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.r_[top_c, top_c[::-1]], y=np.r_[top_h, np.zeros_like(top_h)],
                  fill="toself", fillcolor=theme.rgba(top_color, 0.85), line=dict(width=0),
                  name=top_title, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=np.r_[bot_c, bot_c[::-1]], y=np.r_[-bot_h, np.zeros_like(bot_h)],
                  fill="toself", fillcolor=theme.rgba(bot_color, 0.85), line=dict(width=0),
                  name=bot_title, hoverinfo="skip"))
    fig.add_hline(y=0, line=dict(color=theme.AXIS, width=0.8))

    step = 10.0 if hi - lo > 60 else 5.0
    xt = _ticks(x0, x1, step)
    fig.update_layout(theme.plotly_template())
    fig.update_layout(
        title=None, height=height, showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
        xaxis=dict(range=[x0, x1], showgrid=False, tickvals=xt,
                   title=dict(text=x_title)),
        yaxis=dict(range=[-1.2, 1.2], tickvals=[-1.0, -0.5, 0.0, 0.5, 1.0],
                   ticktext=["−1", "−0.5", "0", "0.5", "1"],
                   title=dict(text="normalized intensity"), zeroline=False))
    return fig


def profile_figures(result) -> Dict[str, str]:
    """The Profile section's two Plotly figures.

    * ``profile_before_after`` — the input run before (raw RT) and after
      (Cal_RT) calibration, mirror-style on one shared RT axis.
    * ``profile_calibrated`` — the calibrated input on top, the reference
      sample run below, both on the same reference-column RT axis.
    """
    from rt_anchor.viz import theme

    rt = result.rt_minutes().to_numpy(dtype=float)
    cal = result.values("Cal_RT_min").to_numpy(dtype=float)
    inten = _intensity(result)

    before = _pseudo_tic(rt, inten)
    after = _pseudo_tic(cal, inten)
    ref = _reference_tic(result)

    figs = {}
    if before is not None and after is not None:
        figs["profile_before_after"] = _tic_figure(
            before[0], before[1], "before · raw RT", theme.ACCENT,
            after[0], after[1], "after · Cal_RT", theme.PRIMARY,
            "RT (min)").to_json()
    if after is not None:
        if ref is not None:
            figs["profile_calibrated"] = _tic_figure(
                after[0], after[1], "input · Cal_RT", theme.PRIMARY,
                ref[0], ref[1], "reference", theme.PRIMARY_LT,
                "RT on the reference column (min)").to_json()
        else:
            fig = _tic_figure(
                after[0], after[1], "input · Cal_RT", theme.PRIMARY,
                after[0], after[1], "reference", theme.PRIMARY_LT,
                "RT on the reference column (min)")
            fig.update_layout(annotations=[dict(
                text="reference profile unavailable for this run",
                showarrow=False, x=0.5, y=0.5, xref="paper", yref="paper",
                font=dict(size=14, color=theme.TXT2))])
            figs["profile_calibrated"] = fig.to_json()
    return figs


# --------------------------------------------------------------- table preview -

#: Result columns shown in the preview, in order, with their display names.
PREVIEW_COLUMNS = [
    ("Cal_RT_min", "Cal_RT (min)"),
    ("Cal_RT_uncertainty_min", "Cal_RT ± (min)"),
    ("iRT", "iRT"),
    ("iRT_uncertainty", "iRT ±"),
    ("iRT_reliability", "reliability"),
    ("RI_spread", "iRT spread"),
    ("is_extrapolated", "extrapolated"),
]


def table_preview(result, n: int = 120) -> Dict:
    tbl = result.table
    # RI_spread is structurally NaN outside the per-injection tier — a column of
    # dashes is noise, not information, so it only appears when it has values.
    per_sample = result.model.get("calibration_scope") == "sample"
    wanted = [(result.rt_col, "RT"), (result.mz_col, "m/z")]
    wanted += [(result.col(key), disp) for key, disp in PREVIEW_COLUMNS
               if per_sample or key != "RI_spread"]
    cols = [(c, d) for c, d in wanted if c and c in tbl.columns]
    head = tbl[[c for c, _ in cols]].head(n)

    def cell(v):
        if isinstance(v, (bool, np.bool_)):
            return "yes" if bool(v) else "no"
        if isinstance(v, (int, np.integer)):
            return int(v)
        if isinstance(v, (float, np.floating)):
            return _f(round(float(v), 4))
        return None if pd.isna(v) else str(v)

    rows = [[cell(r[c]) for c, _ in cols] for _, r in head.iterrows()]
    return {"columns": [d for _, d in cols], "rows": rows,
            "n_total": int(len(tbl)), "n_shown": int(len(head))}


# --------------------------------------------------------------- bundle --------

def _notes(result) -> Dict:
    from rt_anchor.viz import metrics as _vm
    return {
        "detection": _vm.DETECTION_NOTE,
        "curve": (
            "The stage-1 curve maps <b>your column's RT</b> onto the "
            "<b>reference column's RT</b>. Every point is an anonymous "
            "m/z-matched feature pair — <b>squares</b> from the two standards "
            "runs, <b>diamonds</b> from the two sample runs (serum covers the "
            "sparse early and late ends the mixture cannot). Rings are the "
            "stage-2 plasma-lipid anchors. The strip below is each pair's "
            "residual about the fit."
        ),
        "profile": (
            "Reconstructed from the feature table by binning feature intensities "
            "— a pseudo-chromatogram, not a scan-level TIC. Each panel draws two "
            "profiles on one shared RT axis, with intensity normalized per "
            "profile. <b>Top panel</b>: your run before (raw RT) and after "
            "calibration. <b>Bottom panel</b>: your calibrated run against the "
            "reference sample run."
        ),
    }


def build_viz(result) -> Dict:
    m = result.model
    n = len(result.table)
    curve_d = curve(result)
    qc = (m.get("detection_qc", {}) or {}).get("user_standards_run", {}) or {}
    irt = m.get("irt", {}) or {}

    struct = {}
    # `default_panel` is True for panel="none" too (the engine treats "no
    # manifest" as the built-in one), but a user who told us they have no panel
    # has no standards to hover — and the depictions cost ~340 KB in the bundle.
    if result.default_panel and result.panel_key != "none":
        try:
            from rt_anchor.viz import structures  # molecular depiction only (not a chart)
            struct = structures.render_default_structures()
        except Exception:
            struct = {}

    # Detection and the two Profile panels are Plotly figures; radar and the
    # curve are the app's own SVG charts.
    from rt_anchor.viz import metrics as _vm
    figures = {"detection": _vm.figure_plotly(result).to_json()}
    figures.update(profile_figures(result))

    ref = m.get("reference", {}) or result.reference or {}
    return {
        "kpis": kpis(result),
        "radar": radar(result),
        "figures": figures,
        "curve": curve_d,
        "table": table_preview(result),
        "structures": struct,
        "notes": _notes(result),
        "meta": {
            "method": m.get("method", "cross-column-v2"),
            "scope": m.get("calibration_scope", "project"),
            "polarity": result.polarity, "source_format": result.source_format,
            "n_features": int(n),
            "panel": result.panel_key,
            "reference_label": _s(ref.get("label")),
            "reference_default": bool(ref.get("is_default", False)),
            "n_detected": int(qc.get("n_detected", 0) or 0),
            "n_panel": int(qc.get("n_panel", 0) or 0),
            "n_landmarks": int(irt.get("n_landmarks", 0) or 0),
            "has_curve": bool(not curve_d.get("empty")),
        },
    }
