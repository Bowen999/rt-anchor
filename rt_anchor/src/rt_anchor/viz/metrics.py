"""Module 1 — key metrics: KPI tiles, a quality radar, and the detection ladder.

Everything here reads the v2 model (spec §5.1) rather than the retired v1 warp:
:func:`compute_metrics` pulls from ``model["curve"]``, ``model["anchors"]``,
``model["irt"]`` and ``model["detection_qc"]``, plus the appended ``iRT`` /
``Cal_RT_min`` columns of the result table.

The **detection ladder** changed meaning with the method. Under v1 it compared a
panel standard's RT in the samples against a reference baseline, because the
panel in the sample *was* the calibration. Under v2 the panel drives nothing:
the ladder is QC, and the question it answers is the only one still worth
asking — *did the calibration put your standards where the reference run has
them?* Both markers therefore live on **one** axis, the reference column's, with
your standards run pushed through the fitted curve first. The panel's nominal
manifest RT is deliberately **not** plotted: it comes from whichever method the
manifest was tabulated on (mix15's is an Orbitrap run, the bundled reference is
the Column 25 QTOF), so a bar drawn to it would show a method difference and read
as a calibration error.

Figure titles are set by the HTML report as copyable text, so the plotly figures
carry no baked-in title.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from . import theme

# plain-language caption for the detection ladder (shared by the HTML report and
# the desktop app, so the marker glyphs are always explained the same way)
DETECTION_NOTE = (
    "<b>QC only — this does not drive the calibration.</b> "
    "Each row is a standard of the chosen panel, on the <b>reference column's</b> "
    "time axis: <b>blue circle</b> = where it eluted in the reference standards run "
    "(these are the iRT landmarks), <b>clay square</b> = where it eluted in "
    "<i>your</i> standards run after the cross-column curve was applied. "
    "The bar between them is the disagreement; a row with one marker missing is a "
    "standard that was found in only one of the two runs."
)

#: Sentinel used where a number genuinely does not exist for this run.
DASH = "—"


# --------------------------------------------------------------- helpers ------

def _f(x, default=np.nan) -> float:
    """float() that survives ``None`` (the model writes null, not NaN)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return default
    return v if np.isfinite(v) else default


def _as_bool(s) -> pd.Series:
    return pd.Series(s).astype(str).str.strip().str.lower().isin(("true", "1", "1.0", "yes"))


def _fmt(v, d=2, suffix="") -> str:
    return DASH if not np.isfinite(_f(v)) else f"{_f(v):.{d}f}{suffix}"


def _rng(lo, hi, d=1) -> str:
    lo, hi = _f(lo), _f(hi)
    return DASH if not (np.isfinite(lo) and np.isfinite(hi)) else f"{lo:.{d}f}–{hi:.{d}f}"


def _panel_targets(result) -> pd.DataFrame:
    """The chosen panel's ion targets (name/class/mz/rt_ref_min), or an empty frame.

    ``panel="none"`` claims no identity for the user's mixture, so there are no
    targets — the caller falls back to the landmark set instead of inventing one.
    """
    stored = (result.panel or {}).get("targets")
    if isinstance(stored, pd.DataFrame) and len(stored):
        return stored.reset_index(drop=True)
    return pd.DataFrame(columns=["name", "class", "mz", "rt_ref_min"])


# --------------------------------------------------------------- metrics ------

def compute_metrics(result) -> Dict:
    """Every number the Overview and Detection sections need, in one dict.

    Reads the model blocks rather than recomputing anything: the report must
    agree with ``model.json`` by construction, not by coincidence.
    """
    model = result.model or {}
    curve = model.get("curve", {}) or {}
    anch = model.get("anchors", {}) or {}
    irt = model.get("irt", {}) or {}
    qc = (model.get("detection_qc", {}) or {}).get("user_standards_run", {}) or {}
    ref = model.get("reference", {}) or result.reference or {}
    pan = model.get("panel", {}) or {}
    cfg = model.get("config", {}) or {}

    tbl = result.table
    n = len(tbl)
    iv = result.values("iRT")
    cal_rt = result.values("Cal_RT_min")
    unc = result.values("Cal_RT_uncertainty_min")
    iunc = result.values("iRT_uncertainty")
    extrap = _as_bool(tbl[result.col("is_extrapolated")])
    rt = result.rt_minutes()

    conf = model.get("confidence_counts", {}) or {}
    conf = {k: int(conf.get(k, 0)) for k in ("high", "medium", "low", "none")}

    ladder = _ladder_rows(result)

    return dict(
        # ---- identity ----
        method=model.get("method", "cross-column-v2"),
        scope=model.get("calibration_scope", "project"),
        polarity=result.polarity, source_format=result.source_format,
        panel_key=pan.get("key", result.panel_key), panel_label=pan.get("label", ""),
        n_standards_panel=int(pan.get("n_standards", 0) or 0),
        reference_key=ref.get("key", ""), reference_label=ref.get("label", ""),
        reference_default=bool(ref.get("is_default", False)),
        reference_sample=ref.get("sample", ""), reference_standards=ref.get("standards", ""),

        # ---- stage 1: the curve ----
        n_pairs=int(curve.get("n_pairs", 0) or 0),
        n_pairs_kept=int(curve.get("n_pairs_kept", 0) or 0),
        n_pairs_standards=int(curve.get("n_pairs_standards", 0) or 0),
        n_pairs_sample=int(curve.get("n_pairs_sample", 0) or 0),
        curve_frac=_f(cfg.get("curve_frac")),
        rt_span_src=tuple(curve.get("rt_span_src_min", (np.nan, np.nan)) or (np.nan, np.nan)),
        rt_span_ref=tuple(curve.get("rt_span_ref_min", (np.nan, np.nan)) or (np.nan, np.nan)),
        resid_median=_f((curve.get("residual_min") or {}).get("median_abs")),
        resid_p90=_f((curve.get("residual_min") or {}).get("p90_abs")),

        # ---- stage 2: the anchors and the gate ----
        gate_engaged=bool(anch.get("engaged", False)),
        gate_reduction=_f(anch.get("gate_mse_reduction")),
        gate_threshold=_f(anch.get("gate_threshold"), 0.2),
        gate_reason=str(anch.get("gate_reason", "") or ""),
        class_aware=bool(anch.get("class_aware", False)),
        n_anchors_used=int(anch.get("n_used", 0) or 0),
        n_anchors_validated=int(anch.get("n_validated", 0) or 0),
        lam_g=_f(anch.get("lam_g")), lam_c=_f(anch.get("lam_c")),
        sigma_anchor=_f(anch.get("sigma_anchor_min"), 0.0),
        loo_median=_f((anch.get("loo_residual_min") or {}).get("median")),
        loo_p90=_f((anch.get("loo_residual_min") or {}).get("p90")),
        loo_max=_f((anch.get("loo_residual_min") or {}).get("max")),

        # ---- the iRT ruler ----
        n_landmarks=int(irt.get("n_landmarks", 0) or 0),
        landmark_rt_span=tuple(irt.get("landmark_rt_span_min", (np.nan, np.nan))
                               or (np.nan, np.nan)),
        irt_scale_range=tuple(irt.get("irt_range", (np.nan, np.nan)) or (np.nan, np.nan)),
        irt_definition=irt.get("definition", "affine on the reference-column RT axis"),
        irt_panel_used=irt.get("panel_used", ""), irt_reason=irt.get("reason", ""),

        # ---- what came out ----
        n_features=n,
        rt_range_samples=(_f(rt.min()), _f(rt.max())),
        cal_rt_range=(_f(cal_rt.min()), _f(cal_rt.max())),
        iRT_range=(_f(iv.min()), _f(iv.max())),
        coverage=(float(iv.notna().mean()) if n else 0.0),
        pct_extrap=(100.0 * float(extrap.sum()) / n if n else 0.0),
        unc_median=_f(unc.median()), irt_unc_median=_f(iunc.median()),
        reliability=conf,
        frac_high=(conf["high"] / n if n else 0.0),
        warp_source=("curve+anchors" if anch.get("engaged") else "curve"),

        # ---- detection QC on the USER's standards run ----
        n_panel=int(qc.get("n_panel", 0) or 0),
        n_detected=int(qc.get("n_detected", 0) or 0),
        native_rt=dict(qc.get("native_rt", {}) or {}),
        rt_range_panel=tuple(qc.get("rt_range_min", (np.nan, np.nan)) or (np.nan, np.nan)),
        n_features_std=qc.get("n_features"),
        qc_note=qc.get("note", "QC only — does not drive the calibration"),

        # ---- the ladder ----
        names=ladder["names"], ref_run_rt=ladder["ref_run_rt"],
        user_cal_rt=ladder["user_cal_rt"], user_raw_rt=ladder["user_raw_rt"],
        ladder_note=ladder["note"],
    )


def _ladder_rows(result) -> Dict:
    """Rows of the detection ladder, all on the reference column's time axis.

    ``ref_run_rt`` is the landmark's RT in the *reference* standards run;
    ``user_cal_rt`` is its RT in the *user's* standards run after the fitted
    curve. ``user_raw_rt`` is kept for the hover, because "22.4 min on your
    column" is the number a chromatographer recognises.
    """
    landmarks = result.landmarks if result.landmarks is not None else pd.DataFrame()
    lm_rt: Dict[str, float] = {}
    if len(landmarks) and "name" in landmarks.columns:
        lm_rt = {str(r["name"]): _f(r["rt_ref_run_min"]) for _, r in landmarks.iterrows()}

    native = {str(k): _f(v) for k, v in ((result.panel or {}).get("native_rt", {}) or {}).items()}

    targets = _panel_targets(result)
    order: List[str] = []
    if len(targets):
        order = [str(v) for v in targets.sort_values("rt_ref_min")["name"]]
    for nm in sorted(lm_rt, key=lambda k: lm_rt[k]):
        if nm not in order:
            order.append(nm)
    for nm in sorted(native, key=lambda k: native[k]):
        if nm not in order:
            order.append(nm)

    cal_rt: Dict[str, float] = {}
    if native and result.calibrator is not None:
        names = list(native)
        pred = result.calibrator.predict(np.array([native[k] for k in names], dtype=float))
        cal_rt = {nm: float(v) for nm, v in zip(names, pred)}

    # sort rows by whatever position we actually know, so the ladder reads in
    # elution order even when one of the two runs missed a standard
    def key(nm: str) -> float:
        for src in (lm_rt, cal_rt):
            if nm in src and np.isfinite(src[nm]):
                return src[nm]
        return np.inf

    order = sorted(order, key=key)
    note = ""
    if not native:
        note = ("no user panel selected — detection QC skipped; rows show the "
                "fallback landmark panel located on the reference standards run")
    return dict(names=order, ref_run_rt=lm_rt, user_cal_rt=cal_rt,
                user_raw_rt=native, note=note)


# ------------------------------------------------------------------ tiles -----

def gate_headline(m: Dict) -> Tuple[str, str]:
    """``(value, note)`` for the anchor-gate tile — the same words everywhere."""
    red = m["gate_reduction"]
    pct = f"{100 * red:.0f}%" if np.isfinite(red) else DASH
    if m["gate_engaged"]:
        return (f"on · {pct}",
                f"λg {m['lam_g']:.1f} · λc {m['lam_c']:.1f} · {m['n_anchors_used']} anchors")
    if m["n_anchors_validated"] or m["n_anchors_used"]:
        return (f"off · {pct}",
                f"below the {100 * m['gate_threshold']:.0f}% gate · stage-1 curve only")
    return ("off", "no plasma-lipid anchors · stage-1 curve only")


def kpi_tiles(result) -> List[Tuple[str, object, str]]:
    """Six KPI tiles for the v2 method.

    A tile is ``(label, value, note)``. ``value`` is either a headline string or
    a list of ``(caption, value)`` pairs rendered as stacked rows (used where two
    quantities belong together — the two runs, the two pair sources, the median
    and the P90).
    """
    m = compute_metrics(result)
    gate_v, gate_n = gate_headline(m)
    nf_std = m["n_features_std"]

    return [
        ("Features", [("samples", f"{m['n_features']:,}"),
                      ("standards run", f"{nf_std:,}" if nf_std is not None else DASH)],
         "feature rows in each table"),
        ("Matched pairs", [("standards", f"{m['n_pairs_standards']:,}"),
                           ("sample", f"{m['n_pairs_sample']:,}")],
         f"{m['n_pairs_kept']:,} of {m['n_pairs']:,} kept after MAD trimming"),
        ("Curve residual", [("median", _fmt(m["resid_median"], 3)),
                            ("P90", _fmt(m["resid_p90"], 3))],
         "|reference − fitted| (min)"),
        ("Anchor gate", gate_v, gate_n),
        ("iRT landmarks", f"{m['n_landmarks']}",
         f"RT {_rng(*m['landmark_rt_span'], d=2)} min on the reference run"),
        ("Output", [("iRT range", _rng(*m["iRT_range"], d=0)),
                    ("extrapolated", f"{m['pct_extrap']:.1f}%")],
         f"{100 * m['frac_high']:.0f}% high-reliability"),
    ]


# ---------------------------------------------------------------- radar --------

def radar_axes(result) -> List[Tuple[str, float]]:
    """Normalised [0,1] quality axes (outward = better).

    The scalings are deliberate, not cosmetic:

    * **Pair support** — 400 kept pairs is a well-supported curve; the validated
      columns land at 425–566.
    * **Curve fit** — 0.5 min median residual is where a 35-min gradient's
      calibration stops being useful.
    * **Anchors** — out of the 17-lipid plasma panel, and zero when the gate is
      off, because an anchor that was not used is not support.
    * **Landmarks** — out of the chosen panel's ionisable standards.
    * **Reliability** — the share of features rated ``high`` (per-project), or
      the inter-injection agreement (per-sample).
    """
    m = compute_metrics(result)
    pairs = min(m["n_pairs_kept"] / 400.0, 1.0)
    resid = m["resid_median"]
    fit = 1.0 - min((resid if np.isfinite(resid) else 0.5) / 0.5, 1.0)
    anchors = (min(m["n_anchors_used"] / 17.0, 1.0) if m["gate_engaged"] else 0.0)
    denom = max(m["n_panel"] or m["n_standards_panel"], 1)
    landmarks = min(m["n_landmarks"] / denom, 1.0)

    if m["scope"] == "sample":
        sp = result.values("RI_spread")
        med = float(sp.median()) if sp.notna().any() else np.nan
        consistency = 1.0 - min((med if np.isfinite(med) else 1.0) / 1.0, 1.0)
        clabel = "Repeatability"
    else:
        consistency, clabel = m["frac_high"], "Reliability"
    return [("Pair support", pairs), ("Curve fit", max(fit, 0.0)),
            ("Anchors", anchors), ("Landmarks", landmarks),
            (clabel, max(min(consistency, 1.0), 0.0))]


def radar_plotly(result):
    import plotly.graph_objects as go
    axes = radar_axes(result)
    labels = [a for a, _ in axes]; vals = [round(v, 3) for _, v in axes]
    fig = go.Figure(go.Scatterpolar(
        r=vals + [vals[0]], theta=labels + [labels[0]], fill="toself",
        fillcolor=theme.rgba(theme.PRIMARY, 0.22), line=dict(color=theme.PRIMARY, width=2),
        marker=dict(size=6, color=theme.PRIMARY),
        hovertemplate="%{theta}: %{r:.2f}<extra></extra>"))
    fig.update_layout(
        paper_bgcolor=theme.PAPER, font=dict(family=theme.FONT_STACK, size=13, color=theme.TXT2),
        margin=dict(l=48, r=48, t=28, b=28), height=320, showlegend=False,
        polar=dict(bgcolor=theme.PANEL,
                   radialaxis=dict(range=[0, 1], tickvals=[0.25, 0.5, 0.75, 1.0],
                                   tickfont=dict(size=10, color=theme.FAINT), gridcolor=theme.GRID,
                                   angle=90, tickangle=90),
                   angularaxis=dict(gridcolor=theme.GRID, tickfont=dict(size=13, color=theme.TXT))))
    return fig


def radar_mpl(result, ax=None):
    import matplotlib.pyplot as plt
    axes = radar_axes(result)
    labels = [a for a, _ in axes]; vals = [v for _, v in axes]
    ang = np.linspace(0, 2 * np.pi, len(labels), endpoint=False)
    ang = np.r_[ang, ang[:1]]; vv = np.r_[vals, vals[:1]]
    if ax is None:
        fig = plt.figure(figsize=(4.6, 4.6)); ax = fig.add_subplot(111, polar=True)
    ax.plot(ang, vv, color=theme.PRIMARY, lw=2)
    ax.fill(ang, vv, color=theme.PRIMARY, alpha=0.22)
    ax.set_xticks(ang[:-1]); ax.set_xticklabels(labels, fontsize=theme.FS_SUB, color=theme.TXT)
    ax.set_yticks([0.5, 1.0]); ax.set_yticklabels(["0.5", "1"], fontsize=theme.FS_SUB - 3, color=theme.FAINT)
    ax.set_ylim(0, 1); ax.grid(color=theme.GRID)
    ax.spines["polar"].set_color(theme.GRID)
    return ax.figure


# ---------------------------------------------------------------- ladder -------

C_REF_RUN = theme.PRIMARY       # the reference standards run
C_USER = theme.ACCENT           # the user's standards run, calibrated


def figure_mpl(result):
    import matplotlib.pyplot as plt
    m = compute_metrics(result)
    names = m["names"]
    height = max(3.2, min(0.34 * max(len(names), 1) + 1.9, 9.0))
    fig, ax = plt.subplots(figsize=(9.5, height))
    _ladder_mpl(ax, m)
    fig.subplots_adjust(left=0.24, right=0.93, top=0.94, bottom=0.10 if height > 5 else 0.16)
    return fig


def _ladder_mpl(ax, m: Dict) -> None:
    names = m["names"][::-1]
    if not names:
        ax.text(0.5, 0.5, "no panel standards to show", ha="center", va="center",
                color=theme.TXT2, fontsize=theme.FS_SUB, transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ("top", "right", "left", "bottom"):
            ax.spines[s].set_visible(False)
        return

    for i, nm in enumerate(names):
        ref = m["ref_run_rt"].get(nm)
        usr = m["user_cal_rt"].get(nm)
        xs = [v for v in (ref, usr) if v is not None and np.isfinite(v)]
        if len(xs) > 1:
            ax.plot([min(xs), max(xs)], [i, i], color=theme.GRID, lw=6,
                    solid_capstyle="round", zorder=1)
        if ref is not None and np.isfinite(ref):
            ax.scatter([ref], [i], s=54, color=C_REF_RUN, zorder=5,
                       label="reference standards run" if i == 0 else "")
        if usr is not None and np.isfinite(usr):
            ax.scatter([usr], [i], s=36, color=C_USER, marker="s", zorder=4,
                       label="your standards run (calibrated)" if i == 0 else "")
        if len(xs) > 1:
            ax.annotate(f"{usr - ref:+.2f}", (max(xs), i), xytext=(9, 0),
                        textcoords="offset points", va="center",
                        fontsize=theme.FS_SUB - 3.5, color=theme.TXT2)

    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([n.replace("_", " ") for n in names], fontsize=theme.FS_SUB - 1)
    ax.set_ylim(-0.7, len(names) - 0.3)
    theme.style_ax(ax)
    ax.set_xlabel("Retention time on the reference column (min)")
    ax.legend(loc="lower right", fontsize=theme.FS_SUB - 1, handletextpad=0.3,
              borderaxespad=0.6)


def figure_plotly(result):
    """Detection ladder only (no title / panel letter — the report supplies the title).

    Traces, in order: one grey connector per row, then the reference-run marker
    and the calibrated user-run marker. ``customdata`` carries the standard's
    name so the report's structure-on-hover can key off it.
    """
    import plotly.graph_objects as go
    m = compute_metrics(result)
    fig = go.Figure()
    names = m["names"][::-1]
    for i, nm in enumerate(names):
        ref, usr = m["ref_run_rt"].get(nm), m["user_cal_rt"].get(nm)
        raw = m["user_raw_rt"].get(nm)
        xs = [v for v in (ref, usr) if v is not None and np.isfinite(v)]
        if len(xs) > 1:
            fig.add_trace(go.Scatter(x=[min(xs), max(xs)], y=[i, i], mode="lines",
                          line=dict(color=theme.GRID, width=6), hoverinfo="skip",
                          showlegend=False))
        if ref is not None and np.isfinite(ref):
            fig.add_trace(go.Scatter(
                x=[ref], y=[i], mode="markers", legendgroup="ref",
                marker=dict(size=11, color=C_REF_RUN), name="reference standards run",
                showlegend=(i == 0), customdata=[nm],
                hovertemplate=f"{nm}<br>reference standards run %{{x:.2f}} min<extra></extra>"))
        if usr is not None and np.isfinite(usr):
            extra = f" (raw {raw:.2f} min)" if raw is not None and np.isfinite(raw) else ""
            fig.add_trace(go.Scatter(
                x=[usr], y=[i], mode="markers", legendgroup="usr",
                marker=dict(size=9, color=C_USER, symbol="square"),
                name="your standards run (calibrated)", showlegend=(i == 0),
                customdata=[nm],
                hovertemplate=f"{nm}<br>your run, calibrated %{{x:.2f}} min{extra}<extra></extra>"))
    fig.update_layout(theme.plotly_template())
    fig.update_layout(height=max(320, min(26 * len(names) + 140, 900)),
                      showlegend=True, title=None,
                      margin=dict(l=70, r=24, t=36, b=48),
                      yaxis=dict(tickmode="array", tickvals=list(range(len(names))),
                                 ticktext=[n.replace("_", " ") for n in names],
                                 showgrid=False, automargin=True),
                      xaxis=dict(title=dict(text="Retention time on the reference column (min)"),
                                 automargin=True))
    if not names:
        fig.add_annotation(text="no panel standards to show", showarrow=False,
                           x=0.5, y=0.5, xref="paper", yref="paper",
                           font=dict(size=14, color=theme.TXT2))
    return fig


def _rgba(hex_color: str, a: float) -> str:
    """Backwards-compatible alias of :func:`theme.rgba`."""
    return theme.rgba(hex_color, a)
