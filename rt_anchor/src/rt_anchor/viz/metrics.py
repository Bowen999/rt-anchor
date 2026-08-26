"""Module 1 — key metrics: KPI tiles, a quality radar, and the detection ladder.

KPI tiles (features, standards, RT/iRT range, RT offset, coverage) + a radar
summarising calibration quality across a few normalised axes (outward = better) +
a reference-vs-panel-vs-samples standard-detection ladder. Figure titles are set
by the HTML report as copyable text, so the plotly figures carry no baked-in title.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from . import theme

# plain-language caption for the detection ladder (shared by the HTML report and
# the desktop app, so the marker glyphs are always explained the same way)
DETECTION_NOTE = (
    "Each row is a standard; its retention time is shown from three sources — "
    "<b>hollow circle</b> = reference baseline, "
    "<b>blue square</b> = standards run, "
    "<b>blue circle</b> = samples (the table being calibrated). "
    "The bar links the three; a shorter bar means the observed RT agrees more closely with the reference."
)


def _panel_targets(result):
    stored = (result.panel or {}).get("targets")
    if isinstance(stored, pd.DataFrame) and len(stored):
        return stored.reset_index(drop=True)
    from ..config import CalibrationConfig
    from ..panel import build_panel
    cfg = CalibrationConfig.from_dict(result.model.get("config", {}))
    return build_panel(result.polarity or "positive", cfg).targets


def compute_metrics(result) -> Dict:
    tbl = result.table
    ri = pd.to_numeric(tbl[result.col("RI")], errors="coerce")
    extrap = _as_bool(tbl[result.col("is_extrapolated")])
    rt = result.rt_minutes()

    targets = _panel_targets(result)
    ref_rt = dict(zip(targets["name"], targets["rt_ref_min"]))
    names = list(targets["name"])
    n_panel = len(names)

    amap = (result.anchors.drop_duplicates("name").set_index("name")
            if len(result.anchors) else pd.DataFrame())
    sample_rt = {nm: float(amap.loc[nm, "rt_obs_min"]) for nm in names if nm in amap.index}
    panel_info = result.panel or {}
    panel_rt = panel_info.get("native_rt", {})

    offsets = [abs(sample_rt[nm] - ref_rt[nm]) for nm in sample_rt]
    n = len(tbl)
    loo = result.model.get("loo_residual_irt", {})
    anc_irt = result.anchors["irt"] if len(result.anchors) else pd.Series(dtype=float)

    return dict(
        n_features=n, scope=result.model.get("calibration_scope", "project"),
        names=names, n_panel=n_panel, ref_rt=ref_rt, sample_rt=sample_rt, panel_rt=panel_rt,
        n_sample_detected=len(sample_rt), n_panel_detected=panel_info.get("n_detected", len(panel_rt)),
        n_features_std=panel_info.get("n_features"),
        rt_range_samples=(float(rt.min()), float(rt.max())),
        rt_range_panel=tuple(panel_info.get("rt_range", (np.nan, np.nan))),
        iRT_range=(float(ri.min()), float(ri.max())),
        anchor_iRT_span=(float(anc_irt.max() - anc_irt.min()) if len(anc_irt) else np.nan),
        coverage=float(ri.notna().mean()) if n else 0.0,
        offset_median=(float(np.median(offsets)) if offsets else None),
        offset_max=(float(np.max(offsets)) if offsets else None),
        loo_median=loo.get("median"),
        pct_extrap=100.0 * float(np.sum(extrap)) / n if n else 0.0,
        polarity=result.polarity, source_format=result.source_format,
    )


def _as_bool(s):
    return pd.Series(s).astype(str).str.strip().str.lower().isin(("true", "1", "1.0", "yes"))


def kpi_tiles(result) -> List[Tuple[str, object, str]]:
    """Key data (6 tiles) — deliberately excludes calibration-scope and reliability.

    A tile is ``(label, value, note)``. ``value`` is either a headline string or a
    list of ``(caption, value)`` pairs rendered as stacked rows (used where the
    samples table and the standards run are reported side by side).
    """
    m = compute_metrics(result)

    def rng(r, d=1):
        if r is None or (isinstance(r, tuple) and not np.isfinite(r[0])):
            return "—"
        return f"{r[0]:.{d}f}–{r[1]:.{d}f}"

    offset = "—" if m["offset_median"] is None else f"{m['offset_median']:.2f}"
    nf_std = m.get("n_features_std")
    return [
        ("Features", [("samples", f"{m['n_features']:,}"),
                      ("standards run", f"{nf_std:,}" if nf_std is not None else "—")],
         "feature rows in each table"),
        ("Standards", [("samples", f"{m['n_sample_detected']}/{m['n_panel']}"),
                       ("standards run", f"{m['n_panel_detected']}/{m['n_panel']}")],
         f"detected · panel of {m['n_panel']}"),
        ("RT range", [("samples", rng(m["rt_range_samples"])),
                      ("standards run", rng(m["rt_range_panel"]))], "minutes"),
        ("iRT range", f"{m['iRT_range'][0]:.0f}–{m['iRT_range'][1]:.0f}", "dimensionless"),
        ("RT offset", offset, "median |obs−ref| (min)"),
        ("Coverage", f"{m['coverage']*100:.0f}%", "features with an RI"),
    ]


# ---------------------------------------------------------------- radar --------

def radar_axes(result) -> List[Tuple[str, float]]:
    """Normalised [0,1] quality axes (outward = better)."""
    m = compute_metrics(result)
    detection = m["n_sample_detected"] / max(m["n_panel"], 1)
    coverage = m["coverage"]
    span = (m["anchor_iRT_span"] / 100.0) if np.isfinite(m["anchor_iRT_span"]) else 0.0
    loo = m["loo_median"]
    interp = 1.0 - min((loo or 0.0) / 3.0, 1.0)   # 3 iRT = the low-reliability edge
    if m["scope"] == "sample":
        sp = pd.to_numeric(result.table[result.col("RI_spread")], errors="coerce")
        med = float(sp.median()) if sp.notna().any() else np.nan
        consistency = 1.0 - min((med if np.isfinite(med) else 1.0) / 1.0, 1.0)
        clabel = "Repeatability"
    else:
        # panel↔samples RT agreement
        diffs = [abs(m["panel_rt"][nm] - m["sample_rt"][nm])
                 for nm in m["sample_rt"] if nm in m["panel_rt"]]
        med = float(np.median(diffs)) if diffs else 0.0
        consistency = 1.0 - min(med / 0.5, 1.0)
        clabel = "Panel match"
    return [("Detection", detection), ("Coverage", coverage), ("Anchor span", min(span, 1.0)),
            ("Interpolation", interp), (clabel, max(consistency, 0.0))]


def radar_plotly(result):
    import plotly.graph_objects as go
    axes = radar_axes(result)
    labels = [a for a, _ in axes]; vals = [round(v, 3) for _, v in axes]
    fig = go.Figure(go.Scatterpolar(
        r=vals + [vals[0]], theta=labels + [labels[0]], fill="toself",
        fillcolor=_rgba(theme.PRIMARY, 0.22), line=dict(color=theme.PRIMARY, width=2),
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

def figure_mpl(result):
    import matplotlib.pyplot as plt
    m = compute_metrics(result)
    fig, ax = plt.subplots(figsize=(9.5, 6.4))
    _ladder_mpl(ax, m)
    fig.subplots_adjust(left=0.22, right=0.97, top=0.95, bottom=0.1)
    return fig


def _ladder_mpl(ax, m):
    names = m["names"][::-1]
    b1 = theme.blue_ramp(3)[1]
    for i, nm in enumerate(names):
        ref = m["ref_rt"].get(nm); p = m["panel_rt"].get(nm); s = m["sample_rt"].get(nm)
        xs = [v for v in (ref, p, s) if v is not None and np.isfinite(v)]
        if xs:
            ax.plot([min(xs), max(xs)], [i, i], color=theme.GRID, lw=6, solid_capstyle="round", zorder=1)
        if ref is not None:
            ax.scatter([ref], [i], s=64, facecolor="none", edgecolor=theme.NEUTRAL, lw=1.4,
                       zorder=3, label="reference" if i == 0 else "")
        if p is not None:
            ax.scatter([p], [i], s=30, color=b1, marker="s", zorder=4, label="standards run" if i == 0 else "")
        if s is not None:
            ax.scatter([s], [i], s=52, color=theme.PRIMARY, zorder=5, label="samples" if i == 0 else "")
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=theme.FS_SUB)
    ax.set_ylim(-0.6, len(names) - 0.4); theme.style_ax(ax)
    ax.set_xlabel("Retention time (min)")
    ax.legend(loc="upper right", fontsize=theme.FS_SUB - 1, handletextpad=0.3, borderaxespad=0.6)


def figure_plotly(result):
    """Detection ladder only (no title / panel letter — the report supplies the title)."""
    import plotly.graph_objects as go
    m = compute_metrics(result)
    b1 = theme.blue_ramp(3)[1]
    fig = go.Figure()
    names = m["names"][::-1]
    for i, nm in enumerate(names):
        ref, p, s = m["ref_rt"].get(nm), m["panel_rt"].get(nm), m["sample_rt"].get(nm)
        xs = [v for v in (ref, p, s) if v is not None and np.isfinite(v)]
        if xs:
            fig.add_trace(go.Scatter(x=[min(xs), max(xs)], y=[i, i], mode="lines",
                          line=dict(color=theme.GRID, width=6), hoverinfo="skip", showlegend=False))
        if ref is not None:
            fig.add_trace(go.Scatter(x=[ref], y=[i], mode="markers", legendgroup="ref",
                          marker=dict(size=11, color="rgba(0,0,0,0)", line=dict(color=theme.NEUTRAL, width=1.6)),
                          name="reference", showlegend=(i == 0), customdata=[nm],
                          hovertemplate=f"{nm}<br>reference %{{x:.2f}} min<extra></extra>"))
        if p is not None:
            fig.add_trace(go.Scatter(x=[p], y=[i], mode="markers", legendgroup="panel",
                          marker=dict(size=8, color=b1, symbol="square"), name="standards run", showlegend=(i == 0),
                          customdata=[nm], hovertemplate=f"{nm}<br>standards run %{{x:.2f}} min<extra></extra>"))
        if s is not None:
            fig.add_trace(go.Scatter(x=[s], y=[i], mode="markers", legendgroup="samp",
                          marker=dict(size=11, color=theme.PRIMARY), name="samples", showlegend=(i == 0),
                          customdata=[nm], hovertemplate=f"{nm}<br>samples %{{x:.2f}} min<extra></extra>"))
    fig.update_layout(theme.plotly_template())
    fig.update_layout(height=520, showlegend=True, title=None,
                      margin=dict(l=70, r=24, t=36, b=48),
                      yaxis=dict(tickmode="array", tickvals=list(range(len(names))), ticktext=names,
                                 showgrid=False, automargin=True),
                      xaxis=dict(title=dict(text="Retention time (min)"), automargin=True))
    return fig


def _rgba(hex_color: str, a: float) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{a})"
