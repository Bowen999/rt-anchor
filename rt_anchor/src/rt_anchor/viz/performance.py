"""Module 3 — the cross-column calibration curve and its residuals.

This is the report's evidence page, and it is a direct port of the two-panel
diagnostic the v2 method was validated with
(``results/rt_calibration/rt_calibration_diagnostics.png``):

**upper panel** — every stage-1 matched pair, with the *standards*-run pairs and
the *sample*-run pairs drawn differently because they answer different
questions (the mixture pins the middle of the gradient, serum pins the sparse
early and late ends); the pairs MAD-trimming rejected, drawn as hollow grey so
they are visible as *excluded* rather than quietly deleted; the fitted anchor-free
curve; the anchor-refined curve when the stage-2 gate engaged; and the stage-2
anchors as rings.

**lower panel** — the residuals of those same pairs about the curve, before and
after the anchor correction. It is the panel that shows whether the correction
did anything, which is the only honest way to present a correction that is
sometimes gated off.

When the gate is **off**, the figure says so with the measured LOO gain and the
threshold it fell short of, exactly as the reference diagnostic does — a report
that silently omitted the correction would leave the reader thinking it applied.

No baked-in title — the HTML report supplies a copyable one.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd

from . import theme

C_STD = theme.PRIMARY        # standards-run pairs
C_SAMP = theme.PRIMARY_LT    # sample-run pairs
C_FIT = theme.TXT            # the anchor-free stage-1 curve
C_REFINED = theme.ACCENT     # the anchor-refined curve
C_RING = theme.RING          # stage-2 anchors
C_REJECT = theme.REJECT      # pairs dropped by MAD trimming


def compute_curve(result, ngrid: int = 400) -> Dict:
    """Everything both renderers need, computed once.

    ``empty`` is True when there is no calibrator or too few pairs to draw; the
    renderers then produce a figure that says so rather than an empty axis.
    """
    cal = result.calibrator
    curve = result.curve
    pairs = result.pairs if result.pairs is not None else pd.DataFrame()
    if cal is None or curve is None or not len(pairs):
        return dict(empty=True, reason="no fitted curve for this run")

    rt_src = pd.to_numeric(pairs["rt_src"], errors="coerce").to_numpy(dtype=float)
    rt_ref = pd.to_numeric(pairs["rt_ref"], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(rt_src) & np.isfinite(rt_ref)
    if ok.sum() < 2:
        return dict(empty=True, reason="fewer than two usable matched pairs")

    src = np.asarray(pairs["source"], dtype=object)[ok]
    kept = np.asarray(pairs["kept"]).astype(bool)[ok] if "kept" in pairs.columns \
        else np.ones(int(ok.sum()), dtype=bool)
    rt_src, rt_ref = rt_src[ok], rt_ref[ok]

    lo, hi = float(np.nanmin(rt_src)), float(np.nanmax(rt_src))
    grid = np.linspace(lo, hi, ngrid)
    fit = curve.predict(grid)

    engaged = bool(cal.anchors_used)
    refined = cal.predict(grid) if engaged else None

    # residuals: observed reference RT minus what the model says (same sign
    # convention as ColumnCalibrator.residuals(), so the numbers in the figure
    # and the numbers in model.json are the same numbers)
    r_curve = rt_ref - curve.predict(rt_src)
    r_refined = (rt_ref - cal.predict(rt_src)) if engaged else None

    anchors = result.anchors if result.anchors is not None else pd.DataFrame()
    if len(anchors) and "dropped_by_sanity_filter" in anchors.columns:
        anchors = anchors[~anchors["dropped_by_sanity_filter"].astype(bool)]

    return dict(
        empty=False, grid=grid, fit=fit, refined=refined,
        rt_src=rt_src, rt_ref=rt_ref, source=src, kept=kept,
        is_std=(src == "standards"), is_samp=(src != "standards"),
        r_curve=r_curve, r_refined=r_refined,
        anchors=anchors, engaged=engaged,
        lam_g=float(cal.lam_g), lam_c=float(cal.lam_c),
        class_aware=bool(cal.class_aware),
        gate_reduction=float(cal.gate_mse_reduction),
        gate_threshold=float(cal.gate_threshold),
        gate_reason=str(cal.gate_reason or ""),
        source_label="your column", ref_label="the reference column",
    )


#: Kept so callers written against the v1 module keep working; the v1 warp it
#: used to describe no longer exists, so it returns the v2 curve diagnostic.
compute_warp = compute_curve


def anchor_legend_label(d: Dict) -> str:
    """The anchor legend entry, which is where the gate state gets said out loud."""
    if d.get("engaged"):
        return "stage-2 anchors"
    red = d.get("gate_reduction", np.nan)
    thr = 100 * d.get("gate_threshold", 0.2)
    if np.isfinite(red):
        return f"stage-2 anchors (gated off: LOO gain {100 * red:.0f}% < {thr:.0f}%)"
    return "stage-2 anchors (not used)"


def fit_legend_label(d: Dict) -> str:
    kind = "class-aware" if d.get("class_aware") else "shrunk"
    return (f"anchor-refined fit · {kind} "
            f"(λg={d['lam_g']:.1f}, λc={d['lam_c']:.1f})")


def _resid_limits(r: np.ndarray) -> tuple:
    """Percentile limits with headroom, so the legend never sits on the data."""
    r = r[np.isfinite(r)]
    if r.size == 0:
        return (-1.0, 1.0)
    lo, hi = np.percentile(r, [0.5, 99.5])
    if hi <= lo:
        pad = max(abs(hi), 0.1)
        return (lo - pad, hi + pad)
    return (lo - 0.2 * (hi - lo), hi + 0.6 * (hi - lo))


# ---------------------------------------------------------------- matplotlib --

def figure_mpl(result):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    d = compute_curve(result)
    fig = plt.figure(figsize=(11.0, 8.0))
    if d.get("empty"):
        fig.text(0.5, 0.5, d.get("reason", "no calibration curve available"),
                 ha="center", va="center", color=theme.TXT2, fontsize=theme.FS_AXIS)
        return fig

    gs = GridSpec(2, 1, height_ratios=[2.05, 1.0], hspace=0.16,
                  left=0.085, right=0.975, top=0.93, bottom=0.085)
    ax = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)

    # ---- upper: pairs + curves + anchors ----
    drop = ~d["kept"]
    if drop.any():
        ax.scatter(d["rt_src"][drop], d["rt_ref"][drop], s=22, facecolor="none",
                   edgecolor=C_REJECT, lw=0.9, zorder=2,
                   label=f"rejected by MAD trimming ({int(drop.sum())})")
    for mask, colour, size, lab in (
            (d["is_samp"] & d["kept"], C_SAMP, 9, "sample-run pairs"),
            (d["is_std"] & d["kept"], C_STD, 13, "standards-run pairs")):
        if mask.any():
            ax.scatter(d["rt_src"][mask], d["rt_ref"][mask], s=size, alpha=0.55,
                       color=colour, lw=0, zorder=3,
                       label=f"{lab} ({int(mask.sum())})")
    ax.plot(d["grid"], d["fit"], color=C_FIT, lw=1.7, zorder=5, label="anchor-free fit")
    if d["refined"] is not None:
        ax.plot(d["grid"], d["refined"], color=C_REFINED, lw=1.7, ls=(0, (5, 3)),
                zorder=6, label=fit_legend_label(d))
    anc = d["anchors"]
    if len(anc):
        ax.scatter(pd.to_numeric(anc["rt_src"], errors="coerce"),
                   pd.to_numeric(anc["rt_ref"], errors="coerce"),
                   s=64, facecolor="none", edgecolor=C_RING, lw=1.5, zorder=7,
                   label=anchor_legend_label(d))
    theme.style_ax(ax)
    ax.set_ylabel("RT on the reference column (min)")
    ax.legend(loc="upper left", fontsize=theme.FS_SUB - 1.5, handletextpad=0.4,
              labelspacing=0.35, borderaxespad=0.5)
    ax.tick_params(labelbottom=False)

    # ---- lower: residuals ----
    axr.axhline(0, color=theme.AXIS, lw=0.8, zorder=1)
    axr.scatter(d["rt_src"], d["r_curve"], s=9, alpha=0.45, color=theme.NEUTRAL,
                lw=0, zorder=3, label="about the anchor-free fit")
    if d["r_refined"] is not None:
        axr.scatter(d["rt_src"], d["r_refined"], s=9, alpha=0.55, color=C_REFINED,
                    lw=0, zorder=4, label="about the anchor-refined fit")
    axr.set_ylim(*_resid_limits(d["r_curve"]))
    theme.style_ax(axr)
    axr.set_xlabel("RT on your column (min)")
    axr.set_ylabel("residual (min)")
    axr.legend(loc="upper right", fontsize=theme.FS_SUB - 1.5, handletextpad=0.4,
               labelspacing=0.3, borderaxespad=0.5, ncol=2)

    fig.suptitle("Cross-column calibration curve  ·  your column -> the reference column",
                 x=0.02, y=0.985, ha="left", fontsize=theme.FS_TITLE, color=theme.TXT)
    if d["gate_reason"]:
        fig.text(0.02, 0.947, theme.mpl_safe(d["gate_reason"]), ha="left", va="top",
                 fontsize=theme.FS_SUB - 1, color=theme.TXT2, style="italic")
    fig.subplots_adjust(top=0.90)
    return fig


# ---------------------------------------------------------------- plotly ------

def figure_plotly(result):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    d = compute_curve(result)
    if d.get("empty"):
        fig = go.Figure()
        fig.update_layout(theme.plotly_template())
        fig.add_annotation(text=d.get("reason", "no calibration curve available"),
                           showarrow=False, x=0.5, y=0.5, xref="paper", yref="paper",
                           font=dict(size=15, color=theme.TXT2))
        fig.update_layout(height=560, title=None,
                          xaxis=dict(visible=False), yaxis=dict(visible=False))
        return fig

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.66, 0.34],
                        vertical_spacing=0.06)

    drop = ~d["kept"]
    if drop.any():
        fig.add_trace(go.Scatter(
            x=d["rt_src"][drop], y=d["rt_ref"][drop], mode="markers",
            name=f"rejected by MAD trimming ({int(drop.sum())})",
            marker=dict(size=7, color="rgba(0,0,0,0)",
                        line=dict(color=C_REJECT, width=1.1)),
            hovertemplate="rejected<br>%{x:.2f} → %{y:.2f} min<extra></extra>"), 1, 1)
    for mask, colour, size, lab in (
            (d["is_samp"] & d["kept"], C_SAMP, 5, "sample-run pairs"),
            (d["is_std"] & d["kept"], C_STD, 6.5, "standards-run pairs")):
        if mask.any():
            fig.add_trace(go.Scatter(
                x=d["rt_src"][mask], y=d["rt_ref"][mask], mode="markers",
                name=f"{lab} ({int(mask.sum())})",
                marker=dict(size=size, color=theme.rgba(colour, 0.6)),
                hovertemplate=f"{lab}<br>%{{x:.2f}} → %{{y:.2f}} min<extra></extra>"), 1, 1)
    fig.add_trace(go.Scatter(x=d["grid"], y=d["fit"], mode="lines", name="anchor-free fit",
                             line=dict(color=C_FIT, width=2.2),
                             hovertemplate="RT %{x:.2f} → %{y:.2f} min<extra>anchor-free</extra>"),
                  1, 1)
    if d["refined"] is not None:
        fig.add_trace(go.Scatter(x=d["grid"], y=d["refined"], mode="lines",
                                 name=fit_legend_label(d),
                                 line=dict(color=C_REFINED, width=2.2, dash="dash"),
                                 hovertemplate="RT %{x:.2f} → %{y:.2f} min<extra>refined</extra>"),
                      1, 1)
    anc = d["anchors"]
    if len(anc):
        labels = list(anc["label"]) if "label" in anc.columns else None
        fig.add_trace(go.Scatter(
            x=pd.to_numeric(anc["rt_src"], errors="coerce"),
            y=pd.to_numeric(anc["rt_ref"], errors="coerce"),
            mode="markers", name=anchor_legend_label(d), text=labels,
            marker=dict(size=13, color="rgba(0,0,0,0)",
                        line=dict(color=C_RING, width=1.8)),
            hovertemplate="%{text}<br>%{x:.2f} → %{y:.2f} min<extra></extra>"), 1, 1)

    fig.add_trace(go.Scatter(x=d["rt_src"], y=d["r_curve"], mode="markers",
                             name="residual · anchor-free",
                             marker=dict(size=5, color=theme.rgba(theme.NEUTRAL, 0.55)),
                             hovertemplate="RT %{x:.2f} min<br>residual %{y:+.3f} min<extra></extra>"),
                  2, 1)
    if d["r_refined"] is not None:
        fig.add_trace(go.Scatter(x=d["rt_src"], y=d["r_refined"], mode="markers",
                                 name="residual · anchor-refined",
                                 marker=dict(size=5, color=theme.rgba(C_REFINED, 0.65)),
                                 hovertemplate="RT %{x:.2f} min<br>residual %{y:+.3f} min<extra></extra>"),
                      2, 1)
    fig.add_hline(y=0, line=dict(color=theme.AXIS, width=0.8), row=2, col=1)

    fig.update_layout(theme.plotly_template())
    fig.update_xaxes(title_text="RT on your column (min)", row=2, col=1)
    fig.update_yaxes(title_text="RT on the reference column (min)", row=1, col=1)
    fig.update_yaxes(title_text="residual (min)", row=2, col=1,
                     range=list(_resid_limits(d["r_curve"])))
    fig.update_layout(height=700, title=None, hovermode="closest",
                      margin=dict(l=76, r=28, t=30, b=56),
                      legend=dict(orientation="v", x=0.01, y=0.99,
                                  bgcolor=theme.rgba(theme.PANEL, 0.75),
                                  font=dict(size=12, color=theme.TXT2)))
    return fig
