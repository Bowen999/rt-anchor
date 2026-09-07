"""Module 4 — inter-injection RI repeatability (per-sample tier only).

The headline QC of per-sample self-anchoring: do the different injections agree on
where each feature lands? `RI_spread` = MAD of a feature's RI across contributing
injections. (Per VIM/IUPAC, within-instrument within-method agreement is
*repeatability*; "reproducibility" is reserved for cross-lab/cross-platform conditions.)
Two panels:
  (a) distribution of RI_spread across features (most should be small);
  (b) RI_spread vs iRT — where injections disagree (usually the sparse tail).

Only applicable when `calibration_scope == 'sample'`; NaN under per-project.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd

from . import theme


def is_applicable(result) -> bool:
    if result.model.get("calibration_scope") != "sample":
        return False
    s = result.values("RI_spread")
    return bool(s.notna().any())


def compute_repro(result) -> Dict:
    tbl = result.table
    spread = result.values("RI_spread")
    ri = result.values("iRT")
    nc = pd.to_numeric(tbl[result.col("n_contributing")], errors="coerce")
    ok = spread.notna() & ri.notna() & (nc >= 2)
    s, r = spread[ok].to_numpy(), ri[ok].to_numpy()
    return dict(spread=s, ri=r, n=int(ok.sum()),
                median=(float(np.median(s)) if s.size else np.nan),
                p90=(float(np.percentile(s, 90)) if s.size else np.nan),
                conf_high=result.model.get("config", {}).get("conf_high_irt", 0.6))


def _binned(ri, spread, lo, hi, nbins=40):
    edges = np.linspace(lo, hi, nbins + 1)
    idx = np.clip(np.digitize(ri, edges) - 1, 0, nbins - 1)
    cx, med, q1, q3 = [], [], [], []
    for b in range(nbins):
        v = spread[idx == b]
        if v.size >= 3:
            cx.append(0.5 * (edges[b] + edges[b + 1]))
            med.append(np.median(v)); q1.append(np.percentile(v, 25)); q3.append(np.percentile(v, 75))
    return map(np.array, (cx, med, q1, q3))


# ---------------------------------------------------------------- matplotlib --

def figure_mpl(result):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    d = compute_repro(result)
    fig = plt.figure(figsize=(12.5, 5.4))
    gs = GridSpec(1, 2, width_ratios=[1.0, 1.5], wspace=0.24,
                  left=0.08, right=0.97, top=0.86, bottom=0.15)

    # (a) distribution
    ax = fig.add_subplot(gs[0])
    s = d["spread"]
    hi = np.percentile(s, 99) if s.size else 1.0
    ax.hist(np.clip(s, 0, hi), bins=40, color=theme.PRIMARY, alpha=0.85)
    if np.isfinite(d["median"]):
        ax.axvline(d["median"], color=theme.ACCENT, lw=1.4, ls=(0, (4, 3)))
        ax.text(d["median"], ax.get_ylim()[1] * 0.58, f"   median {d['median']:.2f}",
                color=theme.ACCENT_DK, fontsize=theme.FS_SUB - 1, va="center", ha="left")
    theme.style_ax(ax)
    ax.set_xlabel("RI spread across injections (iRT)"); ax.set_ylabel("features")
    ax.set_title("Repeatability distribution", loc="left",
                 fontsize=theme.FS_SUB + 2, color=theme.TXT, fontweight="bold")

    # (b) spread vs iRT
    axb = fig.add_subplot(gs[1])
    if s.size:
        cx, med, q1, q3 = _binned(d["ri"], s, float(np.nanmin(d["ri"])), float(np.nanmax(d["ri"])))
        if len(cx):
            axb.fill_between(cx, q1, q3, color=theme.PRIMARY, alpha=0.18)
            axb.plot(cx, med, color=theme.PRIMARY, lw=1.8)
        axb.axhline(d["conf_high"], color=theme.NEUTRAL, lw=0.8, ls=(0, (4, 3)))
    theme.style_ax(axb)
    axb.set_xlabel("iRT (dimensionless)"); axb.set_ylabel("RI spread (iRT)")
    axb.set_ylim(0, None)
    axb.set_title("Repeatability across the gradient  ·  median ± IQR", loc="left",
                  fontsize=theme.FS_SUB + 2, color=theme.TXT, fontweight="bold")

    fig.suptitle(f"Inter-injection RI repeatability  ·  {d['n']:,} features, "
                 f"{result.model.get('n_injections_used', '?')} injections",
                 x=0.02, y=0.97, ha="left", fontsize=theme.FS_TITLE, color=theme.TXT)
    return fig


# ---------------------------------------------------------------- plotly ------

def figure_plotly(result):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    d = compute_repro(result)
    s = d["spread"]
    fig = make_subplots(rows=1, cols=2, column_widths=[0.4, 0.6], horizontal_spacing=0.1,
                        subplot_titles=("Repeatability distribution",
                                        "Repeatability across the gradient (median ± IQR)"))
    hi = np.percentile(s, 99) if s.size else 1.0
    fig.add_trace(go.Histogram(x=np.clip(s, 0, hi), nbinsx=40, marker_color=theme.PRIMARY,
                  showlegend=False, hovertemplate="spread %{x:.2f}<br>%{y} features<extra></extra>"), 1, 1)
    if np.isfinite(d["median"]):
        fig.add_vline(x=d["median"], line=dict(color=theme.ACCENT, width=1.4, dash="dash"), row=1, col=1)
        fig.add_annotation(x=d["median"], xref="x", y=0.58, yref="y domain", xshift=8,
                           text=f"median {d['median']:.2f}", showarrow=False, xanchor="left",
                           font=dict(size=12, color=theme.ACCENT_DK), row=1, col=1)
    if s.size:
        cx, med, q1, q3 = _binned(d["ri"], s, float(np.nanmin(d["ri"])), float(np.nanmax(d["ri"])))
        if len(cx):
            fig.add_trace(go.Scatter(x=np.r_[cx, cx[::-1]], y=np.r_[q3, q1[::-1]], fill="toself",
                          fillcolor="rgba(43,93,125,0.16)", line=dict(width=0), showlegend=False,
                          hoverinfo="skip"), 1, 2)
            fig.add_trace(go.Scatter(x=cx, y=med, mode="lines", line=dict(color=theme.PRIMARY, width=2),
                          showlegend=False, hovertemplate="iRT %{x:.0f}<br>spread %{y:.2f}<extra></extra>"), 1, 2)
        fig.add_hline(y=d["conf_high"], line=dict(color=theme.NEUTRAL, width=1, dash="dash"), row=1, col=2)
    fig.update_xaxes(title_text="RI spread (iRT)", row=1, col=1)
    fig.update_yaxes(title_text="features", row=1, col=1)
    fig.update_xaxes(title_text="iRT (dimensionless)", row=1, col=2)
    fig.update_yaxes(title_text="RI spread (iRT)", row=1, col=2, rangemode="tozero")
    fig.update_layout(theme.plotly_template())
    fig.update_layout(title=None, height=460, showlegend=False)
    return fig
