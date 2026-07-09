"""Module 3 — the fitted monotone calibration warp RT → iRT.

Just the warp curve with the anchor points and an anchor rug. (The heuristic
uncertainty landscape was removed from the report; sigma is still computed and
written to the RI_uncertainty column.) No baked-in title — the HTML report
supplies a copyable title.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd

from . import theme


def compute_warp(result, ngrid: int = 400) -> Dict:
    w = result.warp
    if w is None or w.rt.size < 2:
        return dict(empty=True)
    grid = np.linspace(w.rt_min, w.rt_max, ngrid)
    return dict(empty=False, grid=grid, curve=w.predict(grid, extrapolate=False),
                anchor_rt=w.rt, anchor_irt=w.irt,
                anchors=result.anchors.drop_duplicates("name") if len(result.anchors) else result.anchors)


# ---------------------------------------------------------------- matplotlib --

def figure_mpl(result):
    import matplotlib.pyplot as plt
    d = compute_warp(result)
    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    if d.get("empty"):
        fig.text(0.5, 0.5, "no warp available", ha="center", color=theme.TXT2)
        return fig
    ax.plot(d["grid"], d["curve"], color=theme.PRIMARY, lw=2.2, zorder=3)
    ax.scatter(d["anchor_rt"], d["anchor_irt"], s=52, color=theme.PRIMARY_DK,
               edgecolor="white", lw=1.0, zorder=4)
    anc = d["anchors"]
    if len(anc) and "name" in anc:
        for _, r in anc.iterrows():
            ax.annotate(r["name"].replace("_", " "), (r["rt_obs_min"], r["irt"]),
                        xytext=(5, -2), textcoords="offset points",
                        fontsize=theme.FS_SUB - 3, color=theme.TXT2)
    lo = float(np.min(d["anchor_irt"]))
    for x in d["anchor_rt"]:
        ax.plot([x, x], [lo - 3, lo - 6], color=theme.FAINT, lw=1.0, clip_on=False)
    theme.style_ax(ax)
    ax.set_xlabel("Retention time (min)"); ax.set_ylabel("iRT (dimensionless)")
    fig.subplots_adjust(left=0.09, right=0.97, top=0.95, bottom=0.11)
    return fig


# ---------------------------------------------------------------- plotly ------

def figure_plotly(result):
    import plotly.graph_objects as go
    d = compute_warp(result)
    fig = go.Figure()
    if d.get("empty"):
        fig.update_layout(theme.plotly_template()); return fig
    fig.add_trace(go.Scatter(x=d["grid"], y=d["curve"], mode="lines",
                  line=dict(color=theme.PRIMARY, width=2.4), name="warp",
                  hovertemplate="RT %{x:.2f} min → iRT %{y:.1f}<extra></extra>"))
    anc = d["anchors"]
    txt = list(anc["name"]) if len(anc) and "name" in anc else None
    fig.add_trace(go.Scatter(x=d["anchor_rt"], y=d["anchor_irt"], mode="markers",
                  marker=dict(size=11, color=theme.PRIMARY_DK, line=dict(color="white", width=1.2)),
                  text=txt, customdata=txt, name="anchor",
                  hovertemplate="%{text}<br>RT %{x:.2f} → iRT %{y:.1f}<extra></extra>"))
    fig.update_layout(theme.plotly_template())
    fig.update_layout(height=520, showlegend=False, title=None,
                      xaxis=dict(title=dict(text="Retention time (min)")),
                      yaxis=dict(title=dict(text="iRT (dimensionless)")))
    return fig
