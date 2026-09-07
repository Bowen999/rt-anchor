"""Module 2 — reconstructed feature-intensity profile (mirror tanglegram).

Feature tables carry no scan-level signal, so this is NOT a true total-ion
chromatogram: a chromatogram-like profile ("pseudo-chromatogram") is *reconstructed*
by binning feature apex intensities (summed over samples) and lightly smoothing.

The figure is a **mirror**: the before trace (raw RT) points up, the after trace
(iRT) points down. Each axis is normalised to its own range; the two axes are
labelled independently (RT on top, iRT on the bottom). Identified species are
**highlighted and connected by dashed lines** whose slant shows the RT offset the
calibration removes (near-flat early, fanning in the compressed tail). Those
lines come from the stage-2 anchors when the run has them, and otherwise from a
sample of the stage-1 matched pairs — see :func:`tie_points`.

Styles: ``clean`` (honest binned intensity) / ``realistic`` (adds a simulated
baseline + noise). ``clean`` is the report default.
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from . import theme

NBINS = 600
C_BEFORE = theme.ACCENT     # clay
C_AFTER = theme.PRIMARY     # blue


def _reconstruct(x, inten, lo, hi, style="clean", nbins=NBINS, seed=0):
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.array([lo, hi]), np.array([0.0, 0.0])
    edges = np.linspace(lo, hi, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    ok = np.isfinite(x) & np.isfinite(inten)
    peaks = gaussian_filter1d(np.histogram(x[ok], bins=edges, weights=inten[ok])[0].astype(float), 1.5)
    if style != "realistic":
        return centers, peaks
    vmax = peaks.max() or 1.0
    rng = np.random.default_rng(seed)
    broad = gaussian_filter1d(peaks, 40); broad = broad / (broad.max() or 1.0)
    sig = peaks + vmax * (0.020 + 0.060 * broad)
    noise = rng.normal(0, 1, sig.size) * (0.012 * vmax + 0.030 * np.sqrt(np.clip(sig, 0, None) * vmax))
    return centers, np.clip(sig + noise, 0.0, None)


def _intensity(result) -> np.ndarray:
    df = result.table
    if result.sample_cols:
        return df[result.sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1).to_numpy()
    low = {c.lower(): c for c in df.columns}
    for key in ("peak_height", "peak_area", "maxo", "height", "area", "intensity"):
        if key in low:
            return pd.to_numeric(df[low[key]], errors="coerce").fillna(0).to_numpy()
    return np.ones(len(df))


MAX_PAIR_TIES = 15          # tie lines drawn when there are no stage-2 anchors


def tie_points(result) -> pd.DataFrame:
    """Species to connect across the mirror: ``name``, ``rt``, ``irt``, ``kind``.

    Under v2 there is no single set of anchors living on both axes, so the tie
    lines are drawn from whatever the run actually identified:

    * **stage-2 anchors** when they exist — endogenous lipids located in the
      sample run itself (``rt_src``) whose reference-column RT maps straight
      onto the iRT ruler. These are named species, which is the point of the
      tanglegram: *this* lipid moved *this* far.
    * otherwise an evenly-spaced sample of the **stage-1 matched pairs**, drawn
      unnamed. A gated-off or non-plasma run still deserves to show the offset
      the calibration removed, even when nothing in it has a name.

    The iRT ruler may be absent (``result.irt is None``); then no tie line can
    be placed and an empty frame comes back.
    """
    cols = ["name", "rt", "irt", "kind"]
    if result.irt is None:
        return pd.DataFrame(columns=cols)

    anc = result.anchors
    if len(anc) and "dropped_by_sanity_filter" in anc.columns:
        anc = anc[~anc["dropped_by_sanity_filter"].astype(bool)]
    if len(anc):
        rt = pd.to_numeric(anc["rt_src"], errors="coerce")
        ref = pd.to_numeric(anc["rt_ref"], errors="coerce")
        ok = rt.notna() & ref.notna()
        if ok.any():
            return pd.DataFrame({
                "name": anc.loc[ok, "label"].astype(str).to_numpy(),
                "rt": rt[ok].to_numpy(),
                "irt": result.irt.to_irt(ref[ok].to_numpy()),
                "kind": "anchor",
            })

    pairs = result.pairs
    if pairs is None or not len(pairs):
        return pd.DataFrame(columns=cols)
    if "kept" in pairs.columns:
        pairs = pairs[pairs["kept"].astype(bool)]
    pairs = pairs.sort_values("rt_a")
    if not len(pairs):
        return pd.DataFrame(columns=cols)
    idx = np.unique(np.linspace(0, len(pairs) - 1, min(MAX_PAIR_TIES, len(pairs)))
                    .round().astype(int))
    sub = pairs.iloc[idx]
    return pd.DataFrame({
        "name": [""] * len(sub),
        "rt": sub["rt_a"].to_numpy(dtype=float),
        "irt": result.irt.to_irt(sub["rt_b"].to_numpy(dtype=float)),
        "kind": "pair",
    })


def compute_tic(result, style: str = "clean") -> Dict:
    rt = result.rt_minutes().to_numpy()
    ri = result.values("iRT").to_numpy()
    inten = _intensity(result)
    ties = tie_points(result)
    rt_lo, rt_hi = float(np.nanmin(rt)), float(np.nanmax(rt))
    fr = ri[np.isfinite(ri)]
    ir_lo, ir_hi = (float(np.nanmin(fr)), float(np.nanmax(fr))) if fr.size else (0.0, 100.0)
    # extend both axes so a tie marker can never fall outside the plotted range
    if len(ties):
        t_rt = pd.to_numeric(ties["rt"], errors="coerce")
        t_ir = pd.to_numeric(ties["irt"], errors="coerce")
        if t_rt.notna().any():
            rt_lo, rt_hi = min(rt_lo, float(t_rt.min())), max(rt_hi, float(t_rt.max()))
        if t_ir.notna().any():
            ir_lo, ir_hi = min(ir_lo, float(t_ir.min())), max(ir_hi, float(t_ir.max()))
    cb, yb = _reconstruct(rt, inten, rt_lo, rt_hi, style, seed=1)
    ca, ya = _reconstruct(ri, inten, ir_lo, ir_hi, style, seed=2)
    return dict(cb=cb, yb=yb / (yb.max() or 1), ca=ca, ya=ya / (ya.max() or 1),
                rt_range=(rt_lo, rt_hi), ir_range=(ir_lo, ir_hi),
                ties=ties, style=style)


def _norm(v, lo, hi):
    return (np.asarray(v, float) - lo) / (hi - lo) if hi > lo else np.zeros_like(np.asarray(v, float))


def _ticks(lo, hi, step):
    import math
    start = math.ceil(lo / step) * step
    return [v for v in np.arange(start, hi + 1e-9, step)]


# ---------------------------------------------------------------- matplotlib --

def figure_mpl(result, style: str = "clean"):
    import matplotlib.pyplot as plt
    d = compute_tic(result, style)
    rt_lo, rt_hi = d["rt_range"]; ir_lo, ir_hi = d["ir_range"]
    nb = _norm(d["cb"], rt_lo, rt_hi); na = _norm(d["ca"], ir_lo, ir_hi)

    fig, ax = plt.subplots(figsize=(12.5, 6.6))
    ax.fill_between(nb, 0, d["yb"], color=C_BEFORE, alpha=0.85, lw=0)
    ax.plot(nb, d["yb"], color=theme.darker(C_BEFORE, 0.6), lw=1.0)
    ax.fill_between(na, 0, -d["ya"], color=C_AFTER, alpha=0.85, lw=0)
    ax.plot(na, -d["ya"], color=theme.darker(C_AFTER, 0.6), lw=1.0)
    ax.axhline(0, color=theme.AXIS, lw=0.8)

    ties = d["ties"]
    if len(ties):
        for _, r in ties.iterrows():
            xb = float(_norm(r["rt"], rt_lo, rt_hi))
            xa = float(_norm(r["irt"], ir_lo, ir_hi))
            ax.plot([xb, xa], [0.07, -0.07], color=theme.NEUTRAL, lw=0.8, ls=(0, (4, 3)), zorder=2)
            ax.scatter([xb], [0.07], s=26, color=C_BEFORE, zorder=4, edgecolor="white", lw=0.6)
            ax.scatter([xa], [-0.07], s=26, color=C_AFTER, zorder=4, edgecolor="white", lw=0.6)

    ax.text(0.004, 1.06, "before  ·  raw RT", color=theme.darker(C_BEFORE, 0.6),
            fontsize=theme.FS_SUB + 1, fontweight="bold")
    ax.text(0.004, -1.16, "after  ·  iRT", color=theme.darker(C_AFTER, 0.6),
            fontsize=theme.FS_SUB + 1, fontweight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(-1.2, 1.2); ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    # bottom axis = iRT
    it = _ticks(ir_lo, ir_hi, 20)
    ax.set_xticks(_norm(it, ir_lo, ir_hi)); ax.set_xticklabels([f"{v:g}" for v in it])
    ax.set_xlabel("iRT (dimensionless)", color=theme.darker(C_AFTER, 0.6))
    ax.spines["bottom"].set_color(theme.AXIS)
    # top axis = RT (min)
    axt = ax.twiny(); axt.set_xlim(0, 1)
    rtt = _ticks(rt_lo, rt_hi, 5)
    axt.set_xticks(_norm(rtt, rt_lo, rt_hi)); axt.set_xticklabels([f"{v:g}" for v in rtt])
    axt.set_xlabel("Retention time (min)", color=theme.darker(C_BEFORE, 0.6), labelpad=8)
    for s in ("right", "left", "bottom"):
        axt.spines[s].set_visible(False)
    axt.spines["top"].set_color(theme.AXIS); axt.tick_params(colors=theme.AXIS, labelcolor=theme.TXT2)

    note = ("reconstructed from the feature table (binned intensity)" if style == "clean"
            else "reconstructed + SIMULATED baseline/noise")
    fig.suptitle("Reconstructed feature-intensity profile  ·  before vs after calibration", x=0.02, y=0.985,
                 ha="left", fontsize=theme.FS_TITLE, color=theme.TXT)
    kind = ("stage-2 anchors" if (len(d["ties"]) and d["ties"]["kind"].iloc[0] == "anchor")
            else "sampled matched pairs")
    fig.text(0.02, 0.93, f"{kind} connected — line slant = RT offset · {note}",
             ha="left", fontsize=theme.FS_SUB, color=theme.TXT2, style="italic")
    fig.subplots_adjust(left=0.04, right=0.98, top=0.82, bottom=0.11)
    return fig


# ---------------------------------------------------------------- plotly ------

def figure_plotly(result, style: str = "clean"):
    import plotly.graph_objects as go
    d = compute_tic(result, style)
    rt_lo, rt_hi = d["rt_range"]; ir_lo, ir_hi = d["ir_range"]
    nb = _norm(d["cb"], rt_lo, rt_hi); na = _norm(d["ca"], ir_lo, ir_hi)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.r_[nb, nb[::-1]], y=np.r_[d["yb"], np.zeros_like(d["yb"])],
                  fill="toself", fillcolor=_rgba(C_BEFORE, 0.85), line=dict(width=0),
                  name="before · raw RT", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=np.r_[na, na[::-1]], y=np.r_[-d["ya"], np.zeros_like(d["ya"])],
                  fill="toself", fillcolor=_rgba(C_AFTER, 0.85), line=dict(width=0),
                  name="after · iRT", hoverinfo="skip"))
    fig.add_hline(y=0, line=dict(color=theme.AXIS, width=0.8))

    ties = d["ties"]
    if len(ties):
        for _, r in ties.iterrows():
            xb = float(_norm(r["rt"], rt_lo, rt_hi)); xa = float(_norm(r["irt"], ir_lo, ir_hi))
            label = str(r["name"]) or "matched feature"
            fig.add_trace(go.Scatter(x=[xb, xa], y=[0.07, -0.07], mode="lines",
                          line=dict(color=theme.NEUTRAL, width=0.8, dash="dash"),
                          showlegend=False, hoverinfo="skip"))
            fig.add_trace(go.Scatter(x=[xb], y=[0.07], mode="markers", showlegend=False,
                          marker=dict(size=8, color=C_BEFORE, line=dict(color="white", width=0.6)),
                          hovertemplate=f"{label}<br>RT %{{customdata:.2f}} min<extra></extra>",
                          customdata=[r["rt"]]))
            fig.add_trace(go.Scatter(x=[xa], y=[-0.07], mode="markers", showlegend=False,
                          marker=dict(size=8, color=C_AFTER, line=dict(color="white", width=0.6)),
                          hovertemplate=f"{label}<br>iRT %{{customdata:.1f}}<extra></extra>",
                          customdata=[r["irt"]]))

    it = _ticks(ir_lo, ir_hi, 20); rtt = _ticks(rt_lo, rt_hi, 5)
    fig.update_layout(theme.plotly_template())
    fig.update_layout(
        title=None,
        height=560, showlegend=True,
        xaxis=dict(range=[0, 1], showgrid=False, tickvals=list(_norm(it, ir_lo, ir_hi)),
                   ticktext=[f"{v:g}" for v in it], title=dict(text="iRT (dimensionless)")),
        xaxis2=dict(range=[0, 1], overlaying="x", side="top", showgrid=False,
                    tickvals=list(_norm(rtt, rt_lo, rt_hi)), ticktext=[f"{v:g}" for v in rtt],
                    title=dict(text="Retention time (min)")),
        yaxis=dict(range=[-1.2, 1.2], showticklabels=False, showgrid=False, zeroline=False))
    # anchor the top axis to a trace
    fig.add_trace(go.Scatter(x=[0], y=[0], xaxis="x2", mode="markers",
                  marker=dict(opacity=0), showlegend=False, hoverinfo="skip"))
    return fig


def _rgba(hex_color: str, a: float) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{a})"
