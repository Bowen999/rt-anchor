"""Module 4 — the stage-2 sample anchors and the gate's verdict.

Stage 2 is the part of the method a reader is most entitled to be sceptical
about: it moves calibrated RT using a handful of endogenous lipids identified
from MS1 alone. This page is the evidence for (or against) that move.

One row per anchor, ordered by where it elutes on the user's column:

* the **filled dot** is the anchor's residual against the stage-1 curve alone —
  how far the anchor-free calibration missed that lipid;
* the **ring** is its *leave-one-out* residual under the correction — what the
  correction achieves for an anchor it was not allowed to see. Comparing the two
  is the only honest way to judge a correction fitted on the same points it is
  being scored on;
* colour is lipid class, because the whole reason stage 2 is class-aware is that
  the residuals separate by class on long gradients (SM/CE high, PC/TG low). If
  the classes visibly stack on one side of zero, the class-offset term is doing
  real work.

Anchors the sanity filter dropped are drawn hollow and greyed: they were found
and then rejected, which is different from never having been found.

When the gate is **off** the rows are still drawn — they are the evidence for
the decision — under a banner giving the measured LOO gain and the threshold it
missed. A report that showed only the engaged case would be a report that never
explains itself when it matters most.

No baked-in title — the HTML report supplies a copyable one.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd

from . import theme

C_PRE = theme.PRIMARY        # residual against the stage-1 curve
C_POST = theme.ACCENT_DK     # leave-one-out residual under the correction
C_DROP = theme.REJECT        # anchors removed by the sanity filter


def is_applicable(result) -> bool:
    """True when the run has stage-2 anchors worth a page (used or rejected)."""
    a = getattr(result, "anchors", None)
    return a is not None and len(a) > 0


def compute_anchors(result) -> Dict:
    """Everything both renderers need, computed once.

    ``empty`` is True when stage 2 never produced an anchor — a non-plasma
    matrix, or a panel-free run. The renderers then say so instead of drawing
    an empty axis.
    """
    a = getattr(result, "anchors", None)
    if a is None or not len(a):
        return dict(empty=True,
                    reason="no plasma-lipid anchors validated in this run "
                           "(stage 1 curve only)")

    a = a.copy()
    if "dropped_by_sanity_filter" not in a.columns:
        a["dropped_by_sanity_filter"] = False
    a["dropped_by_sanity_filter"] = a["dropped_by_sanity_filter"].astype(bool)
    for col in ("residual_min", "loo_residual_min", "rt_src", "rt_ref"):
        if col not in a.columns:
            a[col] = np.nan
        a[col] = pd.to_numeric(a[col], errors="coerce")
    if "lipid_class" not in a.columns:
        a["lipid_class"] = ""
    a["lipid_class"] = a["lipid_class"].astype(str)
    a = a.sort_values("rt_src").reset_index(drop=True)

    cal = getattr(result, "calibrator", None)
    engaged = bool(getattr(result, "anchors_used", False))
    kept = a[~a["dropped_by_sanity_filter"]]

    # class offsets: the quantity the lam_c term shrinks toward zero
    offsets = (kept.groupby("lipid_class")["residual_min"].median().to_dict()
               if len(kept) else {})

    def _rms(v):
        v = np.asarray(v, dtype=float)
        v = v[np.isfinite(v)]
        return float(np.sqrt(np.mean(v ** 2))) if v.size else np.nan

    return dict(
        empty=False,
        anchors=a,
        n_used=int(len(kept)),
        n_dropped=int(a["dropped_by_sanity_filter"].sum()),
        engaged=engaged,
        class_offsets=offsets,
        rms_pre=_rms(kept["residual_min"]) if len(kept) else np.nan,
        rms_post=_rms(kept["loo_residual_min"]) if len(kept) else np.nan,
        lam_g=float(getattr(cal, "lam_g", np.nan)) if cal is not None else np.nan,
        lam_c=float(getattr(cal, "lam_c", np.nan)) if cal is not None else np.nan,
        class_aware=bool(getattr(cal, "class_aware", False)) if cal is not None else False,
        gate_reduction=float(getattr(cal, "gate_mse_reduction", np.nan)) if cal is not None else np.nan,
        gate_threshold=float(getattr(cal, "gate_threshold", np.nan)) if cal is not None else np.nan,
        gate_reason=str(getattr(cal, "gate_reason", "") or "") if cal is not None else "",
    )


def verdict_line(d: Dict) -> str:
    """One sentence stating what stage 2 did and what it bought."""
    if d.get("empty"):
        return str(d.get("reason", ""))
    if d["engaged"]:
        gain = ""
        if np.isfinite(d["rms_pre"]) and np.isfinite(d["rms_post"]):
            gain = (f" — anchor RMS {d['rms_pre']:.3f} → {d['rms_post']:.3f} min "
                    f"leave-one-out")
        lam = (f"λg={d['lam_g']:.1f}, λc={d['lam_c']:.1f}"
               if np.isfinite(d["lam_g"]) else "")
        return f"correction APPLIED ({lam}){gain}"
    pct = (f"{100 * d['gate_reduction']:.0f}%"
           if np.isfinite(d["gate_reduction"]) else "n/a")
    thr = (f"{100 * d['gate_threshold']:.0f}%"
           if np.isfinite(d["gate_threshold"]) else "n/a")
    return (f"correction NOT applied — leave-one-out gain {pct} did not reach "
            f"the {thr} threshold; these anchors are shown as the evidence for "
            f"that decision")


def _empty_fig_mpl(reason: str):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12.5, 3.2))
    ax.axis("off")
    ax.text(0.5, 0.5, f"Stage 2 not applicable\n{reason}", ha="center", va="center",
            fontsize=theme.FS_AXIS, color=theme.TXT2)
    return fig


# ---------------------------------------------------------------- matplotlib --

def figure_mpl(result):
    import matplotlib.pyplot as plt
    d = compute_anchors(result)
    if d["empty"]:
        return _empty_fig_mpl(d["reason"])

    a = d["anchors"]
    n = len(a)
    fig, ax = plt.subplots(figsize=(12.5, max(3.4, 0.34 * n + 2.0)))
    y = np.arange(n)[::-1]                     # first-eluting anchor at the top

    ax.axvline(0, color=theme.AXIS, lw=1.0, zorder=1)
    for yi, (_, r) in zip(y, a.iterrows()):
        dropped = bool(r["dropped_by_sanity_filter"])
        col = C_DROP if dropped else theme.class_color(r["lipid_class"])
        pre, post = r["residual_min"], r["loo_residual_min"]
        if np.isfinite(pre) and np.isfinite(post) and not dropped:
            ax.plot([pre, post], [yi, yi], color=col, lw=1.0, alpha=0.5, zorder=2)
        if np.isfinite(pre):
            ax.scatter([pre], [yi], s=54, color=col, zorder=4,
                       edgecolor="white", lw=0.6,
                       alpha=0.45 if dropped else 1.0)
        if np.isfinite(post) and not dropped:
            ax.scatter([post], [yi], s=54, facecolor="none", edgecolor=C_POST,
                       lw=1.4, zorder=5)

    labels = [f"{r['label']}" + ("  (dropped)" if r["dropped_by_sanity_filter"] else "")
              for _, r in a.iterrows()]
    ax.set_yticks(y[::-1])
    ax.set_yticklabels(labels[::-1], fontsize=theme.FS_TICK)
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_xlabel("residual against the stage-1 curve (min)   ·   "
                  "filled = anchor-free, ring = leave-one-out under the correction",
                  fontsize=theme.FS_AXIS, color=theme.TXT2)
    theme.style_ax(ax)

    off = "   ".join(f"{k} {v:+.3f}" for k, v in sorted(d["class_offsets"].items()) if k)
    sub = verdict_line(d)
    if off:
        sub += f"\nclass median offsets (min):  {off}"
    if d["n_dropped"]:
        sub += f"\n{d['n_dropped']} anchor(s) removed by the sanity filter"
    fig.suptitle("Stage-2 sample anchors", x=0.02, y=0.985, ha="left",
                 fontsize=theme.FS_TITLE, color=theme.TXT)
    fig.text(0.02, 0.945, theme.mpl_safe(sub), ha="left", va="top", fontsize=theme.FS_SUB,
             color=theme.TXT2, style="italic")
    fig.subplots_adjust(left=0.16, right=0.98, top=0.86 - min(0.06, 0.012 * sub.count("\n")),
                        bottom=0.16)
    return fig


# ---------------------------------------------------------------- plotly ------

def figure_plotly(result):
    import plotly.graph_objects as go
    d = compute_anchors(result)
    fig = go.Figure()
    fig.update_layout(theme.plotly_template())
    if d["empty"]:
        fig.add_annotation(text=f"Stage 2 not applicable — {d['reason']}",
                           showarrow=False, x=0.5, y=0.5, xref="paper", yref="paper",
                           font=dict(size=14, color=theme.TXT2))
        fig.update_layout(height=220, xaxis=dict(visible=False), yaxis=dict(visible=False))
        return fig

    a = d["anchors"]
    labels = [f"{r['label']}" + ("  (dropped)" if r["dropped_by_sanity_filter"] else "")
              for _, r in a.iterrows()]

    for _, (lab, r) in enumerate(zip(labels, [r for _, r in a.iterrows()])):
        dropped = bool(r["dropped_by_sanity_filter"])
        col = C_DROP if dropped else theme.class_color(r["lipid_class"])
        pre, post = r["residual_min"], r["loo_residual_min"]
        if np.isfinite(pre) and np.isfinite(post) and not dropped:
            fig.add_trace(go.Scatter(x=[pre, post], y=[lab, lab], mode="lines",
                                     line=dict(color=col, width=1.4),
                                     opacity=0.5, showlegend=False, hoverinfo="skip"))
        if np.isfinite(pre):
            fig.add_trace(go.Scatter(
                x=[pre], y=[lab], mode="markers", showlegend=False,
                marker=dict(size=11, color=col, line=dict(color="white", width=1),
                            opacity=0.45 if dropped else 1.0),
                hovertemplate=(f"{r['label']} ({r['lipid_class']})<br>"
                               f"RT {r['rt_src']:.2f} → ref {r['rt_ref']:.2f} min<br>"
                               f"anchor-free residual %{{x:.3f}} min<extra></extra>")))
        if np.isfinite(post) and not dropped:
            fig.add_trace(go.Scatter(
                x=[post], y=[lab], mode="markers", showlegend=False,
                marker=dict(size=11, color="rgba(0,0,0,0)",
                            line=dict(color=C_POST, width=2)),
                hovertemplate=(f"{r['label']}<br>leave-one-out residual "
                               f"%{{x:.3f}} min<extra></extra>")))

    fig.add_vline(x=0, line=dict(color=theme.AXIS, width=1))
    fig.update_layout(
        height=max(260, 30 * len(a) + 130),
        xaxis=dict(title=dict(text="residual against the stage-1 curve (min) · "
                                   "filled = anchor-free, ring = leave-one-out")),
        yaxis=dict(autorange="reversed", categoryorder="array",
                   categoryarray=labels, title=None),
        margin=dict(l=140, r=30, t=60, b=60))
    fig.add_annotation(text=verdict_line(d), showarrow=False, xref="paper", yref="paper",
                       x=0, y=1.06, xanchor="left", align="left",
                       font=dict(size=12, color=theme.TXT2))
    return fig
