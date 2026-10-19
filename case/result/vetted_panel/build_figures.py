"""Figures for the vetted-lipid validation panel (RT Anchor 1.2.3 on the reprocessed data).

Reads ONLY the CSVs written by build_vetted_panel.py (in this folder) and writes
figures/<set>/ for the two sets plus one QC figure:

  figures/high_confidence/            high-confidence panel lipids
  figures/high_and_low_confidence/    high-confidence + low-confidence panel lipids
      overview_heatmap.png            median dRT per display group and column, before | after calibration
      rt_error_ecdf.png               cumulative distribution of |dRT| (pooled cells)
      rt_error_violin_linear.png      signed dRT per column, violin pairs, linear minute axis
      rt_error_violin_symlog.png      the same on the symmetric-log axis of the heatmap colour scale
  figures/vetted_rt_check_column_20.png   QC of the vetted-RT -> column_20-RT mapping (set independent)

dRT = RT - reference RT (column_22) before calibration and calibrated RT - reference RT after. All heatmaps
share one symmetric-log colour scale (NORM / CMAP, the one of result_v1.2.3/build_comparison.py); the ECDF
and the violins follow result_v1.2.3/build_error_distribution.py. For the pooled figures (ECDF, violins) each
cell is counted once: key = (column, pick Alignment ID, reference pick Alignment ID), high lipids first. The
overview uses each group's own lipids. The per-class heatmaps are circular and drawn by
build_ring_heatmaps.py (the rectangular per-class heatmaps were dropped on 2026-10-14 at the owner's request).

Besides the two panel CSVs it reads vetted_rt_mapping.csv and vetted_rt_check_groups.csv (the per-group
points of the QC figure: sheet, vetted RT of the matched row, column_20 pick RT, status), which
build_vetted_panel.py writes for exactly this purpose.

Run:  python build_figures.py        (about 10 s; needs build_vetted_panel.py to have run)
"""
import os
import sys
import time
from decimal import ROUND_HALF_UP, Decimal

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.scale import SymmetricalLogTransform
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))           # .../vetted_panel
FIG_DIR = os.path.join(HERE, "figures")

NONREF = ["column_12", "column_16", "column_20", "column_41", "column_85"]
SETS = {"high_confidence": "high_confidence_panel.csv",
        "high_and_low_confidence": "high_and_low_confidence_panel.csv"}
GROUP_ORDER = ["CAR", "LPC", "EtherLPC", "LPE", "EtherLPE", "PC", "EtherPC", "PE", "EtherPE", "PI",
               "SM", "Cer_NS", "HexCer_NS", "DG", "TG", "CE", "Other", "Low confidence"]
DPI = 200
WITHIN = 0.2                                                   # min, the headline accuracy threshold

# ----------------------------------------------------------------------------- style (as in the earlier scripts)
BG = "#ffffff"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
UNUSED = "#d9d9d6"
GRID = "#e1e0d9"
BEFORE_FILL, BEFORE_LINE = "#c3c2b7", "#898781"                # before calibration: gray
AFTER_FILL, AFTER_LINE, AFTER = "#2a78d6", "#1c5cab", "#2a78d6"  # after calibration: blue
RED = "#e34948"                                                # QC figure: excluded by the identity rule

NORM = SymLogNorm(linthresh=0.1, linscale=1.0, vmin=-100, vmax=100, base=10)
CMAP = LinearSegmentedColormap.from_list(
    "rt_delta", ["#104281", "#2a78d6", "#9ec5f4", "#ffffff", "#f3a9a8", "#e34948", "#a01f1f"])
CB_TICKS = [-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100]
# the symmetric-log mapping of the colour scale: linear within +-0.1 min, one decade per 0.1 beyond
SYMLOG = SymmetricalLogTransform(base=10, linthresh=0.1, linscale=1.0)
YTICKS = [-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100]

plt.rcParams["font.sans-serif"] = ["Arial", "Segoe UI", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"

LOG = []
DIFFS = []


def log(*parts):
    line = " ".join(str(p) for p in parts)
    print(line)
    LOG.append(line)


# ----------------------------------------------------------------------------- small helpers
def col_name(c):
    return "Column " + c.split("_")[1]


def _minus(v, _pos=None):
    """Tick label with a real minus sign (U+2212) and no trailing zeros."""
    return f"{v:g}".replace("-", "−")


def _round_half_up(v, nd):
    """v rounded to nd decimals, ties away from zero (the shortest decimal form of the float is
    rounded, so 0.15 -> 0.2 and -0.05 -> -0.1; numpy would round half to even)."""
    return Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP)


def _num(v, nd):
    """Fixed-width number with a real minus sign; a value that rounds to zero is shown without a sign."""
    r = _round_half_up(v, nd)
    return f"{abs(r) if r == 0 else r:.{nd}f}".replace("-", "−")


def _signed(v, nd):
    """Like _num, but with an explicit plus sign for positive values."""
    r = _round_half_up(v, nd)
    if r == 0:
        return f"{abs(r):.{nd}f}"
    return ("+" if r > 0 else "−") + f"{abs(r):.{nd}f}"


def _T(v):
    """The symmetric-log transform of the heatmaps' colour scale, applied to a scalar or array."""
    return SYMLOG.transform(np.asarray(v, dtype=float))


_SCRATCH = plt.figure(figsize=(1, 1), dpi=100)
_RENDERER = _SCRATCH.canvas.get_renderer()


def text_width(s, fs, weight="normal"):
    """Advance width of s in inches, as the Agg renderer lays it out."""
    w, _h, _d = _RENDERER.get_text_width_height_descent(
        s, FontProperties(size=fs, weight=weight), ismath=False)
    return w / _SCRATCH.dpi


def wrap_to_width(text, fs, max_in):
    """Break text at spaces so that no line is wider than max_in inches (existing line breaks are kept)."""
    out = []
    for para in text.split("\n"):
        cur = ""
        for word in para.split(" "):
            trial = (cur + " " + word).strip()
            if cur and text_width(trial, fs) > max_in:
                out.append(cur)
                cur = word
            else:
                cur = trial
        out.append(cur)
    return "\n".join(out)


def save(fig, path):
    """fig.savefig with retries (the output folder is on OneDrive: stale placeholders / brief locks)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for attempt in range(8):
        try:
            if os.path.exists(path):
                os.remove(path)
            fig.savefig(path, dpi=DPI, facecolor=BG)
            return
        except OSError:
            if attempt == 7:
                raise
            time.sleep(1.0)


def write_text(path, text):
    for attempt in range(8):
        try:
            if os.path.exists(path):
                os.remove(path)
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text)
            return
        except OSError:
            if attempt == 7:
                raise
            time.sleep(1.0)


def text_color_on(v):
    """White text on dark cells, ink on light ones (relative luminance of the cell colour)."""
    r, g, b, _ = CMAP(float(NORM(v)))

    def lin(u):
        return u / 12.92 if u <= 0.04045 else ((u + 0.055) / 1.055) ** 2.4

    return "#ffffff" if 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b) < 0.25 else INK


# ----------------------------------------------------------------------------- data
def load_set(csv_name):
    """One row per panel lipid (CSV order) with the used flags and dRT values of the 5 columns."""
    df = pd.read_csv(os.path.join(HERE, csv_name), dtype=str, keep_default_na=False)
    P = pd.DataFrame({"lipid": df.lipid, "display_group": df.display_group, "confidence": df.confidence,
                      "same_ion_as": df.same_ion_as,
                      "ref_rt": pd.to_numeric(df["column_22_RT_min"]),
                      "ref_aid": pd.to_numeric(df["column_22_alignment_id"])})
    assert (df["column_22_used"] == "True").all(), "the reference column must be usable for every panel lipid"
    for c in NONREF:
        used = df[f"{c}_used"] == "True"
        P[f"{c}_used"] = used
        P[f"{c}_aid"] = pd.to_numeric(df[f"{c}_alignment_id"].replace("", np.nan))
        P[f"{c}_before"] = pd.to_numeric(df[f"{c}_dRT_before_min"].replace("", np.nan))
        P[f"{c}_after"] = pd.to_numeric(df[f"{c}_dRT_after_min"].replace("", np.nan))
        assert (used == P[f"{c}_before"].notna()).all() and (used == P[f"{c}_after"].notna()).all(), \
            f"{c}: dRT values must be filled exactly for the used cells"
    return P


def pooled_cells(P):
    """Used non-reference cells, each counted once: key = (column, pick Alignment ID, reference pick
    Alignment ID); high lipids first, then low lipids, each in CSV order."""
    parts = []
    for k, c in enumerate(NONREF):
        m = P[f"{c}_used"]
        parts.append(pd.DataFrame({
            "rank": (P.confidence[m] != "high").astype(int), "row": P.index[m], "k": k, "lipid": P.lipid[m],
            "column": c, "aid": P[f"{c}_aid"][m], "ref_aid": P.ref_aid[m],
            "before": P[f"{c}_before"][m], "after": P[f"{c}_after"][m]}))
    E = pd.concat(parts, ignore_index=True).sort_values(["rank", "row", "k"], kind="stable")
    n_all = len(E)
    E = E.drop_duplicates(["column", "aid", "ref_aid"], keep="first").reset_index(drop=True)
    return E, n_all


def groups_present(P):
    groups = [g for g in GROUP_ORDER if (P.display_group == g).any()]
    assert set(groups) == set(P.display_group), "a display group is missing from GROUP_ORDER"
    return groups


# ----------------------------------------------------------------------------- heatmaps
def draw_heatmap(path, row_labels, B, A, footnote, row_h, values, row_fs=11.5, cell_w=0.95):
    """Two blocks side by side ("Before calibration" | "After calibration"), one row per entry of
    row_labels, five columns; B and A are (rows x 5) arrays of dRT in minutes (NaN = not used).
    Cell colour: CMAP(NORM(value)); unused: UNUSED gray. values=True prints the value in each cell
    (before 1 decimal, after 2 decimals, sign shown, white text on dark cells). One shared horizontal
    colour bar under the blocks, a "not used" swatch next to it, the footnote at the bottom.
    The layout is computed in inches from the measured label widths, so nothing is clipped."""
    n = len(row_labels)
    label_w = max(text_width(s, row_fs) for s in row_labels)
    pad_l, label_gap = 0.25, 0.12
    left = pad_l + max(label_w, text_width("ΔRT (min)", 12.5)) + label_gap     # x of the first block
    bw = 5 * cell_w
    gap, pad_r = 0.55, 0.3
    width = left + bw + gap + bw + pad_r
    cap_fs, title_fs, cell_fs = 11.5, 14.0, 10.5
    bottom_labels = n > 16                                  # tall figure: repeat the column names below
    y_title, y_collab, y_rows = 0.28, 0.66, 0.86            # inches from the top
    y_rows_end = y_rows + n * row_h
    y_bar = y_rows_end + (0.40 if bottom_labels else 0.0) + 0.42
    bar_h = 0.20
    y_foot = y_bar + bar_h + 0.58
    height = y_foot + 0.36

    fig = plt.figure(figsize=(width, height), facecolor=BG)

    def put(x_in, y_in, s, **kw):
        return fig.text(x_in / width, 1 - y_in / height, s, **kw)

    def axes_in(x_in, y_in, w_in, h_in):
        return fig.add_axes([x_in / width, 1 - (y_in + h_in) / height, w_in / width, h_in / height])

    for b, (title, M) in enumerate((("Before calibration", B), ("After calibration", A))):
        x0 = left + b * (bw + gap)
        ax = axes_in(x0, y_rows, bw, n * row_h)
        ax.set_xlim(0, 5)
        ax.set_ylim(n, 0)
        ax.axis("off")
        for i in range(n):
            for j in range(5):
                v = M[i, j]
                ok = np.isfinite(v)
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor=CMAP(NORM(v)) if ok else UNUSED,
                                       edgecolor="#ffffff", linewidth=1.6))
                if values and ok:
                    ax.text(j + 0.5, i + 0.5, _signed(v, 1 if b == 0 else 2), ha="center", va="center",
                            fontsize=cell_fs, color=text_color_on(v))
        put(x0 + bw / 2, y_title, title, ha="center", va="center", fontsize=title_fs, fontweight="bold", color=INK)
        for j, c in enumerate(NONREF):
            put(x0 + (j + 0.5) * cell_w, y_collab, col_name(c), ha="center", va="center", fontsize=cap_fs, color=INK)
            if bottom_labels:
                put(x0 + (j + 0.5) * cell_w, y_rows_end + 0.2, col_name(c), ha="center", va="center",
                    fontsize=cap_fs, color=INK)
    for i, s in enumerate(row_labels):
        put(left - label_gap, y_rows + (i + 0.5) * row_h, s, ha="right", va="center", fontsize=row_fs, color=INK)

    # --- shared colour bar (horizontal), its title in the label column, the "not used" swatch to its right
    bar_w = 6.0
    cax = axes_in(left, y_bar, bar_w, bar_h)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP, norm=NORM), cax=cax, orientation="horizontal",
                      ticks=CB_TICKS, format=_minus)
    cb.ax.tick_params(labelsize=11.5, colors=INK2, length=3.5)
    cb.ax.minorticks_off()
    cb.outline.set_visible(False)
    put(left - label_gap, y_bar + bar_h / 2, "ΔRT (min)", ha="right", va="center", fontsize=12.5, color=INK2)
    sw = axes_in(left + bar_w + 0.45, y_bar, 0.34, bar_h)
    sw.set_facecolor(UNUSED)
    sw.set_xticks([])
    sw.set_yticks([])
    for side in sw.spines.values():
        side.set_visible(False)
    put(left + bar_w + 0.45 + 0.34 + 0.1, y_bar + bar_h / 2, "not used", ha="left", va="center",
        fontsize=11.5, color=INK2)

    put(pad_l, y_foot, footnote, ha="left", va="center", fontsize=10.5, color=INK3)
    save(fig, path)
    plt.close(fig)


def overview_figure(P, path):
    groups = groups_present(P)
    labels, B, A = [], [], []
    for g in groups:
        rows = P.index[P.display_group == g]
        labels.append(f"{g} (n={len(rows)})")
        B.append([np.nanmedian(P.loc[rows, f"{c}_before"].to_numpy(float)) if P.loc[rows, f"{c}_used"].any()
                  else np.nan for c in NONREF])
        A.append([np.nanmedian(P.loc[rows, f"{c}_after"].to_numpy(float)) if P.loc[rows, f"{c}_used"].any()
                  else np.nan for c in NONREF])
    B, A = np.array(B, dtype=float), np.array(A, dtype=float)
    draw_heatmap(path, labels, B, A,
                 "Median ΔRT per class and column. ΔRT = RT − reference RT (Column 22) before calibration; "
                 "calibrated RT − reference RT after.", row_h=0.46, values=True)
    return groups, B, A


# ----------------------------------------------------------------------------- ECDF
def ecdf_figure(E, path, footnote):
    """Cumulative distribution of |dRT| over the pooled cells, before and after calibration, as step
    curves on a log |dRT| axis; the share within 0.2 min is labelled on both curves."""
    fig = plt.figure(figsize=(7.5, 5.9), facecolor=BG)
    ax = fig.add_axes([0.105, 0.20, 0.865, 0.765])
    ax.set_facecolor(BG)

    n = len(E)
    curves = {}
    for stage, color, lw in (("before", BEFORE_LINE, 2.0), ("after", AFTER, 2.5)):
        a = np.clip(np.sort(E[stage].abs().to_numpy()), 0.001, None)           # < 0.001 -> left edge
        share = 100.0 * np.arange(1, n + 1) / n
        ax.step(np.concatenate([[0.001], a, [100.0]]), np.concatenate([[0.0], share, [100.0]]),
                where="post", color=color, lw=lw, zorder=3)
        curves[stage] = a

    ax.set_xscale("log")
    ax.set_xlim(0.001, 100)
    xticks = [0.001, 0.01, 0.1, 1, 10, 100]
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{v:g}" for v in xticks], fontsize=12, color=INK2)
    ax.minorticks_off()
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 20))
    ax.tick_params(axis="y", labelsize=12, colors=INK3, length=4)
    for lab in ax.get_yticklabels():
        lab.set_color(INK2)
    ax.set_xlabel("|ΔRT| (min)", fontsize=13, color=INK)
    ax.set_ylabel("share of lipid-by-column points (%)", fontsize=13, color=INK)
    ax.tick_params(axis="x", colors=INK3, length=4)
    ax.grid(color=GRID, lw=0.8, zorder=1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK2)

    # --- threshold line at 0.2 min with the share each curve has there. The labels sit to the right of
    # the line, the before label just under its marker, the after label well under its marker (the
    # after curve climbs from its marker towards 100 % and would run through a label at that height)
    ax.axvline(WITHIN, color=INK3, lw=1.2, ls=(0, (4, 3)), zorder=2)
    for stage, color, y_off in (("before", BEFORE_LINE, -3.0), ("after", AFTER, -6.5)):
        a = curves[stage]
        share = 100.0 * (np.round(a, 9) <= WITHIN).mean()
        ax.plot(WITHIN, share, "o", ms=5.5, color=color, zorder=4)
        ax.text(WITHIN * 1.22, share + y_off, f"{share:.1f}% within 0.2 min", ha="left", va="center",
                fontsize=11, color=INK2)

    # --- median markers at 50 %
    # (2 / 3 decimals as in the earlier script: the pooled values come from the 4-decimal CSVs, so a 4th
    # decimal of the median is a rounding tie, not information)
    for stage, color, nd in (("before", BEFORE_LINE, 2), ("after", AFTER, 3)):
        med = np.median(curves[stage])
        ax.plot(med, 50, "o", ms=5.5, color=color, zorder=4)
        if stage == "after":
            # right-aligned, ending about 0.08 decade left of the marker
            ax.text(med / 10 ** 0.08, 50, f"median {_num(med, nd)} min", ha="right", va="center",
                    fontsize=11, color=INK2)
        else:
            ax.text(med * 1.3, 50, f"median {_num(med, nd)} min", ha="left", va="top", fontsize=11, color=INK2)

    # --- direct curve labels immediately left of their own curve at about the 72 % level
    for stage, name in (("before", "before calibration"), ("after", "after calibration")):
        x72 = curves[stage][int(np.ceil(0.72 * n)) - 1]
        ax.text(x72 / 10 ** 0.1, 72, name, ha="right", va="center", fontsize=11, color=INK2)
    ax.legend(handles=[Line2D([0], [0], color=BEFORE_LINE, lw=2.0, label="before calibration"),
                       Line2D([0], [0], color=AFTER, lw=2.5, label="after calibration")],
              loc="lower right", frameon=False, fontsize=11, labelcolor=INK, handlelength=1.6, borderaxespad=0.3)

    fig.text(0.01, 0.014, wrap_to_width(footnote, 10, 7.5 - 0.25), ha="left", va="bottom", fontsize=10,
             color=INK3, linespacing=1.3)
    save(fig, path)
    plt.close(fig)


# ----------------------------------------------------------------------------- violins
def violin_figure(E, path, scale, footnote):
    """Signed dRT per column before/after calibration: one violin pair per column plus an 'all columns'
    pair. scale="linear": y axis in minutes, density estimated in minutes. scale="symlog": the
    symmetric-log scale of the heatmaps (ticks -100 ... 100) with the density estimated in the
    transformed coordinate. Every violin is scaled so its own widest point reaches one shared maximum
    half-width (about 40 % of the distance between the two violins of a group)."""
    fig = plt.figure(figsize=(11, 6), facecolor=BG)
    ax = fig.add_axes([0.075, 0.195, 0.895, 0.745])
    ax.set_facecolor(BG)

    groups = NONREF + ["all columns"]
    gx = {g: float(i) for i, g in enumerate(NONREF)}
    gx["all columns"] = 5.5                                  # set apart by a slightly wider gap
    pair_dx = 0.19                                           # side-by-side offset
    w_max = 0.4 * 2.0 * pair_dx                              # 40 % of the pair's inner distance
    rng = np.random.default_rng(42)                          # fixed seed for the point jitter
    tf = _T if scale == "symlog" else (lambda v: np.asarray(v, dtype=float))

    kdes = {}
    for g in groups:
        sub = E if g == "all columns" else E[E.column == g]
        for stage in ("before", "after"):
            vt = tf(sub[stage].to_numpy())
            kde = gaussian_kde(vt)                           # Scott bandwidth (default)
            grid = np.linspace(vt.min(), vt.max(), 256)
            kdes[(g, stage)] = (vt, grid, kde(grid), kde)

    x_lo, x_hi = -0.78, 6.45
    med_hw = 0.5 * 1.6 * 3.5 / 72.0 * (x_hi - x_lo) / (0.895 * 11.0)   # median tick half-width, data units
    pad = dict(boxstyle="round,pad=0.18", fc=BG, ec="none")            # white pad behind median labels

    if scale == "linear":
        lo = min(E.before.min(), E.after.min())
        hi = max(E.before.max(), E.after.max())
        span = hi - lo
        ylim = (lo - 0.05 * span, hi + 0.05 * span)
        label_dy = 0.011 * (ylim[1] - ylim[0])               # label offsets above / below a violin tip

    for g in groups:
        sub = E if g == "all columns" else E[E.column == g]
        n_g = len(kdes[(g, "before")][0])
        for stage, sign, fill, line in (
                ("before", -1, BEFORE_FILL, BEFORE_LINE), ("after", +1, AFTER_FILL, AFTER_LINE)):
            vt, grid, dens, kde = kdes[(g, stage)]
            x0 = gx[g] + sign * pair_dx
            w = dens / dens.max() * w_max                    # each violin's own peak -> w_max
            xs = np.concatenate([x0 - w, (x0 + w)[::-1]])
            ys = np.concatenate([grid, grid[::-1]])
            if stage == "after":
                ax.fill(xs, ys, facecolor=fill, alpha=0.55, edgecolor="none", zorder=2)
            else:
                ax.fill(xs, ys, facecolor=fill, edgecolor="none", zorder=2)
            ax.plot(xs, ys, color=line, lw=1.0, zorder=2.5)

            w_pts = kde(vt) / dens.max() * w_max             # local width at each point
            ax.scatter(x0 + rng.uniform(-0.8, 0.8, vt.size) * w_pts, vt,
                       s=5, color=line, alpha=0.35, linewidths=0, zorder=3)

            q25, med, q75 = np.percentile(vt, [25, 50, 75])
            ax.plot([x0, x0], [q25, q75], color=INK2, lw=3.5, solid_capstyle="butt", zorder=4)
            ax.plot([x0 - med_hw, x0 + med_hw], [med, med], color=INK, lw=1.8, solid_capstyle="butt", zorder=5)
            med_raw = sub[stage].median()
            if scale == "linear" and stage == "after":
                # the "after" violin is a thin lens at 0 with long, hair-thin tails: its median label sits
                # beside the lens, just above the zero line (above the tail tip it would float far from it)
                ax.text(x0 + w_max + 0.03, med + 0.6 * label_dy, _num(med_raw, 2), ha="left", va="bottom",
                        fontsize=10, color=INK2, zorder=6)
            elif scale == "linear" and abs(med_raw) < 3.0 and g != "all columns":
                # a near-zero "before" median has no free spot beside its violin (the pad would meet the
                # neighbouring lens or straddle the zero line), so the label sits below the bottom tip
                ax.text(x0, vt.min() - label_dy, _num(med_raw, 2), ha="center", va="top",
                        fontsize=10, color=INK2, bbox=pad, zorder=6)
            else:
                ax.text(x0 + sign * (w_max + 0.06), med, _num(med_raw, 2),
                        ha="left" if sign > 0 else "right", va="center", fontsize=10,
                        color=INK2, bbox=pad, zorder=6)
        ax.text(gx[g], -0.115, f"n = {n_g}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=10.5, color=INK3)

    # --- axes: for "symlog" a linear axis in transformed units, ticks / labels set by hand
    if scale == "linear":
        t_lo, t_hi = int(np.ceil(ylim[0] / 10.0)) * 10, int(np.floor(ylim[1] / 10.0)) * 10
        yticks = list(range(t_lo, t_hi + 1, 10))
        ax.set_ylim(*ylim)
    else:
        yticks = YTICKS
        ax.set_ylim(_T(-100), _T(100))
    ax.set_yticks(tf(yticks))
    ax.set_yticklabels([_minus(v) for v in yticks], fontsize=12, color=INK2)
    ax.set_ylabel("ΔRT (min)", fontsize=13, color=INK)
    ax.set_xlim(x_lo, x_hi)
    ax.set_xticks([gx[g] for g in groups])
    ax.set_xticklabels([col_name(g) if g != "all columns" else g for g in groups], fontsize=13, color=INK)
    ax.tick_params(axis="x", length=0, pad=6)
    ax.tick_params(axis="y", colors=INK3, length=4)
    for y in tf(yticks):
        ax.axhline(y, color=GRID, lw=0.8, zorder=1)
    ax.axhline(0.0, color=INK3, lw=1.0, zorder=1.6)          # slightly stronger line at 0
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK2)

    # the legend sits where the violins are not: top left on the linear axis, bottom left on the symlog axis
    # (there nothing lies below -10 min, while the Column 12 tail of the "after" violin reaches the top)
    ax.legend(handles=[
        Patch(facecolor=BEFORE_FILL, edgecolor=BEFORE_LINE, label="before calibration"),
        Patch(facecolor=AFTER_FILL, edgecolor=AFTER_LINE, alpha=0.55, label="after calibration")],
        loc="upper left" if scale == "linear" else "lower left", ncol=1 if scale == "linear" else 2,
        columnspacing=1.8, frameon=False, fontsize=12, labelcolor=INK, handlelength=1.4, borderaxespad=0.4)

    wrapped = wrap_to_width(footnote, 10.5, 11 - 0.25)
    assert wrapped.count("\n") <= 1, "violin footnote needs more than two lines"
    fig.text(0.008, 0.018, wrapped, ha="left", va="bottom", fontsize=10.5, color=INK3, linespacing=1.3)
    save(fig, path)
    plt.close(fig)


# ----------------------------------------------------------------------------- QC figure
def qc_figure(Q, M, path):
    """x = vetted RT of the matched row, y = column_20 pick RT in this run, one point per group with a
    column_20 pick; blue = in a panel, red = excluded by the identity rule, gray = excluded otherwise;
    the mapping curve of the sheet (black, flat beyond its ends), a +-0.30 min band around it, y = x."""
    sheets = [("492", "492-lipid sheet (negative ESI)"), ("182", "182-lipid sheet (manuscript)")]
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 6.6), facecolor=BG)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.93, bottom=0.2, wspace=0.17)
    style = [("excluded_other", "#a9a8a0", 4), ("excluded_rule_h", RED, 5), ("in_panel", AFTER, 6)]
    counts = {}
    for ax, (sheet, title) in zip(axs, sheets):
        q = Q[Q.sheet == sheet]
        m = M[M.sheet == sheet].sort_values("vetted_RT_min")
        lim = float(np.ceil(max(q.vetted_RT_min.max(), q.column_20_RT_min.max(), m.column_20_RT_min.max()) + 0.5))
        xs = np.linspace(0, lim, 600)
        ys = np.interp(xs, m.vetted_RT_min, m.column_20_RT_min)               # flat beyond the ends
        ax.fill_between(xs, ys - 0.30, ys + 0.30, color=GRID, lw=0, zorder=1)
        ax.plot([0, lim], [0, lim], color=INK3, lw=1.2, ls=(0, (4, 3)), zorder=2)
        c_lo, c_hi = m.vetted_RT_min.min(), m.vetted_RT_min.max()
        inside = (xs >= c_lo) & (xs <= c_hi)
        ax.plot(xs[inside], ys[inside], color=INK, lw=1.8, zorder=3)
        for seg in (xs <= c_lo, xs >= c_hi):                       # the flat parts beyond the curve's ends
            ax.plot(xs[seg], ys[seg], color=INK, lw=1.2, ls=(0, (1, 2)), zorder=3)
        for status, color, z in style:
            s = q[q.status == status]
            ax.scatter(s.vetted_RT_min, s.column_20_RT_min, s=30, color=color, edgecolors="#ffffff",
                       linewidths=0.5, zorder=z)
            counts[(sheet, status)] = len(s)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal")
        ax.set_xticks(range(0, int(lim) + 1, 2))
        ax.set_yticks(range(0, int(lim) + 1, 2))
        ax.tick_params(labelsize=12, colors=INK3, length=4)
        for lab in ax.get_xticklabels() + ax.get_yticklabels():
            lab.set_color(INK2)
        ax.set_xlabel("Vetted RT (min)", fontsize=13, color=INK)
        ax.set_ylabel("Column 20 RT in this run (min)", fontsize=13, color=INK)
        ax.set_title(title, loc="left", fontsize=13.5, color=INK, pad=10)
        ax.grid(color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(INK2)
        ax.text(0.03, 0.97, f"{len(q)} groups\n{counts[(sheet, 'in_panel')]} in a panel\n"
                f"{counts[(sheet, 'excluded_rule_h')]} excluded, more than 0.30 min off\n"
                f"{counts[(sheet, 'excluded_other')]} excluded otherwise",
                transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=INK3, linespacing=1.35)
    handles = [
        Line2D([0], [0], marker="o", ls="", color=AFTER, markeredgecolor="#ffffff", ms=8, label="in a panel"),
        Line2D([0], [0], marker="o", ls="", color=RED, markeredgecolor="#ffffff", ms=8,
               label="excluded: column_20 pick more than 0.30 min from the mapped vetted RT"),
        Line2D([0], [0], marker="o", ls="", color="#a9a8a0", markeredgecolor="#ffffff", ms=8,
               label="excluded for another reason"),
        Line2D([0], [0], color=INK, lw=1.8, label="mapping curve (dotted: flat beyond its ends)"),
        Patch(facecolor=GRID, label="± 0.30 min around the curve"),
        Line2D([0], [0], color=INK3, lw=1.2, ls=(0, (4, 3)), label="y = x")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=11.5, labelcolor=INK,
               bbox_to_anchor=(0.5, 0.0), columnspacing=2.5, handlelength=1.8)
    save(fig, path)
    plt.close(fig)
    return counts


# ----------------------------------------------------------------------------- main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    t0 = time.time()
    log("== figures for the vetted-lipid validation panel (read from the CSVs of build_vetted_panel.py) ==")
    log("font used:", os.path.basename(findfont(FontProperties(family=list(plt.rcParams["font.sans-serif"])))))
    written = []
    for name, csv_name in SETS.items():
        P = load_set(csv_name)
        E, n_all = pooled_cells(P)
        out = os.path.join(FIG_DIR, name)
        n_lipids = len(P)
        a, b = E.after.abs().to_numpy(), E.before.abs().to_numpy()
        ra, rb = np.round(a, 9), np.round(b, 9)
        log("")
        log(f"-- set {name}: {n_lipids} lipids, {n_all} used non-reference cells, {len(E)} pooled cells "
            f"(expected 830) --")
        med_a, med_b = np.median(a), np.median(b)
        log(f"   pooled: n cells {len(E)}; median |dRT_after| {med_a:.5f} min; "
            f"{100 * (ra <= WITHIN).mean():.1f} % within 0.2 min ({int((ra <= WITHIN).sum())} of {len(E)}); "
            f"median |dRT_before| {med_b:.4f} min; {100 * (rb <= WITHIN).mean():.1f} % within 0.2 min "
            f"({int((rb <= WITHIN).sum())} of {len(E)})")
        log("   pooled cells per column: " + ", ".join(f"{c} {int((E.column == c).sum())}" for c in NONREF)
            + f"; |dRT_after| > 0.5 / 1 / 2 min: {int((ra > 0.5).sum())} / {int((ra > 1).sum())} / {int((ra > 2).sum())}")
        for label, expected, got, tol in (
                ("pooled cells", 830, len(E), 0),
                ("cells within 0.2 min after calibration", 702, int((ra <= WITHIN).sum()), 0),
                ("cells within 0.2 min before calibration", 44, int((rb <= WITHIN).sum()), 0),
                ("cells with |dRT_after| > 0.5 / 1 / 2 min", (28, 15, 7),
                 (int((ra > 0.5).sum()), int((ra > 1).sum()), int((ra > 2).sum())), 0),
                ("% within 0.2 min after calibration", 84.6, 100 * (ra <= WITHIN).mean(), 0.05),
                ("% within 0.2 min before calibration", 5.3, 100 * (rb <= WITHIN).mean(), 0.05),
                ("median |dRT_after| (brief: from unrounded values)", 0.0486, med_a, 1e-4),
                ("median |dRT_before| (brief: from unrounded values)", 6.053, med_b, 1e-3)):
            ok = got == expected if tol == 0 else abs(got - expected) <= tol
            DIFFS.extend([] if ok else [f"{name}: {label}"])
            shown = f"{got:.5g}" if isinstance(got, float) else str(got)
            log(f"   {'ok       ' if ok else 'DIFFERENT'}  {label:<52} expected {str(expected):>12}   obtained {shown:>12}"
                + (f"   (tolerance {tol:g})" if tol else ""))
        log("   note: computed from the 4-decimal values in the CSVs; with an even number of cells the median is the mean of "
            "two of them and can end in a 5 (a rounding tie); build_vetted_panel_console.log has the medians from "
            "unrounded values")

        groups, B, A = overview_figure(P, os.path.join(out, "overview_heatmap.png"))
        written.append(os.path.join(out, "overview_heatmap.png"))
        log(f"   overview heatmap table (median dRT per group and column, minutes; before | after):")
        log("     " + f"{'group':<22}" + "".join(f"{col_name(c):>11}" for c in NONREF) + "   |"
            + "".join(f"{col_name(c):>11}" for c in NONREF))
        for g, rb_, ra_ in zip(groups, B, A):
            n_g = int((P.display_group == g).sum())
            fmt = lambda v, nd: f"{'-':>11}" if not np.isfinite(v) else f"{_num(v, nd):>11}"
            log("     " + f"{g + f' (n={n_g})':<22}" + "".join(fmt(v, 2) for v in rb_) + "   |"
                + "".join(fmt(v, 3) for v in ra_))


        n_lip = f"{n_lipids} lipids"
        if name == "high_confidence":
            foot = f"{len(E)} lipid-by-column points ({n_lip}, five columns); reference column: Column 22."
        else:
            foot = (f"{len(E)} distinct lipid-by-column points ({n_lip}, five columns; identical peaks counted "
                    f"once); reference column: Column 22.")
        ecdf_figure(E, os.path.join(out, "rt_error_ecdf.png"), foot.replace("); reference", ");\nreference")
                    if name != "high_confidence" else foot)
        written.append(os.path.join(out, "rt_error_ecdf.png"))
        share = 100 * (ra <= WITHIN).mean()
        foot_lin = foot + (f"\nAfter calibration {share:.1f} % of the values lie within ±0.2 min; "
                           f"the largest |ΔRT| is {a.max():.1f} min.")
        foot_sym = foot + ("\nSymmetric-log axis, linear within ±0.1 min; stretched around zero, so lengths near 0 "
                           "are not comparable with lengths far from 0.")
        violin_figure(E, os.path.join(out, "rt_error_violin_linear.png"), "linear", foot_lin)
        violin_figure(E, os.path.join(out, "rt_error_violin_symlog.png"), "symlog", foot_sym)
        written += [os.path.join(out, "rt_error_violin_linear.png"), os.path.join(out, "rt_error_violin_symlog.png")]

    log("")
    Q = pd.read_csv(os.path.join(HERE, "vetted_rt_check_groups.csv"), keep_default_na=False)
    M = pd.read_csv(os.path.join(HERE, "vetted_rt_mapping.csv"), keep_default_na=False, dtype={"sheet": str})
    Q["sheet"] = Q["sheet"].astype(str)
    path = os.path.join(FIG_DIR, "vetted_rt_check_column_20.png")
    counts = qc_figure(Q, M, path)
    written.append(path)
    log("-- QC figure: groups with a column_20 pick per sheet and status --")
    for sheet in ("492", "182"):
        log(f"   sheet {sheet}: " + ", ".join(f"{s} {counts[(sheet, s)]}" for s in ("in_panel", "excluded_rule_h", "excluded_other")))
    d = Q.column_20_identity_delta_min.abs()
    p_ok, h_ok = d[Q.status == "in_panel"] <= 0.30, d[Q.status == "excluded_rule_h"] > 0.30
    log(f"   consistency: in-panel points with |delta| <= 0.30 min: {int(p_ok.sum())} of {len(p_ok)}; "
        f"rule-h points with |delta| > 0.30 min: {int(h_ok.sum())} of {len(h_ok)}")
    log("")
    log(f"expected-number check: {len(DIFFS)} expected number(s) differ" + (": " + "; ".join(DIFFS) if DIFFS else ""))
    log(f"wrote {len(written)} figures in {time.time() - t0:.0f} s:")
    for p in written:
        log("  " + os.path.relpath(p, HERE).replace("\\", "/"))
    write_text(os.path.join(HERE, "build_figures_console.log"), "\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
