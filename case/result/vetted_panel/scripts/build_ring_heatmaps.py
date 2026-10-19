"""Circular heatmaps of the vetted-lipid table, one figure per lipid category (RT Anchor 1.2.3, anchor-free).

Reads ../vetted_panel_table.csv (222 ion targets) and writes, for each LIPID MAPS category, one ring for
before and one for after calibration:

  figures/category_ring_heatmaps/<Category>_before.png
  figures/category_ring_heatmaps/<Category>_after.png

Categories (LIPID MAPS main classes), from the `class` column of the table:
  Glycerophospholipids  PC, EtherPC, LPC, EtherLPC, PE, EtherPE, LPE, EtherLPE, DMPE, PI, LPI, LPA
  Glycerolipids         TG, DG
  Others                SM, Cer_NS, Cer_NDS, HexCer_NS, HexCer_HS, Hex2Cer (sphingolipids),
                        CE, BileAcid (sterol lipids), CAR (fatty acyls)

The drawing is the circular heatmap of result_v1.2.3/build_comparison.py: each sector is one row of the table,
each ring one column (inner to outer: column_12, 16, 20, 41, 85), the cell colour is the delta against the
reference column column_22 (not drawn) - RT - reference RT before calibration, calibrated RT - reference RT
after - on one symmetric-log colour scale shared by every figure; an open ring with a 60-degree gap on the
right, lipid names set radially outside, curved ring names, the colour bar and a "not used" patch in the gap,
nothing inside the inner radius. A cell the table leaves empty is neutral gray. Inside a category the rows are
ordered by class (1.5-degree gap between classes) and then by reference RT. A lipid listed with two ions
(TG) carries the ion in its label (+NH4, +Na). Large categories get a larger canvas (10 in up to 45 sectors, then in proportion to the square
root of their number); the 300-degree ring is always filled, so a small category has wide sectors.

Run:  python build_ring_heatmaps.py      (needs make_final_table.py to have run)
"""
import os
import re
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle, Wedge

SCRIPTS = os.path.dirname(os.path.abspath(__file__))           # .../vetted_panel/scripts
HERE = os.path.dirname(SCRIPTS)                                # .../vetted_panel
LOGS = os.path.join(HERE, "logs")                              # console logs
FIG_DIR = os.path.join(HERE, "figures", "category_ring_heatmaps")
TABLE = os.path.join(HERE, "vetted_panel_table.csv")

NONREF = ["column_12", "column_16", "column_20", "column_41", "column_85"]
REF = "column_22"

CATEGORY = {"Glycerophospholipids": ["PC", "EtherPC", "LPC", "EtherLPC", "PE", "EtherPE", "LPE", "EtherLPE",
                                     "DMPE", "MMPE", "PS", "PG", "LPG", "PA", "PI", "LPI", "LPA"],
            "Glycerolipids": ["DG", "TG", "MG"],
            "Others": ["SM", "Cer_NS", "Cer_NDS", "Cer_NP", "HexCer_NS", "HexCer_HS", "HexCer_HDS", "Hex2Cer",
                       "CE", "BileAcid", "SSulfate", "CAR", "FA"]}
CLASS_CATEGORY = {c: cat for cat, cs in CATEGORY.items() for c in cs}
# class order inside a category: lyso, diacyl / ether, then the rest
CLASS_ORDER = ["LPC", "EtherLPC", "LPE", "EtherLPE", "LPI", "LPA", "PC", "EtherPC", "PE", "EtherPE", "DMPE",
               "MMPE", "PS", "PG", "LPG", "PA", "PI", "SM", "Cer_NS", "Cer_NDS", "Cer_NP", "HexCer_NS",
               "HexCer_HS", "HexCer_HDS", "Hex2Cer", "DG", "TG", "MG", "CE", "BileAcid", "SSulfate", "CAR", "FA"]

BG = "#ffffff"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
UNUSED = "#d9d9d6"

NORM = matplotlib.colors.SymLogNorm(linthresh=0.1, linscale=1.0, vmin=-100, vmax=100, base=10)
CMAP = LinearSegmentedColormap.from_list(
    "rt_delta", ["#104281", "#2a78d6", "#9ec5f4", "#ffffff", "#f3a9a8", "#e34948", "#a01f1f"])

CB_TICKS = [-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100]
FOOTNOTE = {
    "before": "ΔRT = RT − reference RT (Column 22). Red: later than the reference; blue: earlier.",
    "after": "ΔRT = calibrated RT − reference RT (Column 22). "
             "Red: later than the reference; blue: earlier.",
}

plt.rcParams["font.sans-serif"] = ["Arial", "Segoe UI", "DejaVu Sans"]


# ----------------------------------------------------------------------------- helpers (as in build_comparison.py)
def _em_width(s):
    """Rough text width in em: spaces count 0.3, other characters 0.55."""
    return sum(0.3 if ch == " " else 0.55 for ch in s)


def _bold_measure(renderer, fs, out_dpi=300.0):
    """measure(s) -> (advance width, ink height) in points of s in bold at fs pt, measured at the size it
    has in the 300-dpi output."""
    m = out_dpi / renderer.dpi
    prop = FontProperties(weight="bold", size=fs * m)
    pt = 72.0 / renderer.dpi / m

    def measure(s):
        w, h, _d = renderer.get_text_width_height_descent(s, prop, ismath=False)
        return w * pt, h * pt
    return measure


def _minus(v, _pos=None):
    """Tick label with a real minus sign (U+2212) and no trailing zeros."""
    return f"{v:g}".replace("-", "−")


def _save(fig, path):
    """fig.savefig with retries (cloud-synced folders: stale placeholders, brief locks after writing)."""
    for attempt in range(8):
        try:
            if os.path.exists(path):
                os.remove(path)
            fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.2, facecolor=BG)
            return
        except OSError:
            if attempt == 7:
                raise
            time.sleep(1.0)


# ----------------------------------------------------------------------------- plotted values
def load_table():
    T = pd.read_csv(TABLE)
    unknown = set(T["class"]) - set(CLASS_CATEGORY)
    assert not unknown, f"classes without a category: {unknown}"
    T["category"] = T["class"].map(CLASS_CATEGORY)
    return T


def short_ion(adduct):
    """'[M+NH4]+' -> '+NH4'."""
    return re.sub(r"^\[M|\]\d*[+-]$", "", adduct)


def delta_table(T, category, stage):
    """The values of one figure: index = labels in plot order (class in CLASS_ORDER, then reference RT),
    columns = the five ring columns, values = dRT in minutes (NaN where the table cell is empty); plus the
    class of each row."""
    S = T[T.category == category].copy()
    S["_cls"] = S["class"].map({c: i for i, c in enumerate(CLASS_ORDER)})
    S = S.sort_values(["_cls", f"{REF}_RT_before_min", "target_adduct"], kind="stable")
    w = "before" if stage == "before" else "after"
    vals = np.array([[r[f"{c}_RT_{w}_min"] - r[f"{REF}_RT_before_min"] for c in NONREF]
                     for _, r in S.iterrows()], dtype=float)
    names = S["lipid"].str.replace(";O2", "", regex=False)
    twice = names.duplicated(keep=False)
    labels = [n + (" " + short_ion(a) if t else "") for n, a, t in zip(names, S["target_adduct"], twice)]
    return pd.DataFrame(vals, index=labels, columns=NONREF), S["class"].tolist()


# ----------------------------------------------------------------------------- circular heatmap
def heatmap(table, classes, stage, path, head=""):
    """One circular heatmap (layout and sizing rules of build_comparison.heatmap). The footnote (a first line
    `head`, then the colour-scale line) sits below the canvas, so it never runs into the bottom lipid names."""
    n = len(table)
    side = 10.0 * max(1.0, np.sqrt(n / 45.0))
    figsize = (side, side)

    gap_cls = 1.5                                    # degrees between consecutive classes
    n_cls = len(dict.fromkeys(classes))
    width = (300.0 - gap_cls * (n_cls - 1)) / n
    t = 50.0
    spans, prev = [], None
    for cls in classes:
        if prev is not None and cls != prev:
            t += gap_cls
        prev = cls
        spans.append((t, t + width))
        t += width

    r0, thick = 1.0, 0.34
    R = r0 + 5 * thick
    r_lab = R + 0.12
    wscale = figsize[0] / 10.0

    labels = list(table.index)
    max_em = max(_em_width(s) for s in labels)
    arc_rad = np.radians(width)
    cap_fs = 20.0 * figsize[0] / 10.0

    def _scale(fs):
        lim_i = r_lab + 1.0
        for _ in range(20):
            lim_i = r_lab + max_em * fs / 72 / (figsize[0] / (2 * lim_i)) + 0.15
        return figsize[0] / (2 * lim_i), lim_i

    lab_fs = cap_fs
    for _ in range(50):
        per_unit, lim = _scale(lab_fs)
        arc_pt = arc_rad * r_lab * per_unit * 72
        fs_new = min(cap_fs, 0.5 * np.floor(0.86 * arc_pt / 0.5))
        if fs_new == lab_fs:
            break
        lab_fs = fs_new
    per_unit, lim = _scale(lab_fs)
    while lab_fs > 0.86 * arc_rad * r_lab * per_unit * 72:
        lab_fs -= 0.5
        per_unit, lim = _scale(lab_fs)

    fig = plt.figure(figsize=figsize, facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    for (label, row), (t1, t2) in zip(table.iterrows(), spans):
        tm = np.radians((t1 + t2) / 2)
        for k, col in enumerate(NONREF):
            v = row[col]
            ax.add_patch(Wedge((0, 0), r0 + (k + 1) * thick, t1, t2, width=thick,
                               facecolor=CMAP(NORM(v)) if np.isfinite(v) else UNUSED,
                               edgecolor="#ffffff", linewidth=0.8))
        right = np.cos(tm) >= 0
        ax.text(r_lab * np.cos(tm), r_lab * np.sin(tm), label,
                rotation=np.degrees(tm) if right else np.degrees(tm) + 180,
                rotation_mode="anchor", ha="left" if right else "right", va="center",
                fontsize=lab_fs, color=INK)

    # ring names, curved along each ring's mid-radius from its 50-degree edge clockwise into the gap
    names = [f"Column {col.split('_')[1]}" for col in NONREF]
    r_mid = [r0 + (k + 0.5) * thick for k in range(len(NONREF))]
    renderer = fig.canvas.get_renderer()
    ring_fs = min(16.0 * wscale, 0.75 * thick * per_unit * 72 / 1.25)
    max_arc = 45.0

    def name_arc(fs, k):
        w_pt = _bold_measure(renderer, fs)(names[k])[0]
        return np.degrees(w_pt / 72 / per_unit / r_mid[k])

    if name_arc(ring_fs, 0) > max_arc:
        ring_fs = 0.5 * np.floor(ring_fs * max_arc / name_arc(ring_fs, 0) / 0.5)
        while name_arc(ring_fs, 0) > max_arc:
            ring_fs -= 0.5
    measure = _bold_measure(renderer, ring_fs)
    gap_pt = 0.32 * ring_fs
    cap_pt = measure("H")[1]
    for k, name in enumerate(names):
        rb = r_mid[k] - cap_pt / 2 / 72 / per_unit
        for i, ch in enumerate(name):
            if ch == " ":
                continue
            centre_pt = measure(name[:i + 1])[0] - measure(ch)[0] / 2
            th = np.radians(50.0) - (gap_pt + centre_pt) / 72 / per_unit / r_mid[k]
            ax.text(rb * np.cos(th), rb * np.sin(th), ch, rotation=np.degrees(th) - 90,
                    rotation_mode="anchor", ha="center", va="baseline",
                    fontsize=ring_fs, fontweight="bold", color=INK)

    # colour bar in the gap, "not used" patch under it
    bar_h = 1.8
    bc = np.array([np.cos(np.radians(12.5)), np.sin(np.radians(12.5))]) * (R + 0.72)
    wf, hf = 0.24 * wscale / figsize[0], bar_h / (2 * lim)
    fx, fy = (bc[0] + lim) / (2 * lim), (bc[1] + lim) / (2 * lim)
    cb_ax = fig.add_axes([fx - wf / 2, fy - hf / 2, wf, hf])
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cb = fig.colorbar(sm, cax=cb_ax, orientation="vertical", ticks=CB_TICKS, format=_minus)
    cb_ax.set_title("ΔRT (min)", fontsize=15.0 * wscale, color=INK2, pad=14.5 * wscale)
    cb.ax.yaxis.set_ticks_position("right")
    cb.ax.tick_params(labelsize=15.0 * wscale, colors=INK2, length=3.5 * wscale)
    cb.ax.minorticks_off()
    cb.outline.set_visible(False)

    pw, ph = 0.24 * wscale / per_unit, 0.13 * wscale / per_unit
    nu_gap = (0.36 * 15.0 + 0.65 * 13.5) * wscale / 72.0 / per_unit
    py = bc[1] - bar_h / 2 - nu_gap - ph
    ax.add_patch(Rectangle((bc[0] - pw / 2, py), pw, ph, facecolor=UNUSED, edgecolor="none"))
    ax.text(bc[0] + pw / 2 + 0.05 * wscale / per_unit, py + ph / 2, "not used",
            ha="left", va="center", fontsize=13.5 * wscale, color=INK2)

    # footnote under the lowest possible label extent (the canvas edge); bbox_inches="tight" keeps it
    ax.text(-0.985 * lim, -lim, "\n".join(x for x in (head, FOOTNOTE[stage]) if x),
            ha="left", va="top", fontsize=13.5 * wscale, color=INK3, linespacing=1.4)

    _save(fig, path)
    plt.close(fig)
    return lab_fs, ring_fs, figsize[0]


# ----------------------------------------------------------------------------- main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(LOGS, exist_ok=True)
    for old in os.listdir(FIG_DIR):                  # no stale figure of an earlier grouping
        if old.endswith(("_before.png", "_after.png")):
            os.remove(os.path.join(FIG_DIR, old))
    T = load_table()
    cats = [c for c in CATEGORY if (T.category == c).any()]
    lines = [f"== circular heatmaps, one per lipid category (plotted dRT vs column_22, minutes) ==",
             f"{len(T)} rows of vetted_panel_table.csv in {len(cats)} categories",
             f"  {'figure':<30} {'rows':>5} {'used':>5} {'not used':>9} {'median |dRT|':>13} {'min':>8} {'max':>8}"
             f"  label pt / ring-name pt / canvas in"]
    drawn = {"before": 0, "after": 0}
    rows_drawn = 0
    for cat in cats:
        n_rows = int((T.category == cat).sum())
        rows_drawn += n_rows
        for stage in ("before", "after"):
            table, classes = delta_table(T, cat, stage)
            # self-check: the plotted values are exactly the table's non-empty cells of this category
            S = T[T.category == cat]
            w = stage
            want = sorted(round(float(r[f"{c}_RT_{w}_min"] - r[f"{REF}_RT_before_min"]), 9)
                          for _, r in S.iterrows() for c in NONREF if pd.notna(r[f"{c}_RT_{w}_min"]))
            got = sorted(round(float(v), 9) for v in table.to_numpy().ravel() if np.isfinite(v))
            assert want == got and len(table) == n_rows, (cat, stage)
            path = os.path.join(FIG_DIR, f"{cat.replace(' ', '_')}_{stage}.png")
            head = f"{cat}: {n_rows} ion targets, one sector each."
            lab_fs, ring_fs, side = heatmap(table, classes, stage, path, head)
            vals = table.to_numpy().ravel()
            vals = vals[np.isfinite(vals)]
            drawn[stage] += vals.size
            lines.append(f"  {cat + ', ' + stage:<30} {len(table):>5} {vals.size:>5} {table.size - vals.size:>9} "
                         f"{np.median(np.abs(vals)):>13.4f} {vals.min():>8.3f} {vals.max():>8.3f}"
                         f"  {lab_fs:.1f} / {ring_fs:.1f} / {side:.1f}")
    n_table = int(sum(T[f"{c}_RT_before_min"].notna().sum() for c in NONREF))
    assert rows_drawn == len(T) * 1 and drawn["before"] == drawn["after"] == n_table, (rows_drawn, drawn, n_table)
    lines.append(f"  rows drawn: {rows_drawn} of {len(T)}; used cells drawn: before {drawn['before']}, "
                 f"after {drawn['after']} (the table has {n_table} non-reference cells)")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(LOGS, "build_ring_heatmaps_console.log"), "w", encoding="utf-8") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
