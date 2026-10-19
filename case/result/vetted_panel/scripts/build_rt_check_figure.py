"""Check of the vetted-RT -> column_20-RT mapping, drawn from the 222 rows of vetted_panel_table.csv.

One panel: x = the RT of the matched list entry (vetted_RT_min), y = its pick on column_20 in this run
(column_20_RT_before_min), one point per ion target (222). The vetted RTs come from two list sheets that
were acquired separately, so each sheet has its own mapping curve (robust LOWESS, flat beyond its ends,
vetted_rt_mapping.csv from build_vetted_panel.py) and a +-0.30 min band around it, the window inside which a
column_20 pick is accepted. The sheet of a row follows from its `adduct`: a negative-mode ion is the 492-row
MS-DIAL sheet (negative ESI); a positive-mode ion is the 182-row manuscript sheet, whose list ion is the
target itself.

Writes figures/vetted_rt_check_column_20.png and logs/build_rt_check_figure_console.log.

Run:  python build_rt_check_figure.py     (needs build_vetted_panel.py and make_final_table.py to have run)
"""
import os
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

SCRIPTS = os.path.dirname(os.path.abspath(__file__))           # .../vetted_panel/scripts
HERE = os.path.dirname(SCRIPTS)                                # .../vetted_panel
WORK = os.path.join(HERE, "work")                              # intermediate tables (removed by run_all.py)
LOGS = os.path.join(HERE, "logs")                              # console logs
TABLE = os.path.join(HERE, "vetted_panel_table.csv")
OUT = os.path.join(HERE, "figures", "vetted_rt_check_column_20.png")

DPI = 200
BAND = 0.30                                                    # min, the accepted window around the curve
BG = "#ffffff"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
GRID = "#e1e0d9"
SHEETS = [("492", "492-lipid sheet (negative ESI)", "#2a78d6", "#1c5cab", "#cfe0f6"),
          ("182", "182-lipid sheet (manuscript)", "#e0861c", "#9a560d", "#f6dfc2")]

plt.rcParams["font.sans-serif"] = ["Arial", "Segoe UI", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    os.makedirs(LOGS, exist_ok=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    T = pd.read_csv(TABLE)
    M = pd.read_csv(os.path.join(WORK, "vetted_rt_mapping.csv"), dtype={"sheet": str})
    T["sheet"] = np.where(T["adduct"].str.endswith("-"), "492", "182")
    pos = T[T.sheet == "182"]
    assert (pos.adduct == pos.target_adduct).all() and np.allclose(pos.mz, pos.target_mz), \
        "a positive-mode list ion differs from its target"
    T = T.rename(columns={"column_20_RT_before_min": "y"})

    lines = [f"== vetted RT -> column_20 RT check, {len(T)} rows of vetted_panel_table.csv =="]
    fig = plt.figure(figsize=(7.6, 8.6), facecolor=BG)
    ax = fig.add_axes([0.11, 0.25, 0.86, 0.69])
    lim = float(np.ceil(max(T.vetted_RT_min.max(), T.y.max(), M.column_20_RT_min.max()) + 0.5))
    xs = np.linspace(0, lim, 600)
    ax.plot([0, lim], [0, lim], color=INK3, lw=1.2, ls=(0, (4, 3)), zorder=2)
    for sheet, _title, fill, dark, tint in SHEETS:
        m = M[M.sheet == sheet].sort_values("vetted_RT_min")
        ys = np.interp(xs, m.vetted_RT_min, m.column_20_RT_min)               # flat beyond the ends
        ax.fill_between(xs, ys - BAND, ys + BAND, color=tint, lw=0, alpha=0.75, zorder=1)
        c_lo, c_hi = m.vetted_RT_min.min(), m.vetted_RT_min.max()
        inside = (xs >= c_lo) & (xs <= c_hi)
        ax.plot(xs[inside], ys[inside], color=dark, lw=1.8, zorder=3)
        for seg in (xs <= c_lo, xs >= c_hi):                       # the flat parts beyond the curve's ends
            ax.plot(xs[seg], ys[seg], color=dark, lw=1.2, ls=(0, (1, 2)), zorder=3)
        s = T[T.sheet == sheet]
        ax.scatter(s.vetted_RT_min, s.y, s=30, color=fill, edgecolors="#ffffff", linewidths=0.5,
                   zorder=5 if sheet == "492" else 6)
        # self-check: every point of the table lies inside the accepted window around its sheet's curve
        d = (s.y - np.interp(s.vetted_RT_min, m.vetted_RT_min, m.column_20_RT_min)).abs()
        assert d.max() <= BAND + 1e-6, (sheet, d.max())
        lines.append(f"  sheet {sheet}: {len(s)} points, max |column_20 RT - mapped vetted RT| = {d.max():.3f} min "
                     f"(window {BAND:.2f}), curve from {len(m)} knots")
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
    ax.grid(color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK2)
    n492, n182 = int((T.sheet == "492").sum()), int((T.sheet == "182").sum())
    ax.text(0.03, 0.97, f"{len(T)} ion targets\n{n492} from the 492-lipid sheet\n{n182} from the 182-lipid sheet",
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=INK3, linespacing=1.35)
    handles = [Line2D([0], [0], marker="o", ls="", color=fill, markeredgecolor="#ffffff", ms=8,
                      label=f"{title}, {int((T.sheet == sheet).sum())} ion targets")
               for sheet, title, fill, _d, _t in SHEETS]
    handles += [Line2D([0], [0], color=INK2, lw=1.8, label="mapping curve of the sheet (dotted: flat beyond its ends)"),
                Patch(facecolor=SHEETS[0][4], label=f"± {BAND:.2f} min around the curve"),
                Line2D([0], [0], color=INK3, lw=1.2, ls=(0, (4, 3)), label="y = x")]
    fig.legend(handles=handles, loc="lower center", ncol=1, frameon=False, fontsize=11.5, labelcolor=INK,
               bbox_to_anchor=(0.5, 0.0), handlelength=1.8)
    for attempt in range(8):                                       # cloud-synced folders: brief locks
        try:
            if os.path.exists(OUT):
                os.remove(OUT)
            fig.savefig(OUT, dpi=DPI, facecolor=BG)
            break
        except OSError:
            if attempt == 7:
                raise
            time.sleep(1.0)
    plt.close(fig)
    lines.append(f"  wrote {os.path.relpath(OUT, HERE)}")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(LOGS, "build_rt_check_figure_console.log"), "w", encoding="utf-8") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
