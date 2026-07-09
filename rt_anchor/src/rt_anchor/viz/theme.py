"""Shared visual theme — one blue system, applied to both matplotlib and Plotly.

Aesthetic follows the Nature-family "quietly expensive" look (muted blues, gray
structure, open frame, generous whitespace) but with **deliberately larger fonts**
than print sizing, because these are on-screen / report figures where readability
beats column-density. Blue is the theme colour; a single warm clay accent marks
the one thing that should stand out (QC sample / highlight).
"""

from __future__ import annotations

from typing import List

# ---- neutrals (structure is gray, data is colour) ----
PAPER = "#FCFBF9"      # near-white warm paper
PANEL = "#FFFFFF"
TXT = "#2B2B2B"        # primary text
TXT2 = "#4D4D4D"       # secondary text
AXIS = "#7A7A7A"       # axis lines / ticks
FAINT = "#B0B0B0"
GRID = "#E9ECEA"       # whisper gridline
NEUTRAL = "#9AA3A8"    # "other" / unweighted data

# ---- blue system ----
ACCENT = "#C08552"     # muted clay/ochre — CVD-safe partner to blue
ACCENT_DK = "#8A5223"
PRIMARY = "#2B5D7D"    # deep slate blue (single-series highlight)
PRIMARY_DK = "#1B4A6B"

# sequential blue ramp (dark -> light); sample traces sample this
_BLUE_RAMP = ["#132C39", "#1B4A6B", "#2B5D7D", "#3D7CA8", "#4E8FA6", "#7FB2C4", "#A9CDD8"]

# confidence categories -> colour
CONF_COLORS = {
    "high": "#2B5D7D",
    "medium": "#7FB2C4",
    "low": "#C6D3D9",
    "none": "#D9D3CB",
}

FONT_STACK = "Helvetica, Arial, 'Helvetica Neue', 'Liberation Sans', sans-serif"
MPL_FONTS = ["Helvetica", "Arial", "Helvetica Neue", "DejaVu Sans"]


def blue_ramp(n: int, lo: float = 0.08, hi: float = 0.92) -> List[str]:
    """n colours evenly spaced along the blue ramp (dark->light)."""
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("rta_blue", _BLUE_RAMP)
    if n <= 1:
        return [PRIMARY]
    xs = np.linspace(lo, hi, n)
    return [_to_hex(cmap(x)) for x in xs]


def _to_hex(rgba) -> str:
    r, g, b = (int(round(255 * c)) for c in rgba[:3])
    return f"#{r:02X}{g:02X}{b:02X}"


def darker(hex_color: str, f: float = 0.6) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"#{int(r*f):02X}{int(g*f):02X}{int(b*f):02X}"


# ---- matplotlib ----
# base sizes (pt) for report PDFs — larger than journal print for readability
FS_BASE = 12
FS_TICK = 11
FS_AXIS = 13
FS_TITLE = 16
FS_SUB = 11.5


def apply_mpl_theme() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": MPL_FONTS,
        "font.size": FS_BASE,
        "axes.unicode_minus": True,
        "text.color": TXT,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.labelcolor": TXT2,
        "axes.labelsize": FS_AXIS,
        "axes.titlesize": FS_TITLE,
        "axes.titlecolor": TXT,
        "axes.facecolor": PANEL,
        "figure.facecolor": PAPER,
        "savefig.facecolor": PAPER,
        "xtick.color": AXIS,
        "ytick.color": AXIS,
        "xtick.labelcolor": TXT2,
        "ytick.labelcolor": TXT2,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.frameon": False,
        "legend.fontsize": FS_SUB,
        "figure.dpi": 110,
    })


def style_ax(ax) -> None:
    """Open frame (no top/right spine), sparse outward ticks, gray structure."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=AXIS, labelcolor=TXT2)


# ---- plotly ----
def plotly_template() -> dict:
    """A layout dict applied to every Plotly figure (blue theme, large fonts)."""
    return dict(
        paper_bgcolor=PAPER,
        plot_bgcolor=PANEL,
        font=dict(family=FONT_STACK, size=15, color=TXT),
        title=dict(font=dict(size=21, color=TXT), x=0.02, xanchor="left"),
        margin=dict(l=70, r=30, t=64, b=56),
        xaxis=dict(showgrid=False, zeroline=False, linecolor=AXIS, ticks="outside",
                   ticklen=4, tickcolor=AXIS, tickfont=dict(size=13, color=TXT2),
                   title=dict(font=dict(size=16, color=TXT2))),
        yaxis=dict(showgrid=True, gridcolor=GRID, gridwidth=1, zeroline=False,
                   linecolor=AXIS, ticks="outside", ticklen=4, tickcolor=AXIS,
                   tickfont=dict(size=13, color=TXT2),
                   title=dict(font=dict(size=16, color=TXT2))),
        legend=dict(font=dict(size=13, color=TXT2), bgcolor="rgba(0,0,0,0)"),
        colorway=blue_ramp(6),
        hoverlabel=dict(font=dict(family=FONT_STACK, size=13)),
    )
