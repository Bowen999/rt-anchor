"""Dark-themed figure + KPI bundling for the desktop app.

Reuses rt_anchor's plotly figure builders, then patches the layout to a dark,
blue-accented, brutalist palette. Figures are returned as JSON strings (Plotly's
own serialiser handles numpy) for the JS front-end to `Plotly.newPlot`.
"""

from __future__ import annotations

from typing import Dict

# dark palette (mirrors web/styles.css)
BG = "#12161D"
GRID = "#242B36"
AXIS = "#3A434F"
TEXT = "#C8D2DE"
TEXT_DIM = "#8A94A3"
BLUE = "#5AA0E0"
BLUE_DK = "#3D7CA8"


def _dark(fig, polar=False):
    fig.update_layout(
        paper_bgcolor=BG, plot_bgcolor=BG,
        font=dict(color=TEXT),
        legend=dict(font=dict(color=TEXT), bgcolor="rgba(0,0,0,0)"),
        hoverlabel=dict(bgcolor="#1A2029", bordercolor=AXIS, font=dict(color=TEXT)),
        title=None,
    )
    fig.update_xaxes(gridcolor=GRID, linecolor=AXIS, tickcolor=AXIS, zerolinecolor=GRID,
                     title_font=dict(color=TEXT_DIM), tickfont=dict(color=TEXT))
    fig.update_yaxes(gridcolor=GRID, linecolor=AXIS, tickcolor=AXIS, zerolinecolor=GRID,
                     title_font=dict(color=TEXT_DIM), tickfont=dict(color=TEXT))
    if polar:
        fig.update_layout(polar=dict(
            bgcolor=BG,
            radialaxis=dict(gridcolor=GRID, tickfont=dict(color=TEXT_DIM), linecolor=GRID),
            angularaxis=dict(gridcolor=GRID, tickfont=dict(color=TEXT), linecolor=AXIS)))
    # recolour for dark: brighten primary blue, darken light structural grays
    for tr in fig.data:
        ln = getattr(tr, "line", None)
        if ln is not None and getattr(ln, "color", None):
            ln.color = _remap(ln.color)
        mk = getattr(tr, "marker", None)
        if mk is not None:
            if getattr(mk, "color", None) and isinstance(mk.color, str):
                mk.color = _remap(mk.color)
            mkl = getattr(mk, "line", None)
            if mkl is not None and getattr(mkl, "color", None):
                mkl.color = _remap(mkl.color)
        fc = getattr(tr, "fillcolor", None)
        if fc:
            tr.fillcolor = _remap_fill(fc)
    for sh in (fig.layout.shapes or []):
        if sh.line and getattr(sh.line, "color", None):
            sh.line.color = _remap(sh.line.color)
    # annotation text (e.g. the repeatability "median" label) — keep it legible on dark
    for an in (fig.layout.annotations or []):
        f = getattr(an, "font", None)
        if f is not None and getattr(f, "color", None):
            f.color = _remap(f.color)
    return fig


# light (report) colours -> dark-mode equivalents
_CMAP = {
    "#2B5D7D": BLUE, "#1B4A6B": BLUE, "#3D7CA8": BLUE,   # primary blue -> bright
    "#E9ECEA": "#39424E", "#E6E6E6": "#39424E", "#E4E9EC": "#39424E",  # light grid/connector
    "#B0B0B0": "#5C6675",                                # faint gray -> dark (white edges kept)
    "#8A5223": "#C08552",                                # dark clay label -> bright clay (legible on dark)
}

# rgba fills that encode the primary slate blue -> bright blue (alpha preserved);
# the warm clay fill (192,133,82) is deliberately left untouched.
_FILLMAP = {"43,93,125": "90,160,224"}


def _remap(c):
    return _CMAP.get(c, c)


def _remap_fill(c):
    if isinstance(c, str):
        for src, dst in _FILLMAP.items():
            if src in c:
                return c.replace(src, dst)
    return c


def build_bundle(result) -> Dict:
    from rt_anchor.viz import metrics, performance, repeatability, structures, tic

    # LIGHT ("day mode") canvas: the plots keep their native light report theme
    # (white paper, blue data) and sit on light cards inside the dark brutalist
    # chrome. The dark-remap (`_dark`) is retained only for the selftest.
    figures = {
        "detection": metrics.figure_plotly(result).to_json(),
        "profile": tic.figure_plotly(result, "clean").to_json(),
        "warp": performance.figure_plotly(result).to_json(),
    }
    if repeatability.is_applicable(result):
        figures["repeatability"] = repeatability.figure_plotly(result).to_json()

    struct = structures.render_default_structures() if result.default_panel else {}
    m = metrics.compute_metrics(result)
    return {
        "kpis": [{"label": l, "value": str(v), "sub": s} for l, v, s in metrics.kpi_tiles(result)],
        "radar": metrics.radar_plotly(result).to_json(),
        "figures": figures,
        "notes": {"detection": metrics.DETECTION_NOTE},
        "structures": struct,
        "meta": {
            "scope": m["scope"], "polarity": m["polarity"], "source_format": m["source_format"],
            "n_features": m["n_features"], "n_detected": m["n_sample_detected"], "n_panel": m["n_panel"],
            "has_repeatability": "repeatability" in figures,
        },
    }
