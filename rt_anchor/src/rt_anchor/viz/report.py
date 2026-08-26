"""Assemble the calibration report (interactive HTML + static PDF).

Layout: a header, a **key-data band** (6 KPI tiles + a quality radar side by side),
then one figure per section — each with a copyable HTML `<h2>` title (the plotly
figures carry no baked-in title and no panel numbering). For the default panel,
hovering a standard reveals its molecular structure.

matplotlib renders the PDF, Plotly the HTML; no kaleido/chromium (PyInstaller-safe).
"""

from __future__ import annotations

import html
import json
from typing import Dict, List

from . import metrics, performance, repeatability, structures, theme, tic


def _tic_styles(style: str) -> List[str]:
    return ["clean", "realistic"] if style == "both" else [style]


# each section = (short consistent copyable HTML title, plotly figure, div id)
def _sections(result, tic_style: str):
    secs = [("Standard detection", metrics.figure_plotly(result), "fig_detection")]
    for st in _tic_styles(tic_style):
        label = "Feature-intensity profile"
        if tic_style == "both":
            label += f" ({st})"
        secs.append((label, tic.figure_plotly(result, style=st), f"fig_profile_{st}"))
    secs.append(("Calibration warp", performance.figure_plotly(result), "fig_warp"))
    if repeatability.is_applicable(result):
        secs.append(("Injection repeatability", repeatability.figure_plotly(result), "fig_repeatability"))
    return secs


# ---------------------------------------------------------------- HTML --------

def _report_date() -> str:
    from datetime import datetime
    d = datetime.now()
    return f"{d.day} {d:%B %Y}"


def build_html(result, tic_style: str = "clean") -> str:
    import plotly.io as pio
    from plotly.offline import get_plotlyjs

    def _tile(label, value, note):
        if isinstance(value, str):
            body = f'<div class="kpi-v">{html.escape(value)}</div>'
        else:
            rows = "".join(
                f'<div class="kpi-r"><span class="kpi-k">{html.escape(k)}</span>'
                f'<span class="kpi-n">{html.escape(v)}</span></div>'
                for k, v in value)
            body = f'<div class="kpi-rows">{rows}</div>'
        n = f'<div class="kpi-s">{html.escape(note)}</div>' if note else ""
        return (f'<div class="kpi"><div class="kpi-l">{html.escape(label)}</div>'
                f'{body}{n}</div>')

    tiles = "".join(_tile(l, v, s) for l, v, s in metrics.kpi_tiles(result))

    radar_div = pio.to_html(metrics.radar_plotly(result), full_html=False,
                            include_plotlyjs=False, div_id="fig_radar",
                            config={"displayModeBar": False})

    body = []
    for title, fig, div_id in _sections(result, tic_style):
        div = pio.to_html(fig, full_html=False, include_plotlyjs=False, div_id=div_id,
                          config={"displaylogo": False, "toImageButtonOptions": {"format": "svg"}})
        note = f'<p class="fignote">{metrics.DETECTION_NOTE}</p>' if div_id == "fig_detection" else ""
        body.append(f'<section class="fig"><h2 class="figtitle">{html.escape(title)}</h2>{div}{note}</section>')
    sections = "\n".join(body)

    # structure-on-hover (default panel only)
    struct_js = ""
    if result.default_panel:
        struct = structures.render_default_structures()
        if struct:
            struct_js = _STRUCT_JS.format(data=json.dumps(struct))

    js = get_plotlyjs()
    return _HTML_SHELL.format(css=_CSS, plotlyjs=js, tiles=tiles, radar=radar_div,
                              sections=sections, struct_js=struct_js,
                              date=html.escape(_report_date()))


_CSS = f"""
:root {{
  --paper:{theme.PAPER}; --surface:#FFFFFF;
  --ink:#1F2A33; --txt:{theme.TXT}; --txt2:{theme.TXT2}; --muted:#6E767C;
  --line:#E4E9EC; --line-2:#D6DEE3;
  --primary:{theme.PRIMARY}; --primary-2:#3D7CA8; --accent:{theme.ACCENT};
  --r:14px; --r-sm:10px;
  --shadow:0 1px 2px rgba(31,42,51,.045), 0 6px 20px rgba(31,42,51,.05);
}}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--paper); color:var(--txt);
        font-family:{theme.FONT_STACK}; font-size:16px; line-height:1.55;
        -webkit-font-smoothing:antialiased; text-rendering:optimizeLegibility; }}
.wrap {{ max-width:1200px; margin:0 auto; padding:48px 36px 84px; }}

/* masthead */
.masthead {{ display:flex; justify-content:space-between; align-items:flex-end; gap:28px;
             padding-bottom:22px; margin-bottom:30px; border-bottom:1px solid var(--line-2); }}
.eyebrow {{ font-size:12px; letter-spacing:.2em; font-weight:600; color:var(--primary);
            margin-bottom:11px; }}
.masthead h1 {{ font-size:32px; font-weight:600; margin:0; letter-spacing:-0.02em;
                color:var(--ink); line-height:1.08; }}
.masthead-date {{ font-size:13px; color:var(--muted); white-space:nowrap; padding-bottom:2px; }}

/* key-data band: scorecard (hairline-divided) + quality fingerprint */
.keydata {{ display:grid; grid-template-columns:1.5fr 1fr; gap:20px; align-items:stretch;
            margin-bottom:6px; }}
.kpis {{ display:grid; grid-template-columns:repeat(2,1fr); grid-auto-rows:1fr; gap:1px;
         background:var(--line); border:1px solid var(--line); border-radius:var(--r);
         overflow:hidden; box-shadow:var(--shadow); }}
.kpi {{ background:var(--surface); padding:18px 22px; display:flex; flex-direction:column;
        gap:4px; justify-content:center; }}
.kpi-l {{ font-size:11.5px; letter-spacing:.09em; text-transform:uppercase; font-weight:600;
          color:var(--muted); }}
.kpi-v {{ font-size:30px; font-weight:600; color:var(--primary); letter-spacing:-0.02em;
          line-height:1.05; font-variant-numeric:tabular-nums; }}
.kpi-s {{ font-size:12.5px; color:var(--muted); }}
.kpi-rows {{ display:flex; flex-direction:column; gap:6px; margin-top:3px; }}
.kpi-r {{ display:flex; justify-content:space-between; align-items:baseline; gap:14px; }}
.kpi-k {{ font-size:10.5px; letter-spacing:.08em; text-transform:uppercase;
          font-weight:600; color:var(--muted); white-space:nowrap; }}
.kpi-n {{ font-size:17px; font-weight:600; color:var(--primary); letter-spacing:-0.01em;
          line-height:1.1; font-variant-numeric:tabular-nums; }}
.radar {{ background:var(--surface); border:1px solid var(--line); border-radius:var(--r);
          box-shadow:var(--shadow); padding:20px 22px 10px; display:flex; flex-direction:column; }}

/* section divider */
.figs-label {{ display:flex; align-items:center; gap:14px; margin:36px 2px 4px; }}
.figs-label span {{ font-size:12px; letter-spacing:.14em; text-transform:uppercase;
                    font-weight:600; color:var(--muted); }}
.figs-label::after {{ content:""; flex:1; height:1px; background:var(--line-2); }}

/* figure cards */
section.fig {{ background:var(--surface); border:1px solid var(--line); border-radius:var(--r);
              box-shadow:var(--shadow); padding:20px 22px 14px; margin:18px 0; }}
h2.figtitle {{ font-size:20px; font-weight:600; margin:0 2px 12px; color:var(--ink);
               letter-spacing:-0.01em; user-select:text; }}
p.fignote {{ margin:14px 2px 2px; padding-left:14px; border-left:3px solid var(--primary);
             font-size:13.5px; line-height:1.6; color:var(--muted); }}
p.fignote b {{ color:var(--ink); font-weight:600; }}

/* footer */
footer {{ display:flex; justify-content:space-between; flex-wrap:wrap; gap:8px 24px;
          margin-top:40px; padding-top:16px; border-top:1px solid var(--line);
          color:var(--muted); font-size:12.5px; }}

/* structure-on-hover popup */
#rta-struct-pop {{ position:fixed; display:none; z-index:1000; background:var(--surface);
   border:1px solid var(--line-2); border-radius:var(--r-sm); padding:9px 11px;
   box-shadow:0 8px 28px rgba(31,42,51,.16); pointer-events:none; }}
#rta-struct-pop .nm {{ font-size:13px; color:var(--ink); margin-bottom:3px; font-weight:600; }}

@media (max-width:900px) {{
  .keydata {{ grid-template-columns:1fr; }}
  .kpis {{ grid-template-columns:repeat(3,1fr); }}
  .masthead {{ flex-direction:column; align-items:flex-start; gap:10px; }}
}}
@media (max-width:560px) {{
  .wrap {{ padding:32px 20px 64px; }}
  .masthead h1 {{ font-size:26px; }}
  .kpis {{ grid-template-columns:repeat(2,1fr); }}
  .kpi {{ padding:16px 16px; }}
  .kpi-v {{ font-size:26px; }}
  section.fig {{ padding:16px 14px 12px; }}
}}
"""

_STRUCT_JS = """
<div id="rta-struct-pop"></div>
<script>
var RTA_STRUCT = {data};
window.addEventListener('load', function() {{
  var pop = document.getElementById('rta-struct-pop');
  var mx = 0, my = 0;
  document.addEventListener('mousemove', function(e) {{ mx = e.clientX; my = e.clientY; }});
  function attach(id) {{
    var gd = document.getElementById(id);
    if (!gd || !gd.on) return;
    gd.on('plotly_hover', function(d) {{
      var p = d.points && d.points[0];
      var nm = p && p.customdata;
      if (Array.isArray(nm)) nm = nm[0];
      if (nm && RTA_STRUCT[nm]) {{
        pop.innerHTML = '<div class="nm">' + nm + '</div>' + RTA_STRUCT[nm];
        pop.style.display = 'block';
        pop.style.left = Math.min(mx + 14, window.innerWidth - 260) + 'px';
        pop.style.top = Math.min(my + 14, window.innerHeight - 200) + 'px';
      }}
    }});
    gd.on('plotly_unhover', function() {{ pop.style.display = 'none'; }});
  }}
  ['fig_detection', 'fig_warp'].forEach(attach);
}});
</script>
"""

_HTML_SHELL = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>rt_anchor calibration report</title>
<style>{css}</style>
<script>{plotlyjs}</script></head>
<body><div class="wrap">
<header class="masthead">
  <div class="masthead-id">
    <div class="eyebrow">rt_anchor</div>
    <h1>Retention-index calibration report</h1>
  </div>
  <div class="masthead-date">{date}</div>
</header>
<div class="keydata">
  <div class="kpis">{tiles}</div>
  <div class="radar">
    <h2 class="figtitle">Quality fingerprint</h2>{radar}
  </div>
</div>
<div class="figs-label"><span>Detailed figures</span></div>
{sections}
<footer>
  <span>Generated by rt_anchor · dimensionless iRT scale</span>
  <span>Interactive HTML — hover, zoom, pan · a static PDF companion accompanies this file</span>
</footer>
</div>{struct_js}</body></html>"""


# ---------------------------------------------------------------- PDF ---------

def build_pdf(result, path: str, tic_style: str = "clean") -> None:
    theme.apply_mpl_theme()
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    def titled(fig, title):
        fig.suptitle(title, x=0.02, y=0.985, ha="left", fontsize=theme.FS_TITLE, color=theme.TXT)
        return fig

    with PdfPages(path) as pdf:
        pdf.savefig(_cover_page(result)); plt.close("all")
        pdf.savefig(titled(metrics.figure_mpl(result),
                    "Standard detection — reference vs panel vs samples")); plt.close("all")
        for st in _tic_styles(tic_style):
            pdf.savefig(tic.figure_mpl(result, style=st)); plt.close("all")   # tic has its own suptitle
        pdf.savefig(titled(performance.figure_mpl(result), r"Calibration warp — RT $\rightarrow$ iRT")); plt.close("all")
        if repeatability.is_applicable(result):
            pdf.savefig(repeatability.figure_mpl(result)); plt.close("all")   # has its own suptitle


def _cover_page(result):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(11.5, 8.0))
    fig.text(0.06, 0.92, "Retention-index calibration report", fontsize=24,
             color=theme.TXT, weight="bold")
    for i, (lab, val, sub) in enumerate(metrics.kpi_tiles(result)):
        r, c = divmod(i, 3)
        x = 0.06 + c * 0.30
        y = 0.78 - r * 0.15
        fig.text(x, y, lab.upper(), fontsize=10.5, color=theme.TXT2, weight="bold")
        y -= 0.036
        if isinstance(val, str):
            fig.text(x, y, val, fontsize=21, color=theme.PRIMARY, weight="bold")
            y -= 0.034
        else:
            for k, v in val:
                fig.text(x, y, k, fontsize=9.5, color=theme.TXT2)
                fig.text(x + 0.12, y, v, fontsize=14.5, color=theme.PRIMARY, weight="bold")
                y -= 0.030
        if sub:
            fig.text(x, y, sub, fontsize=10, color=theme.TXT2)
    # radar on the lower half
    fig.text(0.06, 0.48, "Quality fingerprint", fontsize=13.5, color=theme.TXT, weight="bold")
    axr = fig.add_axes([0.30, 0.04, 0.4, 0.34], polar=True)
    metrics.radar_mpl(result, ax=axr)
    return fig


# ---------------------------------------------------------------- entry -------

def write_report(result, out_prefix: str, tic_style: str = "clean",
                 formats=("html", "pdf")) -> Dict[str, str]:
    import os
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)
    paths: Dict[str, str] = {}
    if "html" in formats:
        p = f"{out_prefix}_report.html"
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(build_html(result, tic_style=tic_style))
        paths["report_html"] = p
    if "pdf" in formats:
        p = f"{out_prefix}_report.pdf"
        build_pdf(result, p, tic_style=tic_style)
        paths["report_pdf"] = p
    return paths
