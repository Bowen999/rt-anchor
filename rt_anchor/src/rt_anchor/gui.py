"""RT Anchor desktop client — a pywebview app over the public API.

A native window hosts an HTML/CSS/JS UI (``rt_anchor/webui``); this module is the
Python side: it assembles the page (inlining the stylesheet, the app script and
plotly.js) and exposes an :class:`Api` bridge the front-end calls to browse files,
run a calibration, fetch the section figures (reused verbatim from the HTML
report's plotly builders), page the calibrated table, and export the outputs.

Launch:  ``python -m rt_anchor.gui``  ·  ``rt-anchor-gui``  ·  the frozen .exe.
Requires the optional ``pywebview`` dependency (``pip install rt-anchor[app]``).
"""

from __future__ import annotations

import glob
import json
import os
import sys
import webbrowser
from typing import Dict, List, Optional

from .config import CalibrationConfig
from .errors import RtAnchorError
from .helpers import write_results
from .panel import load_manifest_csv, load_reference_csv
from .pipeline import calibrate

GITHUB_URL = "https://github.com/Bowen999/rt-anchor"
_TABLE_ROWS = 300


def _webui(name: str) -> str:
    """Read a bundled webui asset (works from source and from PyInstaller)."""
    base = getattr(sys, "_MEIPASS", None)
    if base:
        path = os.path.join(base, "rt_anchor", "webui", name)
        if not os.path.exists(path):
            path = os.path.join(base, "webui", name)
    else:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "webui", name)
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _build_html() -> str:
    from plotly.offline import get_plotlyjs
    html = _webui("index.html")
    html = html.replace("/*__CSS__*/", _webui("app.css"))
    html = html.replace("/*__PLOTLYJS__*/", get_plotlyjs())
    html = html.replace("/*__APPJS__*/", _webui("app.js"))
    return html


def _fig_to_dict(fig) -> Dict:
    import plotly.io as pio
    return json.loads(pio.to_json(fig))


def _collect_single(spec: Optional[str]) -> Optional[List[str]]:
    if not spec:
        return None
    if os.path.isdir(spec):
        files = sorted(glob.glob(os.path.join(spec, "*.txt")) +
                       glob.glob(os.path.join(spec, "*.csv")))
    else:
        files = sorted(glob.glob(spec))
    return files or None


def _opt(value) -> Optional[str]:
    v = (value or "").strip() if isinstance(value, str) else value
    return None if v in ("", "auto", None) else v


class Api:
    """Bridge exposed to the web UI as ``pywebview.api.*``."""

    def __init__(self) -> None:
        self._window = None
        self._result = None
        self._out_dir: str = ""
        self._prefix: str = "calibrated"
        self._last_paths: Dict[str, str] = {}

    def set_window(self, window) -> None:
        self._window = window

    # ---- file pickers ----------------------------------------------------
    def pick_file(self):
        return self._dialog(folder=False)

    def pick_folder(self):
        return self._dialog(folder=True)

    def _dialog(self, folder: bool):
        import webview
        if self._window is None:
            return None
        kind = webview.FOLDER_DIALOG if folder else webview.OPEN_DIALOG
        types = () if folder else ("Feature tables (*.csv;*.txt;*.tsv)", "All files (*.*)")
        res = self._window.create_file_dialog(kind, allow_multiple=False, file_types=types)
        if not res:
            return None
        return res[0] if isinstance(res, (list, tuple)) else res

    def peek_meta(self, sample_path: str) -> Dict:
        """Best-effort sibling meta.json hints for the form."""
        out: Dict = {"dir": os.path.dirname(sample_path) if sample_path else ""}
        try:
            with open(os.path.join(os.path.dirname(sample_path), "meta.json"), encoding="utf-8") as fh:
                meta = json.load(fh)
        except Exception:
            return out
        pol = str(meta.get("polarity", "")).strip().lower()
        if pol.startswith("pos"):
            out["polarity"] = "positive"
        elif pol.startswith("neg"):
            out["polarity"] = "negative"
        inst = str(meta.get("instrument", "")).lower()
        if "orbi" in inst:
            out["mz"], out["rtw"] = "8.0", "0.3"
        elif "tof" in inst:
            out["mz"], out["rtw"] = "15.0", "0.5"
        return out

    # ---- calibration -----------------------------------------------------
    def run_calibration(self, params: Dict) -> Dict:
        log: List[str] = []
        try:
            cfg = CalibrationConfig()
            cfg.mz_tol_ppm = float(params.get("mz_tol") or 15.0)
            cfg.rt_window_min = float(params.get("rt_window") or 0.5)
            cfg.min_anchors = int(params.get("min_anchors") or 6)
            cfg.extrapolate = bool(params.get("extrapolate"))
            sample = params.get("sample")
            standards = params.get("standards")
            if not sample or not standards:
                return {"ok": False, "error": "Sample table and standards run are required.", "log": log}
            manifest = load_manifest_csv(params["manifest"]) if params.get("manifest") else None
            reference = load_reference_csv(params["reference"]) if params.get("reference") else None
            single = _collect_single(params.get("single"))
            if single:
                log.append(f"per-sample tier: {len(single)} injection tables")
            result = calibrate(
                sample_table=sample, polarity=params.get("polarity", "positive"),
                standards_table=standards, single_files=single, config=cfg,
                manifest=manifest, reference=reference,
                source_format=_opt(params.get("source_format")), rt_unit=_opt(params.get("rt_unit")),
            )
            self._result = result
            log.extend(result.log)
            self._out_dir = os.path.dirname(os.path.abspath(sample))
            self._prefix = os.path.splitext(os.path.basename(sample))[0] + "_irt"
            return {"ok": True, "log": log, **self._state_payload()}
        except RtAnchorError as e:
            return {"ok": False, "error": str(e), "log": log}
        except Exception as e:  # FileNotFoundError etc.
            return {"ok": False, "error": f"{type(e).__name__}: {e}", "log": log}

    def get_state(self) -> Dict:
        if self._result is None:
            return {"has_result": False}
        return {"has_result": True, **self._state_payload()}

    def _state_payload(self) -> Dict:
        from .viz import repeatability
        return {"repeatability": bool(repeatability.is_applicable(self._result)),
                "out_dir": self._out_dir, "prefix": self._prefix}

    # ---- figures (reuse the report's plotly builders) --------------------
    def get_overview(self) -> Dict:
        from .viz import metrics
        return {"tiles": metrics.kpi_tiles(self._result),
                "radar": _fig_to_dict(metrics.radar_plotly(self._result))}

    def get_figure(self, section: str) -> Dict:
        from .viz import metrics, performance, repeatability, tic
        if self._result is None:
            return {}
        if section == "detection":
            return {**_fig_to_dict(metrics.figure_plotly(self._result)), "note": metrics.DETECTION_NOTE}
        if section == "profile":
            return _fig_to_dict(tic.figure_plotly(self._result, style="clean"))
        if section == "warp":
            return _fig_to_dict(performance.figure_plotly(self._result))
        if section == "repeatability":
            if not repeatability.is_applicable(self._result):
                return {}
            return _fig_to_dict(repeatability.figure_plotly(self._result))
        return {}

    # ---- table -----------------------------------------------------------
    def get_table(self) -> Dict:
        import numpy as np
        import pandas as pd
        df = self._result.table
        head = df.head(_TABLE_ROWS)
        cols = [str(c) for c in df.columns]

        def cell(v):
            if isinstance(v, float):
                if not np.isfinite(v):
                    return ""
                return f"{v:.4g}"
            return "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v)

        rows = [[cell(v) for v in rec] for rec in head.to_numpy().tolist()]
        return {"columns": cols, "rows": rows, "total": int(len(df))}

    # ---- export ----------------------------------------------------------
    def export_outputs(self, params: Dict) -> Dict:
        try:
            out_dir = params.get("out_dir") or self._out_dir
            prefix = params.get("prefix") or self._prefix
            if not out_dir:
                return {"ok": False, "error": "Choose an output folder."}
            fmts = tuple(f for f, on in (("html", params.get("html")), ("pdf", params.get("pdf"))) if on)
            self._out_dir = out_dir
            paths = write_results(self._result, os.path.join(out_dir, prefix),
                                  report=bool(fmts), tic_style=params.get("tic", "clean"),
                                  report_formats=fmts)
            self._last_paths = paths
            return {"ok": True, "paths": paths}
        except RtAnchorError as e:
            return {"ok": False, "error": str(e)}
        except Exception as e:
            return {"ok": False, "error": f"{type(e).__name__}: {e}"}

    # ---- open actions ----------------------------------------------------
    def open_report(self):
        p = self._last_paths.get("report_html")
        if p and os.path.exists(p):
            webbrowser.open("file://" + os.path.abspath(p))

    def open_folder(self):
        d = self._out_dir
        if d and os.path.isdir(d):
            try:
                os.startfile(d)  # type: ignore[attr-defined]
            except AttributeError:
                webbrowser.open("file://" + os.path.abspath(d))

    def open_url(self, url: str):
        if str(url).startswith(("http://", "https://")):
            webbrowser.open(url)


def _error_dialog(message: str) -> None:
    """Show a native error box (better than the raw PyInstaller traceback)."""
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(0, message, "RT Anchor", 0x10)  # MB_ICONERROR
    except Exception:
        sys.stderr.write(message + "\n")


def main() -> int:
    try:
        import webview
    except ImportError:
        sys.stderr.write("pywebview is required for the desktop client: "
                         "pip install rt-anchor[app]\n")
        return 1

    api = Api()

    # optional demo autoload (used for smoke tests / screenshots)
    demo_s, demo_std = os.environ.get("RT_ANCHOR_DEMO_SAMPLE"), os.environ.get("RT_ANCHOR_DEMO_STANDARDS")
    if demo_s and demo_std:
        api.run_calibration({"sample": demo_s, "standards": demo_std,
                             "polarity": os.environ.get("RT_ANCHOR_DEMO_POLARITY", "positive")})

    # WebView2's NavigateToString caps HTML at ~2MB and our page inlines plotly.js
    # (~3.6MB), so serve it from a temp file loaded via file:// instead.
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix="rtanchor_")
    html_path = os.path.join(tmpdir, "index.html")
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(_build_html())

    try:
        window = webview.create_window("RT Anchor", url=html_path, js_api=api,
                                       width=1320, height=880, min_size=(1040, 720),
                                       background_color="#0F141A")
        api.set_window(window)

        def _startup():
            view = os.environ.get("RT_ANCHOR_DEMO_VIEW")   # demo/screenshot helper only
            if view:
                import time
                time.sleep(3.0)
                window.evaluate_js(
                    'var b=document.querySelector(\'.nav[data-sec="%s"]\');if(b)b.click();' % view)

        webview.start(_startup)
    except Exception as e:
        _error_dialog(
            "RT Anchor couldn't start its display engine.\n\n"
            "If you downloaded this as a .zip, right-click the .zip → Properties "
            "→ tick “Unblock” → OK, then extract it again and re-run.\n\n"
            "It also needs the Microsoft Edge WebView2 Runtime (preinstalled on "
            "Windows 11) and .NET Framework 4.8.\n\n"
            f"Details: {type(e).__name__}: {e}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
