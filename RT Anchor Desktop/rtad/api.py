"""JS-exposed API bridge for the desktop app (pywebview `js_api`)."""

from __future__ import annotations

import glob
import os
import sys
from typing import Dict, List, Optional


def _resource(rel: str) -> str:
    """Path to a bundled resource, frozen (sys._MEIPASS) or from source."""
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return os.path.join(base, rel)


class Api:
    def __init__(self):
        self.window = None
        self.result = None

    # ---- bundled example dataset ----
    def load_example(self) -> Dict:
        """Resolve the bundled demo dataset so the user can try the app with one click."""
        samples = _resource(os.path.join("example", "samples.txt"))
        standards = _resource(os.path.join("example", "standards.txt"))
        if not (os.path.exists(samples) and os.path.exists(standards)):
            return {"ok": False, "error": "Example dataset not found in this build."}
        return {"ok": True, "samples": samples, "standards": standards, "polarity": "positive"}

    # ---- native file dialogs ----
    def pick_file(self, title: str = "Choose a file") -> Optional[str]:
        import webview
        r = self.window.create_file_dialog(webview.OPEN_DIALOG, allow_multiple=False)
        return r[0] if r else None

    def pick_files(self) -> List[str]:
        import webview
        r = self.window.create_file_dialog(webview.OPEN_DIALOG, allow_multiple=True)
        return list(r) if r else []

    def pick_folder(self) -> Optional[str]:
        import webview
        r = self.window.create_file_dialog(webview.FOLDER_DIALOG)
        return r[0] if r else None

    # ---- calibration ----
    def run_calibration(self, params: Dict) -> Dict:
        try:
            from rt_anchor import calibrate
            from .darkfigures import build_bundle

            samples = (params.get("samples") or "").strip()
            if not samples or not os.path.exists(samples):
                return {"ok": False, "error": "Choose a valid sample feature table."}
            standards = (params.get("standards") or "").strip()
            if not standards or not os.path.exists(standards):
                return {"ok": False, "error": "Choose a valid standards run — it is required."}
            polarity = params.get("polarity") or "positive"
            cfg = self._build_config(params)

            single = self._resolve_single(params.get("single_files"))
            manifest = self._load_manifest(params.get("manifest"))

            res = calibrate(samples, polarity, standards_table=standards, single_files=single,
                            config=cfg, manifest=manifest)
            self.result = res
            bundle = build_bundle(res)
            bundle["ok"] = True
            return bundle
        except Exception as e:  # surface cleanly to the UI
            return {"ok": False, "error": f"{type(e).__name__}: {e}"}

    def export_csv(self) -> Dict:
        import webview
        if self.result is None:
            return {"ok": False, "error": "Nothing to export yet."}
        r = self.window.create_file_dialog(webview.SAVE_DIALOG, save_filename="calibrated.csv")
        if not r:
            return {"ok": False, "error": "cancelled"}
        path = r if isinstance(r, str) else r[0]
        if not path.lower().endswith(".csv"):
            path += ".csv"
        self.result.table.to_csv(path, index=False)
        return {"ok": True, "path": path}

    def export_report(self) -> Dict:
        import webview
        if self.result is None:
            return {"ok": False, "error": "Nothing to export yet."}
        r = self.window.create_file_dialog(webview.SAVE_DIALOG, save_filename="report")
        if not r:
            return {"ok": False, "error": "cancelled"}
        prefix = r if isinstance(r, str) else r[0]
        from rt_anchor.viz.report import write_report
        paths = write_report(self.result, prefix, formats=("html", "pdf"))
        return {"ok": True, "paths": paths}

    # ---- helpers ----
    @staticmethod
    def _build_config(params: Dict):
        """Build a CalibrationConfig from the Advanced matching parameters.

        Instrument type is intentionally NOT an input; the user tunes the actual
        knobs (m/z tolerance, RT window, min anchors, extrapolation). Blank/invalid
        fields fall back to the package defaults.
        """
        from rt_anchor import CalibrationConfig

        def _num(key, cast):
            v = params.get(key)
            if v is None or v == "":
                return None
            try:
                return cast(v)
            except (TypeError, ValueError):
                return None

        overrides = {}
        mz = _num("mz_tol_ppm", float)
        if mz is not None:
            overrides["mz_tol_ppm"] = mz
        rtw = _num("rt_window_min", float)
        if rtw is not None:
            overrides["rt_window_min"] = rtw
        mina = _num("min_anchors", int)
        if mina is not None:
            overrides["min_anchors"] = mina
        if params.get("extrapolate"):
            overrides["extrapolate"] = True
        return CalibrationConfig(**overrides)

    @staticmethod
    def _resolve_single(spec) -> Optional[List[str]]:
        if not spec:
            return None
        if isinstance(spec, list):
            files = spec
        elif os.path.isdir(spec):
            files = sorted(glob.glob(os.path.join(spec, "*.txt")) + glob.glob(os.path.join(spec, "*.csv")))
        else:
            files = sorted(glob.glob(spec))
        return files or None

    @staticmethod
    def _load_manifest(path):
        if not path:
            return None
        from rt_anchor.panel import load_manifest_csv
        return load_manifest_csv(path)
