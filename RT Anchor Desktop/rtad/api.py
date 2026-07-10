"""JS-exposed API bridge for the desktop app (pywebview `js_api`)."""

from __future__ import annotations

import glob
import json
import math
import os
import sys
from typing import Dict, List, Optional


def _resource(rel: str) -> str:
    """Path to a bundled resource, frozen (sys._MEIPASS) or from source."""
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return os.path.join(base, rel)


def _pkg_version() -> str:
    try:
        from importlib.metadata import version
        return version("rt-anchor")
    except Exception:
        try:
            import rt_anchor
            return getattr(rt_anchor, "__version__", "unknown")
        except Exception:
            return "unknown"


def _json_safe(obj):
    """Recursively replace NaN/inf with None so json.dump stays valid."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


class Api:
    def __init__(self):
        self.window = None
        self.result = None
        self.run_info = None

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

            import time
            from datetime import datetime
            t0 = time.time()
            res = calibrate(samples, polarity, standards_table=standards, single_files=single,
                            config=cfg, manifest=manifest)
            elapsed = time.time() - t0
            self.result = res
            bundle = build_bundle(res)

            m = res.model
            self.run_info = _json_safe({
                "app": "RT Anchor",
                "rt_anchor_version": _pkg_version(),
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "duration_seconds": round(elapsed, 3),
                "inputs": {
                    "samples": samples,
                    "standards": standards,
                    "single_files": single or [],
                    "n_single_files": len(single) if single else 0,
                },
                "polarity": polarity,
                "parameters": cfg.to_dict(),
                "results": {
                    "calibration_scope": m.get("calibration_scope"),
                    "n_features": m.get("n_features"),
                    "n_features_extrapolated": m.get("n_features_extrapolated"),
                    "n_anchors_used": m.get("n_anchors_used"),
                    "anchor_span_min": m.get("anchor_span_min"),
                    "loo_residual_irt": m.get("loo_residual_irt"),
                },
            })
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

    def export_run_info(self) -> Dict:
        """Save the run parameters, timing, versions and result summary as JSON."""
        if self.run_info is None:
            return {"ok": False, "error": "Nothing to export yet."}
        path = self._save_dialog("run_info.json")
        if not path:
            return {"ok": False, "error": "cancelled"}
        if not path.lower().endswith(".json"):
            path += ".json"
        with open(path, "w") as fh:
            json.dump(self.run_info, fh, indent=2)
        return {"ok": True, "path": path}

    def export_all(self) -> Dict:
        """Write the full result bundle (CSV + model + anchors + log + report +
        run_info) into a folder the user picks."""
        import webview
        if self.result is None:
            return {"ok": False, "error": "Nothing to export yet."}
        r = self.window.create_file_dialog(webview.FOLDER_DIALOG)
        if not r:
            return {"ok": False, "error": "cancelled"}
        folder = r[0] if isinstance(r, (list, tuple)) else r
        prefix = os.path.join(folder, "rt_anchor_run")
        from rt_anchor.helpers import write_results
        paths = write_results(self.result, prefix, report=True, report_formats=("html", "pdf"))
        if self.run_info is not None:
            info_path = prefix + "_run_info.json"
            with open(info_path, "w") as fh:
                json.dump(self.run_info, fh, indent=2)
            paths["run_info_json"] = info_path
        return {"ok": True, "dir": folder, "n_files": len(paths), "paths": paths}

    # ---- helpers ----
    def _save_dialog(self, default_name: str) -> Optional[str]:
        import webview
        r = self.window.create_file_dialog(webview.SAVE_DIALOG, save_filename=default_name)
        if not r:
            return None
        return r if isinstance(r, str) else r[0]

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
