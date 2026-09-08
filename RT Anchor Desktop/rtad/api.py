"""JS-exposed API bridge for the desktop app (pywebview `js_api`).

Every method here is callable from ``web/app.js`` as ``pywebview.api.<name>()``
and must return a JSON-serialisable dict. The v2 (cross-column) contract:
``run_calibration`` needs **two runs from the same column and gradient** — the
sample and a standards run — plus a panel key, and optionally a replacement
reference pair. Everything else has a default that reproduces the validated
method.
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys
import threading
import time
from typing import Dict, List, Optional


def _resource(rel: str) -> str:
    """Path to a bundled resource, frozen (sys._MEIPASS) or from source."""
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return os.path.join(base, rel)


def _pkg_version() -> str:
    """Engine version, preferring the module's own ``__version__``.

    ``importlib.metadata`` reads the *installed distribution's* version, which
    in an editable install is whatever ``pyproject.toml`` said at install time
    and in a frozen bundle may not exist at all. ``rt_anchor.__version__`` is
    the number that ships with the code actually running.
    """
    try:
        import rt_anchor
        v = getattr(rt_anchor, "__version__", None)
        if v:
            return str(v)
    except Exception:
        pass
    try:
        from importlib.metadata import version
        return version("rt-anchor")
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


#: What the UI shows in the reference pickers until the user overrides them.
BUNDLED_REFERENCE_LABEL = "bundled Column 25"

#: Run stages the UI names while it waits. The engine reports no intermediate
#: progress, so rather than inventing sub-steps we report the boundaries we
#: genuinely control — and the front-end eases the bar *within* a stage instead
#: of pretending to know how far through it is.
RUN_STAGES = [
    "Preparing inputs",
    "Matching features and fitting the curve",
    "Building the result views",
    "Done",
]


class Api:
    def __init__(self):
        # NOTE: the window reference must stay underscore-private. pywebview's
        # bridge generator walks every PUBLIC attribute of this object
        # (dir() + getattr, recursing into anything non-callable), and a
        # public ``window`` makes it descend into the WinForms/WebView2 COM
        # objects — from a worker thread, which WebView2 forbids ("can only be
        # accessed from the UI thread") and which can wedge the message pump
        # (window shows "Not Responding" while the page still paints).
        self._window = None
        self.result = None
        self.run_info = None
        self._running = False
        self._cal_bundle = None
        self._cal_error = None
        self._cancel_event = None
        self._cal_thread = None
        # ---- engine warm-up (see warm_engine) ----
        self._engine_ready = False
        self._engine_error = None
        self._engine_t0 = None
        self._engine_seconds = None
        self._engine_thread = None
        # ---- run progress (see progress) ----
        self._stage = 0
        self._stage_t0 = None

    # ---- engine warm-up -------------------------------------------------
    def warm_engine(self) -> None:
        """Import the calibration engine on a daemon thread, at app start.

        ``import rt_anchor`` is not cheap: the package ``__init__`` pulls in
        ``crosscolumn``, which imports scipy, scikit-learn and statsmodels at
        module level. Measured cold on a mid-range machine that chain costs
        ~10 s; on the low-spec Windows boxes this app has to run on, with a
        spinning disk and a cold file cache, it is far worse.

        It used to be paid *lazily*, by whichever bridge call needed the engine
        first — which is ``mixture_previews()``, fired the moment the input
        screen loads. The window appeared, and then the standards-panel cards
        simply were not there for half a minute, with nothing on screen to say
        why. Doing it here instead overlaps the import with WebView2 start-up
        and gives the front-end something to wait on and animate.
        """
        if self._engine_thread is not None:
            return
        self._engine_t0 = time.time()

        def _warm():
            try:
                import rt_anchor  # noqa: F401
                from rt_anchor import mixtures  # noqa: F401  (the first call needs it)
            except Exception as e:
                self._engine_error = f"{type(e).__name__}: {e}"
            finally:
                self._engine_seconds = round(time.time() - self._engine_t0, 2)
                self._engine_ready = self._engine_error is None

        self._engine_thread = threading.Thread(target=_warm, daemon=True,
                                               name="rtad-engine-warmup")
        self._engine_thread.start()

    def engine_status(self) -> Dict:
        """Is the engine importable yet? Polled by the boot overlay.

        Deliberately touches nothing that could block: the front-end calls this
        every few hundred ms while the splash is up, and a bridge call that
        blocked would defeat the point.
        """
        return {"ok": True, "ready": bool(self._engine_ready),
                "error": self._engine_error,
                "elapsed": round(time.time() - self._engine_t0, 1) if self._engine_t0 else 0.0,
                "seconds": self._engine_seconds}

    # ---- bundled standard mixtures (card choice on the input screen) ----
    def mixture_previews(self) -> Dict:
        """Card meta + per-standard rows for the offered panels (15 / 21 / Others).

        Sourced from ``rt_anchor.mixtures`` so the cards can never describe a
        different panel from the one the engine actually uses.

        The front-end waits on :meth:`engine_status` before calling this, so by
        the time it runs the import is already paid for. A failed import is
        reported rather than raised: the bridge would otherwise turn it into an
        opaque JS error and the cards would just silently never appear.
        """
        try:
            from .mixtures import DEFAULT_MIXTURE, preview_payload
            return {"ok": True, "default": DEFAULT_MIXTURE, "mixtures": preview_payload()}
        except Exception as e:
            return {"ok": False, "error": f"{type(e).__name__}: {e}"}

    # ---- bundled reference datasets (the column results are expressed on) ----
    def reference_info(self) -> Dict:
        """The bundled reference pairs, and the label the pickers show by default.

        ``Cal_RT_min`` is only meaningful relative to a stated reference column,
        so the UI names it rather than leaving it implicit.
        """
        try:
            from rt_anchor.reference import DEFAULT_REFERENCE_KEY, available_references
            refs = []
            for r in available_references():
                refs.append({
                    "key": r.get("key"), "label": r.get("label"), "sub": r.get("sub"),
                    "sample": os.path.basename(r.get("sample_file", "")),
                    "standards": os.path.basename(r.get("standards_file", "")),
                    "directory": r.get("directory_path"),
                    "available": bool(r.get("available")),
                })
            return {"ok": True, "default_key": DEFAULT_REFERENCE_KEY,
                    "placeholder": BUNDLED_REFERENCE_LABEL, "references": refs}
        except Exception as e:
            return {"ok": False, "error": f"{type(e).__name__}: {e}",
                    "placeholder": BUNDLED_REFERENCE_LABEL}

    def reset_reference(self) -> Dict:
        """Clear a user-supplied reference pair, back to the bundled default.

        The engine takes ``reference_sample``/``reference_standards`` as *both
        or neither*, so the reset is a single action rather than two — clearing
        one of them alone would be refused.
        """
        return {"ok": True, "reference_sample": "", "reference_standards": "",
                "placeholder": BUNDLED_REFERENCE_LABEL}

    # ---- bundled example dataset ----
    def load_example(self) -> Dict:
        """Resolve the bundled demo dataset so the user can try the app with one click."""
        samples = _resource(os.path.join("example", "samples.txt"))
        standards = _resource(os.path.join("example", "standards.txt"))
        if not (os.path.exists(samples) and os.path.exists(standards)):
            return {"ok": False, "error": "Example dataset not found in this build."}
        return {"ok": True, "samples": samples, "standards": standards, "polarity": "positive",
                "mixture": "mix15"}   # the bundled demo run is the 15-standard mix

    # ---- open the project page in the system browser ----
    def open_github(self) -> Dict:
        import webbrowser
        webbrowser.open("https://github.com/Bowen999/rt-anchor")
        return {"ok": True}

    # ---- native file dialogs ----
    def pick_file(self, title: str = "Choose a file") -> Optional[str]:
        try:
            import webview
            r = self._window.create_file_dialog(webview.OPEN_DIALOG, allow_multiple=False)
            return r[0] if r else None
        except Exception as e:
            return None

    def pick_files(self) -> List[str]:
        try:
            import webview
            r = self._window.create_file_dialog(webview.OPEN_DIALOG, allow_multiple=True)
            return list(r) if r else []
        except Exception:
            return []

    def pick_folder(self) -> Optional[str]:
        try:
            import webview
            r = self._window.create_file_dialog(webview.FOLDER_DIALOG)
            return r[0] if r else None
        except Exception:
            return None

    # ---- calibration ----
    def run_calibration(self, params: Dict) -> Dict:
        """Launch calibration in a background thread; UI stays responsive.

        Returns ``{"ok": True, "status": "running"}`` immediately.  The front-end
        must then listen for ``window.__onCalibrationDone()`` and call
        ``get_calibration_result()`` to collect the bundle.
        """
        if self._running:
            return {"ok": False, "error": "A calibration is already in progress."}

        # ---- validate inputs (fast, keep on main thread) ----
        try:
            from . import mixtures
            samples = (params.get("samples") or "").strip()
            if not samples or not os.path.exists(samples):
                return {"ok": False, "error": "Choose a valid sample feature table."}
            standards = (params.get("standards") or "").strip()
            if not standards or not os.path.exists(standards):
                return {"ok": False,
                        "error": "Choose a valid standards mixture table — it is required, "
                                 "and it must come from the same column and gradient as the sample."}
            polarity = params.get("polarity") or "positive"
            mix_key = params.get("mixture") or mixtures.DEFAULT_MIXTURE
            try:
                mix_meta = mixtures.get_meta(mix_key)
            except KeyError:
                return {"ok": False,
                        "error": f"Unknown mixture '{mix_key}'. Known: {mixtures.mixture_keys()}."}

            ref_sample, ref_standards, ref_err = self._resolve_reference(params)
            if ref_err:
                return {"ok": False, "error": ref_err}

            cfg = self._build_config(params)
            single = self._resolve_single(params.get("single_files"))
        except Exception as e:
            return {"ok": False, "error": f"{type(e).__name__}: {e}"}

        # ---- reset state and launch background thread ----
        self._cal_bundle = None
        self._cal_error = None
        self._running = True
        self._cancel_event = threading.Event()
        self._stage = 0
        self._stage_t0 = time.time()

        def _worker():
            try:
                from rt_anchor import calibrate
                from .appviz import build_viz

                from datetime import datetime
                t0 = time.time()
                self._set_stage(1)
                res = calibrate(samples, polarity, standards_table=standards,
                                panel=mix_key,
                                reference_sample=ref_sample,
                                reference_standards=ref_standards,
                                single_files=single, config=cfg)
                if self._cancel_event.is_set():
                    self._cal_error = "cancelled"
                    return
                elapsed = time.time() - t0
                self.result = res
                self._set_stage(2)
                bundle = build_viz(res)

                m = res.model
                curve = m.get("curve", {}) or {}
                anch = m.get("anchors", {}) or {}
                irt = m.get("irt", {}) or {}
                self.run_info = _json_safe({
                    "app": "RT Anchor",
                    "method": m.get("method", "cross-column-v2"),
                    "rt_anchor_version": _pkg_version(),
                    "timestamp": datetime.now().isoformat(timespec="seconds"),
                    "duration_seconds": round(elapsed, 3),
                    "inputs": {
                        "samples": samples,
                        "standards": standards,
                        "mixture": {"key": mix_key, "label": mix_meta["label"],
                                    "n_standards": mixtures.n_standards(mix_key)},
                        "reference": m.get("reference", {}),
                        "single_files": single or [],
                        "n_single_files": len(single) if single else 0,
                    },
                    "polarity": polarity,
                    "parameters": cfg.to_dict(),
                    "results": {
                        "calibration_scope": m.get("calibration_scope"),
                        "n_features": m.get("n_features"),
                        "n_features_extrapolated": m.get("n_features_extrapolated"),
                        "confidence_counts": m.get("confidence_counts"),
                        "curve": {k: curve.get(k) for k in
                                  ("n_pairs", "n_pairs_kept", "n_pairs_standards",
                                   "n_pairs_sample", "residual_min")},
                        "anchors": {k: anch.get(k) for k in
                                    ("engaged", "n_validated", "n_used", "lam_g", "lam_c",
                                     "gate_mse_reduction", "gate_threshold", "gate_reason")},
                        "irt": {k: irt.get(k) for k in
                                ("n_landmarks", "landmark_rt_span_min", "panel_used")},
                    },
                })
                bundle["ok"] = True
                self._cal_bundle = bundle
            except Exception as e:
                self._cal_error = f"{type(e).__name__}: {e}"
            finally:
                self._set_stage(3)
                self._running = False
                # Best-effort nudge. `evaluate_js` blocks up to 20 s waiting for
                # the page and then raises, and a WebView2 that is busy laying
                # out a big result can miss it — so this is an optimisation, not
                # the delivery mechanism. The front-end also polls
                # `get_calibration_result()`, which is what actually guarantees
                # a finished run is never left hanging behind a spinner.
                try:
                    if self._window:
                        self._window.evaluate_js("window.__onCalibrationDone()")
                except Exception:
                    pass

        self._cal_thread = threading.Thread(target=_worker, daemon=True,
                                            name="rtad-calibration")
        self._cal_thread.start()
        return {"ok": True, "status": "running", "stages": RUN_STAGES}

    # ---- run progress ---------------------------------------------------
    def _set_stage(self, i: int) -> None:
        self._stage = int(i)
        self._stage_t0 = time.time()

    def progress(self) -> Dict:
        """Which stage the run is in, and how long it has been there.

        Polled by the run card. The stage index is *real* — it is set at the
        boundaries the app controls — while the seconds let the front-end ease
        the bar within a stage without ever claiming the stage is finished.
        """
        return {"ok": True, "running": bool(self._running),
                "stage": int(self._stage), "n_stages": len(RUN_STAGES),
                "label": RUN_STAGES[min(self._stage, len(RUN_STAGES) - 1)],
                "stage_seconds": round(time.time() - self._stage_t0, 1)
                                 if self._stage_t0 else 0.0,
                "done": self._cal_bundle is not None or self._cal_error is not None}

    def get_calibration_result(self) -> Dict:
        """Return the result if calibration is complete, or status."""
        if self._cal_error is not None:
            return {"ok": False, "error": self._cal_error}
        if self._cal_bundle is not None:
            return self._cal_bundle
        if self._running:
            return {"ok": False, "status": "running"}
        return {"ok": False, "error": "No calibration in progress."}

    def cancel_calibration(self) -> Dict:
        """Ask the running calibration to stop.

        ``calibrate()`` is a single uninterruptible call, so this cannot abort
        the fit itself — the flag is read at the next boundary, and the result
        is then discarded rather than rendered. Says so plainly instead of
        letting the UI imply an instant stop.
        """
        if self._cancel_event:
            self._cancel_event.set()
        return {"ok": True, "immediate": False}

    def export_csv(self) -> Dict:
        import webview
        if self.result is None:
            return {"ok": False, "error": "Nothing to export yet."}
        r = self._window.create_file_dialog(webview.SAVE_DIALOG, save_filename="calibrated.csv")
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
        r = self._window.create_file_dialog(webview.SAVE_DIALOG, save_filename="report")
        if not r:
            return {"ok": False, "error": "cancelled"}
        prefix = r if isinstance(r, str) else r[0]
        from rt_anchor.viz.report import write_report
        paths = write_report(self.result, prefix, formats=("html",))
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

    #: The kinds ``write_results`` produces, in the order the UI lists them.
    #: ``pairs_csv`` / ``landmarks_csv`` are the two v2 additions (spec §5).
    EXPORT_KIND_LABELS = {
        "calibrated_csv": "calibrated table",
        "model_json": "model.json",
        "anchors_csv": "stage-2 anchors",
        "pairs_csv": "stage-1 matched pairs",
        "landmarks_csv": "iRT landmarks",
        "log_txt": "run log",
        "report_html": "HTML report",
        "report_pdf": "PDF report",
        "run_info_json": "run info",
    }

    def export_all(self) -> Dict:
        """Write the full result bundle into a folder the user picks.

        Returns the individual files as well as the count, so the UI can name
        the two v2 outputs (``_pairs.csv`` / ``_landmarks.csv``) rather than
        just claiming a number.
        """
        import webview
        if self.result is None:
            return {"ok": False, "error": "Nothing to export yet."}
        r = self._window.create_file_dialog(webview.FOLDER_DIALOG)
        if not r:
            return {"ok": False, "error": "cancelled"}
        folder = r[0] if isinstance(r, (list, tuple)) else r
        prefix = os.path.join(folder, "rt_anchor_run")
        from rt_anchor.helpers import write_results
        paths = write_results(self.result, prefix, report=True, report_formats=("html",))
        if self.run_info is not None:
            info_path = prefix + "_run_info.json"
            with open(info_path, "w") as fh:
                json.dump(self.run_info, fh, indent=2)
            paths["run_info_json"] = info_path
        files = [{"kind": k,
                  "label": self.EXPORT_KIND_LABELS.get(k, k),
                  "name": os.path.basename(p)}
                 for k, p in paths.items()]
        # a report that failed to render is logged by write_results, not raised
        report_ok = any(k.startswith("report_") for k in paths)
        return {"ok": True, "dir": folder, "n_files": len(paths), "paths": paths,
                "files": files, "report": report_ok}

    # ---- helpers ----
    def _save_dialog(self, default_name: str) -> Optional[str]:
        import webview
        r = self._window.create_file_dialog(webview.SAVE_DIALOG, save_filename=default_name)
        if not r:
            return None
        return r if isinstance(r, str) else r[0]

    @staticmethod
    def _resolve_reference(params: Dict):
        """``(sample, standards, error)`` for the optional reference override.

        Both or neither: a run calibrated against one lab's serum and another
        lab's standards produces perfectly plausible-looking nonsense, so the
        half-specified case is caught here with a UI-shaped message rather than
        left to the engine's ``ConfigError``.
        """
        s = (params.get("reference_sample") or "").strip() or None
        t = (params.get("reference_standards") or "").strip() or None
        if s is None and t is None:
            return None, None, None
        if (s is None) != (t is None):
            missing = "reference standards run" if t is None else "reference sample run"
            return None, None, (
                f"A custom reference needs BOTH runs — the {missing} is still set to the "
                f"{BUNDLED_REFERENCE_LABEL}. Choose it too, or use “Reset to bundled "
                f"reference”.")
        for label, path in (("Reference sample", s), ("Reference standards", t)):
            if not os.path.exists(path):
                return None, None, f"{label} file not found: {path}"
        return s, t, None

    @staticmethod
    def _build_config(params: Dict):
        """Build a CalibrationConfig from the Advanced parameters.

        Instrument type is intentionally NOT an input; the user tunes the actual
        knobs. Blank/invalid fields fall back to the package defaults, which are
        the values the v2 method was validated with.
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

        def _flag(key):
            v = params.get(key)
            return None if v is None else bool(v)

        overrides = {}
        for key, cast in (("mz_tol_ppm", float), ("match_mz_tol_ppm", float),
                          ("rt_window_min", float), ("min_anchors", int),
                          ("curve_frac", float)):
            v = _num(key, cast)
            if v is not None:
                overrides[key] = v
        for key in ("extrapolate", "use_sample_anchors", "use_sample_pairs"):
            v = _flag(key)
            if v is not None:
                overrides[key] = v
        em = params.get("extrapolate_mode")
        if em in ("clamp", "linear"):
            overrides["extrapolate_mode"] = em
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
