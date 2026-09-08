"""RT Anchor — desktop app entry (pywebview).

Wraps the rt_anchor calibration engine in a native macOS window rendering a
local HTML/CSS/JS front-end. Works both frozen (PyInstaller, resources under
sys._MEIPASS) and from source.
"""

from __future__ import annotations

import os
import sys
import threading

import webview

from rtad.api import Api


def _resource(rel: str) -> str:
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, rel)


def _selftest() -> int:
    """Prove the frozen bundle can actually calibrate. Invoked with `--selftest`.

    Importing the dependencies is not enough: under the v2 method the stage-1
    curve is ``statsmodels`` LOWESS -> ``sklearn`` isotonic -> ``scipy`` pchip,
    all imported lazily inside the engine, and the whole method needs the
    bundled ``reference_data`` tree to exist inside the bundle. Both are exactly
    the kind of thing PyInstaller drops silently, so the test runs a **real
    two-run calibration** — the bundled reference pair against itself, which is
    the cheapest honest end-to-end check available offline and has a known
    answer: calibrating a run onto its own column must return ``Cal_RT_min ==
    RT`` to machine precision. Anything else means the pipeline ran but is
    wrong, which an import check would never catch.
    """
    try:
        from rt_anchor import CalibrationConfig, calibrate  # noqa: F401
        from rt_anchor.reference import resolve_reference
        from rtad import appviz
        from rtad.api import Api  # noqa: F401  (the JS bridge itself must import)
        # report EXPORT uses the package's plotly/matplotlib figure builders:
        from rt_anchor.viz import report, structures  # noqa: F401
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages  # noqa: F401
        import plotly.graph_objects as go  # noqa: F401  (report html)
        plt.close(plt.figure())

        pair = resolve_reference()          # raises if reference_data is missing
        res = calibrate(pair.sample_path, "positive",
                        standards_table=pair.standards_path, panel="mix21")
        rt = res.rt_minutes()
        cal = res.values("Cal_RT_min")
        err = float((cal - rt).abs().max())
        bundle = appviz.build_viz(res)      # the exact payload the UI receives
        import json
        json.dumps(bundle, allow_nan=False)  # must be valid JSON for the webview
        n_struct = len(structures.render_default_structures())

        if not (err < 1e-6):
            print(f"SELFTEST_FAIL self-calibration off by {err:.3g} min "
                  f"(expected < 1e-6)")
            return 1
        print(f"SELFTEST_OK reference={pair.key} features={len(res.table)} "
              f"pairs={res.model['curve']['n_pairs']} "
              f"landmarks={res.model['irt']['n_landmarks']} "
              f"self_error={err:.2e}min structures={n_struct}")
        return 0
    except Exception as e:  # pragma: no cover
        import traceback
        traceback.print_exc()
        print(f"SELFTEST_FAIL {type(e).__name__}: {e}")
        return 1


# ---- persistent matplotlib font cache (avoid re-indexing on every launch) ----
def _cache_dir(*parts: str) -> str:
    """Per-platform cache root. Used for the matplotlib font cache.

    This used to hardcode ``~/Library/Caches``, which on Windows quietly created
    a stray ``C:\\Users\\<you>\\Library\\Caches`` tree that matplotlib was never
    pointed at anyway (``app_win`` sets ``MPLCONFIGDIR`` first, so the
    ``setdefault`` below is a no-op there).
    """
    if sys.platform == "darwin":
        base = os.path.join(os.path.expanduser("~"), "Library", "Caches")
    elif os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
    else:
        base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return os.path.join(base, *parts)


_mpl_dir = _cache_dir("RTAnchor", "mpl")
try:
    os.makedirs(_mpl_dir, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", _mpl_dir)
except OSError:
    pass


def _close_splash() -> None:
    """Dismiss the PyInstaller splash screen, if this build has one.

    Only Windows builds carry one (``Splash()`` is unsupported on macOS), and
    only frozen ones — ``pyi_splash`` does not exist when running from source,
    so the import failing is the normal case, not an error. The env-var check
    comes first: when the bootloader never started a splash, ``import
    pyi_splash`` dies noisily on the missing ``_PYI_SPLASH_IPC`` (it prints a
    traceback before raising), which a console-mode ``--selftest`` would show.
    """
    if "_PYI_SPLASH_IPC" not in os.environ:
        return
    try:
        import pyi_splash  # type: ignore[import-not-found]
    except Exception:
        return
    try:
        if pyi_splash.is_alive():
            pyi_splash.close()
    except Exception:
        pass


def main() -> None:
    if "--selftest" in sys.argv:
        _close_splash()
        raise SystemExit(_selftest())
    api = Api()
    window = webview.create_window(
        "RT Anchor",
        _resource("web/index.html"),
        js_api=api,
        width=1200, height=820, min_size=(960, 660),
        # the chrome behind the page, seen while WebView2/WKWebView paints.
        # Matches --bg in web/styles.css; it was still the pre-redesign near-black,
        # which on a slow Windows box meant several seconds of dark rectangle
        # before a warm-paper UI faded in.
        background_color="#FCFBF9",
    )
    api._window = window

    # Start importing the engine now, in parallel with WebView2 start-up,
    # instead of leaving the bill for whichever bridge call needs it first.
    api.warm_engine()

    # Hand over from the native splash the moment the page is up. Both events
    # are subscribed because either can be the one that fires first on a given
    # backend, and _close_splash() is idempotent.
    try:
        window.events.loaded += _close_splash
        window.events.shown += _close_splash
    except Exception:
        pass
    # ...and a watchdog, so a page that never loads cannot strand an
    # always-on-top splash over the user's desktop.
    _t = threading.Timer(30.0, _close_splash)
    _t.daemon = True
    _t.start()

    webview.start(debug=("--debug" in sys.argv))


if __name__ == "__main__":
    main()
