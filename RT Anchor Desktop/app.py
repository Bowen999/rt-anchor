"""RT Anchor — desktop app entry (pywebview).

Wraps the rt_anchor calibration engine in a native macOS window rendering a
local HTML/CSS/JS front-end. Works both frozen (PyInstaller, resources under
sys._MEIPASS) and from source.
"""

from __future__ import annotations

import os
import sys

import webview

from rtad.api import Api


def _resource(rel: str) -> str:
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, rel)


def _selftest() -> int:
    """Exercise the Run + Export code paths headlessly, to prove the frozen
    bundle has every (lazily-imported) dependency. Invoked with `--selftest`."""
    try:
        from rt_anchor import CalibrationConfig, calibrate  # noqa: F401
        from rtad import appviz  # noqa: F401  (app-native viz-data extraction)
        # report EXPORT still uses the package's plotly/matplotlib figure builders:
        from rt_anchor.viz import report, structures  # noqa: F401
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages  # noqa: F401
        import plotly.graph_objects as go  # noqa: F401  (report html)
        plt.close(plt.figure())
        n = len(structures.render_default_structures())
        print(f"SELFTEST_OK structures={n}")
        return 0
    except Exception as e:  # pragma: no cover
        import traceback
        traceback.print_exc()
        print(f"SELFTEST_FAIL {type(e).__name__}: {e}")
        return 1


def main() -> None:
    if "--selftest" in sys.argv:
        raise SystemExit(_selftest())
    api = Api()
    window = webview.create_window(
        "RT Anchor",
        _resource("web/index.html"),
        js_api=api,
        width=1200, height=820, min_size=(960, 660),
        background_color="#0C0F14",
    )
    api.window = window
    webview.start(debug=("--debug" in sys.argv))


if __name__ == "__main__":
    main()
