# -*- mode: python ; coding: utf-8 -*-
# RT Anchor desktop — PyInstaller ONEDIR, windowed build for Windows 10/11.
#
#   Build:  build_win.bat   (or:)
#           pyinstaller --noconfirm --clean build\RTAnchor.win.spec
#
# Uses the rt_anchor build venv (`rt_anchor\build_venv`) where the engine is
# pip-installed, so `import rt_anchor` resolves — same layout as
# rt_anchor\rt_anchor_gui.spec (the proven Windows webview build in this repo).
#
# Backend: pywebview → WinForms/WebView2 via pythonnet (.NET Framework).
import os
from PyInstaller.utils.hooks import (collect_all, collect_data_files,
                                     collect_submodules, copy_metadata,
                                     collect_dynamic_libs)

APP = os.path.dirname(os.path.dirname(SPECPATH))   # build/windows/ -> app root, two up
_icon = os.path.join(SPECPATH, "icon.ico")          # sits next to this spec
icon = _icon if os.path.exists(_icon) else None

datas = [(os.path.join(APP, "web"), "web"),
         (os.path.join(APP, "example"), "example")]   # bundled demo dataset

# The v2 method needs the bundled reference pair — without it every calibration
# fails. Resolve it from the sibling engine checkout (same layout the macOS spec
# and build_win.bat's build_venv path assume), and fail loudly if absent.
RT_SRC = os.path.join(os.path.dirname(APP), "rt_anchor", "src")
REF_DATA = os.path.join(RT_SRC, "rt_anchor", "reference_data")
if not os.path.isdir(REF_DATA):
    raise SystemExit(f"reference_data not found at {REF_DATA} — the frozen app "
                     f"cannot calibrate without it")
datas += [(REF_DATA, os.path.join("rt_anchor", "reference_data"))]
binaries = []
hiddenimports = []

# ---- Plotly -----------------------------------------------------------------
# figure_plotly().to_json() needs plotly's package data + template JSON and
# resolves its version via importlib.metadata; graph_objs is imported
# dynamically, so collect the whole tree.
datas += collect_data_files("plotly")
datas += copy_metadata("plotly")
hiddenimports += collect_submodules("plotly")
hiddenimports += collect_submodules("_plotly_utils")

# ---- Matplotlib -------------------------------------------------------------
# Only Agg + the pdf backend are used (report PDF); mpl-data carries the
# bundled DejaVu fonts used by the report figures.
datas += collect_data_files("matplotlib")
hiddenimports += ["matplotlib.backends.backend_agg",
                  "matplotlib.backends.backend_pdf"]

# ---- SciPy ------------------------------------------------------------------
# Native C-extensions static analysis commonly misses. Do NOT
# collect_submodules('scipy') — drags in the whole test suite.
hiddenimports += ["scipy.interpolate",
                  "scipy.ndimage",
                  "scipy._lib.messagestream",
                  "scipy.special.cython_special",
                  "scipy.sparse.csgraph._validation"]

# ---- statsmodels + sklearn (v2 curve fitting) ----
hiddenimports += ["statsmodels.api", "statsmodels.nonparametric.smoothers_lowess",
                  "sklearn.isotonic"]
for _pkg in ("statsmodels", "sklearn"):
    try:
        datas += collect_data_files(_pkg)
        datas += copy_metadata(_pkg)
    except Exception:
        pass

# ---- pywebview (WebView2 desktop shell) -------------------------------------
# Bundled JS + WebView2Loader DLLs + the pythonnet-backed WinForms backend;
# copy_metadata keeps importlib.metadata.version('pywebview') happy.
datas += collect_data_files("webview")
datas += copy_metadata("pywebview")
hiddenimports += collect_submodules("webview")
hiddenimports += ["clr", "webview.platforms.winforms",
                  "webview.platforms.edgechromium"]
# pythonnet / clr_loader supply the .NET bridge the WinForms backend needs.
for _pkg in ("pythonnet", "clr_loader"):
    try:
        datas += collect_data_files(_pkg)
        datas += copy_metadata(_pkg)
    except Exception:
        pass
try:
    binaries += collect_dynamic_libs("clr_loader")
except Exception:
    pass

# ---- rdkit ------------------------------------------------------------------
# Molecular-structure hover for the built-in panel (~90 MB). If rdkit is not
# installed in the build env, drop this block — the app degrades gracefully
# (structures={} and no hover cards).
try:
    import rdkit  # noqa: F401
    _d, _b, _h = collect_all("rdkit")
    datas += _d; binaries += _b; hiddenimports += _h
    HAS_RDKIT = True
except Exception:
    HAS_RDKIT = False

# ---- the calibration engine -------------------------------------------------
hiddenimports += collect_submodules("rt_anchor")
try:
    datas += copy_metadata("rt-anchor")
except Exception:
    pass

a = Analysis(
    [os.path.join(APP, "app_win.py")],
    pathex=[APP],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={"matplotlib": {"backends": ["agg", "pdf"]}},
    runtime_hooks=[],
    excludes=[
        # GUI toolkits we don't use (tkinter included — pywebview doesn't need it)
        "tkinter", "PyQt5", "PyQt6", "PySide2", "PySide6", "shiboken2", "wx",
        # kaleido/chromium — report.py deliberately avoids headless export
        "kaleido",
        # dev / notebook
        "IPython", "ipykernel", "jupyter", "jupyter_client", "notebook",
        "nbconvert", "nbformat", "pytest", "_pytest", "sphinx", "docutils",
        # heavy ML / unrelated stacks that must never sneak into the bundle
        "torch", "tensorflow", "tensorflow_probability", "keras", "onnx",
        "numba", "llvmlite", "transformers", "datasets",
        "sympy", "cv2", "xgboost", "lightgbm",
        "dask", "distributed", "streamlit", "gradio", "seaborn",
        "zmq", "tornado",
        # unused optional deps verified against the full Run+Export pipeline
        "pyarrow", "sqlalchemy", "openpyxl", "xlrd", "bs4", "html5lib",
        "h5py", "tables", "zarr", "pywt",
        # unused viz stacks
        "panel", "vtk", "holoviews", "hvplot", "datashader", "bokeh",
        "altair", "pydeck", "bqplot", "ipywidgets",
        # matplotlib backends we never use
        "matplotlib.backends.backend_qt5agg", "matplotlib.backends.backend_qtagg",
        "matplotlib.backends.backend_tkagg", "matplotlib.backends.backend_wx",
        "matplotlib.backends.backend_wxagg", "matplotlib.backends.backend_gtk3agg",
        "matplotlib.backends.backend_webagg",
        "matplotlib.tests", "numpy.tests", "scipy.tests", "pandas.tests",
    ],
    noarchive=False,
)

pyz = PYZ(a.pure)

# ONEDIR: exe + _internal\ of dependencies -> dist\RT Anchor\. Distributed as a
# .zip — onedir triggers far fewer antivirus false-positives than onefile and
# starts faster (no unpacking).
#
# NOTE: deliberately NO PyInstaller Splash() here. The splash needs Tcl/Tk
# DLLs + script libraries collected into the bundle, and that collection
# silently no-ops for a Microsoft Store Python / venv build host (the tcl tree
# lives under the ACL-restricted WindowsApps dir) — the result is an app that
# pops "failed to load Tcl DLL" / "SPLASH: failed to load tcl/tk shared
# libraries" dialogs on every launch. The HTML boot overlay in web/index.html
# is the loading mechanism instead: it paints with the first frame and narrates
# the bridge + engine wait, which is where the real delay lives on a low-spec
# box. app.py's _close_splash() stays as a harmless no-op guard.
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="RT Anchor",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    icon=icon,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="RT Anchor",                 # -> dist\RT Anchor\"RT Anchor.exe"
)
