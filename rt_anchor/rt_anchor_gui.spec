# rt_anchor_gui.spec  --  PyInstaller ONEFILE, windowed (console=False) build for Windows 11.
#
#   Build:  pyinstaller --noconfirm --clean rt_anchor_gui.spec   (see build.bat)
#
# Assumes rt-anchor is pip-installed into the build venv (so `import rt_anchor`
# resolves) — build.bat does `pip install .[report]` before invoking this.
# This spec is the single source of truth for the build.

import os
from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

datas = []
binaries = []
hiddenimports = []

# ---- Plotly -----------------------------------------------------------------
# get_plotlyjs() (viz/report.py) reads plotly/package_data/plotly.min.js at
# runtime and pio.to_html pulls template JSON — these MUST be bundled, or the
# HTML report crashes with FileNotFoundError. plotly also resolves its own
# version via importlib.metadata, hence copy_metadata. Figures are built through
# plotly's dynamically-imported validator/graph_objs tree (collect_submodules).
datas += collect_data_files("plotly")
datas += copy_metadata("plotly")
hiddenimports += collect_submodules("plotly")
hiddenimports += collect_submodules("_plotly_utils")

# ---- Matplotlib -------------------------------------------------------------
# mpl-data (matplotlibrc + bundled DejaVu Sans, the guaranteed font fallback).
# Only the Agg rasteriser and the pdf backend are ever used.
datas += collect_data_files("matplotlib")
hiddenimports += ["matplotlib.backends.backend_agg",
                  "matplotlib.backends.backend_pdf"]

# ---- SciPy ------------------------------------------------------------------
# Native C-extensions PyInstaller's static analysis commonly misses. numpy and
# pandas need nothing manual (their bundled hooks handle _libs / .libs). Do NOT
# collect_submodules('scipy') — it drags in the whole test suite and bloats.
hiddenimports += ["scipy.interpolate",
                  "scipy.ndimage",
                  "scipy._lib.messagestream",
                  "scipy.special.cython_special",
                  "scipy.sparse.csgraph._validation"]

# ---- web UI assets + pywebview (WebView2 desktop shell) ---------------------
# The HTML/CSS/JS live in rt_anchor/webui and are read at runtime from _MEIPASS.
datas += [(os.path.join("src", "rt_anchor", "webui"), os.path.join("rt_anchor", "webui"))]
# bundled example dataset for the "Run example" one-click quick test
datas += [(os.path.join("src", "rt_anchor", "example"), os.path.join("rt_anchor", "example"))]
# pywebview: its bundled JS + the WebView2 loader DLLs, plus the pythonnet-backed
# WinForms backend. copy_metadata keeps importlib.metadata.version('pywebview') happy.
datas += collect_data_files("webview")
datas += copy_metadata("pywebview")
hiddenimports += collect_submodules("webview")
hiddenimports += ["clr", "webview.platforms.winforms"]
# pythonnet / clr_loader supply the .NET bridge the WinForms backend needs.
for _pkg in ("pythonnet", "clr_loader"):
    try:
        datas += collect_data_files(_pkg)
        datas += copy_metadata(_pkg)
    except Exception:
        pass
try:
    from PyInstaller.utils.hooks import collect_dynamic_libs
    binaries += collect_dynamic_libs("clr_loader")
except Exception:
    pass

# ---- optional icon (missing file must NOT fail the build) -------------------
_icon = os.path.join("src", "rt_anchor", "assets", "rt_anchor.ico")
icon = _icon if os.path.exists(_icon) else None

a = Analysis(
    ["rt_anchor_gui.py"],
    pathex=["src"],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={"matplotlib": {"backends": ["Agg", "pdf"]}},
    runtime_hooks=[],
    excludes=[
        "rdkit", "rdkit.Chem", "rdkit.Chem.Draw",   # structure hover only; heavy; degrades gracefully
        "kaleido",                                   # report.py deliberately avoids kaleido/chromium
        "PyQt5", "PyQt6", "PySide2", "PySide6", "shiboken2", "shiboken6", "wx",
        "IPython", "ipykernel", "jupyter", "jupyter_client",
        "notebook", "nbconvert", "nbformat",
        "pytest", "_pytest", "sphinx", "docutils",
        "tornado", "zmq", "numba", "torch", "tensorflow", "cv2", "sqlalchemy", "pyarrow",
        "matplotlib.backends.backend_qt5agg", "matplotlib.backends.backend_qtagg",
        "matplotlib.backends.backend_qt5", "matplotlib.backends.backend_wx",
        "matplotlib.backends.backend_wxagg", "matplotlib.backends.backend_gtk3agg",
        "matplotlib.backends.backend_webagg",
        "matplotlib.tests", "numpy.tests", "scipy.tests", "pandas.tests",
    ],
    noarchive=False,
)

pyz = PYZ(a.pure)   # PyInstaller 6.x: Analysis.zipped_data was removed

# ONEDIR: the exe + a _internal/ folder of dependencies, collected into
# dist/rt-anchor/. Distributed as a .zip — onedir triggers far fewer antivirus
# false-positives than a onefile bootloader, and starts faster (no unpacking).
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,            # binaries/datas go to COLLECT, not into the exe
    name="rt-anchor",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,                      # no `strip` on Windows
    upx=False,                        # UPX trips antivirus and can corrupt OpenBLAS DLLs
    console=False,                    # windowed GUI, no console window
    disable_windowed_traceback=False, # keep the crash dialog so users can report errors
    icon=icon,                        # None when no .ico is present
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="rt-anchor",                 # -> dist/rt-anchor/rt-anchor.exe (+ _internal/)
)
