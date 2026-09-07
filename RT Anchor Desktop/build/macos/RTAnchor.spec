# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for RT Anchor (pywebview + rt_anchor), macOS arm64.
# Location-independent — resolve the app root from this spec's own directory
# (build/macos/) and the engine's editable src checkout from the sibling repo:
#   pyinstaller build/macos/RTAnchor.spec --distpath /private/tmp/rtad_dist --workpath /private/tmp/rtad_build --noconfirm
import os
from PyInstaller.utils.hooks import collect_all, collect_submodules

APP = os.path.dirname(os.path.dirname(SPECPATH))   # build/macos/ -> app root, two up
RT_SRC = os.path.join(os.path.dirname(APP), "rt_anchor", "src")
ICON = os.path.join(SPECPATH, "icon.icns")

datas = [(os.path.join(APP, "web"), "web"),
         (os.path.join(APP, "example"), "example")]   # bundled demo dataset (Load example button)

# The v2 method calibrates onto a REFERENCE COLUMN, and the bundled reference
# pair is package *data* — PyInstaller will not pick it up from the package on
# its own, and without it every calibration fails at resolve_reference() with
# "Bundled reference 'col35' is incomplete". Ship the whole tree, keeping the
# rt_anchor/reference_data/... layout that reference.reference_root() looks for
# first under sys._MEIPASS.
REF_DATA = os.path.join(RT_SRC, "rt_anchor", "reference_data")
if not os.path.isdir(REF_DATA):
    raise SystemExit(f"reference_data not found at {REF_DATA} — the frozen app "
                     f"cannot calibrate without it")
datas += [(REF_DATA, os.path.join("rt_anchor", "reference_data"))]

binaries = []
hiddenimports = [
    "bottle", "proxy_tools",
    "webview.platforms.cocoa", "objc", "Foundation", "AppKit", "WebKit", "Quartz",
    "scipy.interpolate", "scipy.optimize", "scipy.ndimage",
    # stage 1 of the v2 curve: statsmodels LOWESS -> sklearn isotonic -> pchip.
    # Both are lazily imported inside crosscolumn.py, so the analysis cannot see
    # them without help.
    "statsmodels.api", "statsmodels.nonparametric.smoothers_lowess",
    "sklearn.isotonic",
]

# pywebview (cocoa backend) + rdkit (structure hover) + plotly (figure JSON) +
# the two fitting libraries the cross-column curve is built from
for pkg in ("webview", "rdkit", "plotly", "statsmodels", "sklearn"):
    d, b, h = collect_all(pkg)
    datas += d; binaries += b; hiddenimports += h

# the calibration engine (editable src layout)
hiddenimports += collect_submodules("rt_anchor")

a = Analysis(
    [os.path.join(APP, "app.py")],
    pathex=[APP, RT_SRC],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=[
        # GUI toolkits we don't use
        "tkinter", "PyQt5", "PyQt6", "PySide2", "PySide6",
        # dev/notebook
        "IPython", "notebook", "jupyter", "pytest", "sphinx", "black", "mypy",
        # heavy ML / unrelated libraries pulled from anaconda base (NOT used by rt_anchor)
        # NOTE: sklearn and statsmodels used to be excluded here. They are now
        # REQUIRED — the v2 stage-1 curve is statsmodels LOWESS -> sklearn
        # IsotonicRegression -> scipy PchipInterpolator. Do not re-add them.
        "tensorflow", "tensorboard", "tensorflow_probability", "keras", "torch",
        "torchvision", "torchaudio", "onnxruntime", "onnx", "jax", "jaxlib", "flax",
        "numba", "llvmlite", "bokeh", "transformers", "datasets",
        "sympy", "cv2", "gensim", "spacy", "xgboost", "lightgbm",
        "grpc", "google", "dask", "distributed", "streamlit", "gradio", "seaborn",
        "numexpr",
        # holoviz / viz stack (unused)
        "panel", "vtk", "vtkmodules", "holoviews", "hvplot", "datashader",
        "param", "pyviz_comms", "colorcet", "pyct", "bqplot", "ipywidgets",
        "altair", "pydeck", "plotly_resampler",
        # anaconda-base extras verified unused by the full Run+Export pipeline
        # (calibrate -> build_bundle/structures -> CSV -> HTML+PDF report). pandas
        # soft-imports pyarrow and degrades gracefully when it's absent.
        "pyarrow",                                  # ~96 MB — pandas optional Arrow backend
        "lxml",                                     # ~4 MB  — unused XML/HTML parser
        "skimage",                                  # ~13 MB — scikit-image
        "h5py", "tables", "numcodecs", "zarr", "pywt",   # HDF5 / zarr / wavelet stacks
        "sqlalchemy", "openpyxl", "xlrd", "bs4", "html5lib",
    ],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="RT Anchor",
    console=False,
    argv_emulation=False,
    target_arch=None,
    icon=ICON,
)
coll = COLLECT(exe, a.binaries, a.datas, name="RT Anchor")
app = BUNDLE(
    coll,
    name="RT Anchor.app",
    icon=ICON,
    bundle_identifier="com.caleyrt.rtanchor",
    info_plist={
        "NSHighResolutionCapable": True,
        "LSMinimumSystemVersion": "11.0",
        "CFBundleName": "RT Anchor",
        "CFBundleDisplayName": "RT Anchor",
        "CFBundleShortVersionString": "0.2.0",
        "CFBundleVersion": "0.2.0",
    },
)
