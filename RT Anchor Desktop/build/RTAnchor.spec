# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for RT Anchor (pywebview + rt_anchor). Run from the app dir:
#   pyinstaller build/RTAnchor.spec --distpath /private/tmp/rtad_dist --workpath /private/tmp/rtad_build --noconfirm
import os
from PyInstaller.utils.hooks import collect_all, collect_submodules

APP = "/Users/bowen/Desktop/Caley_RT/RT Anchor Desktop"
RT_SRC = "/Users/bowen/Desktop/Caley_RT/rt_anchor/src"
ICON = os.path.join(APP, "build", "icon.icns")

datas = [(os.path.join(APP, "web"), "web"),
         (os.path.join(APP, "example"), "example")]   # bundled demo dataset (Load example button)
binaries = []
hiddenimports = [
    "bottle", "proxy_tools",
    "webview.platforms.cocoa", "objc", "Foundation", "AppKit", "WebKit", "Quartz",
    "scipy.interpolate", "scipy.optimize", "scipy.ndimage",
]

# pywebview (cocoa backend) + rdkit (structure hover) + plotly (figure JSON)
for pkg in ("webview", "rdkit", "plotly"):
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
        "tensorflow", "tensorboard", "tensorflow_probability", "keras", "torch",
        "torchvision", "torchaudio", "onnxruntime", "onnx", "jax", "jaxlib", "flax",
        "numba", "llvmlite", "bokeh", "transformers", "datasets", "sklearn",
        "scikit_learn", "sympy", "cv2", "gensim", "spacy", "xgboost", "lightgbm",
        "grpc", "google", "dask", "distributed", "streamlit", "gradio", "seaborn",
        "statsmodels", "numexpr",
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
        "CFBundleShortVersionString": "1.1.1",
        "CFBundleVersion": "1.1.1",
    },
)
