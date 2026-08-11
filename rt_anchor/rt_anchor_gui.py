"""PyInstaller entry launcher for the rt-anchor Windows GUI (onefile, windowed).

These runtime guards MUST execute before matplotlib is imported and before the
Tkinter app starts, so they live at module top. The GUI itself is in
``rt_anchor.gui``; this file only sets up the frozen-app environment and calls it.
"""

import os
import sys

# 1. A windowed build (console=False) has sys.stdout / sys.stderr == None; any
#    print() or warning would then raise AttributeError('NoneType' has no 'write').
if sys.stdout is None:
    sys.stdout = open(os.devnull, "w")
if sys.stderr is None:
    sys.stderr = open(os.devnull, "w")

# 2. Pin matplotlib to its non-interactive backend and give it a PERSISTENT
#    config/cache dir. onefile extracts to a fresh %TEMP%\_MEIxxxx every launch,
#    so without this matplotlib rebuilds its font cache on each cold start.
os.environ.setdefault("MPLBACKEND", "Agg")
_mpl_dir = os.path.join(os.environ.get("LOCALAPPDATA", os.path.expanduser("~")),
                        "rt-anchor", "mpl")
try:
    os.makedirs(_mpl_dir, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", _mpl_dir)
except OSError:
    pass

# 3. Force pythonnet onto the system .NET Framework (always present on Win 10/11),
#    so the pywebview WebView2 backend never depends on a .NET Core runtime.
os.environ.setdefault("PYTHONNET_RUNTIME", "netfx")


def _strip_motw() -> None:
    """Remove the 'downloaded-from-the-internet' mark (Mark-of-the-Web) from every
    bundled binary.

    When a user downloads the .zip and extracts it on NTFS, Windows tags each file
    with a Zone.Identifier (remote zone). .NET Framework then refuses to load the
    bundled managed assemblies — pythonnet's ``Python.Runtime.dll`` AND WebView2's
    ``Microsoft.Web.WebView2.*.dll`` — with a NotSupportedException / "Failed to
    resolve ... Loader.Initialize". Deleting that stream (exactly what a file's
    "Unblock" checkbox does) fixes it. We strip every .dll/.exe/.pyd under the
    bundle so no .NET assembly is left blocked. Best-effort; frozen build only.
    """
    base = getattr(sys, "_MEIPASS", None)
    if not base:
        return
    for root, _dirs, files in os.walk(base):
        for name in files:
            if name.lower().endswith((".dll", ".exe", ".pyd")):
                try:
                    os.remove(os.path.join(root, name) + ":Zone.Identifier")
                except OSError:
                    pass


def _run() -> int:
    _strip_motw()
    from rt_anchor.gui import main
    return main()


if __name__ == "__main__":
    import multiprocessing

    # onefile + windowed convention: without this, a worker process spawned by
    # numpy/scipy/joblib would re-launch the whole GUI. Current code uses no
    # multiprocessing, so this is purely defensive.
    multiprocessing.freeze_support()
    raise SystemExit(_run())
