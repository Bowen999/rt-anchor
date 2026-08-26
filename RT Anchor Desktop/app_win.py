"""RT Anchor desktop — Windows frozen entry (PyInstaller onedir, windowed).

These guards MUST run before matplotlib is imported and before the pywebview
window starts, so they live at module top; the actual app is `app.main()`.

  1. A windowed build has sys.stdout / sys.stderr == None; any print() or
     warning would then raise AttributeError('NoneType' has no 'write').
  2. Pin matplotlib to its non-interactive backend and give it a PERSISTENT
     config/cache dir under %LOCALAPPDATA%, so the bundled DejaVu fonts are
     not re-indexed into a fresh temp dir on every launch.
  3. Force pythonnet onto the system .NET Framework (always present on
     Win 10/11), so pywebview's WebView2 backend never depends on a
     .NET Core runtime.
  4. Strip the Mark-of-the-Web (Zone.Identifier) from every bundled binary —
     after a zip download+extract, .NET refuses to load "remote-zone"
     assemblies (Python.Runtime.dll, WebView2Loader.dll) with a
     NotSupportedException. Best-effort, frozen builds only.
"""

import os
import sys

if sys.stdout is None:
    sys.stdout = open(os.devnull, "w")
if sys.stderr is None:
    sys.stderr = open(os.devnull, "w")

os.environ.setdefault("MPLBACKEND", "Agg")
_mpl_dir = os.path.join(os.environ.get("LOCALAPPDATA", os.path.expanduser("~")),
                        "RTAnchor", "mpl")
try:
    os.makedirs(_mpl_dir, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", _mpl_dir)
except OSError:
    pass

os.environ.setdefault("PYTHONNET_RUNTIME", "netfx")


def _strip_motw() -> None:
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
    from app import main
    main()
    return 0


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    raise SystemExit(_run())
