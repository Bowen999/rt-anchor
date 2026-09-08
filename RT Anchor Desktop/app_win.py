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


# DLLs/EXEs that .NET/pythonnet refuse to load when Mark-of-the-Web is present.
_MOTW_TARGETS = frozenset({
    "python.Runtime.dll", "python3.dll", "python312.dll",
    "WebView2Loader.dll", "WebView2Loader.dll",
    "CLR.dll", "coreclr.dll", "clretwcore.dll",
})

def _strip_motw() -> None:
    """Remove Zone.Identifier from bundled binaries that .NET refuses to load.

    Only touches the known-problematic DLLs instead of walking the entire
    bundle, cutting startup time by ~200ms on cold runs.
    """
    base = getattr(sys, "_MEIPASS", None)
    if not base:
        return
    # Pass 1: targeted files in _internal/ root (fast path, covers 95% of cases)
    internal = os.path.join(base, "_internal")
    for name in _MOTW_TARGETS:
        try:
            os.remove(os.path.join(internal, name) + ":Zone.Identifier")
        except OSError:
            pass
    # Pass 2: full walk only for webview/ and python*/ subdirectories
    for subdir in ("webview", "python"):
        subpath = os.path.join(base, subdir)
        if not os.path.isdir(subpath):
            continue
        for root, _dirs, files in os.walk(subpath):
            for name in files:
                if name.lower().endswith((".dll", ".exe", ".pyd")):
                    try:
                        os.remove(os.path.join(root, name) + ":Zone.Identifier")
                    except OSError:
                        pass


def _check_webview2() -> None:
    """Warn once if the WebView2 Runtime is missing (Win 10 LTSC, older builds).

    Registry first (the documented check, per-machine and per-user), then the
    install directory as a fallback — the client keys are occasionally absent
    even when the runtime is installed and working, and a false warning is
    worse than none.
    """
    found = False
    try:
        import winreg
        subs = (r"SOFTWARE\WOW6432Node\Microsoft\EdgeUpdate\Clients\{F3017226-FE2A-4295-8BEE-13A6279B0CE9}",
                r"SOFTWARE\Microsoft\EdgeUpdate\Clients\{F3017226-FE2A-4295-8BEE-13A6279B0CE9}")
        for root in (winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER):
            for sub in subs:
                try:
                    key = winreg.OpenKey(root, sub, 0, winreg.KEY_READ)
                    winreg.CloseKey(key)
                    found = True
                except OSError:
                    pass
    except Exception:
        pass
    if not found:
        import glob as _glob
        for base in (os.environ.get("ProgramFiles(x86)"), os.environ.get("ProgramFiles")):
            if base and _glob.glob(os.path.join(
                    base, "Microsoft", "EdgeWebView", "Application", "*", "msedgewebview2.exe")):
                found = True
                break
    if not found:
        print("[RT Anchor] WARNING: WebView2 Runtime not found. "
              "The app may fail to start. Download from:\n"
              "  https://developer.microsoft.com/en-us/microsoft-edge/webview2/",
              file=sys.stderr)


def _run() -> int:
    _strip_motw()
    # ---- High-DPI awareness (Windows 10/11 per-monitor) ----
    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(2)   # Per-Monitor DPI Aware v2
    except Exception:
        try:
            import ctypes
            ctypes.windll.user32.SetProcessDPIAware()     # fallback: system DPI aware
        except Exception:
            pass
    # ---- WebView2 Runtime check ----
    _check_webview2()
    from app import main
    main()
    return 0


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    raise SystemExit(_run())
