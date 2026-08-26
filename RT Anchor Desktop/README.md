# RT Anchor — desktop app

A native desktop app (pywebview + the `rt_anchor` engine) for standard‑panel
retention‑index (iRT) calibration. Dark, blue, refined‑brutalist UI.
macOS build (`.app` / `.dmg`, arm64) and Windows build (`RT Anchor.exe`,
onedir zip) ship from the same shared source.

## Deliverables
- macOS: `dist/macos/RT Anchor.dmg` — drag‑to‑Applications installer
  (ad‑hoc signed) · `dist/macos/RT Anchor.app` — the app bundle.
- Windows: `dist/win/RT Anchor-win64.zip` → unzip anywhere, run
  `RT Anchor.exe` (needs the WebView2 Runtime — preinstalled on Win 11 /
  most Win 10).

## What it does
- **Input screen** — the **Required** block: a sample feature table (MS‑DIAL /
  MZmine / MassCube / LipidScreener, auto‑detected), a **standards run** (builds
  the native anchor template), and polarity. A collapsed **Advanced** block holds
  optional per‑injection files and the matching parameters (m/z tolerance, RT
  window, min anchors, extrapolation) — there is no "instrument type" to choose.
  Then **Run calibration**.
- **Output screen** — a left sidebar switches sections: **Overview** (KPI
  scorecard + quality radar), **Detection**, **Profile**, **Warp**,
  **Repeatability** (per‑sample only), **Export** (calibrated CSV / HTML+PDF report).
  For the built‑in panel, hovering a standard shows its molecular structure.

## Run from source (dev)
```bash
/opt/anaconda3/bin/python "app.py"        # macOS, needs rt_anchor installed
build_venv python app.py                  # Windows (rt_anchor build venv)
```

## Layout
Shared source + per-platform build folders:
```
app.py                shared entry (pywebview window + Api bridge)
app_win.py            Windows frozen entry (stdout/MPL/pythonnet guards, MOTW strip)
rtad/api.py           JS-exposed API (file dialogs, run_calibration, export)
rtad/appviz.py        native chart-data extraction (no rt_anchor.viz / Plotly)
web/                  index.html · styles.css · app.js · charts.js (hand-authored SVG)
example/              bundled demo dataset ("Load example" button)

build/
  macos/              RTAnchor.spec · icon.icns · build_mac.sh   (macOS arm64)
  windows/            RTAnchor.win.spec · icon.ico · build_win.bat ·
                      RT Anchor.exe.config                       (Windows x64)

dist/
  macos/              RT Anchor.app · RT Anchor.dmg
  win/                "RT Anchor"/RT Anchor.exe · RT Anchor-win64.zip
```

## Rebuild — macOS (.app / .dmg)
One script (PyInstaller runs to `/private/tmp`, avoiding the iCloud‑Desktop
codesign/TCC issues; artefacts land in `dist/macos/`):
```bash
build/macos/build_mac.sh
# or manually:
pyinstaller build/macos/RTAnchor.spec --distpath /private/tmp/rtad_dist --workpath /private/tmp/rtad_build --noconfirm
codesign --force --deep --sign - "/private/tmp/rtad_dist/RT Anchor.app"
# stage + hdiutil create ... -format UDZO  -> dist/macos/RT Anchor.dmg
```
Notes: `pathlib` backport must be uninstalled (breaks PyInstaller); the spec
excludes the heavy unused anaconda libs (tensorflow/torch/onnx/holoviz/…) to keep
the bundle ~650 MB / dmg ~270 MB. arm64‑only build.

## Rebuild — Windows (exe / zip)
One command (uses the shared `../rt_anchor/build_venv`, where rt-anchor,
pywebview, PyInstaller and rdkit are pip‑installed — see `rt_anchor/build.bat`):
```bat
build\windows\build_win.bat     :: -> dist\win\"RT Anchor.exe" + RT Anchor-win64.zip
```
Smoke-test a windowed (no-console) build by exit code:
```bat
"dist\win\RT Anchor\RT Anchor.exe" --selftest && echo OK
```
Notes: the Windows spec (`build/windows/RTAnchor.win.spec`) targets pywebview's
WinForms/WebView2 backend via pythonnet (.NET Framework); `app_win.py` is the
frozen entry — it silences stdout/stderr, pins MPLBACKEND=Agg with a persistent
font-cache dir under `%LOCALAPPDATA%\RTAnchor\mpl`, forces `PYTHONNET_RUNTIME=netfx`,
and strips Mark-of-the-Web from bundled binaries (zip downloads get tagged).
`build/windows/RT Anchor.exe.config` (`loadFromRemoteSources`) is the fallback
for read-only install locations. Current bundle: engine **rt-anchor 1.1.1**
(install it into `../rt_anchor/build_venv` first; the UI reports the engine
version at runtime). ~260 MB / zip ~116 MB; rdkit is optional — without it
everything works except molecular-structure hover.
