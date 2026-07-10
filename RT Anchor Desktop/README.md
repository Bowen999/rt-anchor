# RT Anchor — desktop app

A native macOS app (pywebview + the `rt_anchor` engine) for standard‑panel
retention‑index (iRT) calibration. Dark, blue, refined‑brutalist UI.

## Deliverables
- `dist/RT Anchor.dmg` — drag‑to‑Applications installer (ad‑hoc signed).
- `dist/RT Anchor.app` — the app bundle.

Unsigned/ad‑hoc → on first open, right‑click the app → **Open** (or
`xattr -cr "/Applications/RT Anchor.app"`) to clear Gatekeeper quarantine.

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
/opt/anaconda3/bin/python "app.py"        # needs rt_anchor installed
```

## Rebuild the app / dmg
The build runs to `/private/tmp` (avoids the iCloud‑Desktop codesign/TCC issues):
```bash
cd "RT Anchor Desktop"
pyinstaller build/RTAnchor.spec --distpath /private/tmp/rtad_dist --workpath /private/tmp/rtad_build --noconfirm
codesign --force --deep --sign - "/private/tmp/rtad_dist/RT Anchor.app"
# stage + hdiutil create ... -format UDZO  -> RT Anchor.dmg
```
Notes: `pathlib` backport must be uninstalled (breaks PyInstaller); the spec
excludes the heavy unused anaconda libs (tensorflow/torch/onnx/holoviz/…) to keep
the bundle ~650 MB / dmg ~270 MB. arm64‑only build.

## Layout
```
app.py                entry (pywebview window + Api bridge)
rtad/api.py           JS-exposed API (file dialogs, run_calibration, export)
rtad/appviz.py        native chart-data extraction (no rt_anchor.viz / Plotly)
web/                  index.html · styles.css · app.js · charts.js (hand-authored SVG)
build/RTAnchor.spec   PyInstaller spec · icon.icns · example/ (bundled demo)
```
