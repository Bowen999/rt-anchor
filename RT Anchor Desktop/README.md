# RT Anchor — desktop app

A native desktop app (pywebview + the `rt_anchor` engine) for **cross‑column
retention‑time calibration**: it puts a sample's retention times onto a stated
**reference column's** time axis (`Cal_RT_min`) and onto a dimensionless 1–100
**iRT** index. Dark, blue, refined‑brutalist UI. macOS build (`.app` / `.dmg`,
arm64) and Windows build (`RT Anchor.exe`, onedir zip) ship from the same
shared source.

> **v2 — the method changed.** Earlier versions warped a sample onto iRT using
> the standard panel found *in that same sample*. The engine now calibrates one
> column *onto another*: your standards run tells it how your column's time axis
> maps onto the reference column's, and the sample inherits that map. The
> outputs, the sections and the advanced inputs all changed accordingly — see
> `RI_CALIBRATION_SPEC_V2.md`.

## Deliverables
- macOS: `dist/macos/RT Anchor.dmg` — drag‑to‑Applications installer
  (ad‑hoc signed) · `dist/macos/RT Anchor.app` — the app bundle.
- Windows: `dist/win/RT Anchor-win64.zip` → unzip anywhere, run
  `RT Anchor.exe` (needs the WebView2 Runtime — preinstalled on Win 11 /
  most Win 10).

## What it does

### Input screen
**Required**

| Input | Notes |
|---|---|
| **Sample feature table** | MS‑DIAL / MZmine / MassCube / LipidScreener, auto‑detected |
| **Standards mixture table** | **must come from the same column and gradient as the sample** — this is the assumption the whole method rests on |
| **Standards panel** | **15** (Caley lipid RT‑calibration mix) · **21** (Mix 4.4 extended panel) · **Others** (unknown / no panel) |
| **Polarity** | positive / negative |

The panel names the iRT landmarks and the detection QC. It does **not** drive
the calibration: stage 1 matches features by m/z alone and claims no
identities, which is why **Others** still calibrates normally (the iRT ruler
then falls back to the 21‑standard landmarks on the reference standards run).

The panel cards come from `rt_anchor.mixtures`, not from a copy in the app —
`rtad/mixtures.py` is a thin re‑export so the two can never disagree about what
"Mix 21" is.

**Advanced** (collapsed)

- **Reference column** — *Reference sample* and *Reference standards* pickers,
  both showing “bundled — Column 25” until overridden, plus
  **Reset to bundled reference**. Both runs must be supplied together; a
  half‑set reference is refused, because a run calibrated against one lab's
  serum and another lab's standards produces plausible‑looking nonsense.
  Substituting your own pair moves the time axis, so those results are not
  comparable with default‑reference runs.
- **Matching & fitting** — the two m/z tolerances (an absolute **Da** window for
  anonymous feature matching, default 0.008; a **ppm** window for identifying
  named panel standards), the LOESS fraction `curve_frac` (CV‑chosen — rarely a
  reason to change it), the RT window, min anchors, both stage toggles
  (**use sample pairs in the curve**, **stage‑2 sample anchors**), and
  extrapolation (on by default; out‑of‑span features are still valued and
  flagged).

### Output screen
A left sidebar switches sections:

| # | Section | Content |
|---|---|---|
| 01 | **Overview** | method/provenance band (reference column, panel, pair counts, gate verdict, iRT definition), KPI metrics, quality radar |
| 02 | **Detection** | the chosen panel located in **your standards run**, on the reference column's axis — **QC only, it does not drive the calibration** |
| 03 | **Profile** | intensity‑weighted feature profile, raw RT vs. iRT |
| 04 | **Curve** | the stage‑1 cross‑column curve: matched pairs (standards = squares, sample = diamonds, MAD‑trimmed = hollow grey), the fit, the anchor‑refined fit when engaged, the stage‑2 anchors as rings, and the residual strip |
| 05 | **Anchors** | the stage‑2 endogenous‑plasma‑lipid anchors: residual against the stage‑1 curve vs. the leave‑one‑out residual under the correction, coloured by lipid class, with the gate's verdict and the class median offsets. Disabled when no anchors validated |
| 06 | **Repeatability** | per‑injection only |
| 07 | **Table** | preview: your RT and m/z beside `Cal_RT_min`, `iRT`, their uncertainties, reliability and the extrapolation flag |
| 08 | **Export** | see below |

Detection and Profile are the engine's own Plotly figures; the radar, curve,
anchors and repeatability charts are hand‑authored SVG (`web/charts.js`). The
numbers behind all of them come from the same engine helpers the exported
report uses (`viz.metrics`, `viz.performance.compute_curve`,
`viz.anchors.compute_anchors`, `viz.report.method_facts`), so the app and the
PDF a reviewer receives cannot disagree.

### Export
**Export all outputs** writes one folder:

`_calibrated.csv` · `_model.json` · `_anchors.csv` (stage‑2) ·
**`_pairs.csv`** (stage‑1 matched pairs) · **`_landmarks.csv`** (iRT landmarks
on the reference run) · `_log.txt` · `_report.html` · `_run_info.json`

The last two CSVs are new in v2 — they are the curve's raw evidence and the
ruler's endpoints — and the app lists every file it wrote by name afterwards.
Individual exports (calibrated CSV, HTML report, run info) are also available.

## Run from source (dev)
```bash
/opt/anaconda3/bin/python "app.py"        # macOS, needs rt_anchor installed
/opt/anaconda3/bin/python "app.py" --selftest   # headless end-to-end check
build_venv python app.py                  # Windows (rt_anchor build venv)
```
`--selftest` runs a **real calibration** — the bundled reference pair against
itself — and asserts `Cal_RT_min == RT` to machine precision, then builds the
full UI payload and checks it is valid JSON. Importing the dependencies is not
enough: `statsmodels`, `sklearn` and the bundled `reference_data` tree are all
things PyInstaller drops silently, and only a run that produces the right
number proves they survived.

```
SELFTEST_OK reference=col35 features=610 pairs=827 landmarks=21 self_error=1.09e-11min structures=15
```

## Layout
Shared source + per-platform build folders:
```
app.py                shared entry (pywebview window + Api bridge + --selftest)
app_win.py            Windows frozen entry (stdout/MPL/pythonnet guards, MOTW strip)
rtad/api.py           JS-exposed API (file dialogs, run_calibration, references, export)
rtad/appviz.py        chart-data extraction for the app's own SVG charts
rtad/mixtures.py      thin re-export of rt_anchor.mixtures (panels live in the engine)
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
"/private/tmp/rtad_dist/RT Anchor.app/Contents/MacOS/RT Anchor" --selftest   # must print SELFTEST_OK
```
Notes: `pathlib` backport must be uninstalled (breaks PyInstaller); the spec
excludes the heavy unused anaconda libs (tensorflow/torch/onnx/holoviz/…).
**`sklearn` and `statsmodels` are no longer excludable** — the stage‑1 curve is
statsmodels LOWESS → sklearn `IsotonicRegression` → scipy `PchipInterpolator` —
and `rt_anchor/reference_data/**` must be listed as data files, or the frozen
app fails every calibration with *“Bundled reference 'col35' is incomplete”*.
Both are handled in `build/macos/RTAnchor.spec`; the `--selftest` above is what
proves it. arm64‑only build.

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
> ⚠️ **`build/windows/RTAnchor.win.spec` has not been updated for v2.** It still
> excludes `sklearn` and `statsmodels` and does not bundle
> `rt_anchor/reference_data/**`, so a Windows build made from it will start and
> then fail every calibration. It needs the same two changes the macOS spec
> received. This port was scoped to macOS source only, so the file was left
> untouched deliberately rather than edited blind.

Notes: the Windows spec targets pywebview's WinForms/WebView2 backend via
pythonnet (.NET Framework); `app_win.py` is the frozen entry — it silences
stdout/stderr, pins MPLBACKEND=Agg with a persistent font-cache dir under
`%LOCALAPPDATA%\RTAnchor\mpl`, forces `PYTHONNET_RUNTIME=netfx`, and strips
Mark-of-the-Web from bundled binaries (zip downloads get tagged).
`build/windows/RT Anchor.exe.config` (`loadFromRemoteSources`) is the fallback
for read-only install locations. rdkit is optional — without it everything works
except molecular-structure hover.
