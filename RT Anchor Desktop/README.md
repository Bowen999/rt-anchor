# RT Anchor — desktop app

A native desktop app (pywebview + the `rt_anchor` engine) for **cross‑column
retention‑time calibration**: it puts a sample's retention times onto a stated
**reference column's** time axis (`Cal_RT_min`) and onto a dimensionless 1–100
**iRT** index. Light "scientific" UI (warm paper, slate‑blue accents, Inter).
macOS build (`.app` / `.dmg`, arm64) and Windows build (`RT Anchor.exe`,
onedir zip) ship from the same shared source.

> **v2 — the method changed.** Earlier versions warped a sample onto iRT using
> the standard panel found *in that same sample*. The engine now calibrates one
> column *onto another*: your standards run tells it how your column's time axis
> maps onto the reference column's, and the sample inherits that map. The
> outputs, the sections and the advanced inputs all changed accordingly — see
> [`docs/RI_CALIBRATION_SPEC_V2.md`](../docs/RI_CALIBRATION_SPEC_V2.md).

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
- **Matching & fitting** — the two m/z tolerances, both **15 ppm** by default
  (the *feature‑match window* for anonymous feature matching, and the *panel
  m/z tolerance* for identifying named panel standards), the LOESS fraction
  `curve_frac` (CV‑chosen — rarely a reason to change it), the RT window, min
  anchors, the stage toggles (**use sample pairs in the curve** — on;
  **series term** — the stage‑1b homologous‑series correction, on by default;
  **lattice term** — the stage‑1c cross‑ladder fallback behind it, on by
  default since 1.2.3; **stage‑2 sample anchors (opt‑in)** — off by default
  since 1.2.2), and extrapolation (on by default; out‑of‑span features are still
  valued and flagged). A hint under the toggles says what the series term does:
  it corrects a feature with the curve residuals of its Kendrick homologues —
  recognised from m/z alone, not identifications — and only where a
  leave‑own‑out gate shows it helps. One more sentence says what the lattice
  term does: it joins homologous series one double bond apart into a class‑wide
  family and corrects features the series term cannot reach, gated the same way.

### Output screen
A left sidebar switches sections:

| # | Section | Content |
|---|---|---|
| 01 | **Overview** | KPI metrics (features per run, curve residual median / P90, iRT range and extrapolated share) and the quality radar |
| 02 | **Detection** | the chosen panel located in **your standards run**, on the reference column's axis — **QC only, it does not drive the calibration** |
| 03 | **Profile** | intensity‑weighted feature profile: before vs. after calibration (raw RT vs. `Cal_RT`), then the calibrated profile against the reference |
| 04 | **Curve** | the stage‑1 cross‑column curve: matched pairs (standards = squares, sample = diamonds, MAD‑trimmed = hollow grey), the fit, the anchor‑refined fit when engaged, the stage‑2 anchors as rings, and the residual strip — with the series‑ and lattice‑term decisions (and the anchor gate's, when stage 2 ran) printed under the figure; the series and lattice terms, when engaged, are per‑feature corrections that the curve does not show |
| 05 | **Table** | preview: your RT and m/z beside `Cal_RT_min`, `iRT`, their uncertainties, reliability, the extrapolation flag, the series‑term columns (`series corr (min)`, `series n`) and the lattice‑term columns (`lattice corr (min)`, `lattice n`) |
| 06 | **Export** | see below |

The stage‑2 anchor table (residuals, class offsets, gate verdict) is not a
section of its own: it is in the exported `_anchors.csv`, `_model.json` and the
report. The app has no per‑injection input, so the repeatability view of the
engine's per‑sample tier is not shown here either.

Detection and Profile are the engine's own Plotly figures; the radar and curve
charts are hand‑authored SVG (`web/charts.js`). The
numbers behind all of them come from the same engine helpers the exported
report uses (`viz.metrics`, `viz.performance.compute_curve`,
`viz.anchors.compute_anchors`, `viz.report.method_facts`), so the app and the
PDF a reviewer receives cannot disagree.

The bundled example dataset (**Load example**) is a run on which the series gate
declines and the lattice gate engages: 12 lattice families holding 319 pairs,
leave‑own‑out MSE 24% below the curve alone on 321 pairs, 5,035 of 22,559
features corrected.

### Export
**Export all outputs** writes one folder:

`_calibrated.csv` · `_model.json` · `_anchors.csv` (stage‑2) ·
**`_pairs.csv`** (stage‑1 matched pairs) · **`_landmarks.csv`** (iRT landmarks
on the reference run) · `_log.txt` · `_report.html` · `_run_info.json`

The last two CSVs are new in v2 — they are the curve's raw evidence and the
ruler's endpoints — and the app lists every file it wrote by name afterwards.
Individual exports (calibrated CSV, HTML report, run info) are also available.
`_run_info.json` records the run's parameters and, since 1.2.2, a
`results.series` block (the series term's gate decision and counts:
`enabled`, `engaged`, `gate_mse_reduction`, `gate_threshold`, `gate_reason`,
`n_series`, `n_pairs_covered`, `n_features_corrected`) and, since 1.2.3, a
`results.lattice` block beside it (the lattice term's: `enabled`, `engaged`,
`gate_mse_reduction`, `gate_threshold`, `gate_reason`, `n_families`,
`n_pairs_covered`, `n_features_corrected`).

## Run from source (dev)
```bash
/opt/anaconda3/bin/python "app.py"        # macOS, needs rt_anchor installed
/opt/anaconda3/bin/python "app.py" --selftest   # headless end-to-end check
..\rt_anchor\build_venv\Scripts\python app.py   # Windows (the build venv)
```
`--selftest` runs a **real calibration** — the bundled reference pair against
itself — and asserts `Cal_RT_min == RT` to machine precision, then builds the
full UI payload and checks it is valid JSON. It also asserts that the result's
model carries a `series` and a `lattice` block, each with an `engaged` key, and
prints both gate reasons. Importing the dependencies is not enough:
`statsmodels`, `sklearn` and the bundled `reference_data` tree are all things
PyInstaller drops silently, and only a run that produces the right number
proves they survived.

```
SELFTEST_OK reference=col35 features=610 pairs=825 landmarks=21 self_error=3.92e-11min structures=15 series=no correction needed: the stage-1 curve already reproduces the series pairs lattice=no correction needed: the stage-1 curve already reproduces the lattice pairs
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
One command. On first use it creates a pip‑only build venv at
`../rt_anchor/build_venv` (rt-anchor, pywebview, PyInstaller, and rdkit when it
installs); delete that folder to rebuild it after an engine change:
```bat
build\windows\build_win.bat     :: -> dist\win\"RT Anchor.exe" + RT Anchor-win64.zip
```
The release CI (`.github/workflows/build-windows.yml`) runs the same script on
every `v*` tag, self‑tests the frozen exe, and attaches the zip to the GitHub
Release as `RT-Anchor-<version>-Windows-x64.zip`.
Smoke-test a windowed (no-console) build by exit code:
```bat
"dist\win\RT Anchor\RT Anchor.exe" --selftest && echo OK
```
The Windows spec is v2‑current: it bundles `rt_anchor/reference_data/**` from
the sibling engine checkout (failing loudly if absent) and pulls in
`statsmodels` + `sklearn` (the stage‑1 LOWESS → isotonic chain) as
hiddenimports + data files — the same two changes the macOS spec received.
Startup feedback on slow machines is the **HTML boot overlay**
(`web/index.html`): it paints with the first frame and narrates the bridge +
engine‑import wait. There is intentionally no PyInstaller `Splash()` — its
Tcl/Tk payload fails to collect from a Microsoft Store Python build host and
the frozen app then pops “failed to load Tcl DLL” dialogs on every launch.

Notes: the Windows spec targets pywebview's WinForms/WebView2 backend via
pythonnet (.NET Framework); `app_win.py` is the frozen entry — it silences
stdout/stderr, pins MPLBACKEND=Agg with a persistent font-cache dir under
`%LOCALAPPDATA%\RTAnchor\mpl`, forces `PYTHONNET_RUNTIME=netfx`, and strips
Mark-of-the-Web from bundled binaries (zip downloads get tagged).
`build/windows/RT Anchor.exe.config` (`loadFromRemoteSources`) is the fallback
for read-only install locations. rdkit is optional — without it everything works
except molecular-structure hover.
