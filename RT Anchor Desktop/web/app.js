/* RT Anchor desktop — front-end logic (pywebview bridge + browser preview) */
"use strict";

const REF_PLACEHOLDER = "bundled Column 25";

const state = {
  params: { samples: "", standards: "", polarity: "positive", mixture: "mix21",
            reference_sample: "", reference_standards: "" },
  bundle: null,
  section: "overview",
};

const SECTIONS = [
  { id: "overview", label: "Overview", idx: "01", title: "Overview" },
  { id: "detection", label: "Detection", idx: "02", title: "Standard detection (QC)" },
  { id: "profile", label: "Profile", idx: "03", title: "Feature-intensity profile" },
  { id: "curve", label: "Curve", idx: "04", title: "Cross-column calibration curve" },
  { id: "table", label: "Table", idx: "05", title: "Calibrated table" },
  { id: "export", label: "Export", idx: "06", title: "Export" },
];

const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));

/* ---- API bridge (real pywebview, or a preview stub in a browser) ---- */
function api() {
  if (window.pywebview && window.pywebview.api) return window.pywebview.api;
  return PREVIEW_API;
}
const PREVIEW_API = {
  pick_file: async () => "/preview/example.txt",
  pick_folder: async () => "/preview/single_files",
  load_example: async () => ({ ok: true, samples: "/preview/example/samples.txt",
    standards: "/preview/example/standards.txt", polarity: "positive", mixture: "mix15" }),
  reference_info: async () => ({ ok: true, default_key: "col35", placeholder: REF_PLACEHOLDER,
    references: [{ key: "col35", label: "Column 25",
      sub: "Human serum + Mix 4.4 standards, QTOF positive",
      sample: "reference_sample.csv", standards: "reference_standards.csv", available: true }] }),
  reset_reference: async () => ({ ok: true, reference_sample: "", reference_standards: "",
    placeholder: REF_PLACEHOLDER }),
  mixture_previews: async () => ({ ok: true, default: "mix21", mixtures: [
    { key: "mix15", label: "15", sub: "Li lab 15 lipids mixture", n_standards: 15,
      rows: [{ name: "PC 6:0_6:0", class: "PC", mz_pos: 454.2569, mz_neg: 498.2474, rt_ref_min: 1.06 },
             { name: "PC 12:0_12:0", class: "PC", mz_pos: 622.4448, mz_neg: 666.4362, rt_ref_min: 3.6 },
             { name: "PC 14:0_14:0", class: "PC", mz_pos: 678.5069, mz_neg: 722.4983, rt_ref_min: 5.63 },
             { name: "PC 16:0_18:1", class: "PC", mz_pos: 760.5856, mz_neg: 804.5766, rt_ref_min: 8.4 },
             { name: "DG 16:0_18:0", class: "DG", mz_pos: 614.5723, mz_neg: 641.5367, rt_ref_min: 12.65 },
             { name: "CE 18:1", class: "CE", mz_pos: 668.6345, mz_neg: null, rt_ref_min: 18.39 }] },
    { key: "mix21", label: "21", sub: "Li lab 21 lipids mixture", n_standards: 21,
      rows: [{ name: "PC 6:0_6:0", class: "PC", mz_pos: 454.2561, mz_neg: 498.2474, rt_ref_min: 1.5 },
             { name: "PC 8:0_8:0", class: "PC", mz_pos: 510.319, mz_neg: 554.31, rt_ref_min: 3.2 },
             { name: "PC 12:0_12:0", class: "PC", mz_pos: 622.444, mz_neg: 666.4352, rt_ref_min: 6.8 },
             { name: "PC 14:0_16:0", class: "PC", mz_pos: 706.5381, mz_neg: 750.5291, rt_ref_min: 11.1 },
             { name: "PC 16:0_18:1", class: "PC", mz_pos: 760.5846, mz_neg: 804.576, rt_ref_min: 12.9 },
             { name: "DG 18:0_18:0", class: "DG", mz_pos: 642.6031, mz_neg: 669.5675, rt_ref_min: 18.0 },
             { name: "CE 20:5", class: "CE", mz_pos: 688.6027, mz_neg: null, rt_ref_min: 20.5 },
             { name: "TG 18:0_18:2_18:0", class: "TG", mz_pos: 904.8328, mz_neg: null, rt_ref_min: 21.8 }] },
    { key: "none", label: "Others", sub: "other mixtures", n_standards: 0, rows: [] },
  ] }),
  run_calibration: async () => window.__PREVIEW__ || { ok: false, error: "no preview data" },
  export_csv: async () => ({ ok: true, path: "(preview)/calibrated.csv" }),
  export_report: async () => ({ ok: true, paths: {} }),
  export_run_info: async () => ({ ok: true, path: "(preview)/run_info.json" }),
  export_all: async () => ({ ok: true, dir: "(preview)/outputs", n_files: 9, report: true,
    files: [{ kind: "calibrated_csv", label: "calibrated table", name: "rt_anchor_run_calibrated.csv" },
            { kind: "model_json", label: "model.json", name: "rt_anchor_run_model.json" },
            { kind: "anchors_csv", label: "stage-2 anchors", name: "rt_anchor_run_anchors.csv" },
            { kind: "pairs_csv", label: "stage-1 matched pairs", name: "rt_anchor_run_pairs.csv" },
            { kind: "landmarks_csv", label: "iRT landmarks", name: "rt_anchor_run_landmarks.csv" },
            { kind: "log_txt", label: "run log", name: "rt_anchor_run_log.txt" },
            { kind: "report_html", label: "HTML report", name: "rt_anchor_run_report.html" },
            { kind: "run_info_json", label: "run info", name: "rt_anchor_run_run_info.json" }] }),
  open_github: async () => ({ ok: true }),
};

/* ---- toast ---- */
let toastT;
function toast(msg, err) {
  const t = $("#toast"); t.textContent = msg; t.className = "toast show" + (err ? " err" : "");
  clearTimeout(toastT); toastT = setTimeout(() => (t.className = "toast"), 2600);
}
function hideToast() { const t = $("#toast"); clearTimeout(toastT); t.className = "toast"; }

/* ===================== INPUT VIEW ===================== */
/* last two path components, on either separator. Windows hands us
   "C:\Users\x\data\samples.txt" — splitting on "/" alone left the whole path
   in a field sized for two components, so it overflowed its box. */
function shortPath(p) {
  if (!p) return "";
  const parts = String(p).split(/[\\/]+/).filter(Boolean);
  return parts.slice(-2).join("/");
}

/* file/folder fields, and what their empty state says */
const FILE_FIELDS = [
  ["samples", "no file selected"],
  ["standards", "no file selected"],
  ["reference_sample", REF_PLACEHOLDER],
  ["reference_standards", REF_PLACEHOLDER],
];
const PLACEHOLDER = Object.fromEntries(FILE_FIELDS);

/* set a file/folder field + flip it to the "filled" state */
function markField(target, path) {
  state.params[target] = path;
  const el = $("#path-" + target);
  if (el) { el.textContent = shortPath(path); el.title = path; el.classList.add("set"); }
  const field = document.querySelector(`.field[data-key="${target}"]`);
  if (field) field.classList.add("filled");
}
function clearField(target) {
  state.params[target] = "";
  const el = $("#path-" + target);
  if (el) { el.textContent = PLACEHOLDER[target] || ""; el.title = ""; el.classList.remove("set"); }
  const field = document.querySelector(`.field[data-key="${target}"]`);
  if (field) field.classList.remove("filled");
}

$$(".btn-file").forEach(btn => btn.addEventListener("click", async () => {
  const target = btn.dataset.target, kind = btn.dataset.pick;
  const path = kind === "folder" ? await api().pick_folder() : await api().pick_file();
  if (!path) return;
  markField(target, path);
  if (target.startsWith("reference_")) refHint();
}));

/* the reference pair is both-or-neither — say so the moment one is set */
function refHint() {
  const s = !!state.params.reference_sample, t = !!state.params.reference_standards;
  const el = $("#ref-hint");
  if (!el) return;
  el.classList.toggle("warn", s !== t);
  if (s !== t) {
    el.innerHTML = "<b>Both reference runs are needed.</b> A run calibrated against one lab's " +
      "serum and another lab's standards produces plausible-looking nonsense, so a half-set " +
      "reference is refused. Choose the other run, or reset to the bundled reference.";
  } else if (s && t) {
    el.innerHTML = "<b>Custom reference pair.</b> Results will be on <i>this</i> column's time " +
      "axis and are not comparable with default-reference runs. Both runs must come from the " +
      "same column and gradient.";
  } else {
    el.innerHTML = "Results are expressed on this column's time axis. Leave both on the bundled " +
      "Column 25 unless you want everything on a different column — substituting your own " +
      "pair moves the axis, so those numbers are no longer comparable with default-reference " +
      "runs. Both runs must be supplied together.";
  }
}
const refResetBtn = $("#ref-reset");
if (refResetBtn) refResetBtn.addEventListener("click", async () => {
  try { await api().reset_reference(); } catch (e) { /* local reset is the real one */ }
  clearField("reference_sample"); clearField("reference_standards");
  refHint();
  toast("Reference reset to the bundled Column 25");
});

$$(".seg").forEach(seg => seg.addEventListener("click", e => {
  const b = e.target.closest("button"); if (!b) return;
  $$("button", seg).forEach(x => x.classList.remove("active"));
  b.classList.add("active"); state.params[seg.dataset.seg] = b.dataset.val;
}));

/* load the bundled demo dataset into the Required fields */
const loadExBtn = $("#load-example");
if (loadExBtn) loadExBtn.addEventListener("click", async () => {
  const r = await api().load_example();
  if (!r || !r.ok) { showError((r && r.error) || "Could not load the example dataset."); return; }
  markField("samples", r.samples);
  markField("standards", r.standards);
  if (r.mixture) selectMixture(r.mixture);
  state.params.polarity = r.polarity || "positive";
  const seg = $('[data-seg="polarity"]');
  if (seg) $$("button", seg).forEach(b => b.classList.toggle("active", b.dataset.val === state.params.polarity));
  setStatus("Example dataset loaded — press Run calibration.", "ok");
});

/* ---- standards-panel cards + preview (either/or choice) ---- */
let MIXTURES = [];                       // payload of api().mixture_previews()
let MIX_DEFAULT = "mix21";               // backend's pre-selected panel
function renderMixCards() {
  const wrap = $("#mix-cards"); wrap.innerHTML = "";
  MIXTURES.forEach(m => {
    const b = document.createElement("button");
    b.type = "button"; b.className = "mix-card"; b.dataset.key = m.key;
    if (state.params.mixture === m.key) b.classList.add("active");
    b.innerHTML = `<span class="mc-l">${esc(m.label)}</span>` +
      `<span class="mc-s">${esc(m.sub || "")}</span>`;
    b.addEventListener("click", () => selectMixture(m.key));
    wrap.appendChild(b);
  });
}
function renderMixPreview() {
  const tbl = $("#mix-table");
  const m = MIXTURES.find(x => x.key === state.params.mixture);
  if (!m) { tbl.innerHTML = ""; return; }
  if (!m.rows || !m.rows.length) {
    tbl.innerHTML = `<tbody><tr><td class="mix-empty">` +
      `No panel identity is claimed for your mixture. The calibration is unaffected` +
      `</td></tr></tbody>`;
    return;
  }
  const thead = "<thead><tr><th>#</th><th>Standard</th><th>m/z (+)</th><th>m/z (−)</th><th>Ref RT <span class='unit'>min</span></th></tr></thead>";
  const rows = m.rows.map((r, i) => "<tr>" +
    `<td class="num">${i + 1}</td>` +
    `<td>${esc(r.name)}</td>` +
    `<td class="num">${r.mz_pos != null ? r.mz_pos.toFixed(4) : '<span class="na">—</span>'}</td>` +
    `<td class="num">${r.mz_neg != null ? r.mz_neg.toFixed(4) : '<span class="na">ND</span>'}</td>` +
    `<td class="num">${Number(r.rt_ref_min).toFixed(2)}</td></tr>`).join("");
  tbl.innerHTML = thead + `<tbody>${rows}</tbody>`;
}
function selectMixture(key) {
  state.params.mixture = key;
  $$("#mix-cards .mix-card").forEach(c => c.classList.toggle("active", c.dataset.key === key));
  renderMixPreview();
}
/* wheel-scroll fallback: pywebview's WKWebView can swallow wheel events over
   nested scroll areas — drive the preview table's scrollTop manually. */
function initMixTableScroll() {
  const wrap = $(".mix-table-wrap");
  if (!wrap) return;
  wrap.addEventListener("wheel", e => {
    if (wrap.scrollHeight <= wrap.clientHeight + 1) return;
    e.preventDefault();
    wrap.scrollTop += e.deltaY;
  }, { passive: false });
}
async function initMixtures() {
  try {
    const r = await api().mixture_previews();
    if (r && r.ok && Array.isArray(r.mixtures)) {
      MIXTURES = r.mixtures;
      if (r.default) MIX_DEFAULT = r.default;
      if (!MIXTURES.some(m => m.key === state.params.mixture)) state.params.mixture = MIX_DEFAULT;
    }
  } catch (e) { /* fall through: UI stays usable without previews */ }
  renderMixCards(); renderMixPreview();
  initMixTableScroll();
}
async function initReference() {
  try {
    const r = await api().reference_info();
    if (r && r.references && r.references.length) {
      const d = r.references.find(x => x.key === r.default_key) || r.references[0];
      const label = d.label ? `bundled — ${d.label}` : (r.placeholder || REF_PLACEHOLDER);
      ["reference_sample", "reference_standards"].forEach(k => {
        PLACEHOLDER[k] = label;
        const el = $("#path-" + k);
        if (el && !state.params[k]) el.textContent = label;
      });
    }
  } catch (e) { /* the pickers keep their static placeholder */ }
  refHint();
}

/* Advanced section collapse/expand */
const advToggle = $("#adv-toggle");
if (advToggle) advToggle.addEventListener("click", () => {
  const open = $("#adv").classList.toggle("open");
  advToggle.setAttribute("aria-expanded", open ? "true" : "false");
});

/* pull the Advanced parameters into params (blank -> package default) */
function collectAdvanced() {
  const v = id => { const el = $(id); return el ? el.value.trim() : ""; };
  const chk = (id, dflt) => { const el = $(id); return el ? el.checked : dflt; };
  state.params.match_mz_tol_ppm = v("#adv-mzppm");
  state.params.mz_tol_ppm = v("#adv-mz");
  state.params.curve_frac = v("#adv-curvefrac");
  state.params.rt_window_min = v("#adv-rtw");
  state.params.min_anchors = v("#adv-minanchors");
  state.params.use_sample_pairs = chk("#adv-samplepairs", true);
  state.params.use_sample_anchors = chk("#adv-sampleanchors", true);
  state.params.extrapolate = chk("#adv-extrapolate", true);
  const em = $("#adv-extramode");
  state.params.extrapolate_mode = em ? em.value : "linear";
}

$("#run").addEventListener("click", run);
async function run() {
  if (!state.params.samples) { setStatus("Choose a sample feature table first.", "err"); return; }
  if (!state.params.standards) { setStatus("Choose a standards mixture table — it is required.", "err"); return; }
  const s = !!state.params.reference_sample, t = !!state.params.reference_standards;
  if (s !== t) {
    $("#adv").classList.add("open");
    $("#adv-toggle").setAttribute("aria-expanded", "true");
    refHint();
    setStatus("A custom reference needs both runs — see Advanced.", "err");
    return;
  }
  collectAdvanced();
  const btn = $("#run"); btn.disabled = true;
  setStatus("", "");
  showProgress(true);
  /* warm Plotly while the engine calibrates, so it is ready the moment
     the results need it */
  window._ensurePlotly().catch(() => {});
  let res;
  try { res = await api().run_calibration(state.params); }
  catch (e) { res = { ok: false, error: String(e) }; }
  if (!res || !res.ok) {
    btn.disabled = false; showProgress(false);
    setStatus("", ""); showError((res && res.error) || "Calibration failed."); return;
  }
  /* Async: the backend returned {"ok": true, "status": "running"}.
     The UI stays responsive — __onCalibrationDone() will fire when ready. */
}
function setStatus(msg, cls) { const s = $("#status"); s.textContent = msg; s.className = "status" + (cls ? " " + cls : ""); }

/* ---- staged run progress ------------------------------------------------
   Two sources of truth, blended. The BACKEND knows the real stage boundaries
   (api.progress(): preparing / calibrating / building views) — on a slow
   machine a run can sit in the fit for a minute, and only the backend can say
   so. Where that endpoint is missing (browser preview, older backend) the
   time-based guesses below are all we have: the bar eases toward staged caps
   (always moving, never claiming to be done) and only completes when the
   result actually arrives. */
const RUN_STAGES = [
  { at: 0.0, cap: 0.10 },   // reading & validating inputs
  { at: 1.2, cap: 0.38 },   // matching features
  { at: 6.0, cap: 0.72 },   // fitting the curve
  { at: 14.0, cap: 0.95 },  // rendering results
];
/* backend stage -> front-end stage + bar cap. The backend's stage 1 ("matching
   features and fitting the curve") covers front stages 1 AND 2, so it maps to
   1; its stage 2 ("building the result views") maps to front stage 3. Caps are
   per backend stage: the bar may not pass 75% while the fit is still running,
   however long that takes. */
const BE_TO_FRONT = [0, 1, 3, 3];
const BE_CAPS = [0.10, 0.75, 0.95];
let _runTimer = null, _pollTimer = null, _beStage = null;
function showProgress(on) {
  const card = $("#runcard");
  clearInterval(_runTimer); _runTimer = null;
  clearInterval(_pollTimer); _pollTimer = null; _beStage = null;
  if (!on) { card.classList.add("hidden"); return; }
  card.classList.remove("hidden");
  $$("#rc-stages li").forEach(li => { li.className = ""; });
  const fill = $("#rc-fill"), t0 = performance.now();
  let cur = 0, stage = -1;
  fill.style.width = "0%";
  if (api().progress) {
    _pollTimer = setInterval(async () => {
      try {
        const p = await api().progress();
        if (p && p.ok && typeof p.stage === "number") _beStage = p.stage;
      } catch (e) { /* keep the time-based fallback */ }
    }, 600);
  }
  _runTimer = setInterval(() => {
    const t = (performance.now() - t0) / 1000;
    $("#rc-elapsed").textContent = t.toFixed(1) + "s";
    let s, cap;
    if (_beStage != null) {
      s = BE_TO_FRONT[Math.min(_beStage, BE_TO_FRONT.length - 1)];
      cap = BE_CAPS[Math.min(_beStage, BE_CAPS.length - 1)] * 100;
    } else {
      s = 0;
      RUN_STAGES.forEach((st, i) => { if (t >= st.at) s = i; });
      cap = RUN_STAGES[s].cap * 100;
    }
    if (s !== stage) {
      stage = s;
      $$("#rc-stages li").forEach((li, i) => {
        li.classList.toggle("done", i < stage);
        li.classList.toggle("active", i === stage);
      });
    }
    cur += (cap - cur) * 0.055;   // asymptotic crawl
    fill.style.width = cur.toFixed(1) + "%";
  }, 90);
}
/* result arrived: run the bar to 100%, tick every stage, then hand over */
function finishProgress(cb) {
  clearInterval(_runTimer); _runTimer = null;
  clearInterval(_pollTimer); _pollTimer = null;
  $$("#rc-stages li").forEach(li => { li.classList.add("done"); li.classList.remove("active"); });
  $("#rc-fill").style.width = "100%";
  setTimeout(cb, 420);
}

/* ---- error modal (clearer than an inline line) ---- */
function showModal(tag, title, message, detail) {
  $("#modal-tag").textContent = tag || "Error";
  $("#modal-title").textContent = title || "Something went wrong";
  $("#modal-msg").textContent = message || "";
  const d = $("#modal-detail");
  d.textContent = detail || ""; d.style.display = detail ? "block" : "none";
  $("#modal").classList.remove("hidden");
}
function hideModal() { $("#modal").classList.add("hidden"); }
$("#modal-x").addEventListener("click", hideModal);
$("#modal-ok").addEventListener("click", hideModal);
$("#modal-backdrop").addEventListener("click", hideModal);
document.addEventListener("keydown", e => { if (e.key === "Escape") hideModal(); });

/* Turn a raw backend error into a plain-language title + hint.
   Ordered most-specific first; these are the v2 failure modes, in roughly the
   order they actually happen to people. */
const ERROR_HINTS = [
  [/must be given together|reference_sample|reference_standards/i, "Reference pair incomplete",
    "A custom reference needs BOTH runs — the reference sample and the reference standards, from " +
    "the same column and gradient. Set the other one under Advanced, or use “Reset to bundled " +
    "reference” to go back to the Column 25 reference."],
  [/share (fewer|less) than|<\s*5 matched pairs|too few matched pairs between/i,
    "Your standards run and the reference barely match",
    "Almost always one of two things: the polarity is wrong (a negative-mode run matched against " +
    "a positive-mode reference shares almost no m/z), or the wrong file was chosen for the " +
    "standards run. Check Polarity, then check that the standards file really is the standards " +
    "mixture run."],
  [/CalibrationError/i, "Not enough matched pairs to fit the curve",
    "The curve is built from features matched by m/z between your runs and the reference runs. " +
    "Too few survived. Check that the polarity matches how the data were acquired, that the " +
    "sample and standards runs come from the SAME column and gradient, and that the m/z values " +
    "are on the same scale (Da, not ppm-shifted). Widening the feature-match window under " +
    "Advanced is the last resort, not the first."],
  [/PanelError/i, "Panel or reference-data problem",
    "The chosen panel could not be built for this polarity, or the bundled reference dataset is " +
    "missing from this build. Use “positive” or “negative” to match the acquisition mode; if the " +
    "message mentions reference_data, the app bundle is incomplete."],
  [/ColumnResolutionError|InputFormatError/i, "Unrecognised input table",
    "The file format, or its m/z / retention-time columns, couldn't be read. Confirm it's an " +
    "MS-DIAL, MZmine, MassCube or LipidScreener export — and that you picked the feature table, " +
    "not a summary or a peak list."],
  [/RTUnitError/i, "Retention-time unit problem",
    "The retention-time unit couldn't be inferred from the table. Retention times in seconds and " +
    "in minutes look identical to a parser; set the unit explicitly if your export is unusual."],
  [/A standards run is required/i, "Standards run missing",
    "The cross-column method needs two runs from your column: the sample and a standards run. If " +
    "the standards are spiked into the sample itself, choose that same file for both."],
  [/ConfigError/i, "Missing or invalid input", null],
];
function showError(raw) {
  const err = String(raw || "Calibration failed.");
  const technical = err.includes(": ") ? err.slice(err.indexOf(": ") + 2) : err;
  for (const [re, title, hint] of ERROR_HINTS) {
    if (re.test(err)) { showModal("Error", title, hint || technical, hint ? technical : ""); return; }
  }
  showModal("Error", "Calibration failed", technical, "");
}

/* ===================== OUTPUT VIEW ===================== */
/* cross-fade + slight rise between views, so a finished run never "snaps" */
function _switchView(show, hide) {
  hide.classList.add("pre");
  setTimeout(() => { hide.classList.add("hidden"); hide.classList.remove("pre"); }, 300);
  show.classList.remove("hidden");
  show.classList.add("pre");
  void show.offsetHeight;
  show.classList.remove("pre");
}
function openOutput() { _switchView($("#view-output"), $("#view-input")); }
function openInput() { _switchView($("#view-input"), $("#view-output")); }
/* + New calibration — reset to a completely fresh input screen */
function resetInput() {
  state.params = { samples: "", standards: "", polarity: "positive",
                   mixture: MIX_DEFAULT, reference_sample: "", reference_standards: "" };
  if (MIXTURES.length) selectMixture(MIX_DEFAULT);
  state.bundle = null;
  Object.keys(_secCache).forEach(k => delete _secCache[k]);
  $("#main-body").innerHTML = "";
  FILE_FIELDS.forEach(([key]) => clearField(key));
  const seg = $('[data-seg="polarity"]');
  if (seg) $$("button", seg).forEach(b => b.classList.toggle("active", b.dataset.val === "positive"));
  $$("#adv-body input").forEach(inp => {                 // restore HTML defaults
    if (inp.type === "checkbox") inp.checked = inp.defaultChecked;
    else inp.value = inp.defaultValue;
  });
  $$("#adv-body select").forEach(sel => {
    Array.from(sel.options).forEach(o => { o.selected = o.defaultSelected; });
  });
  const adv = $("#adv");
  adv.classList.remove("open");
  $("#adv-toggle").setAttribute("aria-expanded", "false");
  refHint();
  setStatus("", ""); showProgress(false);
}
$("#new-run").addEventListener("click", () => { resetInput(); openInput(); });
const ghBtn = $("#sb-github");
if (ghBtn) ghBtn.addEventListener("click", () => { if (api().open_github) api().open_github(); });

/* which sections have something to show for this run */
function sectionEnabled(id) {
  const m = state.bundle.meta || {};
  if (id === "curve") return !!m.has_curve;
  return true;
}
function buildNav() {
  const nav = $("#nav"); nav.innerHTML = "";
  SECTIONS.forEach(s => {
    const b = document.createElement("button");
    b.className = "nav-item"; b.dataset.id = s.id;
    b.innerHTML = `<span class="idx">${s.idx}</span>${s.label}`;
    if (!sectionEnabled(s.id)) {
      b.disabled = true;
      b.title = "Not available for this run";
    } else {
      b.addEventListener("click", () => selectSection(s.id));
    }
    nav.appendChild(b);
  });
  const eng = $("#sb-engine-d");
  const m = state.bundle.meta || {};
  if (eng && m.reference_label) {
    eng.textContent = `RT Anchor · cross-column onto ${m.reference_label}`;
  }
}

/* section DOM is rendered once and kept alive (plotly instances included) so
   switching sections never rebuilds or re-parses anything */
const _secCache = {};

function selectSection(id) {
  state.section = id;
  const s = SECTIONS.find(x => x.id === id);
  $("#sec-title").textContent = s.title;
  $$(".nav-item").forEach(n => n.classList.toggle("active", n.dataset.id === id));
  CK.hideAll();
  if (structPop) structPop.style.display = "none";
  const body = $("#main-body");
  if (!_secCache[id]) {
    const node = document.createElement("div");
    node.className = "section-body";
    _secCache[id] = node;          /* cached before render: a render error must never
                                      leave an uncached, unhidden section stacked */
    body.appendChild(node);
    try { renderSection(id, node); }
    catch (e) {
      node.innerHTML = `<div class="panel"><div class="panel-h">${esc(s.title)}</div>` +
        `<div class="fig-note">This section could not be rendered: ${esc(String(e))}</div></div>`;
    }
  }
  Object.keys(_secCache).forEach(k => _secCache[k].classList.toggle("offscreen", k !== id));
  const activeNode = _secCache[id];
  if (activeNode) {
    void activeNode.offsetHeight;
    if (window.Plotly) {
      activeNode.querySelectorAll(".js-plotly-plot").forEach(el => {
        try { Plotly.Plots.resize(el); } catch(_){}
      });
    }
  }
  const main = $(".main");
  if (main) main.scrollTop = 0;
}

/* build a light card with a black header strip; returns the plot host div */
function panelInto(body, title) {
  const panel = document.createElement("div"); panel.className = "panel";
  const h = document.createElement("div"); h.className = "panel-h"; h.textContent = title;
  panel.appendChild(h);
  const host = document.createElement("div"); host.className = "plot"; panel.appendChild(host);
  body.appendChild(panel);
  return { panel, host };
}

function renderSection(id, body) {
  if (id === "overview") return renderOverview(body);
  if (id === "table") return renderTable(body);
  if (id === "export") return renderExport(body);
  if (id === "profile") {
    const figs = state.bundle.figures || {};
    if (figs.profile_before_after) {
      const p1 = panelInto(body, "Before vs after calibration");
      drawPlot(p1.host, figs.profile_before_after);
    }
    if (figs.profile_calibrated) {
      const p2 = panelInto(body, "Calibrated profile — input vs reference");
      drawPlot(p2.host, figs.profile_calibrated);
    }
    const note = state.bundle.notes && state.bundle.notes.profile;
    if (note) {
      const n = document.createElement("div"); n.className = "fig-note"; n.innerHTML = note;
      body.appendChild(n);
    }
    return;
  }
  const sec = SECTIONS.find(x => x.id === id);
  const { panel, host } = panelInto(body, sec ? sec.title : "");
  if (id === "detection") {
    // the package's own Plotly figure, so the app and the report agree
    drawPlot(host, state.bundle.figures[id], attachStructures);
  } else {
    CK.render(id, host, state.bundle);   // curve (SVG)
  }
  const note = state.bundle.notes && state.bundle.notes[id];
  if (note) {
    const n = document.createElement("div"); n.className = "fig-note"; n.innerHTML = note;
    panel.appendChild(n);
  }
}

function renderOverview(body) {
  const k = state.bundle.kpis || [];
  const rowTiles = k.filter(t => t.rows && t.rows.length);
  const valTiles = k.filter(t => !(t.rows && t.rows.length));

  /* metric rows: each value carries its own caption, because the v2 tiles pair
     different things (two runs, two pair sources, a median and a P90) */
  if (rowTiles.length) {
    const rows = rowTiles.map(t => {
      const sub = t.sub ? `<span class="kt-m-s">${esc(t.sub)}</span>` : "";
      const cells = t.rows.map(r =>
        `<td class="kt-v"><span class="kt-c">${esc(r.k)}</span>${esc(r.v)}</td>`).join("");
      return `<tr><th class="kt-m"><span class="kt-m-l">${esc(t.label)}</span>${sub}</th>${cells}</tr>`;
    }).join("");
    const tbl = document.createElement("div"); tbl.className = "kpi-table";
    tbl.innerHTML = `<table class="kpi-metrics"><tbody>${rows}</tbody></table>`;
    body.appendChild(tbl);
  }

  /* single-value metrics stay as bold tiles */
  if (valTiles.length) {
    const grid = document.createElement("div"); grid.className = "kpis kpis-" + valTiles.length;
    grid.innerHTML = valTiles.map(t =>
      `<div class="kpi"><div class="kpi-l">${esc(t.label)}</div><div class="kpi-v">${esc(t.value)}</div>` +
      (t.sub ? `<div class="kpi-s">${esc(t.sub)}</div>` : "") + `</div>`
    ).join("");
    body.appendChild(grid);
  }

  const { host } = panelInto(body, "Quality fingerprint — outward is better");
  CK.radar(host, state.bundle.radar);
}

function renderExport(body) {
  body.innerHTML =
    `<div class="export-all">
       <button class="btn-exp big" id="exp-all"><div class="t">Export all outputs →</div>
         <div class="d">Writes everything into one folder: calibrated CSV · model.json ·
           supporting evidence files · log.txt · interactive HTML report ·
           run_info.json</div></button>
     </div>
     <div class="export-row">
       <button class="btn-exp" id="exp-csv"><div class="t">Calibrated table</div>
         <div class="d">CSV — your original columns plus Cal_RT_min, iRT, their uncertainties,
           iRT_reliability and the extrapolation flag</div></button>
       <button class="btn-exp" id="exp-report"><div class="t">Report</div>
         <div class="d">Interactive HTML — detection QC and the cross-column curve</div></button>
       <button class="btn-exp" id="exp-info"><div class="t">Run info</div>
         <div class="d">JSON — parameters, reference pair, run time, versions and a result
           summary</div></button>
     </div>
     <div id="exp-manifest"></div>`;
  const doExport = async (busyMsg, call, okMsg, after) => {
    if (busyMsg) toast(busyMsg);
    let r;
    try { r = await call(); } catch (e) { r = { ok: false, error: String(e) }; }
    if (r && !r.ok && r.error === "cancelled") { hideToast(); return; }   // user dismissed the dialog
    if (!r || !r.ok) { hideToast(); showError((r && r.error) || "Export failed."); return; }
    toast(okMsg(r));
    if (after) after(r);
  };
  $("#exp-all").addEventListener("click", () => doExport("Writing all outputs…",
    () => api().export_all(), r => `Wrote ${r.n_files} files to ${shortPath(r.dir)}`,
    r => renderManifest(r)));
  $("#exp-csv").addEventListener("click", () => doExport(null,
    () => api().export_csv(), r => "Saved: " + shortPath(r.path)));
  $("#exp-report").addEventListener("click", () => doExport("Writing report…",
    () => api().export_report(), () => "Report written"));
  $("#exp-info").addEventListener("click", () => doExport(null,
    () => api().export_run_info(), r => "Saved: " + shortPath(r.path)));
}

/* list exactly what "Export all" wrote, so the two new v2 files are named */
function renderManifest(r) {
  const host = $("#exp-manifest"); if (!host) return;
  const files = r.files || [];
  if (!files.length) { host.innerHTML = ""; return; }
  const rows = files.map(f =>
    `<tr><td>${esc(f.label)}</td><td class="mono">${esc(f.name)}</td></tr>`).join("");
  const warn = r.report === false
    ? `<div class="fig-note">The report could not be rendered for this run; the data files above
        are unaffected — see the run log for the reason.</div>` : "";
  host.innerHTML = `<div class="panel exp-manifest"><div class="panel-h">Written to ${esc(shortPath(r.dir))}</div>` +
    `<div class="table-wrap"><table class="data-table"><thead><tr><th>Output</th><th>File</th></tr></thead>` +
    `<tbody>${rows}</tbody></table></div>${warn}</div>`;
}

/* ---- calibrated-table preview ---- */
const esc = s => String(s).replace(/[&<>]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));
const NUMERIC_COLS = new Set(["RT", "m/z", "Cal_RT (min)", "Cal_RT ± (min)", "iRT", "iRT ±", "iRT spread"]);
const TABLE_PAGE_SIZE = 100;
let _tablePage = 0;
function renderTable(body) {
  const t = state.bundle.table;
  if (!t || !t.columns || !t.rows.length) {
    body.innerHTML = `<div class="panel"><div class="panel-h">No table available</div></div>`; return;
  }
  const totalPages = Math.ceil(t.rows.length / TABLE_PAGE_SIZE);
  _tablePage = 0;

  function drawPage(page) {
    body.innerHTML = "";
    _tablePage = page;
    const start = page * TABLE_PAGE_SIZE;
    const end = Math.min(start + TABLE_PAGE_SIZE, t.rows.length);
    const pageRows = t.rows.slice(start, end);

    const panel = document.createElement("div"); panel.className = "panel";
    const h = document.createElement("div"); h.className = "panel-h";
    h.textContent = `Calibrated table · rows ${start + 1}–${end} of ${Number(t.n_total).toLocaleString()}`;
    panel.appendChild(h);

    const head = "<thead><tr>" + t.columns.map(c => `<th>${esc(c)}</th>`).join("") + "</tr></thead>";
    const bodyHtml = pageRows.map(r => "<tr>" + r.map((v, i) => {
      const col = t.columns[i];
      if (v === null || v === undefined) return `<td class="num"><span class="na">—</span></td>`;
      if (col === "reliability") return `<td><span class="rel rel-${esc(v)}">${esc(v)}</span></td>`;
      const cls = NUMERIC_COLS.has(col) ? "num" : "";
      return `<td class="${cls}">${esc(v)}</td>`;
    }).join("") + "</tr>").join("");
    const wrap = document.createElement("div"); wrap.className = "table-wrap";
    wrap.innerHTML = `<table class="data-table">${head}<tbody>${bodyHtml}</tbody></table>`;
    panel.appendChild(wrap);

    if (totalPages > 1) {
      const nav = document.createElement("div"); nav.className = "table-nav";
      nav.style.cssText = "display:flex;align-items:center;gap:8px;padding:8px 12px;font:13px sans-serif;";
      const prev = document.createElement("button"); prev.className = "btn-sm"; prev.textContent = "← Prev";
      prev.disabled = page === 0;
      prev.addEventListener("click", () => drawPage(page - 1));
      const info = document.createElement("span");
      info.textContent = `Page ${page + 1} of ${totalPages}`;
      const next = document.createElement("button"); next.className = "btn-sm"; next.textContent = "Next →";
      next.disabled = page >= totalPages - 1;
      next.addEventListener("click", () => drawPage(page + 1));
      nav.appendChild(prev); nav.appendChild(info); nav.appendChild(next);
      panel.appendChild(nav);
    }

    const note = document.createElement("div"); note.className = "fig-note";
    note.innerHTML = "Preview of the calibrated feature table — your <b>RT</b> and <b>m/z</b> beside " +
      "<b>Cal_RT</b> (the same feature's retention time on the reference column) and <b>iRT</b> " +
      "(the dimensionless 1–100 index). The full table, with every original column, is under " +
      "<b>Export</b>.";
    panel.appendChild(note);
    body.appendChild(panel);
  }
  drawPage(0);
}

/* ---- Plotly (detection + profile, the package figures) ---- */
/* `after(div)` runs once newPlot has resolved — earlier than that the div is
   not a Plotly graph yet (no .on()), so event hooks must wait for it. */
function drawPlot(div, figStr, after) {
  const f = typeof figStr === "string" ? JSON.parse(figStr) : figStr;
  window._ensurePlotly().then(function() {
    return Plotly.newPlot(div, f.data, f.layout, { displaylogo: false, responsive: true,
      modeBarButtonsToRemove: ["select2d", "lasso2d"], toImageButtonOptions: { format: "svg" } });
  }).then(function() { if (after) after(div); });
}
let structPop;
function attachStructures(div) {
  const S = state.bundle.structures; if (!S || !Object.keys(S).length) return;
  if (!structPop) {
    structPop = document.createElement("div"); structPop.id = "struct-pop";
    structPop.style.cssText = "position:fixed;display:none;z-index:1000;background:#FFFFFF;" +
      "border:1px solid #EAE8E2;border-radius:10px;padding:9px 11px;" +
      "box-shadow:0 1px 2px rgba(43,43,43,.05),0 12px 32px -8px rgba(43,43,43,.13);pointer-events:none;";
    document.body.appendChild(structPop);
    document.addEventListener("mousemove", e => { structPop._x = e.clientX; structPop._y = e.clientY; });
  }
  div.on("plotly_hover", d => {
    let nm = d.points && d.points[0] && d.points[0].customdata;
    if (Array.isArray(nm)) nm = nm[0];
    if (nm && S[nm]) {
      structPop.innerHTML = `<div style="font:650 12px Inter,sans-serif;color:#2B2B2B;margin-bottom:4px">${nm}</div>` +
        `<div style="background:#fff">${S[nm]}</div>`;
      structPop.style.display = "block";
      structPop.style.left = Math.min(structPop._x + 14, window.innerWidth - 260) + "px";
      structPop.style.top = Math.min(structPop._y + 14, window.innerHeight - 200) + "px";
    }
  });
  div.on("plotly_unhover", () => { if (structPop) structPop.style.display = "none"; });
}

/* ---- async calibration handler ---- */
window.__onCalibrationDone = async function() {
  try {
    var res = await api().get_calibration_result();
    var btn = $("#run");
    btn.disabled = false;
    if (!res || !res.ok) { showProgress(false); setStatus("", ""); showError((res && res.error) || "Calibration failed."); return; }
    setStatus("Preparing results…", "busy");
    const figs = res.figures || {};
    Object.keys(figs).forEach(k => {              // parse figure JSON once, up front
      if (typeof figs[k] === "string") {
        try { figs[k] = JSON.parse(figs[k]); } catch (e) {}
      }
    });
    setStatus("");
    finishProgress(() => showResults(res));
  } catch (e) {
    var btn = $("#run");
    if (btn) btn.disabled = false;
    showProgress(false);
    setStatus("", "");
    showError(String(e));
  }
};

/* Switch to the output view and paint the Overview atomically, then warm up
   Plotly and pre-render the remaining sections offscreen in small idle chunks,
   so the first click on any section never waits on anything. */
function showResults(bundle) {
  state.bundle = bundle;
  buildNav();
  openOutput();
  selectSection("overview");
  window._ensurePlotly().then(() => {
    const b = state.bundle;
    let i = 0;
    const step = () => {
      if (state.bundle !== b) return;             // user started a new run
      while (i < SECTIONS.length) {
        const s = SECTIONS[i++];
        if (s.id !== state.section && !_secCache[s.id]) {
          const node = document.createElement("div");
          node.className = "section-body offscreen";
          $("#main-body").appendChild(node);
          try { renderSection(s.id, node); _secCache[s.id] = node; }
          catch (e) { node.remove(); }   // left uncached: re-rendered on first real click
          return setTimeout(step, 60);
        }
      }
    };
    step();
  }).catch(() => {});
}

/* ===================== BOOT ===================== */
/* Two things must be true before the input screen means anything: the
   pywebview bridge has to exist, and the Python side has to have finished
   importing the calibration engine. On a fast Mac both are done before the
   first paint. On a low-spec Windows machine the bridge lands late and the
   engine import (scipy + scikit-learn + statsmodels) can take half a minute,
   and the old code just called the API on `load` — which meant it either
   silently fell through to the browser PREVIEW stub and showed fake panel
   cards, or left the standards-panel section empty with nothing to explain
   why. The overlay stays up and says which of the two we are waiting on. */
const boot = {
  el: () => document.getElementById("boot"),
  say(msg, sub) {
    const m = document.getElementById("boot-msg");
    const s = document.getElementById("boot-sub");
    if (m && msg != null) m.textContent = msg;
    if (s && sub != null) s.textContent = sub;
  },
  fail(msg, sub) {
    const b = this.el(); if (b) b.classList.add("err");
    this.say(msg, sub);
  },
  done() {
    const b = this.el(); if (!b) return;
    b.classList.add("gone");
    setTimeout(() => b.remove(), 400);
  },
};

/* Resolve the bridge, or decide we are in a plain browser. */
function whenBridgeReady() {
  return new Promise(resolve => {
    if (window.pywebview && window.pywebview.api) return resolve(true);
    if (window.__PREVIEW__) return resolve(false);      // explicit browser preview
    let settled = false;
    const go = ok => { if (!settled) { settled = true; resolve(ok); } };
    window.addEventListener("pywebviewready", () => go(true), { once: true });
    /* pywebviewready can fire before this script runs, so poll as well. */
    const t0 = performance.now();
    const tick = setInterval(() => {
      if (window.pywebview && window.pywebview.api) { clearInterval(tick); go(true); }
      else if (performance.now() - t0 > 20000) { clearInterval(tick); go(false); }
      else if (performance.now() - t0 > 2500) boot.say("Connecting to the app…");
    }, 120);
  });
}

/* Wait for the Python-side engine import, narrating the wait. */
async function whenEngineReady() {
  if (!api().engine_status) return true;              // preview stub
  for (;;) {
    let st;
    try { st = await api().engine_status(); }
    catch (e) { return true; }                        // older backend: just proceed
    if (!st) return true;
    if (st.error) {
      boot.fail("The calibration engine could not be loaded.", st.error);
      return false;
    }
    if (st.ready) return true;
    const s = Number(st.elapsed || 0);
    boot.say("Starting the calibration engine…",
      s > 4 ? `${s.toFixed(0)}s — first launch is the slow one; the libraries are being read from disk.` : "");
    await new Promise(r => setTimeout(r, 250));
  }
}

async function boot_() {
  boot.say("Starting…");
  const bridged = await whenBridgeReady();
  if (!bridged && !window.__PREVIEW__) {
    /* Opened as a file in a browser, or the bridge never arrived. The UI still
       works against the preview stub; say so rather than pretending. */
    boot.say("Preview mode — no calibration backend attached.");
  }
  if (bridged && !(await whenEngineReady())) return;   // overlay stays, showing why
  boot.say("Loading standards panels…");
  await Promise.all([initMixtures(), initReference()]);
  boot.done();
  if (!window.pywebview && window.__PREVIEW__ && window.__PREVIEW__.ok) {
    showResults(window.__PREVIEW__);
  }
}
if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", boot_, { once: true });
} else {
  boot_();
}
