/* RT Anchor desktop — front-end logic (pywebview bridge + browser preview) */
"use strict";

const state = {
  params: { samples: "", standards: "", single_files: "", polarity: "positive", mixture: "mix21" },
  bundle: null,
  section: "overview",
};

const SECTIONS = [
  { id: "overview", label: "Overview", idx: "01", title: "Overview" },
  { id: "detection", label: "Detection", idx: "02", title: "Standard detection" },
  { id: "profile", label: "Profile", idx: "03", title: "Feature-intensity profile" },
  { id: "warp", label: "Warp", idx: "04", title: "Calibration warp" },
  { id: "repeatability", label: "Repeatability", idx: "05", title: "Injection repeatability" },
  { id: "table", label: "Table", idx: "06", title: "Calibrated table" },
  { id: "export", label: "Export", idx: "07", title: "Export" },
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
  mixture_previews: async () => ({ ok: true, default: "mix21", mixtures: [
    { key: "mix15", label: "Mix 15",
      rows: [{ name: "PC 6:0_6:0", mz_pos: 454.2569, mz_neg: 498.2474, rt_ref_min: 1.06 },
             { name: "PC 12:0_12:0", mz_pos: 622.4448, mz_neg: 666.4362, rt_ref_min: 3.6 },
             { name: "PC 14:0_14:0", mz_pos: 678.5069, mz_neg: 722.4983, rt_ref_min: 5.63 },
             { name: "PC 16:0_18:1", mz_pos: 760.5856, mz_neg: 804.5766, rt_ref_min: 8.4 },
             { name: "DG 16:0_18:0", mz_pos: 614.5723, mz_neg: 641.5367, rt_ref_min: 12.65 },
             { name: "CE 18:1", mz_pos: 668.6345, mz_neg: null, rt_ref_min: 18.39 }] },
    { key: "mix21", label: "Mix 21",
      rows: [{ name: "PC 6:0_6:0", mz_pos: 454.2561, mz_neg: 498.2474, rt_ref_min: 1.5 },
             { name: "PC 8:0_8:0", mz_pos: 510.319, mz_neg: 554.31, rt_ref_min: 3.2 },
             { name: "PC 12:0_12:0", mz_pos: 622.444, mz_neg: 666.4352, rt_ref_min: 6.8 },
             { name: "PC 14:0_16:0", mz_pos: 706.5381, mz_neg: 750.5291, rt_ref_min: 11.1 },
             { name: "PC 16:0_18:1", mz_pos: 760.5846, mz_neg: 804.576, rt_ref_min: 12.9 },
             { name: "DG 18:0_18:0", mz_pos: 642.6031, mz_neg: 669.5675, rt_ref_min: 18.0 },
             { name: "CE 20:5", mz_pos: 688.6027, mz_neg: null, rt_ref_min: 20.5 },
             { name: "TG 18:0_18:2_18:0", mz_pos: 904.8328, mz_neg: null, rt_ref_min: 21.8 }] },
  ] }),
  run_calibration: async () => window.__PREVIEW__ || { ok: false, error: "no preview data" },
  export_csv: async () => ({ ok: true, path: "(preview)/calibrated.csv" }),
  export_report: async () => ({ ok: true, paths: {} }),
  export_run_info: async () => ({ ok: true, path: "(preview)/run_info.json" }),
  export_all: async () => ({ ok: true, dir: "(preview)/outputs", n_files: 7 }),
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
function shortPath(p) { return p ? p.split("/").slice(-2).join("/") : ""; }

/* set a file/folder field + flip it to the "filled" (green) state */
function markField(target, path) {
  state.params[target] = path;
  const el = $("#path-" + target);
  el.textContent = shortPath(path); el.title = path; el.classList.add("set");
  const field = document.querySelector(`.field[data-key="${target}"]`);
  if (field) field.classList.add("filled");
}

$$(".btn-file").forEach(btn => btn.addEventListener("click", async () => {
  const target = btn.dataset.target, kind = btn.dataset.pick;
  const path = kind === "folder" ? await api().pick_folder() : await api().pick_file();
  if (!path) return;
  markField(target, path);
}));

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

/* ---- standard-mixture cards + preview (either/or choice) ---- */
let MIXTURES = [];                       // payload of api().mixture_previews()
let MIX_DEFAULT = "mix21";               // backend's pre-selected mixture
function renderMixCards() {
  const wrap = $("#mix-cards"); wrap.innerHTML = "";
  MIXTURES.forEach(m => {
    const b = document.createElement("button");
    b.type = "button"; b.className = "mix-card"; b.dataset.key = m.key;
    if (state.params.mixture === m.key) b.classList.add("active");
    b.textContent = m.label;
    b.addEventListener("click", () => selectMixture(m.key));
    wrap.appendChild(b);
  });
}
function renderMixPreview() {
  const tbl = $("#mix-table");
  const m = MIXTURES.find(x => x.key === state.params.mixture);
  if (!m) { tbl.innerHTML = ""; return; }
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
}

/* Advanced section collapse/expand */
const advToggle = $("#adv-toggle");
if (advToggle) advToggle.addEventListener("click", () => {
  const open = $("#adv").classList.toggle("open");
  advToggle.setAttribute("aria-expanded", open ? "true" : "false");
});

/* pull the Advanced matching parameters into params (blank -> package default) */
function collectAdvanced() {
  const v = id => { const el = $(id); return el ? el.value.trim() : ""; };
  state.params.mz_tol_ppm = v("#adv-mz");
  state.params.rt_window_min = v("#adv-rtw");
  state.params.min_anchors = v("#adv-minanchors");
  const ex = $("#adv-extrapolate");
  state.params.extrapolate = ex ? ex.checked : false;
}

$("#run").addEventListener("click", run);
async function run() {
  if (!state.params.samples) { setStatus("Choose a sample feature table first.", "err"); return; }
  if (!state.params.standards) { setStatus("Choose a standards mixture table — it is required.", "err"); return; }
  collectAdvanced();
  const btn = $("#run"); btn.disabled = true;
  showProgress(true);
  setStatus("Calibrating — identifying anchors, fitting the warp…", "busy");
  let res;
  try { res = await api().run_calibration(state.params); }
  catch (e) { res = { ok: false, error: String(e) }; }
  btn.disabled = false; showProgress(false);
  if (!res || !res.ok) { setStatus("", ""); showError((res && res.error) || "Calibration failed."); return; }
  setStatus("");
  state.bundle = res;
  buildNav(); openOutput(); selectSection("overview");
}
function setStatus(msg, cls) { const s = $("#status"); s.textContent = msg; s.className = "status" + (cls ? " " + cls : ""); }
function showProgress(on) { $("#progress").classList.toggle("hidden", !on); }

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

/* turn a raw backend error into a plain-language title + hint (+ technical detail) */
const ERROR_HINTS = [
  [/AnchorIdentificationError/i, "No standards were detected",
    "None of the panel standards were found. Check that the polarity matches how the data were acquired, that the standards run really contains the panel, and try widening the m/z tolerance or RT window under Advanced."],
  [/PanelError/i, "Panel / polarity problem",
    "The standard panel couldn't be built for this polarity. Use “positive” or “negative” to match the acquisition mode."],
  [/CalibrationError/i, "The warp couldn't be fitted",
    "There were too few or non-monotonic anchors to fit a calibration curve. Check the standards run and the matching tolerances."],
  [/ColumnResolutionError|InputFormatError/i, "Unrecognised input table",
    "The file format, or its m/z / retention-time columns, couldn't be read. Confirm it's an MS-DIAL, MZmine, MassCube or LipidScreener export."],
  [/RTUnitError/i, "Retention-time unit problem",
    "The retention-time unit couldn't be inferred from the table."],
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
function openOutput() { $("#view-input").classList.add("hidden"); $("#view-output").classList.remove("hidden"); }
function openInput() { $("#view-output").classList.add("hidden"); $("#view-input").classList.remove("hidden"); }
/* + New calibration — reset to a completely fresh input screen */
const FILE_FIELDS = [
  ["samples", "no file selected"],
  ["standards", "no file selected"],
  ["single_files", "no folder selected"],
];
function resetInput() {
  state.params = { samples: "", standards: "", single_files: "", polarity: "positive", mixture: MIX_DEFAULT };
  if (MIXTURES.length) selectMixture(MIX_DEFAULT);
  state.bundle = null;
  FILE_FIELDS.forEach(([key, placeholder]) => {
    const el = $("#path-" + key);
    if (el) { el.textContent = placeholder; el.title = ""; el.classList.remove("set"); }
    const field = document.querySelector(`.field[data-key="${key}"]`);
    if (field) field.classList.remove("filled");
  });
  const seg = $('[data-seg="polarity"]');
  if (seg) $$("button", seg).forEach(b => b.classList.toggle("active", b.dataset.val === "positive"));
  $$("#adv-body input").forEach(inp => {                 // restore HTML defaults
    if (inp.type === "checkbox") inp.checked = inp.defaultChecked;
    else inp.value = inp.defaultValue;
  });
  const adv = $("#adv");
  adv.classList.remove("open");
  $("#adv-toggle").setAttribute("aria-expanded", "false");
  setStatus("", ""); showProgress(false);
}
$("#new-run").addEventListener("click", () => { resetInput(); openInput(); });
const ghBtn = $("#sb-github");
if (ghBtn) ghBtn.addEventListener("click", () => { if (api().open_github) api().open_github(); });

function buildNav() {
  const nav = $("#nav"); nav.innerHTML = "";
  const hasRep = state.bundle.meta.has_repeatability;
  const hasWarp = state.bundle.warp && !state.bundle.warp.empty;
  SECTIONS.forEach(s => {
    const b = document.createElement("button");
    b.className = "nav-item"; b.dataset.id = s.id;
    b.innerHTML = `<span class="idx">${s.idx}</span>${s.label}`;
    const off = (s.id === "repeatability" && !hasRep) || (s.id === "warp" && !hasWarp);
    if (off) b.disabled = true;
    else b.addEventListener("click", () => selectSection(s.id));
    nav.appendChild(b);
  });
}

function selectSection(id) {
  state.section = id;
  const s = SECTIONS.find(x => x.id === id);
  $("#sec-title").textContent = s.title;
  $$(".nav-item").forEach(n => n.classList.toggle("active", n.dataset.id === id));
  render(id);
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

function render(id) {
  const body = $("#main-body");
  $$(".js-plotly-plot", body).forEach(d => { try { Plotly.purge(d); } catch (e) {} });
  CK.hideAll();
  if (structPop) structPop.style.display = "none";
  body.innerHTML = "";
  if (id === "overview") return renderOverview(body);
  if (id === "table") return renderTable(body);
  if (id === "export") return renderExport(body);
  const sec = SECTIONS.find(x => x.id === id);
  const { panel, host } = panelInto(body, sec ? sec.title : "");
  if (id === "detection" || id === "profile") {
    // reverted to the package's Plotly figures
    drawPlot(host, state.bundle.figures[id]);
    if (id === "detection") attachStructures(host);
  } else {
    CK.render(id, host, state.bundle);   // radar (via renderOverview) + warp + repeatability
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

  /* shared-column table: one "sample" column + one "standards run" column */
  if (rowTiles.length) {
    const rows = rowTiles.map(t => {
      const kv = {};
      t.rows.forEach(r => { kv[r.k] = r.v; });
      const v = key => (kv[key] !== undefined && kv[key] !== null) ? esc(kv[key]) : "<span class=\"na\">—</span>";
      const sub = t.sub ? `<span class="kt-m-s">${esc(t.sub)}</span>` : "";
      return `<tr><th class="kt-m"><span class="kt-m-l">${esc(t.label)}</span>${sub}</th>` +
             `<td class="kt-v">${v("samples")}</td>` +
             `<td class="kt-v">${v("standards run")}</td></tr>`;
    }).join("");
    const tbl = document.createElement("div"); tbl.className = "kpi-table";
    tbl.innerHTML =
      `<table class="kpi-metrics">
         <thead><tr><th class="kt-m">Metric</th><th class="kt-v">Sample</th><th class="kt-v">Standards run</th></tr></thead>
         <tbody>${rows}</tbody>
       </table>`;
    body.appendChild(tbl);
  }

  /* single-value metrics stay as bold tiles */
  if (valTiles.length) {
    const grid = document.createElement("div"); grid.className = "kpis";
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
         <div class="d">Writes everything into one folder: calibrated CSV · model.json · anchors.csv · log.txt · interactive HTML report · run_info.json</div></button>
     </div>
     <div class="export-row">
       <button class="btn-exp" id="exp-csv"><div class="t">Calibrated table</div>
         <div class="d">CSV — original columns + RI, RI_uncertainty, RI_reliability, RI_spread, flags</div></button>
       <button class="btn-exp" id="exp-report"><div class="t">Report</div>
         <div class="d">Interactive HTML — hover, zoom, pan (detection, warp, repeatability)</div></button>
       <button class="btn-exp" id="exp-info"><div class="t">Run info</div>
         <div class="d">JSON — parameters, run time, versions, and a result summary</div></button>
     </div>`;
  const doExport = async (busyMsg, call, okMsg) => {
    if (busyMsg) toast(busyMsg);
    let r;
    try { r = await call(); } catch (e) { r = { ok: false, error: String(e) }; }
    if (r && !r.ok && r.error === "cancelled") { hideToast(); return; }   // user dismissed the dialog
    if (!r || !r.ok) { hideToast(); showError((r && r.error) || "Export failed."); return; }
    toast(okMsg(r));
  };
  $("#exp-all").addEventListener("click", () => doExport("Writing all outputs…",
    () => api().export_all(), r => `Wrote ${r.n_files} files to ${shortPath(r.dir)}`));
  $("#exp-csv").addEventListener("click", () => doExport(null,
    () => api().export_csv(), r => "Saved: " + shortPath(r.path)));
  $("#exp-report").addEventListener("click", () => doExport("Writing report…",
    () => api().export_report(), () => "Report written"));
  $("#exp-info").addEventListener("click", () => doExport(null,
    () => api().export_run_info(), r => "Saved: " + shortPath(r.path)));
}

/* ---- calibrated-table preview ---- */
const esc = s => String(s).replace(/[&<>]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));
const NUMERIC_COLS = new Set(["RT", "m/z", "RI", "RI_uncertainty", "RI_spread"]);
function renderTable(body) {
  const t = state.bundle.table;
  if (!t || !t.columns || !t.rows.length) {
    body.innerHTML = `<div class="panel"><div class="panel-h">No table available</div></div>`; return;
  }
  const panel = document.createElement("div"); panel.className = "panel";
  const h = document.createElement("div"); h.className = "panel-h";
  h.textContent = `Calibrated table · first ${t.n_shown} of ${Number(t.n_total).toLocaleString()} rows`;
  panel.appendChild(h);
  const head = "<thead><tr>" + t.columns.map(c => `<th>${esc(c)}</th>`).join("") + "</tr></thead>";
  const bodyRows = t.rows.map(r => "<tr>" + r.map((v, i) => {
    const col = t.columns[i];
    if (v === null || v === undefined) return `<td class="num"><span class="na">—</span></td>`;
    if (col === "RI_reliability") return `<td><span class="rel rel-${esc(v)}">${esc(v)}</span></td>`;
    const cls = NUMERIC_COLS.has(col) ? "num" : "";
    return `<td class="${cls}">${esc(v)}</td>`;
  }).join("") + "</tr>").join("");
  const wrap = document.createElement("div"); wrap.className = "table-wrap";
  wrap.innerHTML = `<table class="data-table">${head}<tbody>${bodyRows}</tbody></table>`;
  panel.appendChild(wrap);
  const note = document.createElement("div"); note.className = "fig-note";
  note.innerHTML = "Preview of the calibrated feature table — <b>RT</b> and <b>m/z</b> beside the appended " +
    "<b>RI</b> columns. The full table (all rows and every original column) is available under <b>Export</b>.";
  panel.appendChild(note);
  body.appendChild(panel);
}

/* ---- Plotly (detection + profile only, reverted to the package figures) ---- */
function drawPlot(div, figStr) {
  const f = typeof figStr === "string" ? JSON.parse(figStr) : figStr;
  Plotly.newPlot(div, f.data, f.layout, { displaylogo: false, responsive: true,
    modeBarButtonsToRemove: ["select2d", "lasso2d"], toImageButtonOptions: { format: "svg" } });
}
let structPop;
function attachStructures(div) {
  const S = state.bundle.structures; if (!S || !Object.keys(S).length) return;
  if (!structPop) {
    structPop = document.createElement("div"); structPop.id = "struct-pop";
    structPop.style.cssText = "position:fixed;display:none;z-index:1000;background:#FFFFFF;" +
      "border:3px solid #141210;padding:8px 10px;box-shadow:5px 5px 0 rgba(20,18,16,.28);pointer-events:none;";
    document.body.appendChild(structPop);
    document.addEventListener("mousemove", e => { structPop._x = e.clientX; structPop._y = e.clientY; });
  }
  div.on("plotly_hover", d => {
    let nm = d.points && d.points[0] && d.points[0].customdata;
    if (Array.isArray(nm)) nm = nm[0];
    if (nm && S[nm]) {
      structPop.innerHTML = `<div style="font:800 12px sans-serif;color:#141210;margin-bottom:3px">${nm}</div>` +
        `<div style="background:#fff;padding:2px">${S[nm]}</div>`;
      structPop.style.display = "block";
      structPop.style.left = Math.min(structPop._x + 14, window.innerWidth - 260) + "px";
      structPop.style.top = Math.min(structPop._y + 14, window.innerHeight - 200) + "px";
    }
  });
  div.on("plotly_unhover", () => { if (structPop) structPop.style.display = "none"; });
}

/* preview mode: auto-open output if preview data is injected */
window.addEventListener("load", () => {
  initMixtures();
  if (!window.pywebview && window.__PREVIEW__ && window.__PREVIEW__.ok) {
    state.bundle = window.__PREVIEW__; buildNav(); openOutput(); selectSection("overview");
  }
});
