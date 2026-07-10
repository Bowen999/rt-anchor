/* RT Anchor desktop — front-end logic (pywebview bridge + browser preview) */
"use strict";

const state = {
  params: { samples: "", standards: "", single_files: "", polarity: "positive", manifest: "" },
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
    standards: "/preview/example/standards.txt", polarity: "positive" }),
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
  state.params.polarity = r.polarity || "positive";
  const seg = $('[data-seg="polarity"]');
  if (seg) $$("button", seg).forEach(b => b.classList.toggle("active", b.dataset.val === state.params.polarity));
  setStatus("Example dataset loaded — press Run calibration.", "ok");
});

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
  if (!state.params.standards) { setStatus("Choose a standards run — it is required.", "err"); return; }
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
$("#new-run").addEventListener("click", openInput);
const ghBtn = $("#sb-github");
if (ghBtn) ghBtn.addEventListener("click", () => { if (api().open_github) api().open_github(); });

function buildNav() {
  const nav = $("#nav"); nav.innerHTML = "";
  const hasRep = state.bundle.meta.has_repeatability;
  SECTIONS.forEach(s => {
    const b = document.createElement("button");
    b.className = "nav-item"; b.dataset.id = s.id;
    b.innerHTML = `<span class="idx">${s.idx}</span>${s.label}`;
    if (s.id === "repeatability" && !hasRep) b.disabled = true;
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

function render(id) {
  const body = $("#main-body");
  // tear down any existing plots (frees Plotly's window resize listeners) + the hover popup
  $$(".js-plotly-plot", body).forEach(d => { try { Plotly.purge(d); } catch (e) {} });
  if (structPop) structPop.style.display = "none";
  body.innerHTML = "";
  if (id === "overview") return renderOverview(body);
  if (id === "table") return renderTable(body);
  if (id === "export") return renderExport(body);
  // a single figure section
  const fig = state.bundle.figures[id];
  if (!fig) { body.innerHTML = `<div class="panel"><div class="panel-h">Not available</div></div>`; return; }
  const sec = SECTIONS.find(x => x.id === id);
  const div = plotInto(body, sec ? sec.title : "");
  drawPlot(div, fig);
  if (id === "detection" || id === "warp") attachStructures(div);
  const note = state.bundle.notes && state.bundle.notes[id];
  if (note) {
    const n = document.createElement("div");
    n.className = "fig-note"; n.innerHTML = note;
    div.parentElement.appendChild(n);   // inside the light card, under the plot
  }
}

function renderOverview(body) {
  const k = state.bundle.kpis;
  const grid = document.createElement("div"); grid.className = "kpis";
  grid.innerHTML = k.map(t =>
    `<div class="kpi"><div class="kpi-l">${t.label}</div><div class="kpi-v">${t.value}</div>` +
    `<div class="kpi-s">${t.sub}</div></div>`).join("");
  body.appendChild(grid);
  const panel = document.createElement("div"); panel.className = "panel";
  panel.innerHTML = `<div class="panel-h">Quality fingerprint — outward is better</div>`;
  const div = document.createElement("div"); div.className = "plot"; panel.appendChild(div);
  body.appendChild(panel);
  drawPlot(div, state.bundle.radar);
}

function renderExport(body) {
  body.innerHTML =
    `<div class="export-all">
       <button class="btn-exp big" id="exp-all"><div class="t">Export all outputs →</div>
         <div class="d">Writes everything into one folder: calibrated CSV · model.json · anchors.csv · log.txt · report (HTML + PDF) · run_info.json</div></button>
     </div>
     <div class="export-row">
       <button class="btn-exp" id="exp-csv"><div class="t">Calibrated table</div>
         <div class="d">CSV — original columns + RI, RI_uncertainty, RI_reliability, RI_spread, flags</div></button>
       <button class="btn-exp" id="exp-report"><div class="t">Report</div>
         <div class="d">Interactive HTML + static PDF (detection, warp, repeatability)</div></button>
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

/* ---- plotly helpers ---- */
function plotInto(body, title) {
  const panel = document.createElement("div"); panel.className = "panel";
  if (title) {
    const h = document.createElement("div"); h.className = "panel-h"; h.textContent = title;
    panel.appendChild(h);
  }
  const div = document.createElement("div"); div.className = "plot"; panel.appendChild(div);
  body.appendChild(panel); return div;
}
function drawPlot(div, figStr) {
  const f = typeof figStr === "string" ? JSON.parse(figStr) : figStr;
  Plotly.newPlot(div, f.data, f.layout, { displaylogo: false, responsive: true,
    modeBarButtonsToRemove: ["select2d", "lasso2d"], toImageButtonOptions: { format: "svg" } });
}

/* ---- molecular-structure hover ---- */
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
  if (!window.pywebview && window.__PREVIEW__ && window.__PREVIEW__.ok) {
    state.bundle = window.__PREVIEW__; buildNav(); openOutput(); selectSection("overview");
  }
});
