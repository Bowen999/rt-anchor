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
    standards: "/preview/example/standards.txt", polarity: "positive" }),
  run_calibration: async () => window.__PREVIEW__ || { ok: false, error: "no preview data" },
  export_csv: async () => ({ ok: true, path: "(preview)" }),
  export_report: async () => ({ ok: true, paths: {} }),
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

$$(".btn-file").forEach(btn => btn.addEventListener("click", async () => {
  const target = btn.dataset.target, kind = btn.dataset.pick;
  const path = kind === "folder" ? await api().pick_folder() : await api().pick_file();
  if (!path) return;
  state.params[target] = path;
  const el = $("#path-" + target);
  el.textContent = shortPath(path); el.title = path; el.classList.add("set");
}));

$$(".seg").forEach(seg => seg.addEventListener("click", e => {
  const b = e.target.closest("button"); if (!b) return;
  $$("button", seg).forEach(x => x.classList.remove("active"));
  b.classList.add("active"); state.params[seg.dataset.seg] = b.dataset.val;
}));

/* load the bundled demo dataset into the Required fields */
function setPath(target, path) {
  state.params[target] = path;
  const el = $("#path-" + target);
  el.textContent = shortPath(path); el.title = path; el.classList.add("set");
}
const loadExBtn = $("#load-example");
if (loadExBtn) loadExBtn.addEventListener("click", async () => {
  const r = await api().load_example();
  if (!r || !r.ok) { setStatus((r && r.error) || "Could not load the example.", "err"); return; }
  setPath("samples", r.samples);
  setPath("standards", r.standards);
  state.params.polarity = r.polarity || "positive";
  const seg = $('[data-seg="polarity"]');
  if (seg) $$("button", seg).forEach(b => b.classList.toggle("active", b.dataset.val === state.params.polarity));
  setStatus("Example dataset loaded — press Run calibration.", "busy");
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
  setStatus("Calibrating — identifying anchors, fitting warp…", "busy");
  let res;
  try { res = await api().run_calibration(state.params); }
  catch (e) { res = { ok: false, error: String(e) }; }
  btn.disabled = false;
  if (!res || !res.ok) { setStatus((res && res.error) || "Calibration failed.", "err"); return; }
  setStatus("");
  state.bundle = res;
  buildNav(); openOutput(); selectSection("overview");
}
function setStatus(msg, cls) { const s = $("#status"); s.textContent = msg; s.className = "status" + (cls ? " " + cls : ""); }

/* ===================== OUTPUT VIEW ===================== */
function openOutput() { $("#view-input").classList.add("hidden"); $("#view-output").classList.remove("hidden"); }
function openInput() { $("#view-output").classList.add("hidden"); $("#view-input").classList.remove("hidden"); }
$("#new-run").addEventListener("click", openInput);

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
  const m = state.bundle.meta;
  $("#main-meta").innerHTML =
    `${cap(m.source_format)} · ${m.polarity} &nbsp;|&nbsp; <b>${m.scope}</b> calibration ` +
    `&nbsp;|&nbsp; ${m.n_detected}/${m.n_panel} standards &nbsp;|&nbsp; ${Number(m.n_features).toLocaleString()} features`;
  render(id);
}
const cap = s => s ? s[0].toUpperCase() + s.slice(1) : s;

function render(id) {
  const body = $("#main-body");
  // tear down any existing plots (frees Plotly's window resize listeners) + the hover popup
  $$(".js-plotly-plot", body).forEach(d => { try { Plotly.purge(d); } catch (e) {} });
  if (structPop) structPop.style.display = "none";
  body.innerHTML = "";
  if (id === "overview") return renderOverview(body);
  if (id === "export") return renderExport(body);
  // a single figure section
  const fig = state.bundle.figures[id];
  if (!fig) { body.innerHTML = `<div class="panel"><div class="panel-h">Not available</div></div>`; return; }
  const div = plotInto(body, "");
  drawPlot(div, fig);
  if (id === "detection" || id === "warp") attachStructures(div);
  const note = state.bundle.notes && state.bundle.notes[id];
  if (note) {
    const n = document.createElement("div");
    n.className = "fig-note"; n.innerHTML = note; body.appendChild(n);
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
    `<div class="export-row">
       <button class="btn-exp" id="exp-csv"><div class="t">Calibrated table (CSV)</div>
         <div class="d">Original columns + RI, RI_uncertainty, RI_reliability, RI_spread, flags</div></button>
       <button class="btn-exp" id="exp-report"><div class="t">Report (HTML + PDF)</div>
         <div class="d">The full interactive HTML report and a static PDF</div></button>
     </div>`;
  $("#exp-csv").addEventListener("click", async () => {
    const r = await api().export_csv();
    if (!r.ok && r.error === "cancelled") return;               // user dismissed the save dialog
    toast(r.ok ? "Saved: " + shortPath(r.path) : (r.error || "export failed"), !r.ok);
  });
  $("#exp-report").addEventListener("click", async () => {
    toast("Writing report…");
    const r = await api().export_report();
    if (!r.ok && r.error === "cancelled") { hideToast(); return; }
    toast(r.ok ? "Report written" : (r.error || "export failed"), !r.ok);
  });
}

/* ---- plotly helpers ---- */
function plotInto(body) {
  const panel = document.createElement("div"); panel.className = "panel";
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
    structPop.style.cssText = "position:fixed;display:none;z-index:1000;background:#12161D;" +
      "border:2px solid #333D4A;padding:8px 10px;box-shadow:5px 5px 0 rgba(0,0,0,.5);pointer-events:none;";
    document.body.appendChild(structPop);
    document.addEventListener("mousemove", e => { structPop._x = e.clientX; structPop._y = e.clientY; });
  }
  div.on("plotly_hover", d => {
    let nm = d.points && d.points[0] && d.points[0].customdata;
    if (Array.isArray(nm)) nm = nm[0];
    if (nm && S[nm]) {
      structPop.innerHTML = `<div style="font:600 12px sans-serif;color:#EDF1F6;margin-bottom:3px">${nm}</div>` +
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
