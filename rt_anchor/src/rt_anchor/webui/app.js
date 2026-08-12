/* RT Anchor desktop client — front-end.
   Talks to the Python side via pywebview.api.*; renders plots with Plotly. */
(function () {
  "use strict";
  const $ = (id) => document.getElementById(id);
  const api = () => (window.pywebview && window.pywebview.api) || null;

  let hasResult = false;
  const loaded = {};                 // section -> true once its plot/data is drawn
  const PLOT_CFG = { responsive: true, displaylogo: false,
                     toImageButtonOptions: { format: "svg" } };

  // ---------------------------------------------------------------- views ----
  function showView(sec) {
    document.querySelectorAll(".nav").forEach(n =>
      n.classList.toggle("active", n.dataset.sec === sec));
    document.querySelectorAll(".view").forEach(v =>
      v.classList.toggle("on", v.dataset.view === sec));
    if (hasResult) loadSection(sec);
  }

  async function loadSection(sec) {
    if (sec === "overview") { await ensureOverview(); return; }
    if (sec === "table") { await ensureTable(); return; }
    if (sec === "export") return;
    // plot sections: detection / profile / warp / repeatability
    if (loaded[sec]) { resize("fig-" + sec); return; }
    const fig = await api().get_figure(sec);
    if (fig && fig.data) {
      Plotly.newPlot("fig-" + sec, fig.data, fig.layout, PLOT_CFG);
      loaded[sec] = true;
    }
    if (sec === "detection" && fig && fig.note) $("note-detection").innerHTML = fig.note;
  }

  function resize(id) { const el = $(id); if (el && el.data) Plotly.Plots.resize(el); }

  // ---------------------------------------------------------------- overview -
  async function ensureOverview() {
    if (loaded.overview) { resize("fig-radar"); return; }
    const ov = await api().get_overview();
    renderKpis(ov.tiles);
    if (ov.radar && ov.radar.data)
      Plotly.newPlot("fig-radar", ov.radar.data, ov.radar.layout, PLOT_CFG);
    loaded.overview = true;
  }

  function renderKpis(tiles) {
    $("kpis").innerHTML = (tiles || []).map(t =>
      `<div class="kpi"><div class="l">${esc(t[0])}</div>` +
      `<div class="v">${esc(t[1])}</div><div class="s">${esc(t[2])}</div></div>`).join("");
  }

  // ---------------------------------------------------------------- table ----
  async function ensureTable() {
    if (loaded.table) return;
    const t = await api().get_table();
    $("table-caption").textContent =
      `CALIBRATED FEATURES + RI COLUMNS — showing ${t.rows.length} of ${t.total.toLocaleString()}`;
    const head = "<tr>" + t.columns.map(c => `<th>${esc(c)}</th>`).join("") + "</tr>";
    const body = t.rows.map(r => "<tr>" + r.map(c => `<td>${esc(c)}</td>`).join("") + "</tr>").join("");
    $("table-wrap").innerHTML = `<table class="dt"><thead>${head}</thead><tbody>${body}</tbody></table>`;
    loaded.table = true;
  }

  // ---------------------------------------------------------------- modal ----
  function openModal() { $("modal").classList.add("on"); }
  function closeModal() { $("modal").classList.remove("on"); }

  async function browse(inputId, kind) {
    const a = api(); if (!a) return;
    const p = kind === "dir" ? await a.pick_folder() : await a.pick_file();
    if (p) $(inputId).value = p;
    if (inputId === "f-sample" && p) {
      const meta = await a.peek_meta(p);
      if (meta && meta.polarity) $("f-polarity").value = meta.polarity;
      if (meta && meta.mz) $("f-mz").value = meta.mz;
      if (meta && meta.rtw) $("f-rtw").value = meta.rtw;
      if (!$("ex-dir").value && meta && meta.dir) $("ex-dir").value = meta.dir;
    }
  }

  function preset(kind) {
    if (kind === "qtof") { $("f-mz").value = "15.0"; $("f-rtw").value = "0.5"; }
    else { $("f-mz").value = "8.0"; $("f-rtw").value = "0.3"; }
  }

  async function runCalibration() {
    const a = api(); if (!a) return;
    const params = {
      sample: $("f-sample").value, standards: $("f-standards").value,
      polarity: $("f-polarity").value,
      mz_tol: $("f-mz").value, rt_window: $("f-rtw").value, min_anchors: $("f-minanc").value,
      extrapolate: $("f-extrap").checked,
      source_format: $("f-fmt").value, rt_unit: $("f-unit").value,
      single: $("f-single").value, manifest: $("f-manifest").value, reference: $("f-reference").value,
    };
    if (!params.sample || !params.standards) { logRun("Sample table and standards run are required."); return; }
    showLoading("Calibrating…");
    $("run-log").textContent = "";
    const res = await a.run_calibration(params);
    hideLoading();
    (res.log || []).forEach(logRun);
    if (!res.ok) { logRun("ERROR: " + res.error); return; }
    onResult(res);
    closeModal();
    showView("overview");
  }

  function onResult(res) {
    hasResult = true;
    for (const k in loaded) delete loaded[k];
    $("empty").classList.add("hidden");
    // repeatability availability
    const rep = document.querySelector('.nav[data-sec="repeatability"]');
    if (res.repeatability) rep.classList.remove("disabled");
    else rep.classList.add("disabled");
    if (res.out_dir) $("ex-dir").value = res.out_dir;
    if (res.prefix) $("ex-prefix").value = res.prefix;
  }

  async function runExample() {
    const a = api(); if (!a) return;
    showLoading("Running example dataset…");
    const res = await a.run_example();
    hideLoading();
    if (!res.ok) { alert("Example failed:\n" + res.error); return; }
    onResult(res);
    showView("overview");
  }

  // ---------------------------------------------------------------- export ---
  async function runExport() {
    const a = api(); if (!a) return;
    const params = {
      out_dir: $("ex-dir").value, prefix: $("ex-prefix").value,
      html: $("ex-html").checked, pdf: $("ex-pdf").checked, tic: $("ex-tic").value,
    };
    if (!params.out_dir) { logExport("Choose an output folder."); return; }
    showLoading("Writing outputs…");
    const res = await a.export_outputs(params);
    hideLoading();
    $("ex-log").textContent = "";
    if (!res.ok) { logExport("ERROR: " + res.error); return; }
    Object.entries(res.paths).forEach(([k, v]) => logExport(k + ": " + v));
    $("ex-report").disabled = !res.paths.report_html;
    $("ex-folder").disabled = false;
  }

  // ---------------------------------------------------------------- utils ----
  function logRun(m) { $("run-log").textContent += m + "\n"; }
  function logExport(m) { $("ex-log").textContent += m + "\n"; }
  function showLoading(t) { $("loading-txt").textContent = t || "Working…"; $("loading").classList.add("on"); }
  function hideLoading() { $("loading").classList.remove("on"); }
  function esc(s) { return String(s == null ? "" : s).replace(/[&<>]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c])); }

  // ---------------------------------------------------------------- wire ----
  let wired = false;
  function wire() {
    if (wired) return; wired = true;
    $("nav").addEventListener("click", e => {
      const b = e.target.closest(".nav");
      if (!b || b.classList.contains("disabled")) return;
      if (!hasResult) { openModal(); return; }
      showView(b.dataset.sec);
    });
    $("btn-new").addEventListener("click", openModal);
    $("empty-new").addEventListener("click", openModal);
    $("btn-example").addEventListener("click", runExample);
    $("empty-example").addEventListener("click", runExample);
    $("modal-x").addEventListener("click", closeModal);
    $("btn-github").addEventListener("click", () => { const a = api(); if (a) a.open_url("https://github.com/Bowen999/rt-anchor"); });
    $("adv-toggle").addEventListener("click", () => {
      const adv = $("adv"); adv.hidden = !adv.hidden;
      $("adv-toggle").textContent = (adv.hidden ? "▸" : "▾") + " ADVANCED OPTIONS";
    });
    document.querySelectorAll("[data-browse]").forEach(b =>
      b.addEventListener("click", () => browse(b.dataset.browse, b.dataset.kind)));
    document.querySelectorAll("[data-preset]").forEach(b =>
      b.addEventListener("click", () => preset(b.dataset.preset)));
    $("run-btn").addEventListener("click", runCalibration);
    $("ex-run").addEventListener("click", runExport);
    $("ex-report").addEventListener("click", () => { const a = api(); if (a) a.open_report(); });
    $("ex-folder").addEventListener("click", () => { const a = api(); if (a) a.open_folder(); });
    showView("overview");
  }

  // pywebview injects its api asynchronously
  if (window.pywebview && window.pywebview.api) { wire(); bootstrap(); }
  else window.addEventListener("pywebviewready", () => { wire(); bootstrap(); });
  // also wire immediately so the shell is interactive even before the bridge is ready
  window.addEventListener("DOMContentLoaded", () => { try { wire(); } catch (e) {} });

  async function bootstrap() {
    const a = api(); if (!a) return;
    const st = await a.get_state();
    if (st && !st.has_example) {
      ["btn-example", "empty-example"].forEach(id => { const el = $(id); if (el) el.style.display = "none"; });
    }
    if (st && st.has_result) { onResult(st); showView("overview"); }
    // otherwise the empty state (RUN EXAMPLE / NEW CALIBRATION) stays visible
  }
})();
