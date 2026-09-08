/* RT Anchor — chartkit: hand-authored minimal SVG charts (no Plotly).
   Slate-blue data palette mirroring the engine theme; source→shape,
   class→colour; clay is reserved for hover/active. One shared renderer:
     radar · cross-column curve.
   (Detection and Profile are the engine's own Plotly figures — see app.js.) */
"use strict";
const CK = (() => {
  const NS = "http://www.w3.org/2000/svg";
  /* palette mirrors the engine theme (rt_anchor.viz.theme) / plot-nature blue
     system: structure in soft grays, slate-blue data, clay for interaction */
  const C = {
    ink: "#2B2B2B", sub: "#4D4D4D", well: "#FFFFFF", card: "#FCFBF9",
    grid: "#E6E6E6", axis: "#7A7A7A", muted: "#7A7A7A", muted2: "#B0B0B0",
    blue: "#2B5D7D", navy: "#1B4A6B", teal: "#3E7C73", brick: "#A63D40", hot: "#C08552",
    steel: "#3D7CA8", pale: "#9CC3D5",
  };
  /* Windows ships none of the mac-first faces below, and Chromium/WebView2 only
     honours ui-monospace on macOS — without the Cascadia/Consolas steps the
     chart labels fell all the way back to Courier New. */
  const MONO = '"JetBrains Mono","IBM Plex Mono",ui-monospace,"SF Mono",Menlo,' +
               '"Cascadia Mono",Consolas,"Segoe UI Mono",monospace';
  /* Lipid classes -> colour. The stage-2 anchors are endogenous plasma lipids,
     so the palette has to cover more than the mixture's PC/DG/CE; the ordering
     keeps classes that neighbour each other in RT visually apart. */
  const CLASS_COLORS = {
    LPC: "#1B4A6B", PC: "#2B5D7D", SM: "#4E8FA6", DG: "#7FB2C4",
    CE: "#C08552", TG: "#8A5223", PE: "#A63D40", LPE: "#5A7D9A", CER: "#1F7A8C",
  };
  const classColor = c => CLASS_COLORS[String(c || "").toUpperCase()] || C.muted;
  const fmt = (v, d = 2) => (v == null || !isFinite(v)) ? "—" : (+v).toFixed(d);
  const pct = v => (v == null || !isFinite(v)) ? "—" : (100 * v).toFixed(0) + "%";
  const si = v => { v = +v; if (!isFinite(v)) return "—"; const a = Math.abs(v);
    return a >= 1e6 ? (v / 1e6).toFixed(1) + "M" : a >= 1e3 ? (v / 1e3).toFixed(1) + "k" : v.toFixed(0); };
  const esc = s => String(s == null ? "" : s).replace(/[&<>]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));

  const el = (tag, attrs = {}, kids = []) => {
    const n = document.createElementNS(NS, tag);
    for (const k in attrs) if (attrs[k] != null) n.setAttribute(k, attrs[k]);
    (Array.isArray(kids) ? kids : [kids]).forEach(k =>
      k != null && n.appendChild(typeof k === "string" ? document.createTextNode(k) : k));
    return n;
  };
  const scale = (d0, d1, r0, r1) => {
    const f = v => r0 + (v - d0) * (r1 - r0) / ((d1 - d0) || 1);
    f.inv = p => d0 + (p - r0) * ((d1 - d0) || 1) / ((r1 - r0) || 1);
    return f;
  };
  function niceTicks(min, max, n = 6) {
    if (!isFinite(min) || !isFinite(max) || min === max) return [min || 0];
    const step0 = (max - min) / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const nm = step0 / mag, step = (nm < 1.5 ? 1 : nm < 3 ? 2 : nm < 7 ? 5 : 10) * mag;
    const out = []; for (let v = Math.ceil(min / step) * step; v <= max + 1e-9; v += step) out.push(+v.toFixed(6));
    return out;
  }
  const txt = (x, y, s, o = {}) => el("text", {
    x, y, "font-family": MONO, "font-size": o.size || 10, "font-weight": o.w || 500,
    fill: o.fill || C.ink, "text-anchor": o.anchor || "middle", "letter-spacing": o.ls || "0",
    transform: o.rot ? `rotate(${o.rot} ${x} ${y})` : null, "text-rendering": "geometricPrecision",
  }, s);
  const square = (x, y, s, fill, stroke = C.ink, sw = 2) =>
    el("rect", { x: x - s / 2, y: y - s / 2, width: s, height: s, fill, stroke, "stroke-width": sw, "shape-rendering": "crispEdges" });
  const diamond = (x, y, s, fill, stroke = C.ink, sw = 2) =>
    el("rect", { x: x - s / 2, y: y - s / 2, width: s, height: s, fill, stroke, "stroke-width": sw, transform: `rotate(45 ${x} ${y})` });
  const ring = (x, y, r, stroke, sw = 2.4) =>
    el("circle", { cx: x, cy: y, r, fill: "none", stroke, "stroke-width": sw });
  const line = (x1, y1, x2, y2, stroke, sw = 1, dash = null) =>
    el("line", { x1, y1, x2, y2, stroke, "stroke-width": sw, "stroke-dasharray": dash, "shape-rendering": "crispEdges" });

  /* --- plot area: open frame — left + bottom spines only, soft gray --- */
  function well(g, x, y, w, h) {
    g.appendChild(el("rect", { x, y, width: w, height: h, fill: C.well }));
    g.appendChild(line(x, y, x, y + h, C.axis, 1.2));
    g.appendChild(line(x, y + h, x + w, y + h, C.axis, 1.2));
  }
  function xAxis(g, sx, y0, ticks, dp, title) {
    ticks.forEach(t => { const x = sx(t);
      g.appendChild(line(x, y0, x, y0 + 5, C.axis, 1.2));
      g.appendChild(txt(x, y0 + 17, fmt(t, dp), { size: 10, fill: C.axis })); });
    if (title) g.appendChild(txt((sx.x0 + sx.x1) / 2, y0 + 33, title, { size: 11, w: 600, ls: ".02em", fill: C.sub }));
  }
  function yAxis(g, sy, x0, ticks, dp, title) {
    ticks.forEach(t => { const y = sy(t);
      g.appendChild(line(x0 - 5, y, x0, y, C.axis, 1.2));
      g.appendChild(txt(x0 - 9, y + 3.5, fmt(t, dp), { size: 10, anchor: "end", fill: C.axis })); });
    if (title) g.appendChild(txt(x0 - 42, (sy.y0 + sy.y1) / 2, title, { size: 11, w: 600, ls: ".02em", fill: C.sub, rot: -90 }));
  }
  /* whisper gridlines — horizontal only (verticals are chartjunk here) */
  function gridY(g, sy, x0, x1, ticks) {
    ticks.forEach(t => g.appendChild(line(x0, sy(t), x1, sy(t), C.grid, 1)));
  }
  const svgRoot = (w, h) => el("svg", { viewBox: `0 0 ${w} ${h}`, width: "100%",
    style: "display:block", preserveAspectRatio: "xMidYMid meet" });

  /* legend strip / verdict plate rendered as HTML above the SVG */
  function legendStrip(host, items) {
    const d = document.createElement("div"); d.className = "ck-legend";
    d.innerHTML = items.map(it => {
      const g = it.glyph || "square";
      const sw = g === "line"
        ? `<span class="ck-lg-line" style="background:${it.color}"></span>`
        : g === "dash"
        ? `<span class="ck-lg-line dash" style="background:${it.color}"></span>`
        : g === "ring"
        ? `<span class="ck-lg-ring" style="border-color:${it.color}"></span>`
        : g === "diamond"
        ? `<span class="ck-lg-dia" style="background:${it.fill || it.color};border-color:${it.color}"></span>`
        : `<span class="ck-lg-sq" style="background:${it.fill === undefined ? it.color : it.fill};border-color:${it.color}"></span>`;
      return `<span class="ck-lg">${sw}${esc(it.label)}</span>`;
    }).join("");
    host.appendChild(d);
  }
  function plate(host, cls, html) {
    const d = document.createElement("div"); d.className = "ck-verdict " + (cls || "");
    d.innerHTML = html; host.appendChild(d);
    return d;
  }

  /* ---------- floating overlays (tooltip + structure card) ---------- */
  let tip, card;
  function ensureOverlays() {
    if (!tip) {
      tip = document.createElement("div"); tip.className = "ck-tip"; document.body.appendChild(tip);
      card = document.createElement("div"); card.className = "ck-card"; document.body.appendChild(card);
    }
  }
  function place(node, x, y) {
    const w = node.offsetWidth, h = node.offsetHeight;
    node.style.left = Math.min(x + 14, window.innerWidth - w - 8) + "px";
    node.style.top = Math.min(y + 14, window.innerHeight - h - 8) + "px";
  }
  function showTip(title, cls, rows, x, y) {
    ensureOverlays();
    const swatch = cls ? `<span class="ck-sw" style="background:${classColor(cls)}"></span>` : "";
    tip.innerHTML = `<div class="ck-tip-h">${swatch}${esc(title)}</div>` +
      rows.map(r => `<div class="ck-tip-r"><span>${esc(r[0])}</span><b>${esc(r[1])}</b></div>`).join("");
    tip.style.display = "block"; place(tip, x, y);
  }
  const hideTip = () => { if (tip) tip.style.display = "none"; };
  let STRUCT = {};
  function showCard(std, x, y) {
    ensureOverlays();
    const svg = STRUCT[std.name];
    if (!svg) return;
    card.innerHTML = `<div class="ck-card-h">${esc(std.name)}<span class="ck-chip" style="border-color:${classColor(std.class)};color:${classColor(std.class)}">${esc(std.class || "")}</span></div>` +
      `<div class="ck-card-b">${svg}</div>`;
    card.style.display = "block"; place(card, x, y + 120);
  }
  const hideCard = () => { if (card) card.style.display = "none"; };
  const hideAll = () => { hideTip(); hideCard(); };

  /* ============================ CHARTS ============================ */

  // 1 — RADAR (pentagon)
  function radar(host, data) {
    const W = 460, H = 380, cx = W / 2, cy = H / 2 + 6, R = 132;
    const s = svgRoot(W, H), g = el("g"); s.appendChild(g);
    s.setAttribute("style", "display:block;max-width:560px;margin:0 auto");
    const n = data.length, ang = i => -Math.PI / 2 + i * 2 * Math.PI / n;
    const pt = (i, r) => [cx + Math.cos(ang(i)) * R * r, cy + Math.sin(ang(i)) * R * r];
    [0.25, 0.5, 0.75, 1].forEach(rr => {
      const p = data.map((_, i) => pt(i, rr).join(",")).join(" ");
      g.appendChild(el("polygon", { points: p, fill: "none", stroke: rr === 1 ? C.muted2 : C.grid, "stroke-width": rr === 1 ? 1.2 : 1 }));
    });
    data.forEach((_, i) => { const [x, y] = pt(i, 1); g.appendChild(line(cx, cy, x, y, C.grid, 1)); });
    const poly = data.map((d, i) => pt(i, Math.max(0.02, d.value)).join(",")).join(" ");
    g.appendChild(el("polygon", { points: poly, fill: C.blue, "fill-opacity": .12, stroke: C.blue, "stroke-width": 2.2 }));
    g.appendChild(square(cx, cy, 6, C.ink));
    data.forEach((d, i) => {
      const [vx, vy] = pt(i, Math.max(0.02, d.value)), sq = square(vx, vy, 9, C.blue);
      const [lx, ly] = pt(i, 1.16);
      g.appendChild(txt(lx, ly, d.axis.toUpperCase(), { size: 10, w: 600, fill: C.sub }));
      g.appendChild(txt(lx, ly + 13, fmt(d.value, 2), { size: 9, fill: C.muted }));
      sq.style.cursor = "pointer";
      sq.addEventListener("mousemove", e => { sq.setAttribute("fill", C.hot);
        showTip(d.axis.toUpperCase(), null, [["value", fmt(d.value, 2)]], e.clientX, e.clientY); });
      sq.addEventListener("mouseleave", () => { sq.setAttribute("fill", C.blue); hideTip(); });
      g.appendChild(sq);
    });
    host.appendChild(s);
  }

  // 2 — CROSS-COLUMN CURVE (matched pairs + fit + anchors, residual strip below)
  function curve(host, data) {
    if (!data || data.empty) {
      host.innerHTML = `<div class="ck-plate">NO CROSS-COLUMN CURVE · ${esc((data && data.reason) || "not available")}</div>`;
      return;
    }
    legendStrip(host, [
      { label: "standards pairs", glyph: "square", color: C.ink, fill: C.blue },
      { label: "sample pairs", glyph: "diamond", color: C.ink, fill: C.pale },
      { label: "stage-1 curve", glyph: "line", color: C.navy },
    ].concat(data.refined ? [{ label: "anchor-refined", glyph: "dash", color: C.teal }] : [])
     .concat(data.anchors && data.anchors.length ? [{ label: "stage-2 anchors", glyph: "ring", color: C.brick }] : []));

    const W = 680, L = 66, R = 26, T = 20, mainH = 272, gap = 20, resH = 96, B = 46;
    const H = T + mainH + gap + resH + B;
    const s = svgRoot(W, H), g = el("g"); s.appendChild(g);

    const P = data.pairs;
    const xs = P.x.filter(v => v != null), ys = P.y.filter(v => v != null);
    const xlo = Math.min.apply(null, xs), xhi = Math.max.apply(null, xs);
    const ylo = Math.min.apply(null, ys), yhi = Math.max.apply(null, ys);
    const xpad = (xhi - xlo) * 0.03 || 0.5, ypad = (yhi - ylo) * 0.03 || 0.5;
    const sx = scale(xlo - xpad, xhi + xpad, L, W - R); sx.x0 = L; sx.x1 = W - R;
    const sy = scale(ylo - ypad, yhi + ypad, T + mainH, T); sy.y0 = T; sy.y1 = T + mainH;
    well(g, L, T, W - R - L, mainH);
    const xt = niceTicks(xlo, xhi, 6), yt = niceTicks(ylo, yhi, 5);
    gridY(g, sy, L, W - R, yt);
    yAxis(g, sy, L, yt, 1, "RT on the reference column (min)");

    // matched pairs — MAD-trimmed pairs are excluded from the figure
    const gStd = el("g"), gSamp = el("g");
    for (let i = 0; i < P.x.length; i++) {
      const x = P.x[i], y = P.y[i];
      if (x == null || y == null || !P.kept[i]) continue;
      const px = sx(x), py = sy(y);
      if (P.std[i]) gStd.appendChild(square(px, py, 5.5, C.blue, C.navy, 1));
      else gSamp.appendChild(diamond(px, py, 5.5, C.pale, C.steel, 1));
    }
    g.appendChild(gSamp); g.appendChild(gStd);

    const poly = (pts, stroke, sw, dash) => {
      const p = pts.filter(q => q).join(" ");
      if (p) g.appendChild(el("polyline", { points: p, fill: "none", stroke, "stroke-width": sw, "stroke-dasharray": dash }));
    };
    if (data.refined) {
      poly(data.refined.x.map((x, i) => (x == null || data.refined.y[i] == null) ? null : `${sx(x)},${sy(data.refined.y[i])}`),
           C.teal, 2.5, "7 4");
    }
    poly(data.fit.x.map((x, i) => (x == null || data.fit.y[i] == null) ? null : `${sx(x)},${sy(data.fit.y[i])}`),
         C.navy, 3);

    // stage-2 anchors as rings, with a white halo so they read over dense scatter
    (data.anchors || []).forEach(a => {
      if (a.x == null || a.y == null) return;
      const px = sx(a.x), py = sy(a.y);
      g.appendChild(ring(px, py, 7, "#FFFFFF", 5));
      const r = ring(px, py, 7, C.brick, 2.4); r.style.cursor = "pointer";
      g.appendChild(r);
      const hit = el("circle", { cx: px, cy: py, r: 10, fill: "transparent", style: "cursor:pointer" });
      hit.addEventListener("mousemove", e => {
        r.setAttribute("stroke", C.hot);
        showTip(a.label || "anchor", a.class,
          [["your RT", fmt(a.x, 2)], ["reference RT", fmt(a.y, 2)],
           ["resid vs curve", fmt(a.resid, 3)], ["leave-one-out", fmt(a.loo, 3)]],
          e.clientX, e.clientY);
      });
      hit.addEventListener("mouseleave", () => { r.setAttribute("stroke", C.brick); hideAll(); });
      g.appendChild(hit);
    });

    // hover crosshair — live "your RT -> reference RT" probe
    const gx = data.fit.x, gy = data.fit.y;
    const interp = rt => {
      if (rt <= gx[0]) return gy[0];
      if (rt >= gx[gx.length - 1]) return gy[gy.length - 1];
      let i = 1; while (i < gx.length && gx[i] < rt) i++;
      const t = (rt - gx[i - 1]) / ((gx[i] - gx[i - 1]) || 1); return gy[i - 1] + t * (gy[i] - gy[i - 1]);
    };
    const cross = el("g", { style: "display:none", "pointer-events": "none" });
    const vg = line(0, T, 0, T + mainH, C.hot, 1.5, "4 3"), hg = line(L, 0, W - R, 0, C.hot, 1.5, "4 3");
    const dot = square(0, 0, 9, C.hot, C.ink, 1.5);
    const cxr = el("rect", { fill: C.ink, "shape-rendering": "crispEdges" }), cxt = txt(0, 0, "", { size: 9, w: 700, fill: C.card });
    const cyr = el("rect", { fill: C.ink, "shape-rendering": "crispEdges" }), cyt = txt(0, 0, "", { size: 9, w: 700, fill: C.card, anchor: "end" });
    [vg, hg, dot, cxr, cxt, cyr, cyt].forEach(n => cross.appendChild(n));
    g.appendChild(cross);
    const pt = s.createSVGPoint();
    const wellHit = el("rect", { x: L, y: T, width: W - R - L, height: mainH, fill: "transparent", style: "cursor:crosshair" });
    wellHit.addEventListener("mousemove", e => {
      pt.x = e.clientX; pt.y = e.clientY; const p = pt.matrixTransform(s.getScreenCTM().inverse());
      const rt = Math.max(gx[0], Math.min(gx[gx.length - 1], sx.inv(p.x))), iv = interp(rt);
      const px = sx(rt), py = sy(iv);
      cross.style.display = "";
      vg.setAttribute("x1", px); vg.setAttribute("x2", px);
      hg.setAttribute("y1", py); hg.setAttribute("y2", py);
      dot.setAttribute("x", px - 4.5); dot.setAttribute("y", py - 4.5);
      cxt.textContent = fmt(rt, 2); cxt.setAttribute("x", px); cxt.setAttribute("y", T + mainH + 13);
      cxr.setAttribute("x", px - 21); cxr.setAttribute("y", T + mainH + 3); cxr.setAttribute("width", 42); cxr.setAttribute("height", 15);
      cyt.textContent = fmt(iv, 2); cyt.setAttribute("x", L - 12); cyt.setAttribute("y", py + 3.5);
      cyr.setAttribute("x", L - 50); cyr.setAttribute("y", py - 8); cyr.setAttribute("width", 40); cyr.setAttribute("height", 15);
      showTip("PROBE", null, [["your RT", fmt(rt, 2)], ["reference RT", fmt(iv, 2)]], e.clientX, e.clientY);
    });
    wellHit.addEventListener("mouseleave", () => { cross.style.display = "none"; hideTip(); });
    g.appendChild(wellHit);

    // ---- residual strip ----
    const ry0 = T + mainH + gap;
    well(g, L, ry0, W - R - L, resH);
    const rall = P.r.filter(v => v != null && isFinite(v))
      .concat((P.rr || []).filter(v => v != null && isFinite(v)));
    const srt = rall.slice().sort((a, b) => a - b);
    const q = f => srt.length ? srt[Math.min(srt.length - 1, Math.max(0, Math.round(f * (srt.length - 1))))] : 0;
    const rlim = Math.max(Math.abs(q(0.01)), Math.abs(q(0.99)), 0.05) * 1.25;
    const sr = scale(-rlim, rlim, ry0 + resH - 5, ry0 + 5); sr.y0 = ry0; sr.y1 = ry0 + resH;
    const rt2 = niceTicks(-rlim, rlim, 4);
    rt2.forEach(t => g.appendChild(line(L, sr(t), W - R, sr(t), C.grid, 1)));
    g.appendChild(line(L, sr(0), W - R, sr(0), C.axis, 1.2));
    yAxis(g, sr, L, rt2, 2, "Resid (min)");
    const gr1 = el("g"), gr2 = el("g");
    for (let i = 0; i < P.x.length; i++) {
      const x = P.x[i]; if (x == null || !P.kept[i]) continue;
      const px = sx(x);
      const v = P.r[i];
      if (v != null && isFinite(v)) {
        const y = sr(Math.max(-rlim, Math.min(rlim, v)));
        gr1.appendChild(el("rect", { x: px - 1.8, y: y - 1.8, width: 3.6, height: 3.6,
          fill: C.muted2, "shape-rendering": "crispEdges" }));
      }
      if (P.rr) {
        const w = P.rr[i];
        if (w != null && isFinite(w)) {
          const y = sr(Math.max(-rlim, Math.min(rlim, w)));
          gr2.appendChild(el("rect", { x: px - 1.8, y: y - 1.8, width: 3.6, height: 3.6,
            fill: C.blue, "fill-opacity": .8, "shape-rendering": "crispEdges" }));
        }
      }
    }
    g.appendChild(gr1); g.appendChild(gr2);
    g.appendChild(txt(W - R - 6, ry0 + 13, P.rr ? "GREY = STAGE 1 · BLUE = REFINED" : "STAGE-1 RESIDUALS",
      { size: 8.5, w: 700, anchor: "end", fill: C.muted, ls: ".06em" }));
    xAxis(g, sx, ry0 + resH, xt, 1, "RT on your column (min)");
    host.appendChild(s);
  }

  function render(id, host, bundle) {
    STRUCT = bundle.structures || {};
    host.innerHTML = "";
    if (id === "curve") return curve(host, bundle.curve);
  }
  return { render, radar, C, classColor, hideAll };
})();
