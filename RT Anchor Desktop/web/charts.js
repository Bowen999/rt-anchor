/* RT Anchor — chartkit: hand-authored brutalist SVG charts (no Plotly).
   Cool-only palette; every datum is an axis-aligned square; source→shape,
   class→colour; cyan is reserved for hover/active. One shared renderer, 5 charts. */
"use strict";
const CK = (() => {
  const NS = "http://www.w3.org/2000/svg";
  const C = {
    ink: "#141210", well: "#FFFFFF", card: "#FCFBF8", band: "#F1EDE4",
    grid: "#CED8DC", gridMajor: "#A7B4BE", muted: "#6B7C88", muted2: "#A9B6BC",
    blue: "#1668C0", navy: "#0B2E4F", teal: "#0E8F8A", indigo: "#4C63B6", cyan: "#17C9E6",
  };
  const MONO = '"JetBrains Mono","IBM Plex Mono",ui-monospace,"SF Mono",Menlo,monospace';
  const classColor = c => ({ PC: C.blue, DG: C.teal, CE: C.indigo }[c] || C.navy);
  const fmt = (v, d = 2) => (v == null || !isFinite(v)) ? "—" : (+v).toFixed(d);
  const si = v => { v = +v; if (!isFinite(v)) return "—"; const a = Math.abs(v);
    return a >= 1e6 ? (v / 1e6).toFixed(1) + "M" : a >= 1e3 ? (v / 1e3).toFixed(1) + "k" : v.toFixed(0); };

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
  const line = (x1, y1, x2, y2, stroke, sw = 1, dash = null) =>
    el("line", { x1, y1, x2, y2, stroke, "stroke-width": sw, "stroke-dasharray": dash, "shape-rendering": "crispEdges" });

  /* --- plot well: white rect + 3px ink frame + double baseline + corner Ls --- */
  function well(g, x, y, w, h) {
    g.appendChild(el("rect", { x, y, width: w, height: h, fill: C.well, stroke: C.ink, "stroke-width": 3, "shape-rendering": "crispEdges" }));
    g.appendChild(line(x + 4, y + 4, x + 4, y + h - 4, C.ink, 3));          // inner left
    g.appendChild(line(x + 4, y + h - 4, x + w - 4, y + h - 4, C.ink, 3));  // inner bottom
    const L = 8;
    [[x, y, 1, 1], [x + w, y, -1, 1], [x, y + h, 1, -1], [x + w, y + h, -1, -1]].forEach(([cx, cy, sx, sy]) => {
      g.appendChild(line(cx, cy, cx + sx * L, cy, C.ink, 2));
      g.appendChild(line(cx, cy, cx, cy + sy * L, C.ink, 2));
    });
  }
  function xAxis(g, sx, y0, ticks, dp, title) {
    ticks.forEach(t => { const x = sx(t);
      g.appendChild(line(x, y0, x, y0 + 6, C.ink, 2));
      g.appendChild(txt(x, y0 + 18, fmt(t, dp), { size: 10 })); });
    if (title) g.appendChild(txt((sx.x0 + sx.x1) / 2, y0 + 34, title, { size: 11, w: 700, ls: ".08em" }));
  }
  function yAxis(g, sy, x0, ticks, dp, title) {
    ticks.forEach(t => { const y = sy(t);
      g.appendChild(line(x0 - 6, y, x0, y, C.ink, 2));
      g.appendChild(txt(x0 - 10, y + 3.5, fmt(t, dp), { size: 10, anchor: "end" })); });
    if (title) g.appendChild(txt(x0 - 40, (sy.y0 + sy.y1) / 2, title, { size: 11, w: 700, ls: ".08em", rot: -90 }));
  }
  function gridY(g, sy, x0, x1, ticks) {
    ticks.forEach(t => g.appendChild(line(x0, sy(t), x1, sy(t), C.grid, 1)));
  }
  const svgRoot = (w, h) => el("svg", { viewBox: `0 0 ${w} ${h}`, width: "100%",
    style: "display:block", preserveAspectRatio: "xMidYMid meet" });

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
    tip.innerHTML = `<div class="ck-tip-h">${swatch}${title}</div>` +
      rows.map(r => `<div class="ck-tip-r"><span>${r[0]}</span><b>${r[1]}</b></div>`).join("");
    tip.style.display = "block"; place(tip, x, y);
  }
  const hideTip = () => { if (tip) tip.style.display = "none"; };
  let STRUCT = {};
  function showCard(std, x, y) {
    ensureOverlays();
    const svg = STRUCT[std.name];
    const body = svg ? `<div class="ck-card-b">${svg}</div>`
      : `<div class="ck-card-mono">${(std.class || "?")}</div>`;
    card.innerHTML = `<div class="ck-card-h">${std.name}<span class="ck-chip" style="border-color:${classColor(std.class)};color:${classColor(std.class)}">${std.class || ""}</span></div>` +
      body + `<div class="ck-card-f">iRT <b>${fmt(std.irt, 1)}</b></div>`;
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
      g.appendChild(el("polygon", { points: p, fill: "none", stroke: rr === 1 ? C.ink : C.grid, "stroke-width": rr === 1 ? 2 : 1 }));
    });
    data.forEach((_, i) => { const [x, y] = pt(i, 1); g.appendChild(line(cx, cy, x, y, C.grid, 1)); });
    const poly = data.map((d, i) => pt(i, Math.max(0.02, d.value)).join(",")).join(" ");
    g.appendChild(el("polygon", { points: poly, fill: C.blue, "fill-opacity": .18, stroke: C.blue, "stroke-width": 3 }));
    g.appendChild(square(cx, cy, 6, C.ink));
    data.forEach((d, i) => {
      const [vx, vy] = pt(i, Math.max(0.02, d.value)), sq = square(vx, vy, 9, C.blue);
      const [lx, ly] = pt(i, 1.16);
      g.appendChild(txt(lx, ly, d.axis.toUpperCase(), { size: 10, w: 700 }));
      g.appendChild(txt(lx, ly + 13, fmt(d.value, 2), { size: 9, fill: C.muted }));
      sq.style.cursor = "pointer";
      sq.addEventListener("mousemove", e => { sq.setAttribute("fill", C.cyan);
        showTip(d.axis.toUpperCase(), null, [["value", fmt(d.value, 2)]], e.clientX, e.clientY); });
      sq.addEventListener("mouseleave", () => { sq.setAttribute("fill", C.blue); hideTip(); });
      g.appendChild(sq);
    });
    host.appendChild(s);
  }

  // 2 — DETECTION ladder (dumbbell)
  function detection(host, data) {
    const rows = data.rows.slice().sort((a, b) => (a.ref ?? 1e9) - (b.ref ?? 1e9));
    const W = 640, rowH = 22, top = 20, L = 150, R = 26, bottom = 44;
    const H = top + rows.length * rowH + bottom;
    const s = svgRoot(W, H), g = el("g"); s.appendChild(g);
    const [lo, hi] = data.rt_range, sx = scale(lo, hi, L, W - R); sx.x0 = L; sx.x1 = W - R;
    const y0 = top + rows.length * rowH;
    well(g, L, top, W - R - L, rows.length * rowH);
    const ticks = niceTicks(lo, hi, 7);
    ticks.forEach(t => g.appendChild(line(sx(t), top, sx(t), y0, C.grid, 1)));
    xAxis(g, sx, y0, ticks, 1, "RETENTION TIME / MIN");
    rows.forEach((r, i) => {
      const cy = top + i * rowH + rowH / 2;
      if (i % 2) g.appendChild(el("rect", { x: L, y: top + i * rowH, width: W - R - L, height: rowH, fill: C.band }));
      const xs = [r.ref, r.std, r.samp].filter(v => v != null).map(sx);
      if (xs.length > 1) g.appendChild(line(Math.min(...xs), cy, Math.max(...xs), cy, C.muted2, 4));
      if (r.ref != null) g.appendChild(square(sx(r.ref), cy, 9, C.well, C.navy, 2));   // reference = hollow sq
      if (r.std != null) g.appendChild(square(sx(r.std), cy, 9, C.navy, C.ink, 2));     // standards run = filled sq
      if (r.samp != null) g.appendChild(diamond(sx(r.samp), cy, 9, C.navy, C.ink, 2));  // samples = diamond
      g.appendChild(el("rect", { x: 12, y: cy - 5, width: 10, height: 10, fill: classColor(r.class), stroke: C.ink, "stroke-width": 1.5 }));
      g.appendChild(txt(30, cy + 3.5, r.name, { size: 10, anchor: "start", w: 700 }));
      // hover band
      const hit = el("rect", { x: L, y: top + i * rowH, width: W - R - L, height: rowH, fill: "transparent", style: "cursor:pointer" });
      hit.addEventListener("mousemove", e => {
        showTip(r.name, r.class, [["reference", fmt(r.ref, 2)], ["standards run", fmt(r.std, 2)],
          ["samples", fmt(r.samp, 2)], ["iRT", fmt(r.irt, 1)]], e.clientX, e.clientY);
        if (r.samp != null || r.irt != null) showCard(r, e.clientX, e.clientY);
      });
      hit.addEventListener("mouseleave", hideAll);
      g.appendChild(hit);
    });
    host.appendChild(s);
  }

  // 3 — PROFILE (mirror butterfly)
  function profile(host, data) {
    const W = 640, H = 360, L = 56, R = 24, T = 22, gap = 34, bandH = 128;
    const s = svgRoot(W, H), g = el("g"); s.appendChild(g);
    const [rlo, rhi] = data.rt_range, sxB = scale(rlo, rhi, L, W - R); sxB.x0 = L; sxB.x1 = W - R;
    const sxA = scale(0, 100, L, W - R); sxA.x0 = L; sxA.x1 = W - R;
    const topY0 = T, topY1 = T + bandH, botY0 = topY1 + gap, botY1 = botY0 + bandH;
    const syB = scale(0, 1, topY1, topY0), syA = scale(0, 1, botY0, botY1);
    // wells
    well(g, L, topY0, W - R - L, bandH); well(g, L, botY0, W - R - L, bandH);
    // stepped areas
    const area = (cx, hh, sx, base, sy) => {
      let d = `M ${sx(cx[0])} ${base}`; const bw = (sx(cx[1]) - sx(cx[0])) / 2 || 1;
      cx.forEach((c, i) => { const x0 = sx(c) - bw, x1 = sx(c) + bw, y = sy(hh[i]); d += ` L ${x0} ${y} L ${x1} ${y}`; });
      d += ` L ${sx(cx[cx.length - 1])} ${base} Z`; return d;
    };
    g.appendChild(el("path", { d: area(data.before.c, data.before.h, sxB, topY1, syB), fill: C.blue, "fill-opacity": .18, stroke: C.blue, "stroke-width": 2 }));
    g.appendChild(el("path", { d: area(data.after.c, data.after.h, sxA, botY0, syA), fill: C.teal, "fill-opacity": .18, stroke: C.teal, "stroke-width": 2 }));
    // standard ribbons + ticks
    data.standards.forEach(st => {
      if (st.rt == null || st.irt == null) return;
      const xb = sxB(st.rt), xa = sxA(st.irt), cc = classColor(st.class);
      g.appendChild(line(xb, topY1, xa, botY0, C.muted2, 1.5));
      g.appendChild(square(xb, topY0 + 8, 7, C.well, C.navy, 2));      // hollow both ends
      g.appendChild(square(xa, botY1 - 8, 7, C.well, C.navy, 2));
      g.appendChild(txt(xb + 7, topY0 + 11, (st.class || ""), { size: 8, w: 700, anchor: "start", fill: cc }));
      const onHover = e => { showTip(st.name, st.class, [["raw RT", fmt(st.rt, 2)], ["iRT", fmt(st.irt, 1)]], e.clientX, e.clientY); showCard(st, e.clientX, e.clientY); };
      // hit regions: top tick column, bottom tick column, and the ribbon itself
      [[xb, topY0, topY1 - topY0 + 8], [xa, botY0 - 8, botY1 - botY0 + 8]].forEach(([hx, hy, hh]) => {
        const hit = el("rect", { x: hx - 7, y: hy, width: 14, height: hh, fill: "transparent", style: "cursor:pointer" });
        hit.addEventListener("mousemove", onHover); hit.addEventListener("mouseleave", hideAll); g.appendChild(hit);
      });
      const rhit = el("line", { x1: xb, y1: topY1, x2: xa, y2: botY0, stroke: "transparent", "stroke-width": 10, style: "cursor:pointer" });
      rhit.addEventListener("mousemove", onHover); rhit.addEventListener("mouseleave", hideAll); g.appendChild(rhit);
    });
    xAxis(g, sxB, topY1, niceTicks(rlo, rhi, 6), 1, null);
    g.appendChild(txt((L + W - R) / 2, topY0 - 6, "BEFORE · RAW RT (MIN)", { size: 9.5, w: 700, ls: ".08em", fill: C.muted }));
    xAxis(g, sxA, botY1, niceTicks(0, 100, 6), 0, "iRT INDEX");
    g.appendChild(txt((L + W - R) / 2, botY0 - 6, "AFTER · iRT", { size: 9.5, w: 700, ls: ".08em", fill: C.teal }));
    host.appendChild(s);
  }

  // 4 — WARP curve + residual strip
  function warp(host, data) {
    if (data.empty) { host.innerHTML = `<div class="ck-plate">CALIBRATION WARP NOT AVAILABLE</div>`; return; }
    const W = 640, L = 58, R = 22, T = 18, mainH = 250, gap = 16, resH = 70, B = 42;
    const H = T + mainH + gap + resH + B;
    const s = svgRoot(W, H), g = el("g"); s.appendChild(g);
    const rx = data.rt_span, sx = scale(rx[0], rx[1], L, W - R); sx.x0 = L; sx.x1 = W - R;
    const sy = scale(0, 100, T + mainH, T); sy.y0 = T; sy.y1 = T + mainH;
    well(g, L, T, W - R - L, mainH);
    const xt = niceTicks(rx[0], rx[1], 6), yt = niceTicks(0, 100, 5);
    gridY(g, sy, L, W - R, yt);
    // span wash
    g.appendChild(el("rect", { x: sx(data.anchors[0].x), y: T, width: sx(data.anchors[data.anchors.length - 1].x) - sx(data.anchors[0].x), height: mainH, fill: C.blue, "fill-opacity": .06 }));
    // curve
    const pts = data.curve.x.map((x, i) => `${sx(x)},${sy(data.curve.y[i])}`).filter(p => p.indexOf("NaN") < 0).join(" ");
    g.appendChild(el("polyline", { points: pts, fill: "none", stroke: C.navy, "stroke-width": 3 }));
    yAxis(g, sy, L, yt, 0, "iRT INDEX");
    // hover crosshair — live RT -> iRT probe with axis-readout chips
    const gx = data.curve.x, gy = data.curve.y;
    const interp = rt => {
      if (rt <= gx[0]) return gy[0]; if (rt >= gx[gx.length - 1]) return gy[gy.length - 1];
      let i = 1; while (i < gx.length && gx[i] < rt) i++;
      const t = (rt - gx[i - 1]) / ((gx[i] - gx[i - 1]) || 1); return gy[i - 1] + t * (gy[i] - gy[i - 1]);
    };
    const cross = el("g", { style: "display:none", "pointer-events": "none" });
    const vg = line(0, T, 0, T + mainH, C.cyan, 1.5, "4 3"), hg = line(L, 0, W - R, 0, C.cyan, 1.5, "4 3");
    const dot = square(0, 0, 9, C.cyan, C.ink, 1.5);
    const cxr = el("rect", { fill: C.ink, "shape-rendering": "crispEdges" }), cxt = txt(0, 0, "", { size: 9, w: 700, fill: C.card });
    const cyr = el("rect", { fill: C.ink, "shape-rendering": "crispEdges" }), cyt = txt(0, 0, "", { size: 9, w: 700, fill: C.card, anchor: "end" });
    [vg, hg, dot, cxr, cxt, cyr, cyt].forEach(n => cross.appendChild(n));
    g.appendChild(cross);
    const pt = s.createSVGPoint();
    const wellHit = el("rect", { x: L, y: T, width: W - R - L, height: mainH, fill: "transparent", style: "cursor:crosshair" });
    wellHit.addEventListener("mousemove", e => {
      pt.x = e.clientX; pt.y = e.clientY; const p = pt.matrixTransform(s.getScreenCTM().inverse());
      const rt = Math.max(rx[0], Math.min(rx[1], sx.inv(p.x))), iv = interp(rt), px = sx(rt), py = sy(iv);
      cross.style.display = "";
      vg.setAttribute("x1", px); vg.setAttribute("x2", px);
      hg.setAttribute("y1", py); hg.setAttribute("y2", py);
      dot.setAttribute("x", px - 4.5); dot.setAttribute("y", py - 4.5);
      cxt.textContent = fmt(rt, 2); cxt.setAttribute("x", px); cxt.setAttribute("y", T + mainH + 13);
      cxr.setAttribute("x", px - 21); cxr.setAttribute("y", T + mainH + 3); cxr.setAttribute("width", 42); cxr.setAttribute("height", 15);
      cyt.textContent = fmt(iv, 1); cyt.setAttribute("x", L - 12); cyt.setAttribute("y", py + 3.5);
      cyr.setAttribute("x", L - 46); cyr.setAttribute("y", py - 8); cyr.setAttribute("width", 36); cyr.setAttribute("height", 15);
      showTip("PROBE", null, [["RT", fmt(rt, 2)], ["iRT", fmt(iv, 1)]], e.clientX, e.clientY);
    });
    wellHit.addEventListener("mouseleave", () => { cross.style.display = "none"; hideTip(); });
    g.appendChild(wellHit);
    // residual strip
    const ry0 = T + mainH + gap;
    well(g, L, ry0, W - R - L, resH);
    const loos = data.anchors.map(a => a.loo).filter(v => v != null && isFinite(v));
    const lmax = Math.max(0.5, ...loos.map(Math.abs)) * 1.15;
    const sr = scale(-lmax, lmax, ry0 + resH - 6, ry0 + 6), zeroY = sr(0);
    g.appendChild(line(L, zeroY, W - R, zeroY, C.ink, 2));
    g.appendChild(txt(W - R - 2, ry0 + 12, "LOO", { size: 9, w: 700, anchor: "end", fill: C.muted }));
    // anchors + residual stems
    data.anchors.forEach(a => {
      const x = sx(a.x);
      if (a.loo != null && isFinite(a.loo)) {
        g.appendChild(line(x, zeroY, x, sr(a.loo), C.ink, 2));
        g.appendChild(square(x, sr(a.loo), 7, C.blue, C.ink, 1.5));
      }
      const sq = square(x, sy(a.y), 11, C.blue); sq.style.cursor = "pointer"; g.appendChild(sq);
      const hit = el("rect", { x: x - 9, y: T, width: 18, height: mainH + gap + resH, fill: "transparent", style: "cursor:pointer" });
      hit.addEventListener("mousemove", e => {
        sq.setAttribute("fill", C.cyan);
        showTip(a.name || "anchor", a.class, [["obs RT", fmt(a.x, 2)], ["iRT", fmt(a.y, 1)], ["LOO", fmt(a.loo, 2)]], e.clientX, e.clientY);
        if (a.name) showCard({ name: a.name, class: a.class, irt: a.y }, e.clientX, e.clientY);
      });
      hit.addEventListener("mouseleave", () => { sq.setAttribute("fill", C.blue); hideAll(); });
      g.appendChild(hit);
    });
    xAxis(g, sx, ry0 + resH, xt, 1, "OBSERVED RETENTION TIME / MIN");
    host.appendChild(s);
  }

  // 5 — REPEATABILITY (histogram + IQR strip)
  function repeatability(host, data) {
    if (!data.applicable) {
      host.innerHTML = `<div class="ck-plate">SINGLE-FILE MODE · repeatability needs per-injection files</div>`;
      return;
    }
    const W = 640, H = 320, L = 56, R = 24, T = 20, B = 44, midGap = 40;
    const pw = (W - L - R - midGap) / 2;
    const s = svgRoot(W, H), g = el("g"); s.appendChild(g);
    // Panel A: histogram
    const ed = data.hist.edges, ct = data.hist.counts, aL = L, aR = aL + pw;
    const sxA = scale(ed[0], ed[ed.length - 1], aL, aR); sxA.x0 = aL; sxA.x1 = aR;
    const cmax = Math.max(1, ...ct), syA = scale(0, cmax, T + H - T - B, T); syA.y0 = T; syA.y1 = H - B;
    well(g, aL, T, pw, H - B - T);
    for (let i = 0; i < ct.length; i++) {
      const x0 = sxA(ed[i]), x1 = sxA(ed[i + 1]), y = syA(ct[i]);
      g.appendChild(el("rect", { x: x0, y, width: Math.max(1, x1 - x0 - 1), height: (H - B) - y, fill: C.blue, stroke: C.ink, "stroke-width": 1 }));
    }
    yAxis(g, syA, aL, niceTicks(0, cmax, 4), 0, "FEATURES");
    xAxis(g, sxA, H - B, niceTicks(ed[0], ed[ed.length - 1], 5), 2, "RI SPREAD (iRT)");
    if (data.stats.median != null) g.appendChild(line(sxA(data.stats.median), T, sxA(data.stats.median), H - B, C.navy, 2));
    g.appendChild(txt(aL, T - 6, "DISTRIBUTION", { size: 9.5, w: 700, ls: ".06em", anchor: "start", fill: C.muted }));
    g.appendChild(txt(aR, T - 6, "MED " + fmt(data.stats.median, 2), { size: 9, w: 700, anchor: "end", fill: C.navy }));
    // Panel B: IQR vs iRT
    const bL = aR + midGap, bR = W - R, sxB = scale(0, 100, bL, bR); sxB.x0 = bL; sxB.x1 = bR;
    const bn = data.binned, allq = bn.q3.concat(bn.med); const smax = Math.max(0.1, ...allq) * 1.1;
    const syB = scale(0, smax, H - B, T); syB.y0 = T; syB.y1 = H - B;
    well(g, bL, T, bR - bL, H - B - T);
    gridY(g, syB, bL, bR, niceTicks(0, smax, 4));
    bn.c.forEach((c, i) => {
      const x = sxB(c), bw = 9;
      g.appendChild(el("rect", { x: x - bw / 2, y: syB(bn.q3[i]), width: bw, height: Math.max(1, syB(bn.q1[i]) - syB(bn.q3[i])), fill: C.muted2, "fill-opacity": .5, stroke: C.ink, "stroke-width": 1.5 }));
      g.appendChild(line(x - bw / 2, syB(bn.med[i]), x + bw / 2, syB(bn.med[i]), C.blue, 3));
    });
    const stepPts = bn.c.map((c, i) => `${sxB(c)},${syB(bn.med[i])}`).join(" ");
    g.appendChild(el("polyline", { points: stepPts, fill: "none", stroke: C.blue, "stroke-width": 2 }));
    xAxis(g, sxB, H - B, niceTicks(0, 100, 5), 0, "iRT INDEX");
    yAxis(g, syB, bL, niceTicks(0, smax, 4), 2, null);
    g.appendChild(txt(bL, T - 6, "SPREAD vs iRT", { size: 9.5, w: 700, ls: ".06em", anchor: "start", fill: C.muted }));
    host.appendChild(s);
  }

  function render(id, host, bundle) {
    STRUCT = bundle.structures || {};
    host.innerHTML = "";
    if (id === "detection") return detection(host, bundle.detection);
    if (id === "profile") return profile(host, bundle.profile);
    if (id === "warp") return warp(host, bundle.warp);
    if (id === "repeatability") return repeatability(host, bundle.repeatability);
  }
  return { render, radar, C, classColor, hideAll };
})();
