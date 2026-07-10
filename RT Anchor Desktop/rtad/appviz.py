"""App-native visualization data.

Extracts the raw + lightly-derived numeric arrays each in-app chart needs, straight
from a ``CalibrationResult`` — deliberately WITHOUT importing ``rt_anchor.viz`` (the
package's Plotly/matplotlib figure builders). The desktop front-end draws everything
itself as hand-authored SVG. Only the calibration engine (``rt_anchor.panel`` /
``.config`` / the result object) is used here.
"""

from __future__ import annotations

import math
from typing import Dict, List

import numpy as np
import pandas as pd


# --------------------------------------------------------------- json helpers --

def _f(v):
    """A finite float, or None (JSON-safe)."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _arr(a) -> List:
    return [_f(x) for x in np.asarray(a, dtype=float).ravel()]


# --------------------------------------------------------------- shared bits ---

def _intensity(result) -> np.ndarray:
    df = result.table
    if result.sample_cols:
        return (df[result.sample_cols].apply(pd.to_numeric, errors="coerce")
                .fillna(0.0).sum(axis=1).to_numpy())
    cols = {c.lower(): c for c in df.columns}
    for key in ("peak_height", "peak_area", "maxo", "height", "area", "intensity"):
        if key in cols:
            return pd.to_numeric(df[cols[key]], errors="coerce").fillna(0.0).to_numpy()
    return np.ones(len(df), dtype=float)


def _panel_targets(result):
    from rt_anchor.config import CalibrationConfig
    from rt_anchor.panel import build_panel
    cfg = CalibrationConfig.from_dict(result.model.get("config", {}))
    return build_panel(result.polarity or "positive", cfg).targets


def _anchor_lookup(result):
    a = result.anchors
    if not len(a):
        return {}
    a = a.drop_duplicates("name")
    out = {}
    for _, r in a.iterrows():
        out[r["name"]] = (float(r["rt_obs_min"]), r.get("class", ""))
    return out


# --------------------------------------------------------------- KPIs ----------

def kpis(result) -> List[Dict]:
    tbl = result.table
    ri = pd.to_numeric(tbl[result.col("RI")], errors="coerce")
    rt = result.rt_minutes()
    panel = result.panel or {}
    anc = result.anchors.drop_duplicates("name") if len(result.anchors) else result.anchors
    n_panel = int(panel.get("n_panel", len(anc)))
    n_samp = int(anc["name"].nunique()) if len(anc) else 0
    n_panel_det = int(panel.get("n_detected", n_samp))
    pr = panel.get("rt_range", [np.nan, np.nan])
    offs = (np.abs(anc["rt_obs_min"].to_numpy() - anc["rt_ref_min"].to_numpy())
            if len(anc) else np.array([]))
    off = _f(np.median(offs)) if offs.size else None
    panel_rng = "—" if not (pr and math.isfinite(pr[0])) else f"{pr[0]:.1f}–{pr[1]:.1f}"

    def rng(a, b, d):
        a, b = _f(a), _f(b)
        return "—" if (a is None or b is None) else f"{a:.{d}f}–{b:.{d}f}"

    def tile(label, value, sub):
        return {"label": label, "value": value, "sub": sub}

    return [
        tile("Features", f"{len(tbl):,}", "calibrated"),
        tile("Standards", f"{n_samp}/{n_panel}", f"samples · panel {n_panel_det}/{n_panel}"),
        tile("RT range", rng(rt.min(), rt.max(), 1), f"samples min · panel {panel_rng}"),
        tile("iRT range", rng(ri.min(), ri.max(), 0), "dimensionless"),
        tile("RT offset", "—" if off is None else f"{off:.2f}", "median |obs−ref| (min)"),
        tile("Coverage", f"{ri.notna().mean() * 100:.0f}%", "features with an RI"),
    ]


# --------------------------------------------------------------- radar ---------

def radar(result) -> List[Dict]:
    tbl = result.table
    ri = pd.to_numeric(tbl[result.col("RI")], errors="coerce")
    panel = result.panel or {}
    anc = result.anchors.drop_duplicates("name") if len(result.anchors) else result.anchors
    n_panel = max(int(panel.get("n_panel", len(anc))), 1)
    detection = (anc["name"].nunique() / n_panel) if len(anc) else 0.0
    coverage = float(ri.notna().mean()) if len(tbl) else 0.0
    if len(anc):
        span = (float(anc["irt"].max() - anc["irt"].min()) / 100.0)
    else:
        span = 0.0
    loo = (result.model.get("loo_residual_irt", {}) or {}).get("median")
    interp = 1.0 - min((loo or 0.0) / 3.0, 1.0)

    if result.model.get("calibration_scope") == "sample":
        sp = pd.to_numeric(tbl[result.col("RI_spread")], errors="coerce")
        nc = pd.to_numeric(tbl[result.col("n_contributing")], errors="coerce")
        sp = sp[nc >= 2]                       # exclude single-injection (structural spread 0)
        med = float(sp.median()) if sp.notna().any() else 1.0
        last = ("Repeatability", 1.0 - min(med, 1.0))
    else:
        native = panel.get("native_rt", {})
        alu = _anchor_lookup(result)
        diffs = [abs(native[n] - alu[n][0]) for n in alu if n in native]
        med = float(np.median(diffs)) if diffs else 0.0
        last = ("Panel match", 1.0 - min(med / 0.5, 1.0))

    axes = [("Detection", detection), ("Coverage", coverage), ("Anchor span", min(span, 1.0)),
            ("Interpolation", interp), last]
    return [{"axis": a, "value": max(0.0, min(1.0, _f(v) or 0.0))} for a, v in axes]


# --------------------------------------------------------------- detection -----

def detection(result) -> Dict:
    targets = _panel_targets(result)
    native = (result.panel or {}).get("native_rt", {})
    alu = _anchor_lookup(result)
    rows = []
    for _, t in targets.iterrows():
        nm = t["name"]
        rows.append({
            "name": nm, "class": t.get("class", ""),
            "ref": _f(t["rt_ref_min"]), "std": _f(native.get(nm)),
            "samp": _f(alu[nm][0]) if nm in alu else None, "irt": _f(t["irt"]),
        })
    xs = [v for r in rows for v in (r["ref"], r["std"], r["samp"]) if v is not None]
    rng = [min(xs), max(xs)] if xs else [0.0, 1.0]
    return {"rows": rows, "rt_range": rng}


# --------------------------------------------------------------- warp ----------

def warp(result) -> Dict:
    w = result.warp
    if w is None:
        return {"empty": True}
    grid = np.linspace(w.rt_min, w.rt_max, 220)
    curve = w.predict(grid, extrapolate=False)
    # label fitted anchors by iRT (exact & unique per standard; w.irt carries it
    # verbatim) — matching by rt_obs breaks in per-sample mode where w.rt is the
    # per-standard MEDIAN but pooled anchors keep the first injection's rt_obs.
    look = {}
    adf = result.anchors.drop_duplicates("name") if len(result.anchors) else result.anchors
    for _, r in adf.iterrows():
        look[round(float(r["irt"]), 6)] = (r["name"], r.get("class", ""))
    anchors = []
    loo = getattr(w, "loo_resid", np.full(w.rt.shape, np.nan))
    for x, y, lo in zip(w.rt, w.irt, loo):
        nm, cl = look.get(round(float(y), 6), (None, None))
        anchors.append({"x": _f(x), "y": _f(y), "loo": _f(lo), "name": nm, "class": cl})
    return {"empty": False, "curve": {"x": _arr(grid), "y": _arr(curve)},
            "anchors": anchors, "rt_span": [_f(w.rt_min), _f(w.rt_max)], "irt_range": [0.0, 100.0]}


# --------------------------------------------------------------- profile -------

def profile(result, nbins: int = 120) -> Dict:
    rt = result.rt_minutes().to_numpy()
    ri = pd.to_numeric(result.table[result.col("RI")], errors="coerce").to_numpy()
    inten = _intensity(result)

    okb = np.isfinite(rt) & np.isfinite(inten)
    rt_lo, rt_hi = (float(np.nanmin(rt[okb])), float(np.nanmax(rt[okb]))) if okb.any() else (0.0, 1.0)
    eb = np.linspace(rt_lo, rt_hi, nbins + 1)
    hb = np.histogram(rt[okb], bins=eb, weights=inten[okb])[0].astype(float)

    oka = np.isfinite(ri) & np.isfinite(inten)
    ea = np.linspace(0.0, 100.0, nbins + 1)
    ha = np.histogram(ri[oka], bins=ea, weights=inten[oka])[0].astype(float)

    def _norm(h):
        m = h.max()
        return (h / m) if m > 0 else h

    cb = (eb[:-1] + eb[1:]) / 2.0
    ca = (ea[:-1] + ea[1:]) / 2.0
    stds = []
    for nm, (rtobs, cl) in _anchor_lookup(result).items():
        stds.append({"name": nm, "class": cl, "rt": _f(rtobs)})
    # attach each standard's iRT from the panel
    irt_by_name = {t["name"]: _f(t["irt"]) for _, t in _panel_targets(result).iterrows()}
    for s in stds:
        s["irt"] = irt_by_name.get(s["name"])
    return {"before": {"c": _arr(cb), "h": _arr(_norm(hb))},
            "after": {"c": _arr(ca), "h": _arr(_norm(ha))},
            "rt_range": [rt_lo, rt_hi], "irt_range": [0.0, 100.0], "standards": stds}


# --------------------------------------------------------------- repeatability -

def repeatability(result) -> Dict:
    scope = result.model.get("calibration_scope")
    sp = pd.to_numeric(result.table[result.col("RI_spread")], errors="coerce").to_numpy()
    ri = pd.to_numeric(result.table[result.col("RI")], errors="coerce").to_numpy()
    nc = pd.to_numeric(result.table[result.col("n_contributing")], errors="coerce").to_numpy()
    ok = np.isfinite(sp) & (nc >= 2)           # a spread needs >=2 contributing injections
    if scope != "sample" or ok.sum() < 5:
        return {"applicable": False}
    s = sp[ok]
    hi = float(np.percentile(s, 98)) or float(s.max()) or 1.0
    counts, edges = np.histogram(np.clip(s, 0, hi), bins=40)

    okj = np.isfinite(sp) & np.isfinite(ri) & (nc >= 2)
    sj, rj = sp[okj], ri[okj]
    ebins = np.linspace(0.0, 100.0, 21)
    c, med, q1, q3 = [], [], [], []
    for i in range(len(ebins) - 1):
        m = (rj >= ebins[i]) & (rj < ebins[i + 1])
        if m.sum() >= 3:
            c.append((ebins[i] + ebins[i + 1]) / 2.0)
            med.append(float(np.median(sj[m]))); q1.append(float(np.percentile(sj[m], 25)))
            q3.append(float(np.percentile(sj[m], 75)))
    return {"applicable": True,
            "hist": {"edges": _arr(edges), "counts": [int(x) for x in counts]},
            "binned": {"c": _arr(c), "med": _arr(med), "q1": _arr(q1), "q3": _arr(q3)},
            "stats": {"median": _f(np.median(s)), "p90": _f(np.percentile(s, 90))}}


# --------------------------------------------------------------- table preview -

def table_preview(result, n: int = 120) -> Dict:
    tbl = result.table
    wanted = [(result.rt_col, "RT"), (result.mz_col, "m/z"),
              (result.col("RI"), "RI"), (result.col("RI_uncertainty"), "RI_uncertainty"),
              (result.col("RI_reliability"), "RI_reliability"),
              (result.col("RI_spread"), "RI_spread"), (result.col("is_extrapolated"), "extrapolated")]
    cols = [(c, d) for c, d in wanted if c and c in tbl.columns]
    head = tbl[[c for c, _ in cols]].head(n)

    def cell(v):
        if isinstance(v, (bool, np.bool_)):
            return "yes" if bool(v) else "no"
        if isinstance(v, (int, np.integer)):
            return int(v)
        if isinstance(v, (float, np.floating)):
            return _f(round(float(v), 4))
        return None if pd.isna(v) else str(v)

    rows = [[cell(r[c]) for c, _ in cols] for _, r in head.iterrows()]
    return {"columns": [d for _, d in cols], "rows": rows,
            "n_total": int(len(tbl)), "n_shown": int(len(head))}


# --------------------------------------------------------------- bundle --------

DETECTION_NOTE = (
    "Each row is a standard; its retention time is shown from three sources — "
    "<b>hollow square</b> = reference baseline, <b>filled square</b> = standards run, "
    "<b>diamond</b> = samples (the table being calibrated). The left chip colours the lipid "
    "class. A shorter span means the observed RT agrees more closely with the reference."
)


def build_viz(result) -> Dict:
    m = result.model
    n = len(result.table)
    ri = pd.to_numeric(result.table[result.col("RI")], errors="coerce")
    panel = result.panel or {}
    anc = result.anchors.drop_duplicates("name") if len(result.anchors) else result.anchors
    rep = repeatability(result)

    struct = {}
    if result.default_panel:
        try:
            from rt_anchor.viz import structures  # molecular depiction only (not a chart)
            struct = structures.render_default_structures()
        except Exception:
            struct = {}

    return {
        "kpis": kpis(result),
        "radar": radar(result),
        "detection": detection(result),
        "warp": warp(result),
        "profile": profile(result),
        "repeatability": rep,
        "table": table_preview(result),
        "structures": struct,
        "notes": {"detection": DETECTION_NOTE},
        "meta": {
            "scope": m.get("calibration_scope", "project"),
            "polarity": result.polarity, "source_format": result.source_format,
            "n_features": int(n),
            "n_detected": int(anc["name"].nunique()) if len(anc) else 0,
            "n_panel": int(panel.get("n_panel", len(anc))),
            "has_repeatability": bool(rep.get("applicable")),
        },
    }
