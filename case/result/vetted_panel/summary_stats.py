"""Summary numbers for the vetted-lipid validation (RT Anchor 1.2.3, reprocessed data).

1. Funnel, from the vetted list to the panel:
   list entries -> entries with a positive-mode target ion (searchable) -> distinct targets -> effective
   targets (targets whose m/z lie within 10 ppm of each other merged; TG / MG carry two ions, either counts)
   -> effective targets with a feature within 10 ppm (and 15 ppm, the compare tables' window) on >= 3 / 4 /
   5 / 6 of the six columns -> panel lipids (identity-checked usable peak, build_vetted_panel.py rules) by the
   number of columns on which the peak is usable.
2. Error metrics before and after calibration, overall and per column, on the high-confidence panel
   (830 lipid-by-column points): median |dRT| (median absolute deviation from the reference RT), 90th
   percentile |dRT|, range of the signed dRT, share within +-0.2 and +-0.5 min. "After" is given in-sample
   (the delivered calibration, whose stage-1 curve was fitted with 83 % of these points as pairs) and held out
   (holdout_cv.py: 10-fold, the lipid and its co-eluting reference features removed from the reference).
   95 % confidence intervals: bootstrap over lipids (2000 resamples; a lipid's points stay together).
3. Checks of sufficiency: points per column and per class, coverage of the reference gradient, the
   identity-robust subset (one single candidate on the column and on the reference), in-sample share.

Writes summary_funnel.csv, summary_metrics.csv and summary_console.log in this folder.
Run:  python summary_stats.py          (after build_vetted_panel.py and holdout_cv.py)
"""
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
CMP = os.path.join(BASE, "data", "compare")
SERUM = {
    "column_12": "12 min/Serum/12 min Serum Area_2_2026_10_12_13_18_31.csv",
    "column_16": "17 min/Serum/17 Serum Area_2_2026_10_12_13_33_45.csv",
    "column_20": "20 min/Serum/20 min Serum Area_2_2026_10_12_12_12_37.csv",
    "column_22": "22 min/Serum/22 min Serum Area_2_2026_10_12_13_24_24.csv",
    "column_41": "45 min/Serum/45 Serum Area_2_2026_10_12_12_48_51.csv",
    "column_85": "90 min/Serum/90 Serum Area_2_2026_10_12_13_11_13.csv",
}
COLS = list(SERUM)
REF = "column_22"
NONREF = [c for c in COLS if c != REF]
B = 2000
RNG = np.random.default_rng(20261014)
OUT = []


def log(s=""):
    print(s)
    OUT.append(s)


def targets_of(s):
    if not isinstance(s, str) or not s.strip():
        return []
    return [float(p.strip().rsplit(" ", 1)[1]) for p in s.split(";")]


def components(tsets, ppm):
    """Union-find over targets: two targets are linked when any of their m/z lie within ppm."""
    parent = list(range(len(tsets)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    flat = sorted((mz, i) for i, t in enumerate(tsets) for mz in t)
    for (m1, i1), (m2, i2) in zip(flat, flat[1:]):
        if (m2 - m1) / m1 * 1e6 <= ppm:
            parent[find(i1)] = find(i2)
    groups = {}
    for i in range(len(tsets)):
        groups.setdefault(find(i), set()).update(tsets[i])
    return list(groups.values())


def n_columns_found(mzs, feats, ppm):
    n = 0
    for c in COLS:
        f = feats[c]
        if any(np.any(np.abs(f - m) <= ppm * 1e-6 * m) for m in mzs):
            n += 1
    return n


def tally(n):
    n = np.asarray(n)
    return {k: int((n >= k).sum()) for k in (3, 4, 5, 6)}


def metrics(d):
    a = np.abs(d)
    return dict(n=len(d), median_abs=np.median(a), p90_abs=np.percentile(a, 90), min=d.min(), max=d.max(),
                p2_5=np.percentile(d, 2.5), p97_5=np.percentile(d, 97.5),
                within_0_2=100 * (a <= 0.2).mean(), within_0_5=100 * (a <= 0.5).mean())


def boot_ci(df, col):
    """95 % CI of median |d|, p90 |d|, % within 0.2 and 0.5 min; resampling lipids with their points."""
    by = [g[col].to_numpy(float) for _, g in df.groupby("lipid")]
    stats = []
    for _ in range(B):
        idx = RNG.integers(0, len(by), len(by))
        a = np.abs(np.concatenate([by[i] for i in idx]))
        stats.append((np.median(a), np.percentile(a, 90), 100 * (a <= 0.2).mean(), 100 * (a <= 0.5).mean()))
    lo, hi = np.percentile(np.array(stats), [2.5, 97.5], axis=0)
    return dict(zip(("median_abs", "p90_abs", "within_0_2", "within_0_5"), zip(lo, hi)))


VARIANTS = (("before calibration", "dRT_before_min"), ("after, in-sample", "dRT_after_insample_min"),
            ("after, held out (10-fold)", "dRT_after_heldout_min"))


def error_section(C, out_csv):
    """Metrics overall (with bootstrap CIs) and per column for one table of held-out cells."""
    recs = []
    head = (f"  {'':<28}{'n':>5}{'median|d|':>11}{'p90|d|':>9}{'range':>20}{'central 95%':>20}"
            f"{'<=0.2 min':>11}{'<=0.5 min':>11}")
    for scope in ["all columns"] + NONREF:
        D = C if scope == "all columns" else C[C.column == scope]
        log(f"\n  {scope}")
        log(head)
        for name, col in VARIANTS:
            m = metrics(D[col].to_numpy(float))
            rng = f"{m['min']:+.2f} .. {m['max']:+.2f}"
            c95 = f"{m['p2_5']:+.2f} .. {m['p97_5']:+.2f}"
            log(f"  {name:<28}{m['n']:>5}{m['median_abs']:>11.3f}{m['p90_abs']:>9.3f}{rng:>20}{c95:>20}"
                f"{m['within_0_2']:>10.1f}%{m['within_0_5']:>10.1f}%")
            rec = dict(scope=scope, variant=name, **m)
            if scope == "all columns" and name != "before calibration":
                ci = boot_ci(D, col)
                for k, (a, b) in ci.items():
                    rec[f"{k}_ci95"] = f"{a:.3f}-{b:.3f}" if "abs" in k else f"{a:.1f}-{b:.1f}"
                log(f"  {'':<28}95% CI: median {rec['median_abs_ci95']}, p90 {rec['p90_abs_ci95']}, "
                    f"<=0.2 {rec['within_0_2_ci95']} %, <=0.5 {rec['within_0_5_ci95']} %")
            recs.append(rec)
    pd.DataFrame(recs).to_csv(os.path.join(HERE, out_csv), index=False, float_format="%.4f")


def per_ion_section(L, S, feats):
    """The same funnel and metrics with every ion of a lipid counted as its own target (per_ion_panel.py,
    holdout_cv.py --panel per_ion_panel.csv --prefix holdout_cv_per_ion)."""
    log("\n\n######## per ion: each ion of a lipid is its own target (TG, MG list two) ########")
    rows = []

    def add(step, high, all_, note=""):
        rows.append(dict(step=step, high_confidence=high, high_and_low=all_, note=note))
        log(f"  {step:<58} high {str(high):>22}   all {str(all_):>22}   {note}")

    Sh = S[S.conf == "high"]
    ion_t = lambda d: list({(m,) for t in d.t for m in (round(x, 4) for x in t)})
    th_, ta_ = ion_t(Sh), ion_t(S)
    add("ion targets (one per distinct ion m/z)", len(th_), len(ta_))
    ch, ca = components(th_, 10.0), components(ta_, 10.0)
    add("effective ion targets (m/z within 10 ppm merged)", len(ch), len(ca))
    for ppm in (10.0, 15.0):
        th, ta = tally([n_columns_found(c, feats, ppm) for c in ch]), tally([n_columns_found(c, feats, ppm) for c in ca])
        add(f"effective ion targets with a feature within {ppm:g} ppm on >= 3/4/5/6 columns",
            "/".join(str(th[k]) for k in (3, 4, 5, 6)), "/".join(str(ta[k]) for k in (3, 4, 5, 6)),
            "any feature at the m/z; not an identification")
    PI = pd.read_csv(os.path.join(HERE, "per_ion_panel.csv"))
    PA = pd.read_csv(os.path.join(HERE, "high_and_low_confidence_panel.csv"))
    low = PA[PA.confidence == "low"]          # low-confidence lipids have one ion each: unchanged
    tp = tally(PI.n_columns_used)
    tpa = tally(pd.concat([PI.n_columns_used, low.n_columns_used]))
    add("panel: identity-checked usable ion target on >= 3/4/5/6 columns",
        "/".join(str(tp[k]) for k in (3, 4, 5, 6)), "/".join(str(tpa[k]) for k in (3, 4, 5, 6)),
        f"{PI.fold_unit.nunique()} lipids (high); per_ion_panel.py")
    two = PI[PI["class"].isin(["TG", "MG"])]
    log("  TG / MG ion targets in the panel: " + ", ".join(f"{k} {v}" for k, v in two.ion.value_counts().items())
        + f"; TG lipids with both ions {int((two.groupby('fold_unit').size() == 2).sum())}, one ion only "
        f"{int((two.groupby('fold_unit').size() == 1).sum())}")
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "summary_funnel_per_ion.csv"), index=False)
    C = pd.read_csv(os.path.join(HERE, "holdout_cv_per_ion_cells.csv"))
    assert len(C) == sum(PI[f"{c}_used"].astype(str).str.lower().eq("true").sum() for c in NONREF)
    log("\n  error vs the reference RT (min), per-ion high-confidence panel")
    error_section(C, "summary_metrics_per_ion.csv")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    hi = pd.read_csv(os.path.join(CMP, "lipid_detection_high_confidence.csv"))
    lo = pd.read_csv(os.path.join(CMP, "lipid_detection_low_confidence.csv"))
    L = pd.concat([hi.assign(conf="high"), lo.assign(conf="low")], ignore_index=True)
    L["t"] = L.positive_target_ions.apply(targets_of)
    L["species"] = L["name"].str.split("|").str[0].str.strip()
    feats = {c: pd.read_csv(os.path.join(BASE, "data", SERUM[c]), skiprows=4, low_memory=False)["Average Mz"]
             .to_numpy(float) for c in COLS}

    # ------------------------------------------------------------------ 1. funnel
    log("== 1. from the vetted list to the panel ==")
    S = L[L.t.apply(len) > 0]
    rows = []

    def add(step, high, all_, note=""):
        rows.append(dict(step=step, high_confidence=high, high_and_low=all_, note=note))
        log(f"  {step:<58} high {str(high):>22}   all {str(all_):>22}   {note}")

    add("list entries", len(hi), len(L), "492-row MS-DIAL sheet + 182-row manuscript sheet")
    add("distinct species names", hi.name.str.split("|").str[0].str.strip().nunique(), L.species.nunique(),
        "sum composition as written")
    add("entries with a positive-mode target (searchable)", int((hi.positive_target_ions.notna()).sum()), len(S),
        "the 14 FA have no positive ion")
    n_ions = lambda d: int(d.t.apply(len).sum())
    n_dist = lambda d: len({round(m, 4) for t in d.t for m in t})
    Sh = S[S.conf == "high"]
    add("target m/z values (TG / MG list two ions)", n_ions(Sh), n_ions(S))
    add("distinct target m/z values", n_dist(Sh), n_dist(S))
    tsets_h = list({tuple(sorted(round(m, 4) for m in t)) for t in Sh.t})
    tsets_a = list({tuple(sorted(round(m, 4) for m in t)) for t in S.t})
    add("distinct targets (one per lipid ion set)", len(tsets_h), len(tsets_a))
    comp_h, comp_a = components(tsets_h, 10.0), components(tsets_a, 10.0)
    add("effective targets (m/z within 10 ppm merged)", len(comp_h), len(comp_a))
    for ppm in (10.0, 15.0):
        nh = [n_columns_found(c, feats, ppm) for c in comp_h]
        na = [n_columns_found(c, feats, ppm) for c in comp_a]
        th, ta = tally(nh), tally(na)
        add(f"effective targets with a feature within {ppm:g} ppm on >= 3/4/5/6 columns",
            "/".join(str(th[k]) for k in (3, 4, 5, 6)), "/".join(str(ta[k]) for k in (3, 4, 5, 6)),
            "any feature at the m/z; not an identification")
    P = pd.read_csv(os.path.join(HERE, "high_confidence_panel.csv"))
    PA = pd.read_csv(os.path.join(HERE, "high_and_low_confidence_panel.csv"))
    tp, tpa = tally(P.n_columns_used), tally(PA.n_columns_used)
    add("panel: identity-checked usable peak on >= 3/4/5/6 columns",
        "/".join(str(tp[k]) for k in (3, 4, 5, 6)), "/".join(str(tpa[k]) for k in (3, 4, 5, 6)),
        "reference column included; 15 ppm candidates, rules of build_vetted_panel.py")
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "summary_funnel.csv"), index=False)

    # ------------------------------------------------------------------ 2. errors
    C = pd.read_csv(os.path.join(HERE, "holdout_cv_cells.csv"))
    assert len(C) == sum(P[f"{c}_used"].astype(str).str.lower().eq("true").sum() for c in NONREF)
    log("\n== 2. error vs the reference RT (min), high-confidence panel ==")
    error_section(C, "summary_metrics.csv")

    # ------------------------------------------------------------------ 3. sufficiency
    log("\n== 3. sufficiency checks ==")
    log("  points per column: " + ", ".join(f"{c} {int((C.column == c).sum())}" for c in NONREF))
    log("  lipids per display group: " + ", ".join(f"{g} {n}" for g, n in P.display_group.value_counts().items()))
    ref = pd.read_csv(os.path.join(BASE, "data", SERUM[REF]), skiprows=4, low_memory=False)["Average Rt(min)"]
    r = np.sort(P.column_22_RT_min.to_numpy(float))
    gaps = np.diff(r)
    log(f"  reference RT of the panel lipids {r.min():.2f}-{r.max():.2f} min (reference run features "
        f"{ref.min():.2f}-{ref.max():.2f}); largest gap between neighbouring panel lipids {gaps.max():.2f} min "
        f"(between {r[gaps.argmax()]:.2f} and {r[gaps.argmax() + 1]:.2f}); lipids per 2-min bin: "
        + ", ".join(f"{int(a)}-{int(a) + 2}: {int(((r >= a) & (r < a + 2)).sum())}" for a in range(0, 24, 2)))
    uniq = []
    for c in NONREF:
        m = (P[f"{c}_used"].astype(str).str.lower() == "true") & (P[f"{c}_n_candidates"] == 1) & \
            (P.column_22_n_candidates == 1)
        uniq += [(c, l) for l in P.lipid[m]]
    U = C.merge(pd.DataFrame(uniq, columns=["column", "lipid"]), on=["column", "lipid"])
    for name, col in VARIANTS[1:]:
        m = metrics(U[col].to_numpy(float))
        log(f"  identity-robust subset (single candidate on the column and on the reference), {name}: n {m['n']}, "
            f"median |d| {m['median_abs']:.3f}, p90 {m['p90_abs']:.3f}, <=0.2 {m['within_0_2']:.1f} %, "
            f"<=0.5 {m['within_0_5']:.1f} %")
    for name, col in VARIANTS[1:]:
        a = C[col].abs()
        log(f"  points beyond 0.5 / 1 / 2 min, {name}: {(a > 0.5).sum()} / {(a > 1).sum()} / {(a > 2).sum()}")
    if os.path.exists(os.path.join(HERE, "holdout_cv_per_ion_cells.csv")):
        per_ion_section(L, S, feats)
    with open(os.path.join(HERE, "summary_console.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
