"""Panel with each ion of a lipid counted as its own point ("per ion"), for the statistics only.

In the delivered panel (build_vetted_panel.py) a lipid listed with two ions - TG ([M+NH4]+, [M+Na]+) and
MG ([M+Na]+, [M+H-H2O]+) - is one target: the most intense candidate over both ions is its peak. Here every
ion is a target of its own, with the same rules applied to it alone:

  candidates  features within 15 ppm of the ion
  pick        the most intense candidate (mean area over the 4 injections; tie: smallest Alignment ID)
  usable      pick in >= 3 of 4 injections and >= 2x the next candidate (or the only one)
  panel       reference column usable, usable on >= 3 of 6 columns, column_20 usable and its pick within
              0.30 min of the mapped vetted RT (the delivered mapping curves, vetted_rt_mapping.csv)

Every other class has one ion per lipid, so its targets are unchanged. Running the same code with
--mode lipid (both ions of a lipid pooled, the other adduct of the same peak not a rival) must reproduce the
delivered high_confidence_panel.csv exactly; the script checks that and stops if it does not.

Writes ../work/per_ion_panel.csv in the column layout of the delivered panel CSVs (plus `ion` and `fold_unit`), so
holdout_cv.py can run on it. Only high-confidence rows are used, as in the delivered statistics.

Run:  python per_ion_panel.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd

SCRIPTS = os.path.dirname(os.path.abspath(__file__))           # .../vetted_panel/scripts
HERE = os.path.dirname(SCRIPTS)                                # .../vetted_panel
WORK = os.path.join(HERE, "work")                              # intermediate tables (removed by run_all.py)
RUN = os.path.dirname(HERE)
BASE = os.path.dirname(RUN)
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
PPM, DOM, ADDUCT_WIN, ID_TOL = 15.0, 2.0, 0.10, 0.30


def load():
    feats, cal = {}, {}
    for c in COLS:
        t = pd.read_csv(os.path.join(BASE, "data", SERUM[c]), skiprows=4, low_memory=False)
        cols = list(t.columns)
        i = cols.index("MS/MS spectrum")
        s = cols[i + 1:i + 5]
        t["area"] = t[s].astype(float).mean(axis=1)
        t["n_det"] = (t[s].astype(float) > 0).sum(axis=1)
        feats[c] = t
        cal[c] = pd.read_csv(os.path.join(RUN, c, "rt_anchor_run_calibrated.csv"), low_memory=False).set_index("Alignment ID")
    M = pd.read_csv(os.path.join(WORK, "vetted_rt_mapping.csv"))
    maps = {str(s): (g.vetted_RT_min.to_numpy(float), g.column_20_RT_min.to_numpy(float))
            for s, g in M.groupby("sheet")}
    return feats, cal, maps


def units(df, mode):
    """Targets: mode 'lipid' = the row's whole ion set, mode 'ion' = each ion alone."""
    out = {}
    for _, r in df.iterrows():
        ions = [(p.strip().rsplit(" ", 1)[0], round(float(p.strip().rsplit(" ", 1)[1]), 4))
                for p in r.positive_target_ions.split(";")]
        full = tuple(sorted(ions))
        for key in ([full] if mode == "lipid" else [(i,) for i in ions]):
            u = out.setdefault(key, dict(key=key, rows=[], fold_unit=None))
            u["rows"].append(r)
            u["fold_unit"] = u["fold_unit"] or "|".join(f"{a} {b}" for a, b in full)
    return list(out.values())


def pick(t, targets):
    mz = t["Average Mz"].to_numpy(float)
    via = {}
    for ion, v in targets:
        for j in np.where(np.abs(mz - v) <= PPM * 1e-6 * v)[0]:
            via.setdefault(int(j), set()).add(ion)
    if not via:
        return None
    idx = sorted(via)
    sub = t.iloc[idx].assign(_j=idx).sort_values(["area", "Alignment ID"], ascending=[False, True])
    p = sub.iloc[0]
    rivals = [q.area for _, q in sub.iloc[1:].iterrows()
              if not (not (via[int(q._j)] & via[int(p._j)]) and abs(q["Average Rt(min)"] - p["Average Rt(min)"]) <= ADDUCT_WIN + 1e-9)]
    rmax = max(rivals) if rivals else 0.0
    dom = np.inf if rmax <= 0 else p.area / rmax
    return dict(aid=int(p["Alignment ID"]), rt=float(p["Average Rt(min)"]), mz=float(p["Average Mz"]),
                area=float(p.area), n_det=int(p.n_det), n_cand=len(idx), dom=float(dom),
                usable=bool(p.n_det >= 3 and dom >= DOM))


def build(mode, feats, cal, maps, hi):
    recs = []
    for u in units(hi, mode):
        rows = pd.DataFrame(u["rows"])
        P = {c: pick(feats[c], u["key"]) for c in COLS}
        usable = {c: bool(P[c] and P[c]["usable"]) for c in COLS}
        sheet = ["492" if s.startswith("492") else "182" for s in rows.sheet]
        mapped = [float(np.interp(v, *maps[s])) for v, s in zip(rows.vetted_RT_min, sheet)]
        if P["column_20"] is not None:
            d = [abs(P["column_20"]["rt"] - m) for m in mapped]
            k = int(np.argmin(d))
        else:
            k, d = 0, [np.nan]
        ok = (usable[REF] and sum(usable.values()) >= 3 and usable["column_20"] and d[k] <= ID_TOL + 1e-9)
        if not ok:
            continue
        r = rows.iloc[k]
        rec = dict(lipid=re.split(r"\s*\|\s*", r["name"])[0].strip(), ion="; ".join(i for i, _ in u["key"]),
                   **{"class": r["class"]}, display_group=r["class"], fold_unit=u["fold_unit"],
                   vetted_names="; ".join(rows.name), n_columns_used=sum(usable.values()))
        ref_rt = P[REF]["rt"]
        for c in COLS:
            p = P[c]
            rec[f"{c}_used"] = usable[c]
            rec[f"{c}_alignment_id"] = p["aid"] if p else np.nan
            rec[f"{c}_RT_min"] = p["rt"] if p else np.nan
            rec[f"{c}_n_candidates"] = p["n_cand"] if p else 0
        for c in NONREF:
            if usable[c]:
                q = cal[c].loc[P[c]["aid"]]
                rec[f"{c}_cal_RT_min"] = float(q["Cal_RT_min"])
                rec[f"{c}_dRT_before_min"] = P[c]["rt"] - ref_rt
                rec[f"{c}_dRT_after_min"] = float(q["Cal_RT_min"]) - ref_rt
            else:
                rec[f"{c}_cal_RT_min"] = rec[f"{c}_dRT_before_min"] = rec[f"{c}_dRT_after_min"] = np.nan
        recs.append(rec)
    return pd.DataFrame(recs)


def main():
    os.makedirs(WORK, exist_ok=True)
    sys.stdout.reconfigure(encoding="utf-8")
    feats, cal, maps = load()
    hi = pd.read_csv(os.path.join(CMP, "lipid_detection_high_confidence.csv"))

    # self-check: per-lipid mode reproduces the delivered panel, point by point
    L = build("lipid", feats, cal, maps, hi)
    D = pd.read_csv(os.path.join(WORK, "high_confidence_panel.csv"))
    m = L.merge(D, on="vetted_names", suffixes=("", "_d"))
    assert len(L) == len(D) == len(m), (len(L), len(D), len(m))
    for c in NONREF:
        a, b = m[f"{c}_dRT_after_min"].to_numpy(float), m[f"{c}_dRT_after_min_d"].to_numpy(float)
        assert np.array_equal(np.isnan(a), np.isnan(b)) and np.nanmax(np.abs(a - b)) < 6e-5, c
        assert (m[f"{c}_used"] == m[f"{c}_used_d"].astype(str).str.lower().eq("true")).all(), c
    print(f"self-check: per-lipid mode reproduces the delivered panel ({len(L)} lipids, "
          f"{sum(L[f'{c}_used'].sum() for c in NONREF)} points)")

    I = build("ion", feats, cal, maps, hi)
    D = D.assign(display_group=D["class"])
    order = ["CAR", "LPC", "EtherLPC", "LPE", "EtherLPE", "PC", "EtherPC", "PE", "EtherPE", "PI", "SM", "Cer_NS",
             "HexCer_NS", "DG", "TG", "CE", "MG"]
    I["_g"] = I["class"].map({g: i for i, g in enumerate(order)}).fillna(len(order))
    I = I.sort_values(["_g", "column_22_RT_min", "ion"], kind="stable").drop(columns="_g")
    I.to_csv(os.path.join(WORK, "per_ion_panel.csv"), index=False, float_format="%.4f")
    two = I[I["class"].isin(["TG", "MG"])]
    print(f"per-ion panel: {len(I)} ion targets ({I.fold_unit.nunique()} lipids), "
          f"{sum(I[f'{c}_used'].sum() for c in NONREF)} non-reference points; TG / MG ion targets: "
          + ", ".join(f"{k} {v}" for k, v in two.ion.value_counts().items()))


if __name__ == "__main__":
    main()
