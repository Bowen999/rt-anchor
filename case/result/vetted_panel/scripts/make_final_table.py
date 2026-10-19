"""The final single table of the vetted-lipid panel: one row per ion target, RT before / after RT Anchor.

RT Anchor 1.2.3, anchor-free (stage 1 curve + series term + lattice term; the stage-2 sample anchors are
off), six columns calibrated onto the reference column_22. Same targets, picks and rules as
per_ion_panel.py (imported, not copied): an ion target is in the table when its peak is identity-checked
usable on the reference column and on >= 3 of the 6 columns, usable on column_20 and within 0.30 min of the
mapped vetted RT. The 219 high-confidence ion targets of per_ion_panel.csv come first; the low-confidence
list rows that pass the same rules (LPA / DMPE, 3 ion targets) follow, each with its reason.

Writes vetted_panel_table.csv, one row per ion target:

  lipid, class, formula           the matched list entry (the one whose mapped vetted RT is nearest the
                                  column_20 pick) and its neutral formula
  mz, adduct                      the list's own ion. The 492-row MS-DIAL sheet is negative ESI, so this is its
                                  [M-H]- / [M+HCOO]- with m/z recomputed from the formula; the 182-row
                                  manuscript sheet lists the positive precursor itself, so there mz / adduct
                                  equal the target
  target_mz, target_adduct        the positive-mode ion that was searched in the six serum tables
  confidence, confidence_reason   high / low as in data/compare/lipid_detection_*_confidence.csv. Low: the
                                  `reason` of that file for the matched list entry. High: the file's own
                                  definition (the m/z hit can point to the listed lipid)
  also_listed_as                  other list entries with the same target ion but another species name
                                  (isomers / isobars the m/z cannot tell apart)
  vetted_RT_min                   RT in the list (an earlier acquisition of the 20-min method)
  column_XX_RT_before_min         RT of the pick on that column, as in the MS-DIAL table
  column_XX_RT_after_min          Cal_RT_min of that feature (reference column_22: unchanged by design)

A column cell is empty when the lipid is not usable there (pick in < 3 of 4 injections, or not >= 2x
the next candidate, or no candidate within 15 ppm).
Column ids: column_12 = 12 min, column_16 = 17 min, column_20 = 20 min, column_22 = 22 min (reference),
column_41 = 45 min, column_85 = 90 min gradient.

Also writes work/final_table_ids.csv (Alignment ID of each cell, read by build_figures.py).

Run order:  python build_vetted_panel.py     (writes work/vetted_rt_mapping.csv, read here)
            python make_final_table.py       (or python run_all.py for the table, the figures and the logs)
"""
import os
import re
import sys

import numpy as np
import pandas as pd

SCRIPTS = os.path.dirname(os.path.abspath(__file__))           # .../vetted_panel/scripts
HERE = os.path.dirname(SCRIPTS)                                # .../vetted_panel
sys.path.insert(0, SCRIPTS)
import per_ion_panel as P                      # noqa: E402  load / units / pick and the rule constants

sys.path.insert(0, P.CMP)
from lipid_detection import ADDUCT, mass       # noqa: E402  theoretical adduct shifts

OUT = os.path.join(HERE, "vetted_panel_table.csv")
HIGH_REASON = "m/z hit can point to the listed lipid"
ORDER = ["CAR", "LPC", "EtherLPC", "LPE", "EtherLPE", "PC", "EtherPC", "PE", "EtherPE", "PI", "SM", "Cer_NS",
         "HexCer_NS", "DG", "TG", "CE", "MG"]


def species(name):
    return re.split(r"\s*\|\s*", name)[0].strip()


def rows_for(df, confidence, feats, cal, maps):
    """P.build(mode='ion') with the matched list entry kept, so its formula / adduct / RT / reason are known."""
    out = []
    for u in P.units(df, "ion"):
        rows = pd.DataFrame(u["rows"])
        picks = {c: P.pick(feats[c], u["key"]) for c in P.COLS}
        usable = {c: bool(picks[c] and picks[c]["usable"]) for c in P.COLS}
        sheet = ["492" if s.startswith("492") else "182" for s in rows.sheet]
        mapped = [float(np.interp(v, *maps[s])) for v, s in zip(rows.vetted_RT_min, sheet)]
        if picks["column_20"] is not None:
            d = [abs(picks["column_20"]["rt"] - m) for m in mapped]
            k = int(np.argmin(d))
        else:
            k, d = 0, [np.nan]
        if not (usable[P.REF] and sum(usable.values()) >= 3 and usable["column_20"] and d[k] <= P.ID_TOL + 1e-9):
            continue
        r = rows.iloc[k]
        (adduct_t, mz_t), = u["key"]
        if isinstance(r.vetted_neg_adduct, str):
            adduct, mz = r.vetted_neg_adduct, mass(r.formula) + ADDUCT[r.vetted_neg_adduct]
        else:
            adduct, mz = adduct_t, mz_t
        lipid = species(r["name"])
        rec = {"lipid": lipid, "class": r["class"], "formula": r.formula,
               "mz": round(mz, 4), "adduct": adduct, "target_mz": mz_t, "target_adduct": adduct_t,
               "confidence": confidence,
               "confidence_reason": r["reason"] if confidence == "low" else HIGH_REASON,
               "also_listed_as": "; ".join(dict.fromkeys(
                   species(n) for n in rows["name"] if species(n) != lipid)),
               "vetted_RT_min": r.vetted_RT_min, "_n_used": sum(usable.values())}
        for c in P.COLS:
            rec[f"_aid_{c}"] = picks[c]["aid"] if usable[c] else np.nan
            if usable[c]:
                rec[f"{c}_RT_before_min"] = picks[c]["rt"]
                rec[f"{c}_RT_after_min"] = float(cal[c].loc[picks[c]["aid"], "Cal_RT_min"])
            else:
                rec[f"{c}_RT_before_min"] = rec[f"{c}_RT_after_min"] = np.nan
        out.append(rec)
    return pd.DataFrame(out)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    feats, cal, maps = P.load()
    hi = pd.read_csv(os.path.join(P.CMP, "lipid_detection_high_confidence.csv"))
    lo = pd.read_csv(os.path.join(P.CMP, "lipid_detection_low_confidence.csv"))
    lo = lo[lo.positive_target_ions.notna()]            # the 14 FA rows have no positive-mode ion

    H = rows_for(hi, "high", feats, cal, maps)
    L = rows_for(lo, "low", feats, cal, maps)

    # self-check 1: the high rows are exactly the panel of per_ion_panel.py (same ions, same picks, same RTs)
    ref = P.build("ion", feats, cal, maps, hi)
    assert len(H) == len(ref), (len(H), len(ref))
    key = lambda d, ion, name: d.set_index([name, ion]).sort_index()   # noqa: E731
    a = key(H.rename(columns={"target_adduct": "ion"}), "ion", "lipid")
    b = key(ref, "ion", "lipid")
    assert a.index.equals(b.index), "ion targets differ from per_ion_panel.py"
    for c in P.COLS:
        assert np.allclose(a[f"{c}_RT_before_min"], b[f"{c}_RT_min"].where(b[f"{c}_used"].astype(bool)),
                           equal_nan=True, atol=1e-9), c
    for c in P.NONREF:
        assert np.allclose(a[f"{c}_RT_after_min"], b[f"{c}_cal_RT_min"], equal_nan=True, atol=1e-9), c
    # self-check 2: the reference column is unchanged by its own calibration
    d = (H[f"{P.REF}_RT_after_min"] - H[f"{P.REF}_RT_before_min"]).abs().max()
    assert d < 1e-6, d

    T = pd.concat([H, L], ignore_index=True)
    T["_g"] = T["class"].map({g: i for i, g in enumerate(ORDER)}).fillna(len(ORDER))
    T["_c"] = (T.confidence == "low").astype(int)
    T = T.sort_values(["_c", "_g", f"{P.REF}_RT_before_min", "target_adduct"], kind="stable")
    cols = (["lipid", "class", "formula", "mz", "adduct", "target_mz", "target_adduct", "confidence",
             "confidence_reason", "also_listed_as", "vetted_RT_min"]
            + [f"{c}_RT_{w}_min" for c in P.COLS for w in ("before", "after")])
    ids = T[["lipid", "target_adduct"] + [f"_aid_{c}" for c in P.COLS]].rename(
        columns={f"_aid_{c}": f"{c}_alignment_id" for c in P.COLS})
    T = T[cols].copy()
    for c in cols:
        if c.endswith("_before_min") or c == "vetted_RT_min":
            T[c] = T[c].round(3)
        elif c.endswith("_after_min"):
            T[c] = T[c].round(4)
    T.to_csv(OUT, index=False, encoding="utf-8")
    os.makedirs(P.WORK, exist_ok=True)          # the pick of every cell, so that build_figures.py counts a peak once
    ids.to_csv(os.path.join(P.WORK, "final_table_ids.csv"), index=False)
    msg = [f"self-check: the {len(H)} high-confidence rows reproduce per_ion_panel.py "
           f"(ions, picks, RT before / after); reference column unchanged (max |d| {d:.1e} min)",
           f"wrote {os.path.basename(OUT)}: {len(T)} ion targets = {len(H)} high + {len(L)} low  "
           f"(expected 219 + 3: summary_funnel_per_ion.csv, 'panel' row)"]
    print("\n".join(msg))
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    with open(os.path.join(HERE, "logs", "make_final_table_console.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(msg) + "\n")


if __name__ == "__main__":
    main()
