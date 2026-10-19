"""Accuracy table before and after calibration, from the 222 rows of vetted_panel_table.csv.

Writes ../vetted_panel_metrics.csv: one row per statistic, two columns (Before calibration, After
calibration), over the five non-reference columns pooled:

  Median Absolute Deviation (min)         median of |dRT|
  90th Percentile Absolute Error (min)    90th percentile of |dRT| (linear interpolation)
  Deviation Range (min)                   smallest to largest signed dRT
  Within +-0.2 min (%)                    share of points with |dRT| <= 0.2 min
  Within +-0.5 min (%)                    share of points with |dRT| <= 0.5 min

Points: ion-by-column cells, a peak shared by two rows counted once (the same points as the ECDF and violin
figures, via build_figures.py). dRT = RT - reference RT (the column_22 pick of the same row) before
calibration, calibrated RT - reference RT after. "After" is the delivered calibration, so it is in-sample:
most of these points are themselves stage-1 pairs of the RT Anchor runs. The reference column itself is not
a row (its dRT is 0 by design). The console log also lists the same statistics for each column and n.

Run:  python make_metrics_table.py     (needs make_final_table.py to have run; run_all.py does both)
"""
import os
import sys

import numpy as np
import pandas as pd

SCRIPTS = os.path.dirname(os.path.abspath(__file__))           # .../vetted_panel/scripts
HERE = os.path.dirname(SCRIPTS)                                # .../vetted_panel
LOGS = os.path.join(HERE, "logs")
OUT = os.path.join(HERE, "vetted_panel_metrics.csv")
sys.path.insert(0, SCRIPTS)
import build_figures as F                                      # noqa: E402  load_table / pooled_cells, as the figures

SLACK = 1e-9                                                   # RTs have 3 decimals: no float noise at the edges


def metrics(d):
    a = np.abs(d)
    return {"n": len(d),
            "Median Absolute Deviation (min)": round(float(np.median(a)), 3),
            "90th Percentile Absolute Error (min)": round(float(np.percentile(a, 90)), 3),
            "Deviation Range Min (min)": round(float(d.min()), 2),
            "Deviation Range Max (min)": round(float(d.max()), 2),
            "Within ±0.2 min (%)": round(100 * float((a <= 0.2 + SLACK).mean()), 1),
            "Within ±0.5 min (%)": round(100 * float((a <= 0.5 + SLACK).mean()), 1)}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    os.makedirs(LOGS, exist_ok=True)
    P = F.load_table()
    E, n_all = F.pooled_cells(P)
    rows = []                                    # long form, for the log: scope x stage
    scopes = [("All columns", E)] + [(F.col_name(c), E[E.column == c]) for c in F.NONREF]
    for scope, e in scopes:
        for stage, col in (("Before calibration", "before"), ("After calibration", "after")):
            rows.append({"Scope": scope, "Stage": stage, **metrics(e[col].to_numpy(float))})
    L = pd.DataFrame(rows)
    assert L.loc[L.Scope == "All columns", "n"].eq(len(E)).all() and L.groupby("Stage").n.sum().eq(2 * len(E)).all()

    def cell(stage, key):
        return L[(L.Scope == "All columns") & (L.Stage == stage)].iloc[0][key]

    out = {}
    for stage in ("Before calibration", "After calibration"):
        out[stage] = [cell(stage, "Median Absolute Deviation (min)"),
                      cell(stage, "90th Percentile Absolute Error (min)"),
                      f"{cell(stage, 'Deviation Range Min (min)'):.2f} to {cell(stage, 'Deviation Range Max (min)'):.2f}",
                      cell(stage, "Within ±0.2 min (%)"), cell(stage, "Within ±0.5 min (%)")]
    T = pd.DataFrame(out, index=["Median Absolute Deviation (min)", "90th Percentile Absolute Error (min)",
                                 "Deviation Range (min)", "Within ±0.2 min (%)", "Within ±0.5 min (%)"])
    T.index.name = "Statistic"
    T.to_csv(OUT, encoding="utf-8")
    text = (f"{len(P)} ion targets, {n_all} used non-reference cells, {len(E)} pooled points "
            f"(a peak shared by two rows counted once)\nwrote {os.path.basename(OUT)}:\n" + T.to_string()
            + "\n\nthe same statistics per scope and stage:\n" + L.to_string(index=False))
    print(text)
    with open(os.path.join(LOGS, "make_metrics_table_console.log"), "w", encoding="utf-8") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
