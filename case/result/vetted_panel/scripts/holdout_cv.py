"""Held-out (10-fold cross-validated) calibration error for the high-confidence vetted panel.

Why: 83 % of the panel's lipid-by-column points are themselves kept stage-1 pairs of the RT Anchor runs
(the same column feature paired with the same reference feature), so the delivered "after calibration"
errors are mostly in-sample. Here every panel lipid is predicted by a calibration that never saw it.

How: the 182 high-confidence lipids are split into 10 folds (sorted by display group and reference RT, then
dealt round-robin, so every fold spans the classes and the gradient). For each fold the reference serum table
(column_22) is rewritten without the fold's reference peaks and without every reference feature eluting
within 0.04 min of them (the lipid's adducts, isotopes and in-source fragments - and, conservatively, anything
else co-eluting). The five non-reference columns are then calibrated against this reduced reference with
RT Anchor 1.2.3, default settings, exactly as in ../run_all.sh (only the reference sample table differs).
The column tables are untouched, so each held-out lipid keeps its pick and gets a calibrated RT from a
curve (and series / lattice terms) fitted without any pair that involves it; its error is measured against
the reference RT of its original pick.

Writes (../work/): holdout_cv_cells.csv (one row per used lipid-by-column point, in-sample and held-out
values side by side) and holdout_cv_folds.csv (per fold and column: features removed, pairs, gate outcomes).
The 50 engine runs go to a work folder that is deleted afterwards (each run's calibrated table is read, the
needed rows are kept). About 3-5 minutes.

Run:  python holdout_cv.py
      python holdout_cv.py --panel per_ion_panel.csv --prefix holdout_cv_per_ion
      (any panel CSV in the delivered column layout; with a `fold_unit` column, the rows of one unit - the
      two ions of one TG - are held out together)
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

SCRIPTS = os.path.dirname(os.path.abspath(__file__))           # .../vetted_panel/scripts
HERE = os.path.dirname(SCRIPTS)                                # .../vetted_panel
WORK = os.path.join(HERE, "work")                              # intermediate tables (removed by run_all.py)
RUN = os.path.dirname(HERE)                                           # result_v1.2.3_reprocessed
BASE = os.path.dirname(RUN)
SRC = os.environ.get("RT_ANCHOR_SRC", os.path.join(BASE, "software", "rt-anchor", "rt_anchor", "src"))
D = os.path.join(BASE, "data")
REF_S = os.path.join(D, "22 min", "Serum", "22 min Serum Area_2_2026_10_12_13_24_24.csv")
REF_T = os.path.join(D, "22 min", "STDs", "Res", "22 min stds Area_1_2026_08_16_23_14_28.csv")
SAMPLES = {
    "column_12": os.path.join(D, "12 min", "Serum", "12 min Serum Area_2_2026_10_12_13_18_31.csv"),
    "column_16": os.path.join(D, "17 min", "Serum", "17 Serum Area_2_2026_10_12_13_33_45.csv"),
    "column_20": os.path.join(D, "20 min", "Serum", "20 min Serum Area_2_2026_10_12_12_12_37.csv"),
    "column_41": os.path.join(D, "45 min", "Serum", "45 Serum Area_2_2026_10_12_12_48_51.csv"),
    "column_85": os.path.join(D, "90 min", "Serum", "90 Serum Area_2_2026_10_12_13_11_13.csv"),
}
STANDARDS = {
    "column_12": os.path.join(D, "12 min", "stds", "12 min STDs Area_1_2026_08_28_02_03_07.csv"),
    "column_16": os.path.join(D, "17 min", "STDs", "17 STDs Area_1_2026_08_29_19_34_08.csv"),
    "column_20": os.path.join(D, "20 min", "STDs", "Results", "20 min STDs Area_1_2026_08_17_00_20_03.csv"),
    "column_41": os.path.join(D, "45 min", "STDS", "45 stds Area_1_2026_08_29_17_12_33.csv"),
    "column_85": os.path.join(D, "90 min", "STDs", "90 STDS Area_1_2026_08_27_23_47_07.csv"),
}
NONREF = list(SAMPLES)
N_FOLDS = 10
COELUTE = 0.04          # min; reference features this close to a held-out reference peak are removed too


def write_reduced_reference(drop_ids, path):
    """Copy the MS-DIAL reference table byte for byte, leaving out the data lines whose Alignment ID is in
    drop_ids (4 meta lines + the header are kept)."""
    with open(REF_S, "rb") as fh:
        lines = fh.read().split(b"\n")
    out, dropped = [], 0
    for i, line in enumerate(lines):
        if i >= 5 and line.strip():
            aid = int(line.split(b",", 1)[0])
            if aid in drop_ids:
                dropped += 1
                continue
        out.append(line)
    with open(path, "wb") as fh:
        fh.write(b"\n".join(out))
    return dropped


def gate_status(log_text, term):
    m = re.search(rf"{term} term: ([^;\n]*)", log_text)
    s = m.group(1) if m else "?"
    return "engaged" if "engaged" in s else ("off" if ("gated off" in s or "no correction" in s) else s)


def main():
    os.makedirs(WORK, exist_ok=True)
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default="high_confidence_panel.csv")
    ap.add_argument("--prefix", default="holdout_cv")
    args = ap.parse_args()
    P = pd.read_csv(os.path.join(WORK, args.panel))
    ref = pd.read_csv(REF_S, skiprows=4, low_memory=False)
    if "fold_unit" in P.columns:      # deal whole units (all ions of a lipid) round-robin
        U = (P.groupby("fold_unit", sort=False)
             .agg(display_group=("display_group", "first"), rt=("column_22_RT_min", "min"))
             .sort_values(["display_group", "rt"], kind="stable"))
        fold_of = dict(zip(U.index, np.arange(len(U)) % N_FOLDS))
        P["fold"] = P.fold_unit.map(fold_of)
    else:                             # one row per lipid (the delivered panel)
        order = P.sort_values(["display_group", "column_22_RT_min"], kind="stable").index
        P["fold"] = -1
        P.loc[order, "fold"] = np.arange(len(order)) % N_FOLDS

    work = tempfile.mkdtemp(prefix="rt_anchor_cv_")
    env = dict(os.environ, PYTHONPATH=SRC)
    cells, folds = [], []
    try:
        for k in range(N_FOLDS):
            F = P[P.fold == k]
            picks = F.column_22_RT_min.to_numpy(float)
            near = np.abs(ref["Average Rt(min)"].to_numpy(float)[:, None] - picks[None, :]).min(axis=1) <= COELUTE + 1e-9
            drop = set(ref.loc[near, "Alignment ID"].astype(int)) | set(F.column_22_alignment_id.astype(int))
            red = os.path.join(work, f"reference_fold{k}.csv")
            n_drop = write_reduced_reference(drop, red)
            assert n_drop == len(drop), (k, n_drop, len(drop))
            for c in NONREF:
                out = os.path.join(work, f"fold{k}", c, "rt_anchor_run")
                os.makedirs(os.path.dirname(out), exist_ok=True)
                cmd = [sys.executable, "-m", "rt_anchor.cli", "calibrate", "--samples", SAMPLES[c],
                       "--standards", STANDARDS[c], "--polarity", "positive", "--panel", "mix21",
                       "--reference-sample", red, "--reference-standards", REF_T, "--no-report", "--out", out]
                r = subprocess.run(cmd, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace")
                if r.returncode != 0:
                    raise RuntimeError(f"fold {k} {c}: exit {r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
                cal = pd.read_csv(out + "_calibrated.csv", low_memory=False).set_index("Alignment ID")
                pairs = pd.read_csv(out + "_pairs.csv")
                with open(out + "_log.txt", encoding="utf-8", errors="replace") as fh:
                    log = fh.read()
                folds.append(dict(fold=k, column=c, n_lipids_held_out=len(F), n_reference_features_removed=n_drop,
                                  n_pairs=len(pairs), n_pairs_kept=int(pairs.kept.sum()),
                                  series_term=gate_status(log, "series"), lattice_term=gate_status(log, "lattice")))
                for _, row in F[F[f"{c}_used"].astype(str).str.lower() == "true"].iterrows():
                    aid = int(row[f"{c}_alignment_id"])
                    q = cal.loc[aid]
                    assert abs(float(q["Average Rt(min)"]) - row[f"{c}_RT_min"]) < 1e-3
                    cells.append(dict(
                        fold=k, lipid=row.lipid, **({"ion": row.ion} if "ion" in row.index else {}),
                        display_group=row.display_group, column=c, alignment_id=aid,
                        RT_min=row[f"{c}_RT_min"], reference_RT_min=row.column_22_RT_min,
                        dRT_before_min=row[f"{c}_dRT_before_min"],
                        cal_RT_insample_min=row[f"{c}_cal_RT_min"], dRT_after_insample_min=row[f"{c}_dRT_after_min"],
                        cal_RT_heldout_min=float(q["Cal_RT_min"]),
                        dRT_after_heldout_min=float(q["Cal_RT_min"]) - row.column_22_RT_min,
                        warp_source_heldout=q["warp_source"]))
                shutil.rmtree(os.path.dirname(out))
                print(f"fold {k} {c}: {len(F)} lipids held out, {n_drop} reference features removed, "
                      f"{len(pairs)} pairs, series {folds[-1]['series_term']}, lattice {folds[-1]['lattice_term']}",
                      flush=True)
            os.remove(red)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    C = pd.DataFrame(cells)
    C.to_csv(os.path.join(WORK, f"{args.prefix}_cells.csv"), index=False, float_format="%.4f")
    Fd = pd.DataFrame(folds)
    Fd.to_csv(os.path.join(WORK, f"{args.prefix}_folds.csv"), index=False)
    a_in, a_out = C.dRT_after_insample_min.abs(), C.dRT_after_heldout_min.abs()
    n_exp = sum(P[f"{c}_used"].astype(str).str.lower().eq("true").sum() for c in NONREF)
    print(f"\n{len(C)} held-out points (expected {n_exp}); median |dRT_after| in-sample {a_in.median():.4f} -> held-out "
          f"{a_out.median():.4f} min; within 0.2 min {(a_in <= 0.2).mean():.1%} -> {(a_out <= 0.2).mean():.1%}")
    print("gate outcomes per column (folds engaged / 10):")
    for c in NONREF:
        f = Fd[Fd.column == c]
        print(f"  {c}: series {(f.series_term == 'engaged').sum()}, lattice {(f.lattice_term == 'engaged').sum()}")


if __name__ == "__main__":
    main()
