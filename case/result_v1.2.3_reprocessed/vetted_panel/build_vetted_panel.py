"""Vetted-lipid validation panel for the RT Anchor 1.2.3 run on the reprocessed data.

Six LC methods (one QTOF, positive ESI, human serum) are calibrated onto the reference method
column_22 with RT Anchor 1.2.3 (default settings). The manually vetted lipid list was matched by
precursor m/z only to the six serum tables (data/compare/lipid_detection_*_confidence.csv). This
script turns that list into a panel of lipids whose peak can be pinned down on each column, and
tabulates the retention-time deltas against the reference column before and after calibration.

Inputs  (read only):  data/<n> min/Serum/*Area*.csv           MS-DIAL alignment tables, six columns
                      result_v1.2.3_reprocessed/column_*/rt_anchor_run_calibrated.csv
                      data/compare/lipid_detection_{high,low}_confidence.csv
Outputs (this folder): high_confidence_panel.csv, high_and_low_confidence_panel.csv,
                      excluded_lipids.csv, vetted_rt_mapping.csv, vetted_rt_check_groups.csv (the
                      per-group points of the QC figure, so build_figures.py reads CSVs only),
                      build_vetted_panel_console.log

Rules, in one place (the numbers in brackets are the ones this script prints next to the expected ones):
 1. Targets: the 4-decimal m/z values of `positive_target_ions`, each re-computed from `formula`.
 2. Groups: rows of one file with the identical set of (ion, 4-dp m/z) form one panel candidate.
 3. Candidates: features within 15 ppm of any target of the group; pick = largest mean area (tie:
    smallest Alignment ID); rivals = the other candidates except the other adduct of the same
    compound (matched only through an ion that did not match the pick and within 0.10 min of it);
    dominance = pick area / largest rival area; usable = pick in >= 3 of 4 injections and
    dominance >= 2.
 4. Mapping of the vetted RT (an earlier acquisition of the 20-min method) onto today's column_20
    run: one robust LOWESS curve per vetted-list sheet (492 / 182).
 5. Matched row = the row whose mapped vetted RT is nearest the column_20 pick.
 6. Panel = reference column usable, >= 3 of 6 columns usable, column_20 usable, and the column_20
    pick within 0.30 min of the mapped vetted RT. Everything else goes to excluded_lipids.csv.
 7. Used non-reference cell = a column on which the lipid is usable: dRT_before = RT - reference RT,
    dRT_after = calibrated RT - reference RT.

Conventions where the rules leave a choice (all visible in the CSVs):
 - confidence is "high" / "low"; the <col>_used flags are True / False; "inf" for an infinite dominance.
 - n_columns_used counts the columns (of the 6) on which the lipid is usable, the reference included;
   n_columns_series_corrected / _lattice_corrected count used non-reference cells with a correction != 0.
 - A column without any candidate has _used False, _n_candidates 0 and blank pick fields.
 - "within 0.10 min" (other adduct), "within 0.2 min" and the 15 ppm window are inclusive; RT differences
   are compared with 1e-9 min of slack so that 3-decimal RTs do not fall victim to float noise.
 - excluded_lipids.csv lists the high groups, then the low groups (each in list order of first
   appearance), then the 14 FA rows; the FA rows were not searched, so their counts / dominance / RT are blank.
 - Dominance is printed with 2 decimals, so a printed 2.00 can still be a pick that is not usable
   (1.9993 on column_16 for PE 38:4); the decisions use the unrounded value.

Run:  python build_vetted_panel.py        (a few seconds; writes next to this script)
"""
import collections
import csv
import io
import os
import re
import sys
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

HERE = os.path.dirname(os.path.abspath(__file__))           # .../result_v1.2.3_reprocessed/vetted_panel
RUN_DIR = os.path.dirname(HERE)                             # .../result_v1.2.3_reprocessed
BASE = os.path.dirname(RUN_DIR)                             # project root
DATA = os.path.join(BASE, "data")

COLS = ["column_12", "column_16", "column_20", "column_22", "column_41", "column_85"]
REF = "column_22"
NONREF = [c for c in COLS if c != REF]
SERUM = {
    "column_12": "12 min/Serum/12 min Serum Area_2_2026_10_12_13_18_31.csv",
    "column_16": "17 min/Serum/17 Serum Area_2_2026_10_12_13_33_45.csv",
    "column_20": "20 min/Serum/20 min Serum Area_2_2026_10_12_12_12_37.csv",
    "column_22": "22 min/Serum/22 min Serum Area_2_2026_10_12_13_24_24.csv",
    "column_41": "45 min/Serum/45 Serum Area_2_2026_10_12_12_48_51.csv",
    "column_85": "90 min/Serum/90 Serum Area_2_2026_10_12_13_11_13.csv",
}
COMPARE = {"high": "lipid_detection_high_confidence.csv", "low": "lipid_detection_low_confidence.csv"}

# --- rule constants
TOL = 15e-6                     # m/z tolerance, relative
MIN_DET = 3                     # injections (of 4) the pick must be detected in
MIN_DOM = 2.0                   # dominance for "usable"
OTHER_ADDUCT_DRT = 0.10         # min; a candidate matched only through the other adduct, this close
TRAIN_DOM = 5.0                 # dominance needed for a training pair when there are > 1 candidates
LOWESS_FRAC, LOWESS_IT, FIT_ROUNDS, KEEP_TOL = 0.5, 3, 3, 0.5
MAX_IDENTITY_DELTA = 0.30       # min
WITHIN = 0.2                    # min, the headline accuracy threshold
FLOAT_EDGE = 1e-9               # guard when comparing differences of 3-decimal RTs to a threshold

# target self-check constants (rule 1)
AM = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956,
      "P": 30.97376163, "S": 31.97207100, "Na": 22.9897692809}
ELECTRON = 0.00054857990946
PROTON = AM["H"] - ELECTRON
ION_SHIFT = {"[M+H]+": PROTON,
             "[M+NH4]+": AM["N"] + 4 * AM["H"] - ELECTRON,
             "[M+Na]+": AM["Na"] - ELECTRON,
             "[M+H-H2O]+": PROTON - (2 * AM["H"] + AM["O"])}
TARGET_CHECK_TOL = 6e-5         # Da

GROUP_ORDER = ["CAR", "LPC", "EtherLPC", "LPE", "EtherLPE", "PC", "EtherPC", "PE", "EtherPE", "PI",
               "SM", "Cer_NS", "HexCer_NS", "DG", "TG", "CE", "Other", "Low confidence"]
MIN_CLASS_LIPIDS = 5            # classes with fewer lipids in the high panel -> "Other" (owner's choice for
                                # the ring heatmaps, 2026-10-14; was 3)
SAME_ION_RE = re.compile(r"\[M\+H\]\+ of ((?:PE|LPE) \d+:\d+)")

REASON = {
    "a": "no candidate within 15 ppm on column_22",
    "b": "column_22 pick detected in fewer than 3 of 4 injections",
    "c": "column_22 pick less than 2x the next candidate",
    "d": "usable on fewer than 3 columns",
    "e": "no candidate within 15 ppm on column_20",
    "f": "column_20 pick detected in fewer than 3 of 4 injections",
    "g": "column_20 pick less than 2x the next candidate",
    "h": "column_20 pick more than 0.30 min from the vetted RT",
    "fa": "no positive-mode target ion (not searched)",
}

# ------------------------------------------------------------------------------------- logging
LOG = []
DIFFS = []
# float edge cases met on the way (reported at the end; every kind is listed, also when it was met 0 times)
EDGE = collections.Counter({
    "mz at the 15 ppm edge (relative distance within 1e-6)": 0,
    "picks decided by an area tie (smallest Alignment ID)": 0,
    "candidate pairs |dRT| within 1e-6 of 0.10 min (other-adduct test)": 0,
    "training pair x round with |residual| within 1e-6 of 0.5 min": 0,
    "identity deltas within 1e-6 of 0.30 min": 0,
    "cells with |dRT_before| within 1e-6 of 0.2 min": 0,
    "cells with |dRT_after| within 1e-6 of 0.2 min": 0,
    "used high-panel corrections that are non-zero but print as 0.0000": 0})
EDGE_NOTES = []                  # one line per edge case that really was met (which lipid, which value)


def log(*parts):
    line = " ".join(str(p) for p in parts)
    print(line)
    LOG.append(line)


def expect(label, expected, got, nd=None):
    """Print one 'expected vs obtained' line; remember the difference, never adjust a rule for it."""
    if nd is None:
        same = got == expected
    else:
        same = round(float(got), nd) == expected
    if not same:
        DIFFS.append(label)
    shown = f"{got:.{nd}f}" if nd is not None else str(got)
    exp_s = f"{expected:.{nd}f}" if nd is not None else str(expected)
    log(f"  {'ok       ' if same else 'DIFFERENT'}  {label:<58} expected {exp_s:>9}   obtained {shown:>9}")
    return same


def write_text(path, text):
    """Write a UTF-8 file with LF newlines. The target is a OneDrive folder: a stale target can be a
    sync placeholder that refuses to open, and the sync client briefly locks fresh files, so remove
    the old file and retry a few times."""
    for attempt in range(8):
        try:
            if os.path.exists(path):
                os.remove(path)
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text)
            return
        except OSError:
            if attempt == 7:
                raise
            time.sleep(1.0)


def format_value(v, kind):
    """Write-time rounding: RT / m/z / dRT / corrections 4 decimals, dominance 2 decimals ('inf' when
    infinite), mean area 1 decimal; counts and IDs as integers; blank for missing."""
    if v is None or (isinstance(v, (float, np.floating)) and np.isnan(v)):
        return ""
    if kind == "str":
        return str(v)
    if kind == "bool":
        return "True" if v else "False"
    if kind == "int":
        return str(int(v))
    v = float(v)
    if np.isinf(v):
        return "inf" if v > 0 else "-inf"
    nd = {"rt": 4, "dom": 2, "area": 1}[kind]
    return f"{round(v, nd) + 0.0:.{nd}f}"            # + 0.0 turns -0.0 into 0.0


STR_COLS = {"lipid", "class", "display_group", "confidence", "same_ion_as", "vetted_names",
            "matched_vetted_name", "formula", "positive_target_ions", "reason", "sheet", "status"}


def kind_of(name):
    if name in STR_COLS:
        return "str"
    if (name.startswith("n_") or name.endswith(("_alignment_id", "_n_candidates", "_n_injections_detected"))):
        return "int"                                  # (before the '_used' test: n_columns_used is a count)
    if name.endswith("_used"):
        return "bool"
    if name.endswith("_dominance"):
        return "dom"
    if name.endswith("_mean_area"):
        return "area"
    return "rt"


def write_csv(path, header, records):
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(header)
    for rec in records:
        w.writerow([format_value(rec.get(h), kind_of(h)) for h in header])
    write_text(path, buf.getvalue())


# ------------------------------------------------------------------------------------- inputs
def read_serum_tables():
    """feature arrays per column from the serum table + calibrated table (same rows, same order)."""
    feat = {}
    for c in COLS:
        t = pd.read_csv(os.path.join(DATA, SERUM[c]), skiprows=4, low_memory=False)
        names = list(t.columns)
        i = names.index("MS/MS spectrum")
        inj = names[i + 1:i + 5]
        assert len(inj) == 4 and len(names) - (i + 5) == 2, \
            f"{c}: expected 4 injection columns + 2 more (average, stdev) after 'MS/MS spectrum'"
        area = t[inj].to_numpy(dtype=float)
        assert np.isfinite(area).all() and (area >= 0).all(), f"{c}: bad injection areas"
        mean_area = area.mean(axis=1)
        avg = t[names[i + 5]].to_numpy(dtype=float)
        assert np.abs(avg - mean_area).max() < 1.0, f"{c}: MS-DIAL average != mean of the 4 columns"
        aid = t["Alignment ID"].to_numpy()
        assert len(np.unique(aid)) == len(aid), f"{c}: duplicate Alignment IDs"

        cal = pd.read_csv(os.path.join(RUN_DIR, c, "rt_anchor_run_calibrated.csv"),
                          usecols=["Alignment ID", "Average Rt(min)", "Cal_RT_min",
                                   "series_correction_min", "lattice_correction_min"],
                          low_memory=False)
        assert len(cal) == len(t), f"{c}: calibrated table has a different number of rows"
        assert (cal["Alignment ID"].to_numpy() == aid).all(), f"{c}: Alignment ID differs row by row"
        assert np.array_equal(cal["Average Rt(min)"].to_numpy(), t["Average Rt(min)"].to_numpy()), \
            f"{c}: Average Rt(min) differs between serum and calibrated table"
        feat[c] = dict(
            aid=aid, rt=t["Average Rt(min)"].to_numpy(dtype=float), mz=t["Average Mz"].to_numpy(dtype=float),
            area=mean_area, ndet=(area > 0).sum(axis=1),
            cal_rt=cal["Cal_RT_min"].to_numpy(dtype=float),
            series=cal["series_correction_min"].to_numpy(dtype=float),
            lattice=cal["lattice_correction_min"].to_numpy(dtype=float))
    return feat


def mass(formula):
    counts = collections.Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[el] += int(n) if n else 1
    return sum(AM[el] * n for el, n in counts.items())


def parse_targets(text):
    """'[M+NH4]+ 736.6450; [M+Na]+ 741.6004' -> {'[M+NH4]+': 736.645, '[M+Na]+': 741.6004}."""
    out = {}
    if text.strip() == "":
        return out
    for part in text.split(";"):
        ion, val = part.strip().rsplit(" ", 1)
        assert re.fullmatch(r"\d+\.\d{4}", val), f"target not given with 4 decimals: {part!r}"
        out[ion] = float(val)
    return out


def read_vetted():
    """All 674 list rows (high file first, then low file), in CSV order."""
    rows = []
    for conf, fn in COMPARE.items():
        df = pd.read_csv(os.path.join(DATA, "compare", fn), keep_default_na=False)
        for i, r in enumerate(df.to_dict("records")):
            assert r["sheet"].startswith(("492", "182")), r["sheet"]
            rows.append(dict(
                conf=conf, i=i, sheet=r["sheet"][:3], name=r["name"], cls=r["class"],
                formula=r["formula"], rt=float(r["vetted_RT_min"]), ions_str=r["positive_target_ions"],
                ions=parse_targets(r["positive_target_ions"]), reason=r.get("reason", ""),
                n_exp=[int(x) for x in r["n_candidate_features_per_column"].split(";")]))
    return rows


# ------------------------------------------------------------------------------------- search
def search(feat, targets, col):
    """Candidates, pick, rivals, dominance and the 'usable' flag of one group on one column."""
    f = feat[col]
    hit = {}
    for ion, t in targets.items():
        d = np.abs(f["mz"] - t)
        tol = TOL * t
        hit[ion] = d <= tol                                          # inclusive
        EDGE["mz at the 15 ppm edge (relative distance within 1e-6)"] += int((np.abs(d / tol - 1) < 1e-6).sum())
    mask = np.any(np.vstack(list(hit.values())), axis=0)
    idx = np.flatnonzero(mask)
    n = len(idx)
    if n == 0:
        return dict(n=0, idx=None, aid=np.nan, mz=np.nan, rt=np.nan, ndet=np.nan, area=np.nan,
                    dom=np.nan, usable=False, tie=False, edge_pairs=[])
    order = np.lexsort((f["aid"][idx], -f["area"][idx]))              # largest area, then smallest ID
    p = int(idx[order[0]])
    top = f["area"][idx].max()
    tie = int((f["area"][idx] == top).sum()) > 1
    if tie:
        EDGE["picks decided by an area tie (smallest Alignment ID)"] += 1
    ions_p = {ion for ion in targets if hit[ion][p]}
    rival_area, edge_pairs = [], []
    for k in idx:
        k = int(k)
        if k == p:
            continue
        ions_k = {ion for ion in targets if hit[ion][k]}
        drt = abs(f["rt"][k] - f["rt"][p])
        if ions_k.isdisjoint(ions_p):
            if abs(drt - OTHER_ADDUCT_DRT) < 1e-6:
                EDGE["candidate pairs |dRT| within 1e-6 of 0.10 min (other-adduct test)"] += 1
                edge_pairs.append((float(drt), bool(drt <= OTHER_ADDUCT_DRT), bool(drt <= OTHER_ADDUCT_DRT + FLOAT_EDGE)))
            if drt <= OTHER_ADDUCT_DRT + FLOAT_EDGE:
                continue                                              # the other adduct of the same compound
        rival_area.append(f["area"][k])
    if not rival_area or max(rival_area) == 0:
        dom = np.inf
    else:
        dom = float(f["area"][p] / max(rival_area))
    ndet = int(f["ndet"][p])
    return dict(n=n, idx=p, aid=int(f["aid"][p]), mz=float(f["mz"][p]), rt=float(f["rt"][p]),
                ndet=ndet, area=float(f["area"][p]), dom=dom, usable=bool(ndet >= MIN_DET and dom >= MIN_DOM),
                tie=tie, edge_pairs=edge_pairs)


# ------------------------------------------------------------------------------------- mapping
def fit_curve(x, y):
    """Robust LOWESS curve of the vetted RT onto the column_20 RT: three rounds; after each fit the
    curve (duplicate fitted x averaged) is evaluated at ALL training pairs and keep = |y - curve| <= 0.5.
    Returns the third curve (cx, cy) and the keep mask after the third fit."""
    keep = np.ones(len(x), dtype=bool)
    for _ in range(FIT_ROUNDS):
        fit = lowess(y[keep], x[keep], frac=LOWESS_FRAC, it=LOWESS_IT, return_sorted=True)
        cx, inv = np.unique(fit[:, 0], return_inverse=True)
        cy = np.bincount(inv, weights=fit[:, 1]) / np.bincount(inv)
        res = np.abs(y - np.interp(x, cx, cy))
        EDGE["training pair x round with |residual| within 1e-6 of 0.5 min"] += int((np.abs(res - KEEP_TOL) < 1e-6).sum())
        keep = res <= KEEP_TOL
    return cx, cy, keep


# ------------------------------------------------------------------------------------- main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    t0 = time.time()
    log("== RT Anchor 1.2.3 (reprocessed data): vetted-lipid validation panel ==")
    feat = read_serum_tables()
    log("serum tables (rows) / calibrated tables read, row alignment asserted:")
    for c in COLS:
        log(f"  {c}: {len(feat[c]['aid'])} features, 4 injection columns + 2 summary columns asserted")

    rows = read_vetted()
    n_rows = {k: sum(r["conf"] == k for r in rows) for k in COMPARE}
    log(f"vetted list rows: high {n_rows['high']}, low {n_rows['low']}, total {len(rows)}")

    # ---- rule 1: targets recomputed from formula + ion
    worst = 0.0
    for r in rows:
        for ion, t in r["ions"].items():
            diff = abs(mass(r["formula"]) + ION_SHIFT[ion] - t)
            assert diff <= TARGET_CHECK_TOL, f"target mismatch {r['name']} {ion}: {diff}"
            worst = max(worst, diff)
    no_target = [r for r in rows if not r["ions"]]
    assert all(r["conf"] == "low" and r["cls"] == "FA" for r in no_target)
    log(f"rule 1 self-check: every target equals formula + ion (largest |difference| {worst:.2e} Da, "
        f"limit {TARGET_CHECK_TOL:.0e}); rows without targets: {len(no_target)} (FA, low file)")

    # ---- rule 2: groups
    groups = []
    for conf in COMPARE:
        seen = {}
        for r in rows:
            if r["conf"] != conf or not r["ions"]:
                continue
            key = frozenset((ion, f"{v:.4f}") for ion, v in r["ions"].items())
            if key not in seen:
                seen[key] = SimpleNamespace(conf=conf, rows=[], targets=dict(r["ions"]))
                groups.append(seen[key])
            seen[key].rows.append(r)
    n_groups = {k: sum(g.conf == k for g in groups) for k in COMPARE}
    n_low_rows = sum(len(g.rows) for g in groups if g.conf == "low")
    log("rule 2: groups (identical set of ion and 4-dp m/z within one file)")
    expect("high groups", 352, n_groups["high"])
    expect("low groups", 36, n_groups["low"])
    expect("low rows with targets", 59, n_low_rows)

    # ---- rule 3: candidates, picks; self-check against n_candidate_features_per_column
    for g in groups:
        g.hits = {c: search(feat, g.targets, c) for c in COLS}
        counts = [g.hits[c]["n"] for c in COLS]
        for r in g.rows:
            assert counts == r["n_exp"], f"candidate counts differ from the compare file: {r['name']} {counts} vs {r['n_exp']}"
    log(f"rule 3 self-check: candidate counts of all {sum(len(g.rows) for g in groups)} searched rows over "
        f"the 6 columns equal n_candidate_features_per_column")
    for g in groups:
        for c in COLS:
            h = g.hits[c]
            for drt, plain, guarded in h["edge_pairs"]:
                EDGE_NOTES.append(
                    f"{g.rows[0]['name'].split('|')[0].strip()} ({g.conf}) on {c}: two candidates {drt!r} min apart, "
                    f"other-adduct test (<= 0.10 min): counted as {'other adduct' if guarded else 'rival'} "
                    f"(a plain float comparison decides {'the same' if plain == guarded else 'DIFFERENTLY'})")
            if h["tie"]:
                EDGE_NOTES.append(f"{g.rows[0]['name'].split('|')[0].strip()} ({g.conf}) on {c}: pick decided by an area tie")

    # ---- rule 4: mapping of the vetted RT onto today's column_20 run
    pairs = {"492": [], "182": []}
    for g in groups:
        if g.conf != "high":
            continue
        sheets = {r["sheet"] for r in g.rows}
        rts = {round(r["rt"], 3) for r in g.rows}
        h20 = g.hits["column_20"]
        if len(sheets) == 1 and len(rts) == 1 and h20["idx"] is not None and (h20["n"] == 1 or h20["dom"] >= TRAIN_DOM):
            pairs[next(iter(sheets))].append((next(iter(rts)), h20["rt"]))
    curve, n_train, n_kept = {}, {}, {}
    for sheet in ("492", "182"):
        pts = sorted(pairs[sheet])
        x = np.array([p[0] for p in pts])
        y = np.array([p[1] for p in pts])
        cx, cy, keep = fit_curve(x, y)
        curve[sheet] = (cx, cy)
        n_train[sheet], n_kept[sheet] = len(pts), int(keep.sum())
    log("rule 4: vetted RT -> column_20 RT, one LOWESS curve per sheet (frac 0.5, 3 rounds)")
    expect('"492" training pairs', 72, n_train["492"])
    expect('"492" kept after the 3rd fit', 60, n_kept["492"])
    expect('"182" training pairs', 38, n_train["182"])
    expect('"182" kept after the 3rd fit', 35, n_kept["182"])

    def mapped(sheet, rt):
        return float(np.interp(rt, curve[sheet][0], curve[sheet][1]))

    # ---- rule 5: matched row, label, class
    for g in groups:
        for r in g.rows:
            r["mapped"] = mapped(r["sheet"], r["rt"])
        h20 = g.hits["column_20"]
        if h20["idx"] is not None:
            dist = np.array([abs(r["mapped"] - h20["rt"]) for r in g.rows])
            g.k = int(np.argmin(dist))                                # first in CSV order on a tie
            g.delta = h20["rt"] - g.rows[g.k]["mapped"]
            EDGE["identity deltas within 1e-6 of 0.30 min"] += int(abs(abs(g.delta) - MAX_IDENTITY_DELTA) < 1e-6)
        else:
            g.k, g.delta = 0, np.nan
        m = g.rows[g.k]
        g.label = m["name"].split("|")[0].strip()
        g.cls = m["cls"]
        g.sheet = m["sheet"]
        g.mapped = m["mapped"]
        g.same_ion = ""
        if g.conf == "low":
            hit = SAME_ION_RE.search(m["reason"])
            g.same_ion = hit.group(1) if hit else ""

    # ---- rule 6: panel membership, first failing rule
    for g in groups:
        ref, c20 = g.hits[REF], g.hits["column_20"]
        n_usable = sum(g.hits[c]["usable"] for c in COLS)
        g.n_usable = n_usable
        if ref["n"] == 0:
            g.reason = REASON["a"]
        elif ref["ndet"] < MIN_DET:
            g.reason = REASON["b"]
        elif ref["dom"] < MIN_DOM:
            g.reason = REASON["c"]
        elif n_usable < 3:
            g.reason = REASON["d"]
        elif c20["n"] == 0:
            g.reason = REASON["e"]
        elif c20["ndet"] < MIN_DET:
            g.reason = REASON["f"]
        elif c20["dom"] < MIN_DOM:
            g.reason = REASON["g"]
        elif abs(g.delta) > MAX_IDENTITY_DELTA:
            g.reason = REASON["h"]
        else:
            g.reason = None
    exp_reasons = {
        "high": {REASON["a"]: 58, REASON["c"]: 53, REASON["h"]: 29, REASON["g"]: 17, REASON["d"]: 10, REASON["e"]: 3},
        "low": {REASON["h"]: 13, REASON["c"]: 11, REASON["a"]: 6, REASON["d"]: 2, REASON["e"]: 1}}
    panel = {k: [g for g in groups if g.conf == k and g.reason is None] for k in COMPARE}
    log("rule 6: panel membership (first failing rule)")
    expect("high groups in panel", 182, len(panel["high"]))
    expect("low groups in panel", 3, len(panel["low"]))
    expect("low panel labels (sorted)", "DMPE 34:3, LPA 18:1, LPA 20:3",
           ", ".join(sorted(g.label for g in panel["low"])))
    for conf in COMPARE:
        counts = collections.Counter(g.reason for g in groups if g.conf == conf and g.reason is not None)
        for reason in list(exp_reasons[conf]) + sorted(set(counts) - set(exp_reasons[conf])):
            expect(f"{conf} excluded: {reason}", exp_reasons[conf].get(reason, 0), counts.get(reason, 0))
        expect(f"{conf} excluded, total", {"high": 170, "low": 33}[conf], sum(counts.values()))

    # ---- rule 8: display groups
    high_cls = collections.Counter(g.cls for g in panel["high"])
    for g in groups:
        if g.conf == "low":
            g.display = "Low confidence"
        else:
            g.display = g.cls if high_cls[g.cls] >= MIN_CLASS_LIPIDS else "Other"
    for g in panel["high"] + panel["low"]:
        assert g.display in GROUP_ORDER, f"display group {g.display!r} (class {g.cls}) is not in the group order"
    size = collections.Counter(g.display for g in panel["high"])
    log("rule 8: display groups of the high panel (lipids)")
    exp_sizes = {"TG": 34, "PC": 24, "SM": 21, "EtherPC": 21, "Other": 20, "LPC": 17, "CE": 11, "EtherPE": 9,
                 "PI": 8, "DG": 7, "Cer_NS": 5, "PE": 5}
    for name in list(exp_sizes) + sorted(set(size) - set(exp_sizes)):
        expect(f"group {name}", exp_sizes.get(name, 0), size.get(name, 0))
    other_classes = collections.Counter(g.cls for g in panel["high"] if g.display == "Other")
    log("  classes folded into Other:", ", ".join(f"{k} {v}" for k, v in sorted(other_classes.items())))

    # ---- rule 7: cells of the used non-reference columns
    def cells_of(g):
        out = []
        ref = g.hits[REF]
        for c in NONREF:
            h = g.hits[c]
            if not h["usable"]:
                continue
            f = feat[c]
            cal = float(f["cal_rt"][h["idx"]])
            out.append(dict(col=c, aid=h["aid"], ref_aid=ref["aid"], rt=h["rt"], ref_rt=ref["rt"], cal_rt=cal,
                            before=h["rt"] - ref["rt"], after=cal - ref["rt"],
                            series=float(f["series"][h["idx"]]), lattice=float(f["lattice"][h["idx"]])))
        return out

    for g in groups:
        g.cells = cells_of(g)
    all_cells = {k: [(g, cl) for g in panel[k] for cl in g.cells] for k in COMPARE}

    def stats(cells, count_edges=False):
        b = np.abs([cl["before"] for cl in cells])
        a = np.abs([cl["after"] for cl in cells])
        rb, ra = np.round(b, 9), np.round(a, 9)                        # threshold tests: no float fuzz on 3-decimal RTs
        if count_edges:
            EDGE["cells with |dRT_before| within 1e-6 of 0.2 min"] += int((np.abs(rb - WITHIN) < 1e-6).sum())
            EDGE["cells with |dRT_after| within 1e-6 of 0.2 min"] += int((np.abs(ra - WITHIN) < 1e-6).sum())
        return dict(n=len(cells), med_after=float(np.median(a)), med_before=float(np.median(b)),
                    in_after=int((ra <= WITHIN).sum()), in_before=int((rb <= WITHIN).sum()),
                    gt05=int((ra > 0.5).sum()), gt1=int((ra > 1.0).sum()), gt2=int((ra > 2.0).sum()))

    sh = stats([cl for _, cl in all_cells["high"]], count_edges=True)
    sl = stats([cl for _, cl in all_cells["low"]])
    for g, cl in all_cells["high"]:
        for stage in ("before", "after"):
            v = abs(cl[stage])
            if abs(round(v, 9) - WITHIN) < 1e-6:
                EDGE_NOTES.append(
                    f"{g.label} on {cl['col']}: |dRT_{stage}| = {v!r} min against the 0.2 min threshold (<=): counted as "
                    f"{'within' if round(v, 9) <= WITHIN else 'outside'} (a plain float comparison decides "
                    f"{'the same' if (v <= WITHIN) == (round(v, 9) <= WITHIN) else 'DIFFERENTLY'})")
    log("rule 7: used non-reference cells of the high panel")
    expect("used non-reference cells", 830, sh["n"])
    per_col = collections.Counter(cl["col"] for _, cl in all_cells["high"])
    for c, e in zip(NONREF, (163, 178, 182, 160, 147)):
        expect(f"  cells on {c}", e, per_col[c])
    expect("median |dRT_after| (min)", 0.0486, sh["med_after"], nd=4)
    expect("within 0.2 min after calibration, cells", 702, sh["in_after"])
    expect("within 0.2 min after calibration, %", 84.6, 100 * sh["in_after"] / sh["n"], nd=1)
    expect("median |dRT_before| (min)", 6.053, sh["med_before"], nd=3)
    expect("within 0.2 min before calibration, cells", 44, sh["in_before"])
    expect("within 0.2 min before calibration, %", 5.3, 100 * sh["in_before"] / sh["n"], nd=1)
    expect("cells with |dRT_after| > 0.5 min", 28, sh["gt05"])
    expect("cells with |dRT_after| > 1 min", 15, sh["gt1"])
    expect("cells with |dRT_after| > 2 min", 7, sh["gt2"])
    expect("low panel: used non-reference cells", 9, sl["n"])

    # ---- sort the panels: group order, then column_22 RT
    def sort_key(g):
        return (GROUP_ORDER.index(g.display), g.hits[REF]["rt"])

    high_sorted = sorted(panel["high"], key=sort_key)
    low_sorted = sorted(panel["low"], key=sort_key)
    sets = {"high_confidence": high_sorted, "high_and_low_confidence": high_sorted + low_sorted}

    # ---- rule 9: pooled cells of a set (each cell once; high lipids first)
    log("rule 9: pooled cells per set (key = column, pick Alignment ID, reference pick Alignment ID)")
    pooled = {}
    for name, gl in sets.items():
        seen, cells = set(), []
        for g in gl:
            for cl in g.cells:
                key = (cl["col"], cl["aid"], cl["ref_aid"])
                if key not in seen:
                    seen.add(key)
                    cells.append(cl)
        pooled[name] = cells
        st = stats(cells)
        expect(f"{name}: pooled cells", 830, st["n"])
        log(f"    {name}: median |dRT_after| {st['med_after']:.4f} min, {100 * st['in_after'] / st['n']:.1f} % within "
            f"0.2 min ({st['in_after']} of {st['n']}); median |dRT_before| {st['med_before']:.3f} min, "
            f"{100 * st['in_before'] / st['n']:.1f} % ({st['in_before']} of {st['n']})")
    dupes = len(all_cells["high"]) - len(pooled["high_confidence"])
    log(f"    duplicate cells inside the high panel: {dupes}; the {sl['n']} low-panel cells are the same peaks as: "
        + ", ".join(sorted({f"{g.label} (= {g.same_ion})" for g in panel['low']})))

    # ---- output: panel CSVs
    head = ["lipid", "class", "display_group", "confidence", "same_ion_as", "vetted_names", "n_vetted_rows",
            "matched_vetted_name", "formula", "positive_target_ions", "vetted_RT_min",
            "vetted_RT_mapped_to_column_20_min", "column_20_identity_delta_min", "n_columns_used",
            "n_columns_series_corrected", "n_columns_lattice_corrected"]
    for c in COLS:
        head += [f"{c}_{s}" for s in ("used", "alignment_id", "mz", "RT_min", "n_candidates", "dominance",
                                      "n_injections_detected", "mean_area")]
    for c in NONREF:
        head += [f"{c}_{s}" for s in ("cal_RT_min", "dRT_before_min", "dRT_after_min",
                                      "series_correction_min", "lattice_correction_min")]

    def panel_record(g):
        m = g.rows[g.k]
        rec = {"lipid": g.label, "class": g.cls, "display_group": g.display, "confidence": g.conf,
               "same_ion_as": g.same_ion, "vetted_names": "; ".join(r["name"] for r in g.rows),
               "n_vetted_rows": len(g.rows), "matched_vetted_name": m["name"], "formula": m["formula"],
               "positive_target_ions": m["ions_str"], "vetted_RT_min": m["rt"],
               "vetted_RT_mapped_to_column_20_min": g.mapped, "column_20_identity_delta_min": g.delta,
               "n_columns_used": sum(g.hits[c]["usable"] for c in COLS),
               "n_columns_series_corrected": sum(cl["series"] != 0 for cl in g.cells),
               "n_columns_lattice_corrected": sum(cl["lattice"] != 0 for cl in g.cells)}
        for c in COLS:
            h = g.hits[c]
            rec.update({f"{c}_used": h["usable"], f"{c}_alignment_id": h["aid"], f"{c}_mz": h["mz"],
                        f"{c}_RT_min": h["rt"], f"{c}_n_candidates": h["n"], f"{c}_dominance": h["dom"],
                        f"{c}_n_injections_detected": h["ndet"], f"{c}_mean_area": h["area"]})
        for cl in g.cells:
            c = cl["col"]
            rec.update({f"{c}_cal_RT_min": cl["cal_rt"], f"{c}_dRT_before_min": cl["before"],
                        f"{c}_dRT_after_min": cl["after"], f"{c}_series_correction_min": cl["series"],
                        f"{c}_lattice_correction_min": cl["lattice"]})
        return rec

    written = []
    for name, gl in sets.items():
        path = os.path.join(HERE, f"{name}_panel.csv")
        write_csv(path, head, [panel_record(g) for g in gl])
        written.append(os.path.basename(path))
        labels = [g.label for g in gl]
        dup = [k for k, v in collections.Counter(labels).items() if v > 1]
        log(f"wrote {os.path.basename(path)}: {len(gl)} lipids, {len(head)} columns"
            + (f"; label used more than once: {dup}" if dup else "; labels unique"))
    corrected = [(cl["series"], cl["lattice"]) for g in high_sorted for cl in g.cells]
    tiny = sum((0 < abs(s) < 5e-5) + (0 < abs(l) < 5e-5) for s, l in corrected)
    EDGE["used high-panel corrections that are non-zero but print as 0.0000"] += tiny

    # ---- output: excluded lipids (every group not in a panel + the FA rows)
    ehead = ["confidence", "lipid", "class", "vetted_names", "positive_target_ions", "reason",
             "n_columns_with_candidates", "n_columns_usable"]
    for c in COLS:
        ehead += [f"{c}_{s}" for s in ("n_candidates", "dominance", "RT_min")]
    excluded = []
    for g in groups:
        if g.reason is None:
            continue
        m = g.rows[g.k]
        rec = {"confidence": g.conf, "lipid": g.label, "class": g.cls,
               "vetted_names": "; ".join(r["name"] for r in g.rows), "positive_target_ions": m["ions_str"],
               "reason": g.reason, "n_columns_with_candidates": sum(g.hits[c]["n"] > 0 for c in COLS),
               "n_columns_usable": g.n_usable}
        for c in COLS:
            h = g.hits[c]
            rec.update({f"{c}_n_candidates": h["n"], f"{c}_dominance": h["dom"], f"{c}_RT_min": h["rt"]})
        excluded.append(rec)
    n_group_excl = len(excluded)
    for r in no_target:
        excluded.append({"confidence": r["conf"], "lipid": r["name"].split("|")[0].strip(), "class": r["cls"],
                         "vetted_names": r["name"], "positive_target_ions": "", "reason": REASON["fa"]})
    path = os.path.join(HERE, "excluded_lipids.csv")
    write_csv(path, ehead, excluded)
    written.append("excluded_lipids.csv")
    log(f"wrote excluded_lipids.csv: {n_group_excl} excluded groups + {len(no_target)} FA rows = {len(excluded)} entries")
    expect("excluded entries: high groups", 170, sum(1 for g in groups if g.conf == "high" and g.reason))
    expect("excluded entries: low groups", 33, sum(1 for g in groups if g.conf == "low" and g.reason))
    expect("excluded entries: FA rows", 14, len(no_target))

    # ---- output: the mapping curves and the per-group points of the QC figure
    recs = [{"sheet": s, "vetted_RT_min": x, "column_20_RT_min": y}
            for s in ("492", "182") for x, y in zip(*curve[s])]
    write_csv(os.path.join(HERE, "vetted_rt_mapping.csv"), ["sheet", "vetted_RT_min", "column_20_RT_min"], recs)
    written.append("vetted_rt_mapping.csv")
    log(f"wrote vetted_rt_mapping.csv: curve points 492 = {len(curve['492'][0])}, 182 = {len(curve['182'][0])}")
    recs = []
    for g in groups:
        if g.hits["column_20"]["idx"] is None:
            continue
        status = "in_panel" if g.reason is None else ("excluded_rule_h" if g.reason == REASON["h"] else "excluded_other")
        recs.append({"confidence": g.conf, "sheet": g.sheet, "lipid": g.label, "status": status, "reason": g.reason or "",
                     "vetted_RT_min": g.rows[g.k]["rt"], "vetted_RT_mapped_to_column_20_min": g.mapped,
                     "column_20_RT_min": g.hits["column_20"]["rt"], "column_20_identity_delta_min": g.delta})
    write_csv(os.path.join(HERE, "vetted_rt_check_groups.csv"),
              ["confidence", "sheet", "lipid", "status", "reason", "vetted_RT_min", "vetted_RT_mapped_to_column_20_min",
               "column_20_RT_min", "column_20_identity_delta_min"], recs)
    written.append("vetted_rt_check_groups.csv")
    log(f"wrote vetted_rt_check_groups.csv: {len(recs)} groups with a column_20 pick "
        f"(in panel {sum(r['status'] == 'in_panel' for r in recs)}, rule h {sum(r['status'] == 'excluded_rule_h' for r in recs)}, "
        f"other {sum(r['status'] == 'excluded_other' for r in recs)})")

    # ---- self-check: every list row accounted for exactly once (by row identity, then via the written CSVs)
    placed = collections.Counter()
    for g in sets["high_and_low_confidence"]:                       # the panel lipids (high, then low)
        placed.update((r["conf"], r["i"]) for r in g.rows)
    for g in groups:                                                # the excluded groups
        if g.reason is not None:
            placed.update((r["conf"], r["i"]) for r in g.rows)
    placed.update((r["conf"], r["i"]) for r in no_target)           # the FA rows
    assert set(placed) == {(r["conf"], r["i"]) for r in rows} and set(placed.values()) == {1}         and sum(placed.values()) == len(rows), "a list row is missing or counted twice"
    all_names = collections.Counter(r["name"] for r in rows)
    got = collections.Counter()
    hp = pd.read_csv(os.path.join(HERE, "high_confidence_panel.csv"), keep_default_na=False, dtype=str)
    hl = pd.read_csv(os.path.join(HERE, "high_and_low_confidence_panel.csv"), keep_default_na=False, dtype=str)
    ex = pd.read_csv(os.path.join(HERE, "excluded_lipids.csv"), keep_default_na=False, dtype=str)
    for names in list(hp.vetted_names) + list(hl[hl.confidence == "low"].vetted_names) + list(ex.vetted_names):
        got.update(names.split("; "))
    assert got == all_names, "vetted_names across the CSVs do not reproduce the 674 list rows"
    assert (hl[hl.confidence == "high"].reset_index(drop=True).equals(hp)), \
        "high rows of high_and_low_confidence_panel.csv differ from high_confidence_panel.csv"
    log(f"self-check: all {len(rows)} list rows appear exactly once across high panel ({len(hp)} lipids), "
        f"low panel ({int((hl.confidence == 'low').sum())} lipids) and excluded table ({len(ex)} entries); "
        f"checked by row index and again by the vetted_names in the written CSVs")

    # ---- float edge cases met on the way (they would be where two implementations could differ)
    log("float edge cases met (counts):")
    for k in sorted(EDGE):
        log(f"  {EDGE[k]:>4}  {k}")
    for note in EDGE_NOTES:
        log("        - " + note)

    # ---- verdict
    log("")
    log(f"expected-number check: {len(DIFFS)} expected number(s) differ" + (": " + "; ".join(DIFFS) if DIFFS else ""))
    log(f"done in {time.time() - t0:.1f} s; wrote: " + ", ".join(written) + ", build_vetted_panel_console.log")
    write_text(os.path.join(HERE, "build_vetted_panel_console.log"), "\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
