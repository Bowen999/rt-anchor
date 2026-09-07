"""The endogenous plasma-lipid anchor panel, and how a feature earns anchor status.

Stage 2 of the cross-column method needs anchors that are actually *there*. The
v2 hit-rate table settled the question of what to use: most mixture standards
are absent from serum, while these 17 high-abundance human plasma lipids are
present in every serum run on every column. They are the anchors, and they are
also the pool the manual-check list is drawn from.

Exact masses are computed from sum composition with monoisotopic atomic masses.
The PC / DG / CE mass rules were cross-checked against the lab's own standards
list (20260619_LipidRTCal_MZandRT.xlsx) and agree to under 1 mDa.

Identification is **sum-composition only** — MS1, no MS/MS. No claim is made
about sn-position or double-bond location, and none should be read into the
labels.

Validation, per lipid per run:

1. enumerate every feature within the m/z window of the expected precursor;
2. require the adduct annotation to match the expected adduct — *where the table
   has adducts at all*; where it does not, this filter is skipped, not faked;
3. prefer fully-filled features (``Fill % >= 1``, likewise skipped where absent)
   and among those take the highest S/N;
4. check the within-class RT-order rules below. A monotone column-to-column
   mapping must preserve elution order within a class, so a violation is
   evidence of an isomer mis-pick, not of interesting chemistry;
5. record the remaining significant in-window features as isomer candidates.

**Curve-assisted isomer refinement.** Highest-S/N selection can grab the wrong
isomer — observed for LPC 18:1 on column 15, where a more intense isomer elutes
0.4 min early. Once the anchor-free stage-1 curve exists, every in-window,
adduct-consistent candidate is mapped through it and the candidate whose
calibrated RT lands closest to the reference run's RT wins, with near-ties
broken by Fill% then S/N. The curve is built from anonymous m/z-matched
features only, so no lipid identity leaks into the refinement that is about to
produce the anchors.

**Matrix caveat.** This panel is human plasma/serum. In any other matrix too few
anchors will validate and the engine will fall back to the stage-1 curve. That
is the designed behaviour, not a failure — but it does mean two users with the
same instrument and different matrices get corrections of different quality with
no error raised, so callers should log how many anchors validated.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .crosscolumn import mz_window
from .io.schema import FeatureTable

# monoisotopic atomic masses
_C, _H, _N, _O, _P = 12.0, 1.007825, 14.003074, 15.994915, 30.973762
PROTON = 1.007276   # mass of a proton, for [M+H]+
NH4 = 18.033823     # NH4 adduct mass, for [M+NH4]+

DEFAULT_ISOMER_SN_FRAC = 0.10
DEFAULT_TIE_TOL_MIN = 0.15


def _mass(nC, nH, nN, nO, nP):
    return nC * _C + nH * _H + nN * _N + nO * _O + nP * _P


def plasma_lipid_candidates() -> pd.DataFrame:
    """The 17-lipid prior-knowledge panel.

    Columns: lipid, lipid_class, sum_comp, formula, exact_mass, adduct,
    precursor_mz.
    """
    rows = []

    def add(name, cls, comp, nC, nH, nN, nO, nP, adduct):
        M = _mass(nC, nH, nN, nO, nP)
        mz = M + (PROTON if adduct == "[M+H]+" else NH4)
        rows.append(dict(lipid=name, lipid_class=cls, sum_comp=comp,
                         formula=f"C{nC}H{nH}" + (f"N{nN}" if nN else "")
                                 + f"O{nO}" + ("P" if nP else ""),
                         exact_mass=round(M, 6), adduct=adduct,
                         precursor_mz=round(mz, 6)))

    for a, b in [(16, 0), (18, 2), (18, 1), (18, 0)]:                     # LPC
        nC = a + 8
        add(f"LPC {a}:{b}", "LPC", f"{a}:{b}", nC, 2 * nC - 2 * b + 2, 1, 7, 1, "[M+H]+")
    for a, b in [(34, 2), (34, 1), (36, 4), (36, 2), (36, 1), (38, 4)]:   # PC
        nC = a + 8
        add(f"PC {a}:{b}", "PC", f"{a}:{b}", nC, 2 * nC - 2 * b, 1, 8, 1, "[M+H]+")
    for a, b in [(34, 1), (36, 1)]:                                       # SM d18:1/*
        nC = a + 5
        add(f"SM {a}:{b};O2", "SM", f"{a}:{b};O2", nC, 2 * nC + 3 - 2 * b, 2, 6, 1, "[M+H]+")
    for a, b in [(18, 2), (18, 1)]:                                       # CE
        add(f"CE {a}:{b}", "CE", f"{a}:{b}", 27 + a, 44 + 2 * a - 2 * b, 0, 2, 0, "[M+NH4]+")
    for a, b in [(52, 3), (52, 2), (54, 3)]:                              # TG
        add(f"TG {a}:{b}", "TG", f"{a}:{b}", a + 3, 2 * a - 2 * b + 2, 0, 6, 0, "[M+NH4]+")

    return pd.DataFrame(rows)


# Expected within-class elution order, earlier -> later. Checked against the
# observed RTs of every run: a monotone column-to-column mapping must preserve
# these orders, so a violation flags a likely isomer mis-pick.
RT_ORDER_RULES: List[Tuple[str, List[str]]] = [
    ("LPC", ["LPC 18:2", "LPC 18:1", "LPC 18:0"]),   # C18, unsaturation decreases
    ("LPC", ["LPC 16:0", "LPC 18:1"]),               # 16:0 elutes before 18:1 here
    ("LPC", ["LPC 16:0", "LPC 18:0"]),               # saturated, chain lengthens
    ("PC", ["PC 36:4", "PC 36:2", "PC 36:1"]),       # C36
    ("PC", ["PC 34:2", "PC 34:1"]),                  # C34
    ("PC", ["PC 34:2", "PC 36:2"]),                  # b=2, chain lengthens
    ("PC", ["PC 34:1", "PC 36:1"]),                  # b=1, chain lengthens
    ("SM", ["SM 34:1;O2", "SM 36:1;O2"]),
    ("CE", ["CE 18:2", "CE 18:1"]),
    ("TG", ["TG 52:3", "TG 52:2"]),
]


def feature_frame(ft: FeatureTable) -> pd.DataFrame:
    """Flatten a feature table into the columns this module needs.

    ``adduct`` and ``fill`` are ``None``-valued columns when the source table has
    no such concept; every consumer here checks for that and skips the
    corresponding filter rather than substituting a default.
    """
    n = len(ft.df)
    adducts = ft.adducts()
    fills = ft.fill()
    return pd.DataFrame({
        "row": np.arange(n),
        "mz": ft.mz().to_numpy(dtype=float),
        "rt": ft.rt_minutes().to_numpy(dtype=float),
        "sn": ft.sn().to_numpy(dtype=float),
        "intensity": ft.intensity().to_numpy(dtype=float),
        "fid": ft.feature_id().to_numpy(),
        "adduct": (None if adducts is None else adducts.to_numpy()),
        "fill": (np.full(n, np.nan) if fills is None else fills.to_numpy(dtype=float)),
    }).assign(has_adduct=adducts is not None, has_fill=fills is not None)


def _candidate_pool(frame: pd.DataFrame, precursor_mz: float, adduct: str,
                    tol: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """In-window features, and the adduct-consistent sub-pool preferred among them.

    When no in-window feature carries the expected adduct the whole window is
    kept: the annotation may simply be missing, and an empty pool would lose a
    real anchor. When the table has no adduct column at all the filter never
    runs.
    """
    in_win = frame.loc[(frame["mz"] - precursor_mz).abs() < tol]
    if in_win.empty:
        return in_win, in_win
    if not bool(in_win["has_adduct"].iloc[0]):
        return in_win, in_win
    right = in_win[in_win["adduct"].astype(str) == adduct]
    return in_win, (right if not right.empty else in_win)


def _best_of(pool: pd.DataFrame) -> pd.Series:
    """The pool's pick: highest S/N, or the first in-window row when S/N is absent.

    A bare ``m/z`` + ``RT`` table (a minimal MassCube export, and the shipped
    example inputs) carries neither an S/N column nor any abundance column, so
    ``sn`` is all-NaN and ``idxmax`` has nothing to rank — it raises
    "Encountered all NA values". There is genuinely no evidence to prefer one
    in-window feature over another there, so the pick falls back to the lowest
    row index: deterministic, and honest about the absence of a ranking.
    """
    sn = pool["sn"]
    return pool.loc[sn.idxmax()] if sn.notna().any() else pool.iloc[0]


def _prefer_filled(pool: pd.DataFrame) -> pd.DataFrame:
    """Restrict to fully-filled features when the table reports Fill% at all."""
    if pool.empty or not bool(pool["has_fill"].iloc[0]):
        return pool
    full = pool[pool["fill"] >= 1.0]
    return full if not full.empty else pool


def validate_candidates(tables: Dict[str, FeatureTable],
                        candidates: Optional[pd.DataFrame] = None,
                        config=None,
                        *,
                        ref_key: Optional[str] = None,
                        isomer_sn_frac: Optional[float] = None,
                        mz_tol_ppm: Optional[float] = None,
                        mz_tol_min_da: Optional[float] = None,
                        ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Locate every candidate lipid in every supplied sample run.

    ``tables`` maps a run key (any label — the pipeline uses ``"source"`` and
    ``"reference"``) to its feature table.

    Returns ``(picks, summary)``: one row per (lipid, run) and one row per lipid.
    """
    if candidates is None:
        candidates = plasma_lipid_candidates()
    keys = list(tables)
    ref_key = ref_key if ref_key is not None else (keys[-1] if keys else None)
    frac = float(isomer_sn_frac if isomer_sn_frac is not None
                 else getattr(config, "isomer_sn_frac", None) or DEFAULT_ISOMER_SN_FRAC)
    frames = {k: feature_frame(ft) for k, ft in tables.items()}

    pick_rows = []
    for _, lip in candidates.iterrows():
        tol = float(mz_window(lip["precursor_mz"], config,
                              mz_tol_ppm=mz_tol_ppm, mz_tol_min_da=mz_tol_min_da))
        for key in keys:
            in_win, pool = _candidate_pool(frames[key], lip["precursor_mz"],
                                           lip["adduct"], tol)
            if in_win.empty:
                continue
            pool = _prefer_filled(pool)
            best = _best_of(pool)
            notable = _isomer_census(in_win, best, frac)
            pick_rows.append(dict(
                lipid=lip["lipid"], lipid_class=lip["lipid_class"], column=key,
                rt=float(best["rt"]), mz_obs=float(best["mz"]),
                sn=float(best["sn"]), intensity=float(best["intensity"]),
                row=int(best["row"]), align_id=best["fid"],
                adduct_ok=bool(
                    (not bool(in_win["has_adduct"].iloc[0]))
                    or (in_win["adduct"].astype(str) == lip["adduct"]).any()),
                n_isomer_candidates=int(len(notable)),
                isomer_rts=";".join(f"{r:.2f}" for r in sorted(notable["rt"])),
                pick_refined=False,
            ))
    picks = pd.DataFrame(pick_rows)
    return picks, summarize_picks(picks, candidates, keys, ref_key=ref_key)


def summarize_picks(picks: pd.DataFrame,
                    candidates: pd.DataFrame,
                    all_keys: Sequence[str],
                    ref_key: Optional[str] = None) -> pd.DataFrame:
    """Per-lipid validation audit: where it was found and whether the order holds.

    A lipid is ``validated`` only when it was found in **every** supplied run and
    breaks none of its RT-order rules anywhere. That is deliberately strict: an
    anchor that is missing on one column would silently change the anchor set
    between columns and make their results incomparable.
    """
    all_keys = list(all_keys)
    ref_key = ref_key if ref_key is not None else (all_keys[-1] if all_keys else None)
    if picks is None or picks.empty:
        picks = pd.DataFrame(columns=["lipid", "column", "rt", "sn",
                                      "n_isomer_candidates"])
    summary_rows = []
    for _, lip in candidates.iterrows():
        sub = picks[picks["lipid"] == lip["lipid"]]
        detected = sorted(sub["column"].tolist())
        order_ok = {}
        for key in detected:
            ok = True
            for _cls, chain in RT_ORDER_RULES:
                if lip["lipid"] not in chain:
                    continue
                rts = {}
                for member in chain:
                    r = picks[(picks["lipid"] == member) & (picks["column"] == key)]
                    if not r.empty:
                        rts[member] = r["rt"].iloc[0]
                seq = [rts[m] for m in chain if m in rts]
                if len(seq) >= 2 and any(a >= b for a, b in zip(seq, seq[1:])):
                    ok = False
            order_ok[key] = ok
        ref = sub[sub["column"] == ref_key]
        summary_rows.append(dict(
            lipid=lip["lipid"], lipid_class=lip["lipid_class"],
            sum_comp=lip["sum_comp"], formula=lip["formula"],
            exact_mass=lip["exact_mass"], adduct=lip["adduct"],
            precursor_mz=lip["precursor_mz"],
            detected_in=",".join(detected), n_columns=len(detected),
            rt_ref_min=float(ref["rt"].iloc[0]) if not ref.empty else np.nan,
            sn_ref=float(ref["sn"].iloc[0]) if not ref.empty else np.nan,
            rt_order_ok=all(order_ok.values()) if order_ok else False,
            max_isomer_candidates=int(sub["n_isomer_candidates"].max()) if not sub.empty else 0,
            validated=(len(detected) == len(all_keys)) and bool(order_ok) and all(order_ok.values()),
        ))
    return pd.DataFrame(summary_rows)


def refine_picks_with_curves(picks: pd.DataFrame,
                             tables: Dict[str, FeatureTable],
                             curves: Dict[str, object],
                             candidates: Optional[pd.DataFrame] = None,
                             config=None,
                             *,
                             ref_key: str = "reference",
                             tie_tol_min: Optional[float] = None,
                             isomer_sn_frac: Optional[float] = None,
                             mz_tol_ppm: Optional[float] = None,
                             mz_tol_min_da: Optional[float] = None) -> pd.DataFrame:
    """Re-pick in-window isomer candidates using the anchor-free curve.

    For each lipid that has a reference-run pick, every in-window,
    adduct-consistent candidate in a source run is mapped through that run's
    stage-1 curve, and the candidate whose calibrated RT lands closest to the
    reference RT wins. Candidates within ``tie_tol_min`` of the best distance
    count as tied; among ties, prefer a filled feature, then the higher S/N.
    Reference-run picks are never changed — they define the target.

    Sets ``pick_refined`` where the pick actually moved.
    """
    if candidates is None:
        candidates = plasma_lipid_candidates()
    tie = float(tie_tol_min if tie_tol_min is not None
                else getattr(config, "anchor_tie_tol_min", None) or DEFAULT_TIE_TOL_MIN)
    frac = float(isomer_sn_frac if isomer_sn_frac is not None
                 else getattr(config, "isomer_sn_frac", None) or DEFAULT_ISOMER_SN_FRAC)
    frames = {k: feature_frame(ft) for k, ft in tables.items()}
    cand = candidates.set_index("lipid")
    picks = picks.copy()
    if "pick_refined" not in picks.columns:
        picks["pick_refined"] = False

    for lip_name, lip in cand.iterrows():
        ref = picks[(picks["lipid"] == lip_name) & (picks["column"] == ref_key)]
        if ref.empty:
            continue
        rt_ref = float(ref["rt"].iloc[0])
        tol = float(mz_window(lip["precursor_mz"], config,
                              mz_tol_ppm=mz_tol_ppm, mz_tol_min_da=mz_tol_min_da))
        for key, curve in curves.items():
            if key == ref_key or key not in frames:
                continue
            hit = picks[(picks["lipid"] == lip_name) & (picks["column"] == key)].index
            if len(hit) == 0:
                continue
            idx = hit[0]
            in_win, pool = _candidate_pool(frames[key], lip["precursor_mz"],
                                           lip["adduct"], tol)
            if pool.empty:
                continue
            cal = curve.predict(pool["rt"].to_numpy())
            d = np.abs(cal - rt_ref)
            tied = pool[d <= d.min() + tie].copy()
            tied["_full"] = (tied["fill"] >= 1.0) if bool(tied["has_fill"].iloc[0]) else False
            best = tied.sort_values(["_full", "sn"], ascending=False).iloc[0]
            changed = best["fid"] != picks.loc[idx, "align_id"]
            notable = _isomer_census(in_win, best, frac)
            picks.loc[idx, "rt"] = float(best["rt"])
            picks.loc[idx, "mz_obs"] = float(best["mz"])
            picks.loc[idx, "sn"] = float(best["sn"])
            picks.loc[idx, "intensity"] = float(best["intensity"])
            picks.loc[idx, "row"] = int(best["row"])
            picks.loc[idx, "align_id"] = best["fid"]
            picks.loc[idx, "n_isomer_candidates"] = int(len(notable))
            picks.loc[idx, "isomer_rts"] = ";".join(f"{r:.2f}" for r in sorted(notable["rt"]))
            picks.loc[idx, "pick_refined"] = bool(changed)
    return picks


def build_anchor_table(picks: pd.DataFrame,
                       source_key: str,
                       validated_lipids: Sequence[str],
                       candidates: Optional[pd.DataFrame] = None,
                       ref_key: str = "reference") -> pd.DataFrame:
    """The stage-2 anchor frame for one source run.

    Columns ``label``, ``rt_src``, ``rt_ref``, ``lipid_class``,
    ``n_isomer_candidates``, ``isomer_rts``, ``pick_refined`` — exactly what
    :func:`rt_anchor.crosscolumn.build_calibrator` consumes and what
    ``<prefix>_anchors.csv`` reports.
    """
    if candidates is None:
        candidates = plasma_lipid_candidates()
    if picks is None or picks.empty:
        return pd.DataFrame(columns=["label", "rt_src", "rt_ref", "lipid_class"])
    keep = list(validated_lipids)
    sub = picks[picks["lipid"].isin(keep)]
    extra = [c for c in ("n_isomer_candidates", "isomer_rts", "pick_refined")
             if c in sub.columns]
    src = (sub[sub["column"] == source_key][["lipid", "rt"] + extra]
           .rename(columns={"rt": "rt_src"}))
    ref = (sub[sub["column"] == ref_key][["lipid", "rt"]]
           .rename(columns={"rt": "rt_ref"}))
    a = src.merge(ref, on="lipid").rename(columns={"lipid": "label"})
    a["lipid_class"] = a["label"].map(candidates.set_index("lipid")["lipid_class"])
    return a.sort_values("rt_src").reset_index(drop=True)


# ---------------------------------------------------------------- internals ---

def _isomer_census(in_win: pd.DataFrame, best: pd.Series, sn_frac: float) -> pd.DataFrame:
    """Other in-window features whose S/N is a notable fraction of the pick's."""
    others = in_win[in_win["fid"] != best["fid"]]
    return others[others["sn"] > sn_frac * best["sn"]]
