"""Identify the panel standards (anchors) inside a feature table.

Two steps mirror the spec:

1. :func:`build_native_template` — from the standard-panel run, find where each
   anchor actually elutes *on this method* (the "native" RT). This centres the
   later search windows, which is essential on QTOF where RT is shifted several
   minutes from the Orbitrap reference.
2. :func:`identify_anchors` — in the table being calibrated, find each anchor in
   a tight window around its native RT (or, if no template, a wide window around
   the reference RT), pick the best candidate, and enforce monotonic order.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .config import CalibrationConfig
from .errors import AnchorIdentificationError
from .io.schema import FeatureTable
from .panel import Panel


def _intensity(ft: FeatureTable) -> pd.Series:
    """Row intensity for candidate ranking.

    Preference: sum over declared sample columns; else a single abundance column
    if present (``peak_height``/``peak_area``/``maxo``/``height``/``area``);
    else NaN (ranking then falls back to RT-closeness).
    """
    if ft.sample_cols:
        sub = ft.df[ft.sample_cols].apply(pd.to_numeric, errors="coerce")
        return sub.sum(axis=1, skipna=True)
    cols = {c.lower(): c for c in ft.df.columns}
    for key in ("peak_height", "peak_area", "maxo", "height", "area", "intensity"):
        if key in cols:
            return pd.to_numeric(ft.df[cols[key]], errors="coerce")
    return pd.Series(np.nan, index=ft.df.index)


def _quality_mask(ft: FeatureTable, config: CalibrationConfig) -> pd.Series:
    """Boolean mask keeping rows that pass the optional peak-quality gate.

    Only applies when the columns exist (they are optional / used-if-present).
    """
    mask = pd.Series(True, index=ft.df.index)
    cols = {c.lower(): c for c in ft.df.columns}
    if "is_isotope" in cols:
        mask &= ~_as_bool(ft.df[cols["is_isotope"]])
    if "is_in_source_fragment" in cols:
        mask &= ~_as_bool(ft.df[cols["is_in_source_fragment"]])
    if "charge" in cols:
        ch = pd.to_numeric(ft.df[cols["charge"]], errors="coerce")
        mask &= (ch.isna() | (ch.abs() == 1))
    if config.min_gaussian_similarity > 0 and "gaussian_similarity" in cols:
        g = pd.to_numeric(ft.df[cols["gaussian_similarity"]], errors="coerce")
        mask &= (g.isna() | (g >= config.min_gaussian_similarity))
    lo, hi = config.asymmetry_range
    if (lo > 0 or hi < 1e9) and "asymmetry_factor" in cols:
        a = pd.to_numeric(ft.df[cols["asymmetry_factor"]], errors="coerce")
        mask &= (a.isna() | ((a >= lo) & (a <= hi)))
    return mask


def _as_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "1.0", "yes"))


def build_native_template(std_table: FeatureTable,
                          panel: Panel,
                          config: CalibrationConfig) -> Dict[str, float]:
    """Return {anchor_name: native_rt_min} from the standard-panel run.

    Uses a wide RT window around the reference RT (to tolerate cross-instrument
    shift) and, within it, the most intense (or RT-closest) candidate.
    """
    mz = std_table.mz().to_numpy()
    rt = std_table.rt_minutes().to_numpy()
    inten = _intensity(std_table).to_numpy()
    qmask = _quality_mask(std_table, config).to_numpy()

    template: Dict[str, float] = {}
    for _, t in panel.targets.iterrows():
        tol = max(t["mz"] * config.mz_tol_ppm / 1e6, config.mz_tol_min_da)
        mzmask = qmask & (np.abs(mz - t["mz"]) <= tol)
        if not mzmask.any():
            continue
        # A dedicated standards run has the standard as the dominant peak at its
        # m/z, so when intensity is available search the WHOLE gradient (this is
        # what lets us catch anchors shifted several minutes from the reference,
        # e.g. QTOF). Without intensity, fall back to a reference-seeded window.
        if np.any(np.isfinite(inten[mzmask])):
            sel = np.where(mzmask)[0]
        else:
            sel = np.where(mzmask & (np.abs(rt - t["rt_ref_min"]) <= config.rt_window_seed_min))[0]
        if sel.size == 0:
            continue
        template[t["name"]] = float(_pick(sel, rt, inten, center=t["rt_ref_min"]))
    return _enforce_monotone_template(template, panel)


def identify_anchors(table: FeatureTable,
                     panel: Panel,
                     config: CalibrationConfig,
                     native_template: Optional[Dict[str, float]] = None) -> pd.DataFrame:
    """Find each anchor in ``table``; return one row per anchor found.

    Columns: name, mz, rt_ref_min, irt, rt_obs_min, intensity, n_candidates,
    endogenous, void.
    """
    mz = table.mz().to_numpy()
    rt = table.rt_minutes().to_numpy()
    inten = _intensity(table).to_numpy()
    qmask = _quality_mask(table, config).to_numpy()

    rows: List[dict] = []
    for _, t in panel.targets.iterrows():
        if native_template and t["name"] in native_template:
            center = native_template[t["name"]]
            win = config.rt_window_min
        else:
            center = t["rt_ref_min"]
            win = config.rt_window_seed_min
        tol = max(t["mz"] * config.mz_tol_ppm / 1e6, config.mz_tol_min_da)
        sel = np.where(qmask & (np.abs(mz - t["mz"]) <= tol) & (np.abs(rt - center) <= win))[0]
        if sel.size == 0:
            continue
        i = _pick(sel, rt, inten, center=center, return_index=True)
        rows.append({
            "name": t["name"], "class": t.get("class", ""), "adduct": t["adduct"],
            "mz": t["mz"], "rt_ref_min": t["rt_ref_min"], "irt": t["irt"],
            "rt_obs_min": float(rt[i]),
            "intensity": (float(inten[i]) if np.isfinite(inten[i]) else np.nan),
            "n_candidates": int(sel.size),
            "endogenous": bool(t["endogenous"]), "void": bool(t["void"]),
        })
    if not rows:
        raise AnchorIdentificationError(
            "No standards identified. Check polarity, m/z tolerance (mz_tol_ppm), "
            "the RT window, and that the standards are actually present/spiked."
        )
    anchors = pd.DataFrame(rows).sort_values("rt_ref_min").reset_index(drop=True)
    anchors = _drop_nonmonotone(anchors, config)
    return anchors


# ---------------------------------------------------------------- internals ---

def _pick(sel: np.ndarray, rt: np.ndarray, inten: np.ndarray,
          center: float, return_index: bool = False):
    """Pick one candidate: most intense if intensity known, else RT-closest.

    Deterministic tie-break by lowest original row index.
    """
    if np.any(np.isfinite(inten[sel])):
        order = sorted(sel, key=lambda k: (-(inten[k] if np.isfinite(inten[k]) else -np.inf), k))
    else:
        order = sorted(sel, key=lambda k: (abs(rt[k] - center), k))
    best = order[0]
    return best if return_index else rt[best]


def _enforce_monotone_template(template: Dict[str, float], panel: Panel) -> Dict[str, float]:
    """Keep only anchors whose native RT increases with reference RT (LIS)."""
    ordered = [(n, template[n]) for n in panel.targets["name"] if n in template]
    if not ordered:
        return template
    kept = _longest_increasing([rt for _, rt in ordered])
    return {ordered[i][0]: ordered[i][1] for i in kept}


def _drop_nonmonotone(anchors: pd.DataFrame, config: CalibrationConfig) -> pd.DataFrame:
    """Remove anchors that break strict increase of rt_obs with rt_ref (LIS)."""
    rts = anchors["rt_obs_min"].to_numpy()
    keep = _longest_increasing(rts, eps=config.tie_epsilon_min)
    dropped = len(anchors) - len(keep)
    if dropped > config.max_drop_anchors:
        # too many violations -> keep the LIS anyway but this signals a bad table
        pass
    out = anchors.iloc[keep].reset_index(drop=True)
    out.attrs["n_dropped_nonmonotone"] = dropped
    return out


def _longest_increasing(values, eps: float = 0.0) -> List[int]:
    """Indices of a longest strictly-increasing subsequence (patience sort)."""
    vals = list(values)
    n = len(vals)
    if n == 0:
        return []
    import bisect
    tails_val: List[float] = []
    tails_idx: List[int] = []
    prev = [-1] * n
    for i, v in enumerate(vals):
        pos = bisect.bisect_left(tails_val, v - eps)
        if pos == len(tails_val):
            tails_val.append(v)
            tails_idx.append(i)
        else:
            tails_val[pos] = v
            tails_idx[pos] = i
        prev[i] = tails_idx[pos - 1] if pos > 0 else -1
    # reconstruct
    res = []
    k = tails_idx[-1]
    while k != -1:
        res.append(k)
        k = prev[k]
    return res[::-1]
