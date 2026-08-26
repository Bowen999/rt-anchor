"""Bundled standard mixtures (the two panels this lab actually runs).

The input screen offers an either/or choice instead of a free-form manifest
file. Both manifests are pre-baked here (NOT read from .xlsx at runtime —
openpyxl is excluded from the frozen build) as canonical ``rt_anchor``
manifest frames: name, class, exact_mass, adduct_pos, adduct_neg,
rt_ref_min, endogenous, void.

* ``mix15`` — the Caley lipid RT-calibration mix. Identical to the engine's
  built-in ``rt_anchor.panel.DEFAULT_MANIFEST`` (used verbatim so the two can
  never drift apart). Source: 20260619_LipidRTCal_MZandRT.xlsx (Exp RT
  Orbitrap column).
* ``mix21`` — the extended Mix 4.4 panel. Source: v2/"Mix 4.4 standards
  information for Bowen.xlsx". That sheet lists ion m/z directly (not exact
  masses), so ``exact_mass`` is back-computed here from the POS m/z minus the
  adduct delta ([M+H]+ for PC, [M+NH4]+ for DG/CE/TG) — verified against the
  known masses of the overlapping species (<0.001 Da). Reference RTs are the
  sheet's "Carly 35 min approximate times"; approximate is fine because the
  native template from the actual standards run locates the anchors.

``endogenous``/``void`` are informational only (exported with the anchors):
shared species inherit the built-in panel's flags, new species follow the same
convention (odd-chain surrogates and short-chain PCs -> not endogenous;
natural long-chain species -> endogenous).
"""

from __future__ import annotations

from typing import Dict, List

import pandas as pd

from rt_anchor.panel import ADDUCTS, DEFAULT_MANIFEST, adduct_mz

MIX15 = "mix15"
MIX21 = "mix21"
DEFAULT_MIXTURE = MIX21          # pre-selected on the input screen

_MANIFEST_COLS = list(DEFAULT_MANIFEST.columns)


def _mix21_manifest() -> pd.DataFrame:
    """Build the Mix-4.4 (21-standard) manifest from its sheet values."""
    # name, class, pos_mz, neg_mz (None = not detected), rt_ref_min (Carly 35-min),
    # endogenous, void
    rows = [
        ("PC 6:0_6:0",        "PC", 454.2561, 498.2474,  1.5, False, True),
        ("PC 7:0_7:0",        "PC", 482.2877, 526.2787,  2.2, False, False),
        ("PC 8:0_8:0",        "PC", 510.3190, 554.3100,  3.2, False, False),
        ("PC 9:0_9:0",        "PC", 538.3502, 582.3413,  4.1, False, True),
        ("PC 10:0_10:0",      "PC", 566.3816, 610.3726,  4.9, False, False),
        ("PC 11:0_11:0",      "PC", 594.4128, 638.4039,  5.8, False, False),
        ("PC 12:0_12:0",      "PC", 622.4440, 666.4352,  6.8, False, False),
        ("PC 13:0_13:0",      "PC", 650.4754, 694.4665,  8.1, False, False),
        ("PC 14:0_14:0",      "PC", 678.5065, 722.4978,  9.5, False, False),
        ("PC 14:0_16:0",      "PC", 706.5381, 750.5291, 11.1, True,  False),
        ("PC 16:0_18:2",      "PC", 758.5688, 802.5604, 11.7, True,  False),
        ("PC 16:0_18:1",      "PC", 760.5846, 804.5760, 12.9, True,  False),
        ("PC 18:0_18:1",      "PC", 788.6161, 832.6073, 14.5, True,  False),
        ("PC 18:0_18:0",      "PC", 790.6315, 834.6230, 15.8, True,  False),
        ("DG 16:0_18:0",      "DG", 614.5718, 641.5362, 17.0, True,  False),
        ("DG 18:0_18:0",      "DG", 642.6031, 669.5675, 18.0, True,  False),
        ("CE 20:5",           "CE", 688.6027, None,     20.5, False, False),
        ("CE 18:2",           "CE", 666.6183, None,     21.2, False, False),
        ("TG 16:0_18:1_16:0", "TG", 850.7858, None,     21.3, True,  False),
        ("CE 18:1",           "CE", 668.6345, None,     21.6, False, False),
        ("TG 18:0_18:2_18:0", "TG", 904.8328, None,     21.8, True,  False),
    ]
    out = []
    for name, cls, pos_mz, neg_mz, rt, endo, void in rows:
        adduct_pos = "[M+H]+" if cls == "PC" else "[M+NH4]+"
        out.append({
            "name": name,
            "class": cls,
            # neutral mass back-computed from the sheet's POS m/z
            "exact_mass": round(float(pos_mz) - ADDUCTS[adduct_pos], 6),
            "adduct_pos": adduct_pos,
            "adduct_neg": "[M+HCOO]-" if neg_mz is not None else "ND",
            "rt_ref_min": float(rt),
            "endogenous": bool(endo),
            "void": bool(void),
        })
    return pd.DataFrame(out, columns=_MANIFEST_COLS)


_MANIFESTS: Dict[str, pd.DataFrame] = {
    MIX15: DEFAULT_MANIFEST.loc[:, _MANIFEST_COLS],
    MIX21: _mix21_manifest(),
}

_META: Dict[str, Dict] = {
    MIX15: {
        "key": MIX15,
        "label": "Mix 15",
        "sub": "Caley lipid RT-calibration mix",
        "source": "20260619_LipidRTCal_MZandRT.xlsx",
        "method_note": "Orbitrap reference RT",
    },
    MIX21: {
        "key": MIX21,
        "label": "Mix 21",
        "sub": "Extended panel (Mix 4.4)",
        "source": "Mix 4.4 standards information for Bowen.xlsx",
        "method_note": "Carly 35-min method, approximate RT",
    },
}


def mixture_keys() -> List[str]:
    return list(_MANIFESTS)


def get_meta(key: str) -> Dict:
    if key not in _META:
        raise KeyError(f"Unknown mixture '{key}'. Known: {mixture_keys()}")
    return dict(_META[key])


def get_manifest(key: str) -> pd.DataFrame:
    """A fresh copy of the mixture's manifest frame (engine-ready)."""
    if key not in _MANIFESTS:
        raise KeyError(f"Unknown mixture '{key}'. Known: {mixture_keys()}")
    return _MANIFESTS[key].copy()


def preview_payload() -> List[Dict]:
    """JSON-safe description of every mixture: card meta + one row per standard.

    Row fields: name, class, mz_pos, mz_neg (None when not detected), rt_ref_min.
    """
    out = []
    for key in mixture_keys():
        man = _MANIFESTS[key]
        meta = get_meta(key)
        classes = sorted(set(man["class"]))
        rts = man["rt_ref_min"].astype(float)
        rows = []
        for _, r in man.iterrows():
            mz_pos = adduct_mz(r["exact_mass"], r["adduct_pos"])
            mz_neg = adduct_mz(r["exact_mass"], r["adduct_neg"])
            rows.append({
                "name": r["name"],
                "class": r.get("class", ""),
                "mz_pos": round(mz_pos, 4) if mz_pos is not None else None,
                "mz_neg": round(mz_neg, 4) if mz_neg is not None else None,
                "rt_ref_min": float(r["rt_ref_min"]),
            })
        out.append({
            **meta,
            "n_standards": int(len(man)),
            "classes": classes,
            "rt_range": [round(float(rts.min()), 2), round(float(rts.max()), 2)],
            "rows": rows,
        })
    return out
