"""Bundled standard mixtures — the panels this lab actually runs.

These used to live in the desktop app (``rtad/mixtures.py``), which meant the
engine and the UI each held their own idea of what "Mix 21" is. They live here
now so the two can never drift; ``rtad.mixtures`` is a thin re-export.

The manifests are pre-baked as canonical ``rt_anchor`` manifest frames (name,
class, exact_mass, adduct_pos, adduct_neg, rt_ref_min, endogenous, void) rather
than read from .xlsx at runtime — openpyxl is excluded from the frozen build,
and a spreadsheet that ships separately from the code is a spreadsheet that
eventually disagrees with it.

* ``mix15`` — the Caley lipid RT-calibration mix. Identical to
  :data:`rt_anchor.panel.DEFAULT_MANIFEST` (used verbatim, not copied).
  Source: 20260619_LipidRTCal_MZandRT.xlsx (Exp RT Orbitrap column).
* ``mix21`` — the extended Mix 4.4 panel. Source: "Mix 4.4 standards
  information for Bowen.xlsx". That sheet lists ion m/z directly rather than
  neutral masses, so ``exact_mass`` is back-computed from the POS m/z minus the
  adduct delta ([M+H]+ for PC, [M+NH4]+ for DG/CE/TG) — checked against the
  known masses of the overlapping species (< 0.001 Da). Reference RTs are the
  sheet's "Carly 35 min approximate times" (Column 25); approximate is fine, because the
  panel is only ever *located* in a real standards run, never trusted blind.
* ``none`` — the user's mixture is unknown, or there is none. Under the v2
  cross-column method this costs nothing for the calibration itself (stage 1
  matches features by m/z and claims no identities); it only means the engine
  cannot report detection QC, and must fall back to ``irt_landmark_panel`` for
  the iRT ruler. See the spec, §6.
* ``mix21lpc`` — ``mix21`` plus three LPC anchors (18:2 / 16:0 / 18:0). Kept
  **registered but not offered**, so older results and ``v2/5 chroms/iRT/
  run_iRT.py`` stay reproducible. The LPCs are endogenous serum lipids, not
  part of the mix; their ``rt_ref_min`` values were back-calculated from the
  Column22 (Holpie) STDs full-range self-calibrated model via
  ``rt_ref = 1 + iRT * 17.4 / 100``. That inversion belongs to the retired v1
  panel-warp iRT scale, which is why the panel is no longer offered.

Which panel is chosen moves the iRT scale, because mix15 and mix21 have
different earliest/latest detected landmarks on the reference run: iRT values
are comparable within a panel choice, not across one.

``endogenous`` / ``void`` are informational only (they are exported with the
anchors): shared species inherit the built-in panel's flags, and new species
follow the same convention — odd-chain surrogates and short-chain PCs are not
endogenous, natural long-chain species are.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import pandas as pd

from .panel import ADDUCTS, DEFAULT_MANIFEST, adduct_mz

MIX15 = "mix15"
MIX21 = "mix21"
MIX21LPC = "mix21lpc"
NONE = "none"

#: Panels offered as cards on the input screen, in display order.
OFFERED_MIXTURES: List[str] = [MIX15, MIX21, NONE]

#: Pre-selected panel.
DEFAULT_MIXTURE = MIX21

_MANIFEST_COLS = list(DEFAULT_MANIFEST.columns)


def _mix21_manifest() -> pd.DataFrame:
    """Build the Mix-4.4 (21-standard) manifest from its sheet values."""
    # name, class, pos_mz, neg_mz (None = not detected), rt_ref_min (Column 25),
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


def _mix21lpc_manifest() -> pd.DataFrame:
    """Mix-4.4 panel + three endogenous LPC anchors (legacy, not offered).

    The LPC reference RTs are not measured in any standards run — LPCs are not
    in the mix. They are back-calculated from the Column22 (Holpie) STDs
    full-range self-calibrated model via ``rt_ref = 1 + iRT * 17.4 / 100``,
    which pins them to the v1 panel's linear iRT scale.
    """
    # name, pos_mz ([M+H]+), rt_ref_min (iRT-derived, see docstring)
    rows = [
        ("LPC 18:2", 520.3398, 3.946),  # iRT 16.93
        ("LPC 16:0", 496.3398, 4.483),  # iRT 20.02
        ("LPC 18:0", 524.3711, 5.293),  # iRT 24.67
    ]
    out = []
    for name, pos_mz, rt in rows:
        out.append({
            "name": name,
            "class": "LPC",
            "exact_mass": round(float(pos_mz) - ADDUCTS["[M+H]+"], 6),
            "adduct_pos": "[M+H]+",
            "adduct_neg": "ND",
            "rt_ref_min": float(rt),
            "endogenous": True,
            "void": False,
        })
    lpc = pd.DataFrame(out, columns=_MANIFEST_COLS)
    man = pd.concat([_mix21_manifest(), lpc], ignore_index=True)
    return man.sort_values("rt_ref_min").reset_index(drop=True)


_MANIFESTS: Dict[str, Optional[pd.DataFrame]] = {
    MIX15: DEFAULT_MANIFEST.loc[:, _MANIFEST_COLS],
    MIX21: _mix21_manifest(),
    MIX21LPC: _mix21lpc_manifest(),
    NONE: None,
}

_META: Dict[str, Dict] = {
    MIX15: {
        "key": MIX15,
        "label": "15",
        "sub": "Li lab 15 lipids mixture",
        "source": "20260619_LipidRTCal_MZandRT.xlsx",
        "method_note": "Orbitrap reference RT",
        "offered": True,
    },
    MIX21: {
        "key": MIX21,
        "label": "21",
        "sub": "Li lab 21 lipids mixture",
        "source": "Mix 4.4 standards information for Bowen.xlsx",
        "method_note": "Column 25 method, approximate RT",
        "offered": True,
    },
    NONE: {
        "key": NONE,
        "label": "Others",
        "sub": "other mixtures",
        "source": "",
        "method_note": "No identity is claimed for the mixture; "
                       "stage-1 calibration is unaffected",
        "offered": True,
    },
    MIX21LPC: {
        "key": MIX21LPC,
        "label": "Mix 21 + LPC",
        "sub": "Extended panel + LPC anchors (18:2/16:0/18:0)",
        "source": "Mix 4.4 + LPC anchors derived from Column22 STDs full-range model",
        "method_note": "Column 25 method, approximate RT (legacy, not offered)",
        "offered": False,
    },
}


def mixture_keys() -> List[str]:
    """Every registered key, including the ones not offered as cards."""
    return list(_MANIFESTS)


def offered_keys() -> List[str]:
    """The keys the input screen offers, in display order."""
    return list(OFFERED_MIXTURES)


def get_meta(key: str) -> Dict:
    if key not in _META:
        raise KeyError(f"Unknown mixture '{key}'. Known: {mixture_keys()}")
    return dict(_META[key])


def get_manifest(key: str) -> Optional[pd.DataFrame]:
    """A fresh copy of the mixture's manifest frame, or ``None`` for ``none``.

    ``None`` is a real answer, not a failure: the user told us they have no
    identifiable panel, and every caller must handle that (skip detection QC,
    fall back to the configured landmark panel for iRT).
    """
    if key not in _MANIFESTS:
        raise KeyError(f"Unknown mixture '{key}'. Known: {mixture_keys()}")
    man = _MANIFESTS[key]
    return None if man is None else man.copy()


def n_standards(key: str) -> int:
    man = _MANIFESTS[key] if key in _MANIFESTS else None
    return 0 if man is None else int(len(man))


def preview_payload(keys: Optional[List[str]] = None) -> List[Dict]:
    """JSON-safe description of each mixture: card meta + one row per standard.

    Row fields: name, class, mz_pos, mz_neg (None when not ionisable in that
    polarity), rt_ref_min. Defaults to the offered panels — pass ``keys`` to
    include a registered-but-not-offered one such as ``mix21lpc``.
    """
    out = []
    for key in (offered_keys() if keys is None else keys):
        man = _MANIFESTS[key]
        meta = get_meta(key)
        if man is None:
            out.append({**meta, "n_standards": 0, "classes": [],
                        "rt_range": None, "rows": []})
            continue
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
