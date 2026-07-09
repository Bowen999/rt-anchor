"""Standard panel: manifest (what the standards are) + reference baseline (the
frozen iRT ruler).

Two things live here:

* ``DEFAULT_MANIFEST`` — the 15-lipid panel: name, class, exact (neutral
  monoisotopic) mass, and the expected adduct per polarity. Target ion m/z is
  computed from ``exact_mass`` + adduct, so no molecular formula is needed.
* ``DEFAULT_REFERENCE`` — per-anchor Orbitrap reference RT (minutes) → this
  defines the dimensionless iRT scale via the fixed affine map (see config).

Both can be overridden by the user (Mode B): ``load_manifest_csv`` /
``load_reference_csv`` read a CSV with the same columns.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .config import CalibrationConfig
from .errors import PanelError

# adduct m/z deltas (monoisotopic; electron mass included for the ion)
_ELECTRON = 0.000548579909
ADDUCTS: Dict[str, float] = {
    "[M+H]+": 1.0072764,           # proton
    "[M+NH4]+": 18.0338254,        # NH4+  (NH4 - e)
    "[M+Na]+": 22.9892214,
    "[M-H]-": -1.0072764,
    "[M+HCOO]-": 44.9982010,       # formate
    "[M+CH3COO]-": 59.0138510,
    "[M+Cl]-": 34.9694011,
}


def adduct_mz(exact_mass: float, adduct: str) -> Optional[float]:
    """Compute ion m/z from neutral monoisotopic mass and an adduct label."""
    if adduct is None or str(adduct).upper() in ("ND", "NAN", "NONE", ""):
        return None
    if adduct not in ADDUCTS:
        raise PanelError(
            f"Unknown adduct '{adduct}'. Known: {sorted(ADDUCTS)}. "
            f"Add it to rt_anchor.panel.ADDUCTS or use a supported adduct."
        )
    return float(exact_mass) + ADDUCTS[adduct]


# ---- default 15-standard panel (Caley_RT lipid RT-calibration mix) ----
# columns: name, class, exact_mass, adduct_pos, adduct_neg, rt_ref_min, endogenous, void
_DEFAULT_ROWS = [
    ("PC 6:0_6:0",   "PC", 453.249150, "[M+H]+",   "[M+HCOO]-",  1.06, False, True),
    ("PC 9:0_9:0",   "PC", 537.343055, "[M+H]+",   "[M+HCOO]-",  1.81, False, True),
    ("PC 11:0_11:0", "PC", 593.405655, "[M+H]+",   "[M+HCOO]-",  2.88, False, False),
    ("PC 12:0_12:0", "PC", 621.436955, "[M+H]+",   "[M+HCOO]-",  3.60, False, False),
    ("PC 13:0_13:0", "PC", 649.468256, "[M+H]+",   "[M+HCOO]-",  4.52, False, False),
    ("PC 14:0_14:0", "PC", 677.499556, "[M+H]+",   "[M+HCOO]-",  5.63, False, False),
    ("PC 16:0_18:3", "PC", 755.546506, "[M+H]+",   "[M+HCOO]-",  6.63, True,  False),
    ("PC 16:0_18:2", "PC", 757.562156, "[M+H]+",   "[M+HCOO]-",  7.37, True,  False),
    ("PC 16:0_18:1", "PC", 759.577806, "[M+H]+",   "[M+HCOO]-",  8.40, True,  False),
    ("PC 18:0_18:1", "PC", 787.609106, "[M+H]+",   "[M+HCOO]-",  9.79, True,  False),
    ("PC 18:0_18:0", "PC", 789.624756, "[M+H]+",   "[M+HCOO]-", 11.07, True,  False),
    ("DG 16:0_18:0", "DG", 596.537976, "[M+NH4]+", "[M+HCOO]-", 12.65, True,  False),
    ("DG 18:0_18:0", "DG", 624.569276, "[M+NH4]+", "[M+HCOO]-", 14.99, True,  False),
    ("CE 20:5",      "CE", 670.568881, "[M+NH4]+", "ND",        17.37, False, False),
    ("CE 18:1",      "CE", 650.600181, "[M+NH4]+", "ND",        18.39, False, False),
]
_MANIFEST_COLS = ["name", "class", "exact_mass", "adduct_pos", "adduct_neg",
                  "rt_ref_min", "endogenous", "void"]

DEFAULT_MANIFEST = pd.DataFrame(_DEFAULT_ROWS, columns=_MANIFEST_COLS)

# Canonical SMILES for the default panel (validated: rdkit exact mass matches the
# manifest to <0.05 mmu). Used only for the optional structure-on-hover preview.
_FA = {"6:0": "CCCCC", "9:0": "CCCCCCCC", "11:0": "CCCCCCCCCC", "12:0": "CCCCCCCCCCC",
       "13:0": "CCCCCCCCCCCC", "14:0": "CCCCCCCCCCCCC", "16:0": "CCCCCCCCCCCCCCC",
       "18:0": "CCCCCCCCCCCCCCCCC", "18:1": r"CCCCCCC/C=C\CCCCCCCC",
       "18:2": r"CCCCCCC/C=C\C/C=C\CCCCC", "18:3": r"CCCCCCC/C=C\C/C=C\C/C=C\CC",
       "20:5": r"CCC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CC"}


def _pc(a, b):
    return f"[C@H](COC(=O){_FA[a]})(OC(=O){_FA[b]})COP(=O)([O-])OCC[N+](C)(C)C"


def _dg(a, b):
    return f"[C@H](COC(=O){_FA[a]})(OC(=O){_FA[b]})CO"


def _ce(a):
    return f"CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](OC(=O){_FA[a]})CC[C@]4(C)[C@H]3CC[C@]12C"


DEFAULT_SMILES = {
    "PC 6:0_6:0": _pc("6:0", "6:0"), "PC 9:0_9:0": _pc("9:0", "9:0"),
    "PC 11:0_11:0": _pc("11:0", "11:0"), "PC 12:0_12:0": _pc("12:0", "12:0"),
    "PC 13:0_13:0": _pc("13:0", "13:0"), "PC 14:0_14:0": _pc("14:0", "14:0"),
    "PC 16:0_18:3": _pc("16:0", "18:3"), "PC 16:0_18:2": _pc("16:0", "18:2"),
    "PC 16:0_18:1": _pc("16:0", "18:1"), "PC 18:0_18:1": _pc("18:0", "18:1"),
    "PC 18:0_18:0": _pc("18:0", "18:0"), "DG 16:0_18:0": _dg("16:0", "18:0"),
    "DG 18:0_18:0": _dg("18:0", "18:0"), "CE 20:5": _ce("20:5"), "CE 18:1": _ce("18:1"),
}


def is_default_manifest(manifest) -> bool:
    """True when the manifest is the built-in 15-standard panel (enables structure hover)."""
    if manifest is None:
        return True
    try:
        return set(manifest["name"]) == set(DEFAULT_MANIFEST["name"])
    except Exception:
        return False


@dataclass
class Panel:
    """Resolved panel for one polarity: targets (name, m/z, iRT, flags)."""
    polarity: str
    targets: pd.DataFrame  # name, mz, rt_ref_min, irt, endogenous, void, adduct

    def __len__(self) -> int:
        return len(self.targets)


def build_panel(polarity: str,
                config: CalibrationConfig,
                manifest: Optional[pd.DataFrame] = None,
                reference: Optional[pd.DataFrame] = None) -> Panel:
    """Resolve the manifest + reference into per-polarity ion targets with iRT.

    ``manifest`` overrides identities (Mode B); ``reference`` overrides the RT
    ruler. Standards with adduct 'ND' in this polarity are dropped.
    """
    pol = _norm_polarity(polarity)
    man = DEFAULT_MANIFEST if manifest is None else manifest
    _validate_manifest(man)

    adduct_col = "adduct_pos" if pol == "positive" else "adduct_neg"

    # reference RT (minutes) per anchor: from override if given, else manifest rt_ref_min
    if reference is not None:
        _validate_reference(reference)
        ref_map = dict(zip(reference["anchor"], reference["RT_ref_min"]))
    else:
        ref_map = dict(zip(man["name"], man["rt_ref_min"]))

    rows = []
    for _, r in man.iterrows():
        adduct = r[adduct_col]
        mz = adduct_mz(r["exact_mass"], adduct)
        if mz is None:
            continue  # not ionisable in this polarity
        rt_ref = ref_map.get(r["name"])
        if rt_ref is None or (isinstance(rt_ref, float) and np.isnan(rt_ref)):
            continue
        rows.append({
            "name": r["name"],
            "class": r.get("class", ""),
            "adduct": adduct,
            "mz": mz,
            "rt_ref_min": float(rt_ref),
            "irt": config.irt_from_rt(float(rt_ref)),
            "endogenous": bool(r.get("endogenous", False)),
            "void": bool(r.get("void", False)),
        })
    if not rows:
        raise PanelError(f"No standards are ionisable in polarity '{pol}'.")
    targets = pd.DataFrame(rows).sort_values("rt_ref_min").reset_index(drop=True)
    return Panel(polarity=pol, targets=targets)


def _norm_polarity(p: str) -> str:
    s = str(p).strip().lower()
    if s in ("positive", "pos", "+", "p"):
        return "positive"
    if s in ("negative", "neg", "-", "n"):
        return "negative"
    raise PanelError(f"Unrecognised polarity '{p}'. Use 'positive' or 'negative'.")


def _validate_manifest(man: pd.DataFrame) -> None:
    need = {"name", "exact_mass", "adduct_pos", "adduct_neg"}
    missing = need - set(man.columns)
    if missing:
        raise PanelError(
            f"Manifest is missing required columns: {sorted(missing)}. "
            f"Required: name, class, exact_mass, adduct_pos, adduct_neg, rt_ref_min."
        )


def _validate_reference(ref: pd.DataFrame) -> None:
    need = {"anchor", "RT_ref_min"}
    missing = need - set(ref.columns)
    if missing:
        raise PanelError(
            f"Reference baseline is missing columns: {sorted(missing)}. "
            f"Required: anchor, RT_ref_min (iRT is recomputed from the scale constants)."
        )


def load_manifest_csv(path: str) -> pd.DataFrame:
    man = pd.read_csv(path)
    _validate_manifest(man)
    return man


def load_reference_csv(path: str) -> pd.DataFrame:
    ref = pd.read_csv(path, comment="#")
    _validate_reference(ref)
    return ref
