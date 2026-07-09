"""Canonical in-memory representation of a feature table.

A :class:`FeatureTable` wraps the **full, unmodified** source DataFrame (no
columns are ever dropped) plus a small amount of resolved metadata telling the
rest of the package *where* m/z and RT live and *which* columns are per-sample
intensities. Calibration only ever **appends** columns to ``df``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from ..errors import RTUnitError

RESULT_COLUMNS = [
    "RI", "RI_uncertainty", "RI_reliability", "RI_spread", "n_contributing",
    "is_extrapolated", "calibration_scope", "warp_source",
]


@dataclass
class FeatureTable:
    df: pd.DataFrame                 # full original table, all columns preserved
    mz_col: str
    rt_col: str
    rt_unit: str                    # 'min' or 'sec'
    sample_cols: List[str] = field(default_factory=list)
    source_format: str = "unknown"
    polarity: Optional[str] = None
    meta: Dict = field(default_factory=dict)   # e.g. LipidScreener group row, MS-DIAL class header

    # ---- accessors ----
    def mz(self) -> pd.Series:
        return pd.to_numeric(self.df[self.mz_col], errors="coerce")

    def rt_raw(self) -> pd.Series:
        return pd.to_numeric(self.df[self.rt_col], errors="coerce")

    def rt_minutes(self) -> pd.Series:
        rt = self.rt_raw()
        if self.rt_unit == "min":
            return rt
        if self.rt_unit == "sec":
            return rt / 60.0
        raise RTUnitError(f"Unknown rt_unit '{self.rt_unit}' (expected 'min' or 'sec').")

    def n_features(self) -> int:
        return len(self.df)

    def summary(self) -> Dict:
        rt = self.rt_minutes()
        return {
            "source_format": self.source_format,
            "polarity": self.polarity,
            "n_features": int(self.n_features()),
            "mz_col": self.mz_col,
            "rt_col": self.rt_col,
            "rt_unit": self.rt_unit,
            "rt_min_minutes": None if rt.empty else float(np.nanmin(rt)),
            "rt_max_minutes": None if rt.empty else float(np.nanmax(rt)),
            "n_sample_columns": len(self.sample_cols),
            "sample_columns": list(self.sample_cols),
            "n_total_columns": len(self.df.columns),
        }
