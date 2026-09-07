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

#: Columns the calibration appends to the user's table, in order (spec §1).
#: ``Cal_RT_min`` is the feature's RT on the reference column's time axis;
#: ``iRT`` is the dimensionless 1-100 index derived from it. The two are
#: independent outputs — a run with no detectable landmark panel still gets a
#: ``Cal_RT_min``, with ``iRT`` NaN and ``iRT_reliability = "none"``.
RESULT_COLUMNS = [
    "Cal_RT_min", "Cal_RT_uncertainty_min",
    "iRT", "iRT_uncertainty", "iRT_reliability",
    "is_extrapolated", "calibration_scope", "warp_source",
    "RI_spread", "n_contributing",
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

    # ---- optional per-row attributes ----
    #
    # The cross-column engine ranks candidates by abundance, prefers fully
    # filled features and checks adduct consistency. Only MS-DIAL exports carry
    # all four of those; the others carry some or none. Each accessor therefore
    # resolves the best available source and returns ``None`` when the concept
    # is genuinely absent, so callers can *skip* the corresponding filter rather
    # than invent a value for it (a fabricated Fill% or adduct would silently
    # change which isomer is picked).

    def intensity(self) -> pd.Series:
        """Per-row abundance used to rank candidates in an m/z window.

        Sum over the declared sample columns when there are any, else a single
        abundance column, else the MS-DIAL ``S/N average`` (the only abundance
        proxy an alignment export without sample columns has), else NaN.
        """
        # local import: identify imports this module, so bind at call time
        from ..identify import _intensity
        val = _intensity(self)
        if val.notna().any():
            return val
        sn = self._numeric_col("s/n average", "s/n", "sn")
        return val if sn is None else sn

    def sn(self) -> pd.Series:
        """Signal-to-noise, falling back to :meth:`intensity`.

        Used where the reference implementation ranked by ``S/N average``:
        the plasma-lipid pick and the isomer census. Formats without an S/N
        column rank by abundance instead, which orders candidates the same way
        in practice (both are dominated by the peak's own height).
        """
        sn = self._numeric_col("s/n average", "s/n", "sn")
        return self.intensity() if sn is None else sn

    def adducts(self) -> Optional[pd.Series]:
        """Adduct annotation as text, or ``None`` when the table has none.

        ``None`` means "skip the adduct-consistency filter" — never "no adduct
        matched".
        """
        col = self._find_col("adduct type", "adduct", "main_adduct")
        return None if col is None else self.df[col].astype(str)

    def fill(self) -> Optional[pd.Series]:
        """Fraction of samples in which the feature was detected (0-1).

        MS-DIAL ``Fill %`` or MassCube ``detection_rate``; ``None`` when the
        table has neither, in which case the "prefer fully filled features"
        step is skipped.
        """
        return self._numeric_col("fill %", "fill%", "detection_rate")

    def feature_id(self) -> pd.Series:
        """Stable per-row identifier used to tell candidates apart.

        MS-DIAL ``Alignment ID``, MZmine ``row ID``, MassCube ``feature_id``;
        otherwise the row index, which is stable because ``df`` is never
        reordered.
        """
        col = self._find_col("alignment id", "feature_id", "row id", "id")
        if col is None:
            return pd.Series(self.df.index, index=self.df.index)
        ids = pd.to_numeric(self.df[col], errors="coerce")
        if ids.isna().all():
            return self.df[col].astype(str)
        return ids

    def n_features(self) -> int:
        return len(self.df)

    # ---- internals ----
    def _find_col(self, *names: str) -> Optional[str]:
        low = {str(c).strip().lower(): c for c in self.df.columns}
        for n in names:
            if n in low:
                return low[n]
        return None

    def _numeric_col(self, *names: str) -> Optional[pd.Series]:
        col = self._find_col(*names)
        return None if col is None else pd.to_numeric(self.df[col], errors="coerce")

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
