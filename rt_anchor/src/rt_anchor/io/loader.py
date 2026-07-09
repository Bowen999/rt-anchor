"""Unified entry point: detect format, dispatch to the adapter, sanity-check RT."""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np

from ..errors import InputFormatError, RTUnitError
from .adapters import _ADAPTERS, detect_format, load_msdial
from .schema import FeatureTable


def load_feature_table(path: str,
                       source_format: Optional[str] = None,
                       polarity: Optional[str] = None,
                       rt_unit: Optional[str] = None) -> FeatureTable:
    """Load any supported feature table into the canonical representation.

    Parameters
    ----------
    path : str
        The exported feature/alignment table.
    source_format : {'masscube','lipidscreener','msdial','mzmine', None}
        Force a format; ``None`` auto-detects.
    polarity : {'positive','negative', None}
        Declared polarity (used to pick the right adduct set downstream).
    rt_unit : {'min','sec', None}
        Override the RT unit; ``None`` uses the format default + a range check.
    """
    if source_format is None:
        fmt, delim, header_row = detect_format(path)
    else:
        fmt = source_format.lower()
        if fmt not in _ADAPTERS:
            raise InputFormatError(f"Unknown source_format '{source_format}'. Known: {sorted(_ADAPTERS)}")
        # still sniff the delimiter / header row
        try:
            _, delim, header_row = detect_format(path)
        except InputFormatError:
            delim, header_row = ("\t" if fmt in ("masscube", "msdial") else ","), 0

    if fmt == "msdial":
        ft = load_msdial(path, delim=delim, header_row=header_row, polarity=polarity)
    else:
        ft = _ADAPTERS[fmt](path, delim=delim, polarity=polarity)

    if rt_unit is not None:
        if rt_unit not in ("min", "sec"):
            raise RTUnitError(f"rt_unit must be 'min' or 'sec', got '{rt_unit}'.")
        ft.rt_unit = rt_unit

    _sanity_check_rt(ft)
    return ft


def _sanity_check_rt(ft: FeatureTable) -> None:
    """Warn if the RT range looks inconsistent with the declared unit.

    A typical LC gradient is minutes-scale (< ~90 min). If a table declared in
    minutes spans thousands, it is almost certainly seconds (and vice-versa).
    We only *warn* — the unit can always be forced via ``rt_unit=``.
    """
    rt = ft.rt_raw().to_numpy(dtype=float)
    rt = rt[np.isfinite(rt)]
    if rt.size == 0:
        return
    hi = float(np.nanmax(rt))
    if ft.rt_unit == "min" and hi > 90:
        warnings.warn(
            f"[rt_anchor] RT declared in minutes but max is {hi:.0f} — this looks "
            f"like seconds. Pass rt_unit='sec' if so (source={ft.source_format}).",
            stacklevel=2,
        )
    if ft.rt_unit == "sec" and hi < 25:
        warnings.warn(
            f"[rt_anchor] RT declared in seconds but max is {hi:.1f} — this looks "
            f"like minutes. Pass rt_unit='min' if so (source={ft.source_format}).",
            stacklevel=2,
        )
