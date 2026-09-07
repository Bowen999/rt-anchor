"""Format detection and per-software adapters.

Supported input formats (auto-detected):

* ``masscube``      — MassCube ``aligned_feature_table.txt`` (tab; RT in min)
* ``lipidscreener`` — LipidScreener ``*_data_cleansing_complete_*.csv`` (comma;
  RT in **seconds**; a "Group" row directly under the header)
* ``msdial``        — MS-DIAL alignment result (tab/comma; several metadata rows
  above the real header; ``Average Rt(min)`` / ``Average Mz``)
* ``mzmine``        — MZmine 2 (``row m/z`` / ``row retention time`` / ``... Peak
  area``) and MZmine 3 (``mz`` / ``rt`` / ``datafile:...:area``)

Each adapter returns a :class:`FeatureTable` whose ``df`` still contains **every
original column** (no columns are dropped). Format-specific side metadata (the
LipidScreener group row, the MS-DIAL class header) is preserved under ``meta``.
"""

from __future__ import annotations

import io
import re
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

from ..errors import ColumnResolutionError, InputFormatError
from .schema import FeatureTable

# ---------------------------------------------------------------- sniffing ----

_MSDIAL_HEADER_TOKENS = ("Average Rt", "Average Mz", "Alignment ID")


def _read_lines(path: str, n: int = 8) -> List[str]:
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        return [next(fh, "") for _ in range(n)]


def _delimiter(sample_lines: List[str]) -> str:
    joined = "\n".join(sample_lines)
    return "\t" if joined.count("\t") >= joined.count(",") else ","


def detect_format(path: str) -> Tuple[str, str, int]:
    """Return ``(format, delimiter, header_row_index)``.

    ``header_row_index`` is the 0-based line where the real column header lives
    (non-zero only for MS-DIAL, which prepends metadata rows).
    """
    head = _read_lines(path, 8)
    if not head or head[0] == "":
        raise InputFormatError(f"'{path}' appears to be empty.")
    delim = _delimiter(head)

    # MS-DIAL: a header token appears within the first few lines (possibly not line 0)
    for i, line in enumerate(head):
        if any(tok in line for tok in _MSDIAL_HEADER_TOKENS):
            # if the token is not on line 0, it's the classic MS-DIAL metadata block
            if i > 0 or ("Average Rt" in line and "Alignment ID" in line):
                return "msdial", delim, i

    header = [c.strip() for c in head[0].rstrip("\n").split(delim)]
    hset = {c.lower() for c in header}

    # LipidScreener: rt + mz + maxo/main_adduct/type
    if {"rt", "mz"} <= hset and ({"maxo", "main_adduct", "type", "neutral_mass"} & hset):
        return "lipidscreener", delim, 0

    # MassCube: m/z + RT (+ its structural columns)
    if "m/z" in hset and "rt" in hset:
        return "masscube", delim, 0

    # MZmine v3: mz + rt + datafile: columns  ;  v2: row m/z + row retention time
    if {"row m/z", "row retention time"} <= hset or {"row m/z", "row rt"} <= hset:
        return "mzmine", delim, 0
    if {"mz", "rt"} <= hset and any(":" in c or ".mzml" in c.lower() for c in header):
        return "mzmine", delim, 0

    raise InputFormatError(
        f"Could not recognise the format of '{path}'. "
        f"Header seen: {header[:12]}{'...' if len(header) > 12 else ''}. "
        f"Supported: MS-DIAL, MZmine, MassCube, LipidScreener. "
        f"You can bypass detection by passing source_format=... to load_feature_table()."
    )


# ---------------------------------------------------------------- helpers -----

def _find(cols: List[str], *candidates: str, contains: Optional[str] = None) -> Optional[str]:
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    if contains:
        for c in cols:
            if contains.lower() in c.lower():
                return c
    return None


def _rt_unit_from_name(name: str, default: str = "min") -> str:
    n = name.lower()
    if "sec" in n or "(s)" in n:
        return "sec"
    if "min" in n:
        return "min"
    return default


# ---------------------------------------------------------------- adapters ----

# MassCube metadata (non-sample) columns; anything after the last of these is a sample
_MASSCUBE_META = {
    "group_id", "feature_id", "m/z", "mz", "rt", "adduct", "isotope_state",
    "is_isotope", "is_in_source_fragment", "scan_scan_cor", "gaussian_similarity",
    "noise_score", "asymmetry_factor", "peak_shape", "detection_rate",
    "detection_rate_gap_filled", "alignment_reference_file", "isotopes",
    "ms2_reference_file", "ms2_scan_id", "ms2", "precursor_ion_fraction",
    "matched_ms2", "search_mode", "annotation", "formula", "similarity",
    "matched_peak_number", "smiles", "inchikey", "database",
}


def load_masscube(path: str, delim: str = "\t", polarity: Optional[str] = None) -> FeatureTable:
    df = pd.read_csv(path, sep=delim, low_memory=False)
    cols = list(df.columns)
    mz = _find(cols, "m/z", "mz")
    rt = _find(cols, "RT", "rt", "retention time")
    if mz is None or rt is None:
        raise ColumnResolutionError(f"MassCube: could not find m/z and RT columns in {cols[:10]}")
    # Sample columns = everything after the LAST recognised metadata column. This
    # is robust to trimmed tables where 'database' has been dropped.
    last_meta = -1
    for i, c in enumerate(cols):
        if c.strip().lower() in _MASSCUBE_META:
            last_meta = i
    sample_cols = cols[last_meta + 1:] if last_meta >= 0 else []
    return FeatureTable(df=df, mz_col=mz, rt_col=rt, rt_unit="min",
                        sample_cols=sample_cols, source_format="masscube", polarity=polarity)


def load_lipidscreener(path: str, delim: str = ",", polarity: Optional[str] = None) -> FeatureTable:
    df = pd.read_csv(path, sep=delim, low_memory=False, dtype=str)
    cols = list(df.columns)
    mz = _find(cols, "mz", "m/z")
    rt = _find(cols, "rt", "retention time")
    if mz is None or rt is None:
        raise ColumnResolutionError(f"LipidScreener: could not find mz/rt in {cols[:10]}")

    meta = {}
    # A leading "Group" row labels the sample columns; capture then drop it (it is
    # not a feature). This removes a ROW, never a column.
    first = df.iloc[0] if len(df) else None
    if first is not None and str(first.iloc[0]).strip().lower() == "group":
        meta["group_row"] = {c: first[c] for c in cols if isinstance(first[c], str) and first[c].strip()}
        df = df.iloc[1:].reset_index(drop=True)

    # sample columns follow 'type'; fall back to group-row keys
    sample_cols: List[str] = []
    if "type" in cols:
        idx = cols.index("type")
        sample_cols = cols[idx + 1:]
    elif "group_row" in meta:
        meta_cols = {"No", "rt", "mz", "maxo", "sn", "polarity", "MS2_No",
                     "adducts", "main_adduct", "neutral_mass", "width", "type"}
        sample_cols = [c for c in cols if c not in meta_cols]

    pol = polarity
    if pol is None and "polarity" in df.columns:
        vals = df["polarity"].dropna().unique()
        # LipidScreener stores 1.0/0.0 or +/-; leave as-is unless clearly labelled
        pol = None
    return FeatureTable(df=df, mz_col=mz, rt_col=rt, rt_unit="sec",
                        sample_cols=sample_cols, source_format="lipidscreener",
                        polarity=pol, meta=meta)


def _is_sample_class_label(label) -> bool:
    """True when an MS-DIAL class-header cell marks a real per-sample column.

    MS-DIAL appends ``Average`` and ``Stdev`` summary columns after the sample
    columns and labels them ``NA`` in the class row. Counting them as samples
    inflates every abundance we derive from the sample columns (and therefore
    which candidate wins inside an m/z window), so empty, ``NA`` and ``NaN``
    labels are all treated as "not a sample".
    """
    s = str(label).strip()
    return bool(s) and s.lower() not in ("na", "nan", "n/a")


def load_msdial(path: str, delim: str = "\t", header_row: int = 0,
                polarity: Optional[str] = None) -> FeatureTable:
    # metadata rows (0..header_row-1) carry the per-sample class labels
    meta = {}
    class_labels = None
    if header_row > 0:
        # keep_default_na=False: MS-DIAL writes a literal "NA" for the trailing
        # Average/Stdev columns, and that string is meaningful here — it marks a
        # column as *not* a sample (see _is_sample_class_label).
        raw = pd.read_csv(path, sep=delim, header=None, nrows=header_row,
                          dtype=str, keep_default_na=False)
        first_row_labels = None
        for _, row in raw.iterrows():
            vals = [str(v) if not pd.isna(v) else "" for v in row.tolist()]
            key = next((v for v in vals if v.strip()), "")
            meta[f"header_{key or 'row'}"] = vals
            if first_row_labels is None:
                first_row_labels = vals
            if key.strip().lower() in ("class", "class name"):
                class_labels = vals
        # The row is always the Class row; its descriptor cell ("Class") sits
        # above a metadata column and can be missing from a trimmed export, so
        # fall back to position rather than losing the per-sample labels.
        if class_labels is None:
            class_labels = first_row_labels
    df = pd.read_csv(path, sep=delim, header=header_row, low_memory=False)
    cols = list(df.columns)
    mz = _find(cols, "Average Mz", "Average m/z", contains="Average Mz")
    rt = _find(cols, "Average Rt(min)", "Average Rt", contains="Average Rt")
    if mz is None or rt is None:
        raise ColumnResolutionError(
            f"MS-DIAL: could not find 'Average Mz'/'Average Rt' in {cols[:15]}")
    rt_unit = _rt_unit_from_name(rt, "min")

    # sample columns: prefer the class-header alignment; else everything after the
    # known MS-DIAL metadata block.
    meta_tokens = {
        "alignment id", "average rt(min)", "average mz", "metabolite name",
        "adduct type", "post curation result", "fill %", "ms/ms assigned",
        "reference rt", "reference m/z", "formula", "ontology", "inchikey",
        "smiles", "annotation tag", "comment", "manually modified",
        "isotope tracking parent id", "isotope tracking weight number",
        "total score", "rt similarity", "dot product", "reverse dot product",
        "fragment presence %", "s/n average", "spectrum reference file name",
        "ms1 isotopic spectrum", "ms/ms spectrum",
    }
    sample_cols: List[str] = []
    if class_labels is not None and len(class_labels) == len(cols):
        # The class-header row carries a leading descriptor cell (e.g. "Class")
        # that sits above the *metadata* columns (Alignment ID, Average Rt, ...),
        # so a non-empty label alone does not make a column a sample — exclude
        # any column that is itself a known MS-DIAL metadata column.
        sample_cols = [c for c, lab in zip(cols, class_labels)
                       if _is_sample_class_label(lab) and c.strip().lower() not in meta_tokens]
    if not sample_cols:
        last_meta = 0
        for i, c in enumerate(cols):
            if c.strip().lower() in meta_tokens:
                last_meta = i
        sample_cols = cols[last_meta + 1:]
    return FeatureTable(df=df, mz_col=mz, rt_col=rt, rt_unit=rt_unit,
                        sample_cols=sample_cols, source_format="msdial",
                        polarity=polarity, meta=meta)


def load_mzmine(path: str, delim: str = ",", polarity: Optional[str] = None) -> FeatureTable:
    df = pd.read_csv(path, sep=delim, low_memory=False)
    cols = list(df.columns)
    mz = _find(cols, "row m/z", "mz", "m/z")
    rt = _find(cols, "row retention time", "rt", "retention time", "row rt")
    if mz is None or rt is None:
        raise ColumnResolutionError(f"MZmine: could not find m/z / RT in {cols[:12]}")
    rt_unit = _rt_unit_from_name(rt, "min")

    # sample columns: MZmine3 'datafile:<name>:area|height'; MZmine2 '<file> Peak area|height'
    sample_cols = [c for c in cols
                   if re.search(r"(?:^|:)\s*(area|height)\s*$", c, re.I)
                   or re.search(r"peak (area|height)", c, re.I)
                   or c.lower().startswith("datafile:")]
    return FeatureTable(df=df, mz_col=mz, rt_col=rt, rt_unit=rt_unit,
                        sample_cols=sample_cols, source_format="mzmine",
                        polarity=polarity, meta={})


_ADAPTERS = {
    "masscube": load_masscube,
    "lipidscreener": load_lipidscreener,
    "msdial": load_msdial,
    "mzmine": load_mzmine,
}
