"""Typed error behaviour and the harder edge cases."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate, load_feature_table
from rt_anchor.config import CalibrationConfig
from rt_anchor.errors import (
    AnchorIdentificationError,
    ColumnResolutionError,
    InputFormatError,
    PanelError,
    RtAnchorError,
    RTUnitError,
)


# ------------------------------------------------- error taxonomy / messages --

@pytest.mark.parametrize("exc", [
    InputFormatError, ColumnResolutionError, RTUnitError, PanelError,
    AnchorIdentificationError,
])
def test_errors_subclass_base(exc):
    assert issubclass(exc, RtAnchorError)


def test_column_resolution_is_input_format_error():
    # ColumnResolutionError is a specialisation of InputFormatError
    assert issubclass(ColumnResolutionError, InputFormatError)


# ------------------------------------------------------- too-few anchors ------

def test_too_few_anchors_raises(tmp_path, pos_targets):
    # keep only 2 of the 15 standards -> below the hard minimum
    keep = pos_targets.iloc[:2]
    lines = ["m/z\tRT"]
    for _, r in keep.iterrows():
        lines.append(f"{r['mz']:.4f}\t{r['rt_ref_min']:.3f}")
    for rt in (2.0, 3.0, 4.0):          # non-anchor noise
        lines.append(f"100.0\t{rt}")
    fp = tmp_path / "few.txt"
    fp.write_text("\n".join(lines) + "\n")
    with pytest.raises(AnchorIdentificationError):
        calibrate(str(fp), "positive", standards_table=str(fp))


def test_no_anchors_header_only_raises(tmp_path):
    p = tmp_path / "hdr.txt"
    p.write_text("m/z\tRT\n")
    with pytest.raises(AnchorIdentificationError):
        calibrate(str(p), "positive", standards_table=str(p))


# ------------------------------------------------------- panel / polarity -----

def test_bad_polarity_raises_panel_error(orbitrap_samples):
    with pytest.raises(PanelError):
        calibrate(orbitrap_samples, "sideways", standards_table=orbitrap_samples)


def test_unknown_adduct_raises_panel_error():
    from rt_anchor.panel import adduct_mz
    with pytest.raises(PanelError):
        adduct_mz(500.0, "[M+K]+")


def test_all_error_messages_nonempty(orbitrap_samples, tmp_path):
    errs = []
    # bad polarity
    try:
        calibrate(orbitrap_samples, "sideways", standards_table=orbitrap_samples)
    except RtAnchorError as e:
        errs.append(e)
    # unrecognised format
    junk = tmp_path / "j.txt"
    junk.write_text("a,b\n1,2\n")
    try:
        load_feature_table(str(junk))
    except RtAnchorError as e:
        errs.append(e)
    # bad rt_unit
    try:
        load_feature_table(orbitrap_samples, rt_unit="furlongs")
    except RtAnchorError as e:
        errs.append(e)
    assert len(errs) == 3
    assert all(str(e).strip() for e in errs)


# ------------------------------------------------------------- edge inputs ----

def test_nan_rt_mz_rows_are_ignored(make_masscube):
    # inject rows with NaN m/z or RT; calibration should still succeed and
    # assign NaN RI to the NaN-RT feature (never a fabricated value).
    from pathlib import Path
    p = make_masscube()
    txt = Path(p).read_text().rstrip("\n")
    txt += "\nnan\t7.0\n100.0\tnan\n"     # NaN m/z row, NaN RT row
    Path(p).write_text(txt + "\n")
    res = calibrate(p, "positive", standards_table=p)
    t = res.table
    # the NaN-RT row must have NaN RI and NOT be flagged extrapolated
    nan_rt = t[pd.to_numeric(t["RT"], errors="coerce").isna()]
    assert len(nan_rt) == 1
    assert pd.to_numeric(nan_rt["RI"]).isna().all()
    assert (~nan_rt["is_extrapolated"].astype(bool)).all()


def test_zero_sample_column_masscube(make_masscube):
    # minimal table (m/z + RT only) -> no sample columns, still calibrates
    p = make_masscube()
    ft = load_feature_table(p)
    assert ft.sample_cols == []
    res = calibrate(p, "positive", standards_table=p)
    assert res.model["n_anchors_used"] >= 6


def test_duplicate_tied_rt_anchors(tmp_path, pos_targets):
    # place two adjacent anchors at an identical RT -> guard drops one, no crash
    lines = ["m/z\tRT"]
    for _, r in pos_targets.iterrows():
        rt = r["rt_ref_min"]
        if r["name"] in ("PC 12:0_12:0", "PC 13:0_13:0"):
            rt = 4.0
        lines.append(f"{r['mz']:.4f}\t{rt:.3f}")
    p = tmp_path / "tie.txt"
    p.write_text("\n".join(lines) + "\n")
    res = calibrate(str(p), "positive", standards_table=str(p))
    # both cannot both survive the strict-increase guard
    names = set(res.anchors["name"])
    assert not ({"PC 12:0_12:0", "PC 13:0_13:0"} <= names)


def test_negative_lipidscreener_loads(lipidscreener_neg):
    ft = load_feature_table(lipidscreener_neg, polarity="negative")
    assert ft.source_format == "lipidscreener"
    assert ft.rt_unit == "sec"
    assert ft.polarity == "negative"
    assert 15 < float(np.nanmax(ft.rt_minutes())) < 25
