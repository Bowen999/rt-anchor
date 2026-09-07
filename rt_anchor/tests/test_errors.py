"""Typed error behaviour and the harder edge cases."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate, load_feature_table
from rt_anchor.config import CalibrationConfig
from rt_anchor.errors import (
    AnchorIdentificationError,
    CalibrationError,
    ColumnResolutionError,
    ConfigError,
    InputFormatError,
    PanelError,
    RtAnchorError,
    RTUnitError,
)


# ------------------------------------------------- error taxonomy / messages --

@pytest.mark.parametrize("exc", [
    InputFormatError, ColumnResolutionError, RTUnitError, PanelError,
    AnchorIdentificationError, CalibrationError, ConfigError,
])
def test_errors_subclass_base(exc):
    assert issubclass(exc, RtAnchorError)


def test_column_resolution_is_input_format_error():
    # ColumnResolutionError is a specialisation of InputFormatError
    assert issubclass(ColumnResolutionError, InputFormatError)


# ------------------------------------------------------ §8 guard rails --------

def test_missing_standards_run_is_a_config_error(orbitrap_samples):
    """The standards run is one of the two runs the curve is built from."""
    with pytest.raises(ConfigError) as ei:
        calibrate(orbitrap_samples, "positive", standards_table="")
    assert "standards" in str(ei.value).lower()


def test_too_few_shared_pairs_names_the_likely_causes(tmp_path, pos_targets):
    # only 2 panel standards -> almost nothing matches the reference standards run
    keep = pos_targets.iloc[:2]
    lines = ["m/z\tRT"]
    for _, r in keep.iterrows():
        lines.append(f"{r['mz']:.4f}\t{r['rt_ref_min']:.3f}")
    for rt in (2.0, 3.0, 4.0):          # non-matching noise
        lines.append(f"100.0\t{rt}")
    fp = tmp_path / "few.txt"
    fp.write_text("\n".join(lines) + "\n")
    with pytest.raises(CalibrationError) as ei:
        calibrate(str(fp), "positive", standards_table=str(fp))
    msg = str(ei.value).lower()
    assert "polarity" in msg and ("reference" in msg or "mz_tol" in msg)


def test_no_features_at_all_raises(tmp_path):
    p = tmp_path / "hdr.txt"
    p.write_text("m/z\tRT\n")
    with pytest.raises(CalibrationError):
        calibrate(str(p), "positive", standards_table=str(p))


def test_negative_run_against_the_positive_reference_loses_the_irt_ruler(
        lipidscreener_neg):
    """A negative-mode run against the positive-mode bundled reference.

    Honest about a real limit: enough m/z collide by chance (196 pairs here) that
    the <5-pairs guard rail does NOT fire, so the mismatch is not refused. What
    does happen is that the negative-mode panel adducts are absent from the
    positive reference standards run, fewer than 2 landmarks are found, and §4
    takes over — iRT is NaN everywhere with reliability "none", while Cal_RT is
    still produced. The evidence is in the model: the curve's median absolute
    pair residual is minutes rather than the ~0.1 min a genuine cross-column fit
    gives.
    """
    res = calibrate(lipidscreener_neg, "negative", standards_table=lipidscreener_neg,
                    panel="mix21")
    assert res.irt is None
    assert res.values("iRT").isna().all()
    assert (res.table[res.col("iRT_reliability")] == "none").all()
    assert res.values("Cal_RT_min").notna().any()      # Cal_RT is independent
    assert "landmark" in res.model["irt"]["reason"]
    assert res.model["curve"]["residual_min"]["median_abs"] > 0.5


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
    # the NaN-RT row must have NaN Cal_RT and NOT be flagged extrapolated
    # (out-of-span is a claim about where the feature eluted; we do not know)
    nan_rt = t[pd.to_numeric(t["RT"], errors="coerce").isna()]
    assert len(nan_rt) == 1
    for c in ("Cal_RT_min", "Cal_RT_uncertainty_min", "iRT", "iRT_uncertainty"):
        assert pd.to_numeric(nan_rt[c]).isna().all(), c
    assert (nan_rt["iRT_reliability"] == "none").all()
    assert (~nan_rt["is_extrapolated"].astype(bool)).all()


def test_zero_sample_column_masscube(make_masscube):
    # minimal table (m/z + RT only) -> no sample columns, still calibrates
    p = make_masscube()
    ft = load_feature_table(p)
    assert ft.sample_cols == []
    res = calibrate(p, "positive", standards_table=p, panel="mix15")
    assert res.model["curve"]["n_pairs"] >= 5
    assert res.values("Cal_RT_min").notna().any()


def test_duplicate_tied_rt_standards_in_detection_qc(tmp_path, pos_targets):
    # two panel standards at an identical RT -> the LIS guard drops one from the
    # detection-QC template, and the calibration itself is unaffected
    lines = ["m/z\tRT"]
    for _, r in pos_targets.iterrows():
        rt = r["rt_ref_min"]
        if r["name"] in ("PC 12:0_12:0", "PC 13:0_13:0"):
            rt = 4.0
        lines.append(f"{r['mz']:.4f}\t{rt:.3f}")
    p = tmp_path / "tie.txt"
    p.write_text("\n".join(lines) + "\n")
    res = calibrate(str(p), "positive", standards_table=str(p), panel="mix15")
    detected = set(res.panel["native_rt"])
    assert not ({"PC 12:0_12:0", "PC 13:0_13:0"} <= detected)
    assert res.values("Cal_RT_min").notna().any()


def test_negative_lipidscreener_loads(lipidscreener_neg):
    ft = load_feature_table(lipidscreener_neg, polarity="negative")
    assert ft.source_format == "lipidscreener"
    assert ft.rt_unit == "sec"
    assert ft.polarity == "negative"
    assert 15 < float(np.nanmax(ft.rt_minutes())) < 25
