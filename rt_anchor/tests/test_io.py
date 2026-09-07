"""I/O layer: format detection, column resolution, RT-unit handling,
column preservation."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate, describe_input, load_feature_table, which_format
from rt_anchor.errors import ColumnResolutionError, InputFormatError, RTUnitError
from rt_anchor.io import detect_format
from rt_anchor.io.schema import RESULT_COLUMNS


# --------------------------------------------------------- format detection ---

def test_detect_masscube(orbitrap_samples):
    assert which_format(orbitrap_samples) == "masscube"
    fmt, delim, hdr = detect_format(orbitrap_samples)
    assert fmt == "masscube" and delim == "\t" and hdr == 0


def test_detect_lipidscreener_real(lipidscreener_pos):
    assert which_format(lipidscreener_pos) == "lipidscreener"


def test_detect_msdial_synthetic(make_msdial):
    p = make_msdial()
    fmt, delim, hdr = detect_format(p)
    assert fmt == "msdial"
    assert hdr > 0, "MS-DIAL header should not be on line 0 (metadata rows above)"


def test_detect_mzmine2_synthetic(make_mzmine2):
    assert which_format(make_mzmine2()) == "mzmine"


def test_detect_mzmine3_synthetic(make_mzmine3):
    assert which_format(make_mzmine3()) == "mzmine"


def test_detect_unrecognisable_raises(tmp_path):
    p = tmp_path / "junk.txt"
    p.write_text("foo,bar,baz\n1,2,3\n")
    with pytest.raises(InputFormatError):
        which_format(str(p))


def test_detect_empty_raises(tmp_path):
    p = tmp_path / "empty.txt"
    p.write_text("")
    with pytest.raises(InputFormatError):
        load_feature_table(str(p))


# ------------------------------------------------------- column resolution ----

def test_missing_mz_rt_raises_column_resolution(tmp_path):
    # A masscube-looking file forced as masscube but without m/z / RT columns.
    p = tmp_path / "nocols.txt"
    p.write_text("alpha\tbeta\n1\t2\n")
    with pytest.raises(ColumnResolutionError):
        load_feature_table(str(p), source_format="masscube")


def test_bad_source_format_raises(orbitrap_samples):
    with pytest.raises(InputFormatError):
        load_feature_table(orbitrap_samples, source_format="waters")


# ----------------------------------------------------------- RT unit logic ----

def test_lipidscreener_parsed_as_seconds(lipidscreener_pos):
    ft = load_feature_table(lipidscreener_pos)
    assert ft.source_format == "lipidscreener"
    assert ft.rt_unit == "sec"
    mx = float(np.nanmax(ft.rt_minutes()))
    assert 15 < mx < 25, f"expected ~21 min after sec->min conversion, got {mx}"
    # raw (seconds) is ~1260, not minutes
    assert float(np.nanmax(ft.rt_raw())) > 1000


def test_lipidscreener_group_row_dropped(lipidscreener_pos):
    ft = load_feature_table(lipidscreener_pos)
    assert "group_row" in ft.meta
    # the "Group" label must not survive as a feature row
    first_no = str(ft.df.iloc[0][ft.df.columns[0]]).strip().lower()
    assert first_no != "group"


def test_rt_unit_override(orbitrap_samples):
    ft_min = load_feature_table(orbitrap_samples)             # native minutes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # forcing 'sec' on a minutes table warns (expected)
        ft_sec = load_feature_table(orbitrap_samples, rt_unit="sec")
    # forcing 'sec' should divide by 60
    assert np.allclose(ft_sec.rt_minutes().dropna().to_numpy(),
                       (ft_min.rt_raw().dropna().to_numpy() / 60.0))


def test_bad_rt_unit_raises(orbitrap_samples):
    with pytest.raises(RTUnitError):
        load_feature_table(orbitrap_samples, rt_unit="hours")


def test_sanity_warning_fires_for_seconds_in_minutes_table(make_masscube):
    # MassCube declared minutes but values are actually seconds (x60) -> warn.
    p = make_masscube(rt_scale=60.0)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        ft = load_feature_table(p)
    assert ft.rt_unit == "min"
    assert any("looks" in str(w.message).lower() for w in rec), \
        "expected a range-sanity warning that a minutes table holds seconds"


def test_no_sanity_warning_for_normal_minutes(make_masscube):
    p = make_masscube(rt_scale=1.0)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        load_feature_table(p)
    assert not any("looks" in str(w.message).lower() for w in rec)


# --------------------------------------------- sample-column classification ---

def test_msdial_sample_cols_exclude_metadata(make_msdial):
    """Alignment ID / Average Rt / Average Mz must NOT be treated as samples."""
    ft = load_feature_table(make_msdial())
    for meta_col in ["Alignment ID", "Average Rt(min)", "Average Mz",
                     "Metabolite name", "Adduct type"]:
        assert meta_col not in ft.sample_cols, \
            f"{meta_col!r} was wrongly classified as a sample column"
    assert set(ft.sample_cols) == {"S1", "S2"}


def test_mzmine2_sample_cols(make_mzmine2):
    ft = load_feature_table(make_mzmine2())
    assert set(ft.sample_cols) == {"SampleA.mzML Peak area", "SampleB.mzML Peak area"}


def test_mzmine3_sample_cols(make_mzmine3):
    ft = load_feature_table(make_mzmine3())
    assert set(ft.sample_cols) == {"datafile:SampleA.mzML:area",
                                   "datafile:SampleB.mzML:area"}


# ----------------------------------------------------- describe_input helper --

def test_describe_input_masscube(orbitrap_samples):
    info = describe_input(orbitrap_samples)
    assert info["source_format"] == "masscube"
    assert info["n_features"] > 0
    assert info["rt_unit"] == "min"
    assert info["path"] == orbitrap_samples


# ------------------------------------------- column preservation (CRITICAL) ---

def test_calibrate_only_appends_columns(qtof_full_samples, qtof_full_standards):
    orig = pd.read_csv(qtof_full_samples, sep="\t", low_memory=False)
    res = calibrate(qtof_full_samples, "positive",
                    standards_table=qtof_full_standards, panel="mix15")
    out = res.table
    # every original column present with identical values, same order preserved
    assert list(out.columns[:len(orig.columns)]) == list(orig.columns)
    assert len(out) == len(orig)
    for c in orig.columns:
        pd.testing.assert_series_equal(
            out[c].reset_index(drop=True), orig[c].reset_index(drop=True),
            check_names=False)
    appended = [c for c in out.columns if c not in orig.columns]
    assert set(appended) == set(RESULT_COLUMNS)


def test_collision_rename_preserves_existing_column(make_masscube):
    # input already carries an 'iRT' column -> the package must keep it and
    # append its own result under 'rtanchor_iRT'.
    p = make_masscube(extra_cols={"iRT": 42})
    res = calibrate(p, "positive", standards_table=p, panel="mix15")
    out = res.table
    assert "iRT" in out.columns and "rtanchor_iRT" in out.columns
    assert res.col("iRT") == "rtanchor_iRT"
    # original values untouched
    assert set(pd.to_numeric(out["iRT"]).unique()) == {42}
    # the computed iRT lives under the renamed column and differs from 42
    assert not np.allclose(pd.to_numeric(out["rtanchor_iRT"]).dropna().to_numpy(), 42)
