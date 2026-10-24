"""The guards against a calibration that looks plausible and is meaningless.

Found on negative-mode runs calibrated against the positive-mode bundled
reference: every m/z pair was two different molecules, and the run still
produced a curve, a backwards line when few pairs survived, and negative
Cal_RT values from extrapolation. Each guard is tested where it lives.
"""

from __future__ import annotations

import shutil

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate
from rt_anchor.config import CalibrationConfig
from rt_anchor.crosscolumn import MonotoneCurve, fit_robust_curve
from rt_anchor.errors import CalibrationError, ConfigError, PolarityMismatchError
from rt_anchor.io import FeatureTable
from rt_anchor.pipeline import _assemble, _same_file


# ------------------------------------------------------------------ polarity --

def test_polarity_mismatch_is_a_config_error():
    assert issubclass(PolarityMismatchError, ConfigError)


@pytest.mark.parametrize("labels,expected", [
    (["[M+H]+"] * 9 + ["[M+Na]+"] * 3, "positive"),
    (["[M-H]-"] * 10 + ["[M+HCOO]-"] * 2 + ["[M+H]+"], "negative"),
    (["[M+2H]2+"] * 12, "positive"),
    (["[M+H]+"] * 6 + ["[M-H]-"] * 6, None),        # no clear majority
    (["[M+H]+"] * 5, None),                          # too few labels to say
    (["Unknown", "null", ""] * 10, None),            # no signed labels at all
])
def test_adduct_polarity(labels, expected):
    df = pd.DataFrame({"mz": 500.0, "rt": 1.0, "Adduct type": labels})
    ft = FeatureTable(df=df, mz_col="mz", rt_col="rt", rt_unit="min")
    assert ft.adduct_polarity() == expected


def test_adduct_polarity_without_an_adduct_column():
    ft = FeatureTable(df=pd.DataFrame({"mz": [500.0], "rt": [1.0]}),
                      mz_col="mz", rt_col="rt", rt_unit="min")
    assert ft.adduct_polarity() is None


def test_negative_run_against_the_bundled_positive_reference_is_refused(
        orbitrap_samples, orbitrap_standards):
    with pytest.raises(PolarityMismatchError) as ei:
        calibrate(orbitrap_samples, "negative", standards_table=orbitrap_standards,
                  panel="mix21")
    msg = str(ei.value)
    assert "positive mode" in msg and "negative" in msg
    assert "Reference sample" in msg and "--reference-sample" in msg
    # the desktop app titles any message holding the keyword names
    # "Reference pair incomplete", which this is not
    assert "reference_sample" not in msg and "reference_standards" not in msg


def test_own_runs_labelled_in_the_other_polarity_are_refused(
        bundled_reference, orbitrap_standards):
    # the bundled reference sample is an MS-DIAL export full of [M+H]+ labels
    with pytest.raises(PolarityMismatchError, match="your own runs"):
        calibrate(bundled_reference.sample_path, "negative",
                  standards_table=orbitrap_standards)


def test_custom_reference_labelled_in_the_other_polarity_is_refused(
        orbitrap_samples, orbitrap_standards, bundled_reference):
    # the bundled files passed as a custom pair: no provenance, but adduct labels
    with pytest.raises(PolarityMismatchError, match="adduct labels"):
        calibrate(orbitrap_samples, "negative", standards_table=orbitrap_standards,
                  panel="mix21",
                  reference_sample=bundled_reference.sample_path,
                  reference_standards=bundled_reference.standards_path)


def _shift_mz(src, dst, ppm):
    df = pd.read_csv(src, sep="\t")
    df["m/z"] = df["m/z"] * (1.0 + ppm * 1e-6)
    df.to_csv(dst, sep="\t", index=False)
    return str(dst)


@pytest.fixture
def reference_without_the_panel(tmp_path, qtof_samples, qtof_standards):
    """A MassCube pair (no adduct labels, no provenance) whose m/z sit 500 ppm off."""
    return (_shift_mz(qtof_samples, tmp_path / "ref_samples.txt", 500.0),
            _shift_mz(qtof_standards, tmp_path / "ref_standards.txt", 500.0))


def test_reference_holding_none_of_the_panel_is_refused(
        orbitrap_samples, orbitrap_standards, reference_without_the_panel):
    ref_s, ref_t = reference_without_the_panel
    with pytest.raises(PolarityMismatchError, match="panel standards"):
        calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                  panel="mix15", config=CalibrationConfig.orbitrap(),
                  reference_sample=ref_s, reference_standards=ref_t)


def test_panel_none_skips_the_panel_check(
        orbitrap_samples, orbitrap_standards, reference_without_the_panel):
    """No panel claimed, nothing to cross-check: the run goes through.

    The known limit of the guards: two large tables collide by chance often
    enough that the pair-count guard rail does not fire either, so a reference
    without labels, provenance or panel cannot be told apart from a good one.
    """
    ref_s, ref_t = reference_without_the_panel
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="none", config=CalibrationConfig.orbitrap(),
                    reference_sample=ref_s, reference_standards=ref_t)
    assert res.values("Cal_RT_min").notna().any()


# ------------------------------------------------------ a curve that runs back --

def test_line_fallback_refuses_a_backwards_map():
    # few pairs, falling: what 18 false pairs on the 12-min Fiehn run produced
    x = np.linspace(3.7, 7.7, 12)
    with pytest.raises(CalibrationError, match="backwards"):
        fit_robust_curve(x, 20.0 - 0.4 * x, min_points=20)


def test_line_fallback_refuses_a_flat_map():
    x = np.linspace(3.7, 7.7, 12)
    with pytest.raises(CalibrationError, match="backwards"):
        fit_robust_curve(x, np.full_like(x, 12.0), min_points=20)


def test_curve_refuses_decreasing_knots():
    with pytest.raises(CalibrationError, match="must not decrease"):
        MonotoneCurve([1.0, 2.0, 3.0], [3.0, 2.0, 1.0])


def test_curve_accepts_tied_knots():
    # isotonic regression pools violators into ties; that is still monotone
    curve = MonotoneCurve([1.0, 2.0, 3.0], [2.0, 2.0, 3.0])
    assert curve.predict(np.array([2.5]))[0] == pytest.approx(2.5, abs=0.5)


# ------------------------------------------------------- a negative Cal_RT ----

def test_assemble_refuses_a_negative_cal_rt():
    rt = np.array([0.1, 1.0, 2.0, np.nan])
    out = _assemble(rt, np.array([-3.0, 0.5, 2.0, 1.0]), np.full(4, 0.1),
                    np.array([True, False, False, False]),
                    np.array(["curve"] * 4, dtype=object), "project",
                    CalibrationConfig(), None)
    assert np.isnan(out["Cal_RT_min"][0])
    assert np.isnan(out["Cal_RT_uncertainty_min"][0])
    assert bool(out["is_extrapolated"][0])
    assert out["iRT_reliability"][0] == "none"
    assert out["Cal_RT_min"][1] == pytest.approx(0.5)
    assert out.attrs["n_negative_cal_rt"] == 1          # the NaN-RT row is not counted


def test_refused_rows_are_logged(qtof_samples, qtof_standards):
    res = calibrate(qtof_samples, "positive", standards_table=qtof_standards,
                    panel="mix15")
    cal = res.values("Cal_RT_min")
    assert (cal.dropna() >= 0).all()
    assert any("negative Cal_RT" in line for line in res.log)


# --------------------------------------------------- a run as its own reference --

def test_same_file(tmp_path, orbitrap_standards):
    copy = tmp_path / "copy.txt"
    shutil.copy(orbitrap_standards, copy)
    other = tmp_path / "other.txt"
    other.write_text("m/z\tRT\n", encoding="utf-8")
    assert _same_file(orbitrap_standards, orbitrap_standards)
    assert _same_file(orbitrap_standards, str(copy))          # a copy is the same run
    assert not _same_file(orbitrap_standards, str(other))
    assert not _same_file(orbitrap_standards, str(tmp_path / "missing.txt"))


def test_a_run_as_its_own_reference_is_warned_about(orbitrap_samples, orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=CalibrationConfig.orbitrap(),
                    reference_sample=orbitrap_samples,
                    reference_standards=orbitrap_standards)
    assert len(res.model["warnings"]) == 2
    assert sum(line.startswith("WARNING:") for line in res.log) == 2
    assert all("itself" in w for w in res.model["warnings"])


def test_an_ordinary_run_carries_no_warning(orbitrap_samples, orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=CalibrationConfig.orbitrap())
    assert res.model["warnings"] == []
