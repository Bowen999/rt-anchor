"""Pipeline behaviour: the two tiers, the output columns, the model JSON, the
panel choices, the extrapolation flags, and the two ways stage 2 can be absent.

Engine-level units live in ``test_crosscolumn.py``; reference-set and panel
resolution in ``test_reference.py``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate
from rt_anchor.config import REMOVED_FIELDS, CalibrationConfig
from rt_anchor.errors import ConfigError
from rt_anchor.identify import _drop_nonmonotone, _longest_increasing
from rt_anchor.io.schema import RESULT_COLUMNS

RELIABILITY = {"high", "medium", "low", "none"}


@pytest.fixture(scope="module")
def proj_result(orbitrap_samples, orbitrap_standards):
    """A real per-project calibration against the bundled 35-min reference."""
    return calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                     panel="mix15", config=CalibrationConfig.orbitrap())


@pytest.fixture(scope="module")
def sample_result(qtof_full_samples, qtof_full_standards, qtof_full_single_files):
    return calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                     panel="mix15", single_files=list(qtof_full_single_files))


# ------------------------------------------------------- longest-increasing ---
# (still used by the detection-QC / landmark path)

@pytest.mark.parametrize("vals,expected_len", [
    ([1, 2, 3, 4, 5], 5),
    ([5, 4, 3, 2, 1], 1),
    ([1, 3, 2, 4, 5], 4),
    ([1, 2, 2, 3], 3),          # strict -> one of the tied 2s dropped
    ([], 0),
])
def test_longest_increasing_length(vals, expected_len):
    idx = _longest_increasing(vals)
    assert len(idx) == expected_len
    assert idx == sorted(idx)
    picked = [vals[i] for i in idx]
    assert all(b > a for a, b in zip(picked, picked[1:]))


def test_drop_nonmonotone_flags_count():
    df = pd.DataFrame({
        "name": list("abcde"),
        "rt_ref_min": [1, 2, 3, 4, 5],
        "rt_obs_min": [1.0, 2.0, 9.0, 3.0, 4.0],   # the 9.0 breaks the run
    })
    out = _drop_nonmonotone(df, CalibrationConfig())
    assert out.attrs["n_dropped_nonmonotone"] == 1
    assert "c" not in set(out["name"])
    assert list(out["rt_obs_min"]) == [1.0, 2.0, 3.0, 4.0]


# ------------------------------------------------------------ the outputs ----

def test_appends_exactly_the_spec_columns(proj_result):
    for c in RESULT_COLUMNS:
        assert proj_result.col(c) in proj_result.table.columns
    assert RESULT_COLUMNS[:3] == ["Cal_RT_min", "Cal_RT_uncertainty_min", "iRT"]


def test_cal_rt_is_monotone_in_rt_and_lands_on_the_reference_axis(proj_result):
    rt = proj_result.rt_minutes().to_numpy()
    # the stage-1 curve itself is the monotone object; a gated-on anchor
    # correction may bend it locally, so monotonicity is asserted on the curve.
    ok = np.isfinite(rt)
    pred = proj_result.curve.predict(np.sort(rt[ok]))
    assert np.all(np.diff(pred) >= -1e-9)
    cal = proj_result.values("Cal_RT_min")
    # the reference column runs ~1.4-24 min: the calibrated RTs must live there,
    # not on the source column's own axis
    assert 0.0 < np.nanmedian(cal) < 30.0


def test_irt_tracks_cal_rt_affinely(proj_result):
    cal = proj_result.values("Cal_RT_min").to_numpy()
    irt = proj_result.values("iRT").to_numpy()
    ok = np.isfinite(cal) & np.isfinite(irt)
    m = proj_result.irt
    assert m is not None
    assert np.allclose(irt[ok], m.to_irt(cal[ok]))
    # affine: a straight line through the landmark endpoints reproduces it
    lo, hi = m.rt_span
    straight = 1.0 + 99.0 * (cal[ok] - lo) / (hi - lo)
    assert np.max(np.abs(irt[ok] - straight)) < 1e-9


def test_uncertainty_and_reliability_populated(proj_result):
    t = proj_result.table
    cal = proj_result.values("Cal_RT_min")
    unc = proj_result.values("Cal_RT_uncertainty_min")
    iunc = proj_result.values("iRT_uncertainty")
    assert unc[cal.notna()].notna().all()
    assert (unc.dropna() >= 0).all()
    # iRT sigma is the RT sigma carried through the (constant) scale slope
    slope = 99.0 / (proj_result.irt.rt_span[1] - proj_result.irt.rt_span[0])
    ok = cal.notna()
    assert np.allclose(iunc[ok], unc[ok] * slope)
    rel = set(t[proj_result.col("iRT_reliability")].astype(str))
    assert rel <= RELIABILITY


def test_extrapolation_is_flagged_and_still_valued(qtof_samples, qtof_standards):
    res = calibrate(qtof_samples, "positive", standards_table=qtof_standards,
                    panel="mix15")
    ex = res.table[res.col("is_extrapolated")].astype(bool).to_numpy()
    cal = res.values("Cal_RT_min").to_numpy()
    rt = res.rt_minutes().to_numpy()
    x0, x1 = res.curve.x0, res.curve.x1
    assert ex.sum() > 0, "the example run must reach outside the matched-pair span"
    # the flag says exactly "outside the pair span" (NaN RT is not out-of-span)
    assert np.array_equal(ex, ((rt < x0) | (rt > x1)) & np.isfinite(rt))
    # v2 assigns a value anyway (extrapolate defaults to True) ...
    assert np.isfinite(cal[ex]).all()
    # ... and never calls it reliable
    rel = res.table[res.col("iRT_reliability")].to_numpy()[ex]
    assert set(rel) <= {"low", "none"}


def test_per_project_scope_and_qc_columns(proj_result):
    t = proj_result.table
    assert (t[proj_result.col("calibration_scope")] == "project").all()
    assert proj_result.values("RI_spread").isna().all()
    assert (proj_result.values("n_contributing") == 1).all()
    assert set(t[proj_result.col("warp_source")]) <= {"curve", "curve+anchors"}


# --------------------------------------------------------------- two tiers ----

def test_per_sample_tier(sample_result):
    t = sample_result.table
    assert sample_result.model["calibration_scope"] == "sample"
    assert (t[sample_result.col("calibration_scope")] == "sample").all()
    assert sample_result.values("RI_spread").notna().any()
    assert sample_result.values("n_contributing").max() > 1
    assert sample_result.model.get("n_injections_used", 0) >= 2


def test_project_tier_when_no_single_files(qtof_full_samples, qtof_full_standards):
    res = calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                    panel="mix15")
    assert res.model["calibration_scope"] == "project"
    assert res.values("RI_spread").isna().all()


# --------------------------------------------------------- model json (§5.1) --

def test_model_json_structure(proj_result):
    m = proj_result.model
    assert m["method"] == "cross-column-v2"
    assert m["calibration_scope"] == "project"
    for key in ("reference", "panel", "curve", "anchors", "irt", "detection_qc",
                "n_features", "n_features_extrapolated", "confidence_counts", "config"):
        assert key in m, f"model.json is missing '{key}'"
    assert m["reference"]["key"] == "col35" and m["reference"]["is_default"]
    c = m["curve"]
    assert c["n_pairs"] == c["n_pairs_standards"] + c["n_pairs_sample"]
    assert 0 < c["n_pairs_kept"] <= c["n_pairs"]
    assert c["rt_span_src_min"][0] < c["rt_span_src_min"][1]
    assert c["residual_min"]["median_abs"] <= c["residual_min"]["p90_abs"]
    a = m["anchors"]
    assert a["gate_threshold"] == 0.2 and isinstance(a["table"], list)
    assert m["irt"]["definition"] == "affine on the reference-column RT axis"
    qc = m["detection_qc"]["user_standards_run"]
    assert qc["n_detected"] <= qc["n_panel"] == 15
    assert "does not drive the calibration" in qc["note"]
    assert sum(m["confidence_counts"].values()) == m["n_features"]


def test_model_records_the_gate_decision(proj_result):
    a = proj_result.model["anchors"]
    assert a["engaged"] == proj_result.anchors_used
    assert a["gate_reason"]
    if a["engaged"]:
        assert a["gate_mse_reduction"] >= a["gate_threshold"]


def test_companion_frames_have_their_spec_columns(proj_result):
    assert list(proj_result.pairs.columns) == ["mz_src", "rt_src", "rt_ref",
                                               "source", "kept"]
    assert set(proj_result.pairs["source"]) <= {"standards", "sample"}
    assert set(proj_result.landmarks.columns) >= {"name", "class", "mz",
                                                  "rt_ref_run_min", "iRT"}
    assert set(proj_result.anchors.columns) >= {
        "label", "lipid_class", "rt_src", "rt_ref", "residual_min",
        "n_isomer_candidates", "isomer_rts", "pick_refined",
        "dropped_by_sanity_filter"}


# ------------------------------------------------------------ panel choices ---

def test_panel_none_still_calibrates_and_falls_back_for_irt(
        orbitrap_samples, orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="none", config=CalibrationConfig.orbitrap())
    assert res.panel_key == "none"
    # calibration itself is unaffected — stage 1 claims no identities
    assert res.values("Cal_RT_min").notna().any()
    # the ruler falls back to irt_landmark_panel, and says so
    assert "fallback" in res.model["irt"]["panel_used"]
    assert res.values("iRT").notna().any()
    # detection QC is skipped, not faked
    assert res.model["detection_qc"]["user_standards_run"]["n_panel"] == 0


def test_panel_choice_moves_the_irt_scale(orbitrap_samples, orbitrap_standards):
    cfg = CalibrationConfig.orbitrap()
    a = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                  panel="mix15", config=cfg)
    b = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                  panel="mix21", config=cfg)
    assert a.irt.rt_span != b.irt.rt_span, (
        "mix15 and mix21 have different landmark endpoints on the reference run, "
        "so iRT is comparable within a panel choice only (spec §10.2)")


def test_unknown_panel_raises_config_error(orbitrap_samples, orbitrap_standards):
    with pytest.raises(ConfigError) as ei:
        calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                  panel="mix99")
    assert "mix99" in str(ei.value)


# ------------------------------------------ stage 2 absent, on purpose --------

def test_non_plasma_matrix_falls_back_to_stage_one_cleanly(non_plasma_table):
    """The archaeal-lipid case: no anchors is a result, not an error (spec §2)."""
    res = calibrate(non_plasma_table, "positive", standards_table=non_plasma_table,
                    panel="mix21")
    assert res.values("Cal_RT_min").notna().any()
    assert not res.anchors_used
    assert (res.table[res.col("warp_source")] == "curve").all()
    assert res.model["anchors"]["engaged"] is False
    assert res.model["anchors"]["gate_reason"]
    assert any("stage 2" in line for line in res.log)


def test_use_sample_anchors_false_skips_stage_two(orbitrap_samples, orbitrap_standards):
    cfg = CalibrationConfig.orbitrap(use_sample_anchors=False)
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=cfg)
    assert not res.anchors_used
    assert res.model["anchors"]["gate_reason"] == "anchors disabled"


def test_use_sample_pairs_false_shrinks_the_curve_evidence(
        orbitrap_samples, orbitrap_standards):
    cfg = CalibrationConfig.orbitrap(use_sample_pairs=False)
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=cfg)
    assert res.model["curve"]["n_pairs_sample"] == 0
    assert res.model["curve"]["n_pairs_standards"] == res.model["curve"]["n_pairs"]


# ----------------------------------------------------------- config presets ---

def test_config_presets_differ():
    q, o = CalibrationConfig.qtof(), CalibrationConfig.orbitrap()
    assert q.mz_tol_ppm == 15.0 and o.mz_tol_ppm == 8.0
    assert q.rt_window_min == 0.5 and o.rt_window_min == 0.3
    # the anonymous window is a method constant, not an instrument setting
    assert q.match_mz_tol_da == o.match_mz_tol_da == 0.008


def test_config_defaults_are_the_validated_ones():
    c = CalibrationConfig()
    assert c.curve_frac == 0.1 and c.curve_iter == 3 and c.curve_mad_k == 3.0
    assert c.anchor_gate_min_mse_reduction == 0.2 and c.min_anchors == 3
    assert c.extrapolate is True and c.extrapolate_mode == "linear"
    assert c.irt_landmark_panel == "mix21"


def test_config_from_dict_rejects_unknown_key():
    with pytest.raises(ConfigError) as ei:
        CalibrationConfig.from_dict({"mz_tol_ppm": 5.0, "bogus_key": 1})
    assert "bogus_key" in str(ei.value)


@pytest.mark.parametrize("removed", sorted(REMOVED_FIELDS))
def test_config_from_dict_names_the_v2_replacement(removed):
    with pytest.raises(ConfigError) as ei:
        CalibrationConfig.from_dict({removed: 1.0})
    msg = str(ei.value)
    assert removed in msg and len(msg) > len(removed) + 40, (
        "a removed key must be refused with the sentence saying what replaced it")


def test_config_roundtrips_through_dict():
    c = CalibrationConfig.orbitrap(curve_frac=0.2, use_sample_anchors=False)
    back = CalibrationConfig.from_dict(c.to_dict())
    assert back == c


def test_config_overrides_propagate_into_model(orbitrap_samples, orbitrap_standards):
    cfg = CalibrationConfig.orbitrap(mz_tol_ppm=6.0, curve_frac=0.2)
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=cfg)
    assert res.model["config"]["mz_tol_ppm"] == 6.0
    assert res.model["config"]["curve_frac"] == 0.2
