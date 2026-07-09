"""Calibration correctness: monotone warp, extrapolation, two tiers,
uncertainty, model json, config presets, and the monotonicity guard."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate
from rt_anchor.calibrate import apply_warp, fit_warp
from rt_anchor.config import CalibrationConfig
from rt_anchor.errors import CalibrationError
from rt_anchor.identify import _drop_nonmonotone, _longest_increasing


# ------------------------------------------------------- longest-increasing ---

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
    # returned indices must be increasing and select a strictly increasing seq
    assert idx == sorted(idx)
    picked = [vals[i] for i in idx]
    assert all(b > a for a, b in zip(picked, picked[1:]))


def test_drop_nonmonotone_flags_count():
    # rt_ref ascending; rt_obs has one out-of-order point that must be dropped
    df = pd.DataFrame({
        "name": list("abcde"),
        "rt_ref_min": [1, 2, 3, 4, 5],
        "rt_obs_min": [1.0, 2.0, 9.0, 3.0, 4.0],   # the 9.0 breaks the run
    })
    out = _drop_nonmonotone(df, CalibrationConfig())
    assert out.attrs["n_dropped_nonmonotone"] == 1
    assert "c" not in set(out["name"])            # the 9.0 anchor (name 'c')
    assert list(out["rt_obs_min"]) == [1.0, 2.0, 3.0, 4.0]


# ----------------------------------------------------- per-project baseline ---

def test_per_project_monotone_and_extrapolation(make_masscube):
    # anchors at their reference RTs + in-span probes + two beyond-span probes
    extra = [(100.0, 0.2), (100.0, 25.0), (100.0, 6.0), (100.0, 9.0)]
    p = make_masscube(extra_rows=extra)
    res = calibrate(p, "positive", standards_table=p)
    t = res.table.copy()
    t["RTn"] = pd.to_numeric(t["RT"])
    lo, hi = res.model["anchor_span_min"]

    # in-span RI increases monotonically with RT
    inspan = t[(t["RTn"] >= lo) & (t["RTn"] <= hi)].sort_values("RTn")
    ri = pd.to_numeric(inspan["RI"]).to_numpy()
    ri = ri[np.isfinite(ri)]
    assert np.all(np.diff(ri) >= -1e-9)

    # beyond span -> NaN RI (never fabricated) + is_extrapolated True
    beyond = t[(t["RTn"] < lo) | (t["RTn"] > hi)]
    assert pd.to_numeric(beyond["RI"]).isna().all()
    assert beyond["is_extrapolated"].astype(bool).all()

    # within span -> not extrapolated
    assert (~inspan["is_extrapolated"].astype(bool)).all()

    # scope + spread columns for per-project
    assert (t["calibration_scope"] == "project").all()
    assert pd.to_numeric(t["RI_spread"]).isna().all()
    assert (pd.to_numeric(t["n_contributing"]) == 1).all()


def test_uncertainty_and_reliability_populated(orbitrap_samples, orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive",
                    standards_table=orbitrap_standards,
                    config=CalibrationConfig.orbitrap())
    t = res.table
    ri = pd.to_numeric(t["RI"])
    rt = pd.to_numeric(t["RT"], errors="coerce")
    unc = pd.to_numeric(t["RI_uncertainty"])
    # uncertainty is populated wherever RI is finite (the meaningful guarantee)
    assert unc[ri.notna()].notna().all()
    # uncertainty is NaN wherever RT itself is NaN (no basis for an estimate)
    if rt.isna().any():
        assert unc[rt.isna()].isna().all()
    # reliability tier vocabulary; high present for the PC-dense core
    rel = set(t["RI_reliability"].astype(str))
    assert rel <= {"high", "medium", "low", "none"}
    assert "high" in rel


def test_extrapolate_flag_fills_beyond_span(make_masscube):
    extra = [(100.0, 25.0)]     # beyond max anchor
    p = make_masscube(extra_rows=extra)
    off = calibrate(p, "positive", standards_table=p, config=CalibrationConfig(extrapolate=False))
    on = calibrate(p, "positive", standards_table=p,
                   config=CalibrationConfig(extrapolate=True, max_extrapolation_min=1.0))
    t_off, t_on = off.table.copy(), on.table.copy()
    t_off["RTn"] = pd.to_numeric(t_off["RT"])
    t_on["RTn"] = pd.to_numeric(t_on["RT"])
    hi = off.model["anchor_span_min"][1]
    row_off = t_off[t_off["RTn"] > hi]
    row_on = t_on[t_on["RTn"] > hi]
    assert pd.to_numeric(row_off["RI"]).isna().all()
    assert pd.to_numeric(row_on["RI"]).notna().all()   # capped extrapolation fills it


# ------------------------------------------------------- model json sanity ----

def test_model_json_structure(orbitrap_samples, orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive",
                    standards_table=orbitrap_standards,
                    config=CalibrationConfig.orbitrap())
    m = res.model
    assert m["calibration_scope"] == "project"
    assert m["n_anchors_used"] >= 6
    lo, hi = m["anchor_span_min"]
    assert lo < hi
    assert m["loo_residual_irt"]["median"] is not None
    assert m["scale_constants"] == {"RT_lo_min": 1.0, "RT_hi_min": 18.4}
    # anchors present + sorted by reference RT
    anchors = pd.DataFrame(m["anchors"])
    assert (anchors["rt_ref_min"].to_numpy() ==
            np.sort(anchors["rt_ref_min"].to_numpy())).all()


# --------------------------------------------------------------- two tiers ----

def test_per_sample_tier(qtof_full_samples, qtof_full_standards,
                         qtof_full_single_files):
    res = calibrate(qtof_full_samples, "positive",
                    standards_table=qtof_full_standards,
                    single_files=qtof_full_single_files)
    t = res.table
    assert res.model["calibration_scope"] == "sample"
    assert (t["calibration_scope"] == "sample").all()
    # per-sample QC columns populated
    assert pd.to_numeric(t["RI_spread"]).notna().any()
    n_contrib = pd.to_numeric(t["n_contributing"])
    assert n_contrib.max() > 1
    assert res.model.get("n_injections_used", 0) >= 2


def test_project_tier_when_no_single_files(qtof_full_samples, qtof_full_standards):
    res = calibrate(qtof_full_samples, "positive",
                    standards_table=qtof_full_standards)
    assert res.model["calibration_scope"] == "project"
    assert pd.to_numeric(res.table["RI_spread"]).isna().all()


# --------------------------------------------------------- warp unit tests ----

def test_fit_warp_extrapolate_false_gives_nan():
    rt = np.linspace(1, 10, 8)
    irt = np.linspace(0, 100, 8)
    cfg = CalibrationConfig()
    w = fit_warp(rt, irt, cfg)
    out = w.predict(np.array([0.5, 5.0, 50.0]), extrapolate=False)
    assert np.isnan(out[0]) and np.isnan(out[2])   # outside span
    assert np.isfinite(out[1])


def test_fit_warp_rejects_degenerate_anchors():
    # all-tied RT -> collapses below hard minimum -> CalibrationError
    rt = np.array([3.0, 3.0, 3.0])
    irt = np.array([10.0, 20.0, 30.0])
    with pytest.raises(CalibrationError):
        fit_warp(rt, irt, CalibrationConfig())


def test_fit_warp_rejects_nonincreasing_irt():
    rt = np.array([1.0, 2.0, 3.0, 4.0])
    irt = np.array([0.0, 50.0, 40.0, 100.0])   # not strictly increasing
    with pytest.raises(CalibrationError):
        fit_warp(rt, irt, CalibrationConfig())


# ----------------------------------------------------------- config presets ---

def test_config_presets_differ():
    q = CalibrationConfig.qtof()
    o = CalibrationConfig.orbitrap()
    assert q.mz_tol_ppm == 15.0
    assert o.mz_tol_ppm == 8.0
    assert q.rt_window_min == 0.5
    assert o.rt_window_min == 0.3


def test_config_from_dict_rejects_unknown_key():
    with pytest.raises(Exception) as ei:
        CalibrationConfig.from_dict({"mz_tol_ppm": 5.0, "bogus_key": 1})
    assert "bogus_key" in str(ei.value)


def test_config_overrides_propagate_into_model(orbitrap_samples, orbitrap_standards):
    cfg = CalibrationConfig.orbitrap(mz_tol_ppm=6.0, extrapolate=True)
    res = calibrate(orbitrap_samples, "positive",
                    standards_table=orbitrap_standards, config=cfg)
    assert res.model["config"]["mz_tol_ppm"] == 6.0
    assert res.model["config"]["extrapolate"] is True
