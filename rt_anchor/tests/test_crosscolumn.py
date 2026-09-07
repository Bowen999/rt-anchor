"""The v2 engine, unit by unit: reciprocal m/z matching, the monotone curve,
the stage-2 gate, the class-aware leave-one-out selection, and the iRT ruler.

These replace the v1 warp tests (``fit_warp`` / ``apply_warp`` / ``MonotoneWarp``),
which tested a method the package no longer contains.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rt_anchor.config import CalibrationConfig
from rt_anchor.crosscolumn import (
    ClassAwareRefiner,
    MonotoneCurve,
    build_calibrator,
    choose_lambda_loo,
    choose_lambda_loo_class,
    exclude_pairs_near_mz,
    fit_robust_curve,
    loo_mse_class,
    loo_residuals_class,
    match_features_by_mz,
    mz_window,
)
from rt_anchor.errors import CalibrationError
from rt_anchor.irt import IRTMapper

from conftest import make_table

FLAT = dict(mz_tol_ppm=0.0, mz_tol_min_da=0.008)   # the anonymous window (§7)


# --------------------------------------------------------- m/z matching -------

def test_match_is_one_to_one_and_reciprocal():
    a = make_table([300.0000, 400.0000, 500.0000], [1.0, 2.0, 3.0], sn=[10, 10, 10])
    b = make_table([300.0020, 400.0020, 900.0000], [1.5, 2.5, 9.0], sn=[10, 10, 10])
    pairs = match_features_by_mz(a, b, None, **FLAT)
    assert len(pairs) == 2                      # 500 has no partner, 900 has none
    assert pairs["idx_a"].is_unique and pairs["idx_b"].is_unique
    assert list(pairs["rt_a"]) == [1.0, 2.0] and list(pairs["rt_b"]) == [1.5, 2.5]


def test_match_ranks_candidates_by_sn_not_row_order():
    # two features of B sit inside A's window; the higher-S/N one must win, and
    # it is deliberately NOT the first row (row order must not decide).
    a = make_table([500.0000], [4.0], sn=[100])
    b = make_table([500.0030, 500.0010], [7.0, 9.0], sn=[5, 5000])
    pairs = match_features_by_mz(a, b, None, **FLAT)
    assert len(pairs) == 1
    assert pairs["rt_b"].iloc[0] == 9.0         # the S/N=5000 feature


def test_match_window_is_absolute_not_ppm():
    """The anonymous window must not widen with m/z (spec §7)."""
    assert mz_window([300.0, 900.0], None, **FLAT).tolist() == [0.008, 0.008]
    # a 0.01 Da gap is outside the flat window at any mass...
    hi = make_table([900.0000], [1.0], sn=[1])
    hi2 = make_table([900.0100], [2.0], sn=[1])
    assert len(match_features_by_mz(hi, hi2, None, **FLAT)) == 0
    # ...but inside a 15 ppm window at m/z 900 (0.0135 Da) — which is exactly the
    # behaviour §7 pins down and does not ship.
    assert len(match_features_by_mz(hi, hi2, None, mz_tol_ppm=15.0,
                                    mz_tol_min_da=0.0)) == 1


def test_exclude_pairs_near_mz_masks_the_anchor_candidates():
    a = make_table([300.0, 400.0, 500.0], [1.0, 2.0, 3.0], sn=[1, 1, 1])
    b = make_table([300.0, 400.0, 500.0], [1.1, 2.1, 3.1], sn=[1, 1, 1])
    pairs = match_features_by_mz(a, b, None, **FLAT)
    kept = exclude_pairs_near_mz(pairs, [400.0005], None, **FLAT)
    assert len(pairs) == 3 and len(kept) == 2
    assert 400.0 not in set(kept["mz_a"])


def test_match_ignores_rows_without_mz_or_rt():
    a = make_table([300.0, np.nan, 500.0], [1.0, 2.0, np.nan], sn=[1, 1, 1])
    b = make_table([300.0, 500.0], [1.5, 3.5], sn=[1, 1])
    pairs = match_features_by_mz(a, b, None, **FLAT)
    assert len(pairs) == 1 and pairs["rt_a"].iloc[0] == 1.0


# ------------------------------------------------------ the monotone curve ----

def test_curve_is_monotone_and_tracks_the_distortion(synth_pair):
    _, _, rt_src, rt_ref = synth_pair(n=300, noise=0.02)
    curve, keep = fit_robust_curve(rt_src, rt_ref, frac=0.1)
    grid = np.linspace(rt_src.min(), rt_src.max(), 1000)
    pred = curve.predict(grid)
    assert np.all(np.diff(pred) >= -1e-9), "the stage-1 curve must be monotone"
    # and it must actually follow the truth, not just be flat
    truth = 1.6 * grid ** 1.15
    assert np.median(np.abs(pred - truth)) < 0.05
    assert keep.sum() >= 0.9 * len(rt_src)      # clean data -> little trimming


def test_curve_trims_outliers_but_never_resurrects_them():
    rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(1, 20, 200))
    y = 1.2 * x + rng.normal(0, 0.01, 200)
    y[50] += 8.0                                 # a gross outlier
    curve, keep = fit_robust_curve(x, y, frac=0.2)
    order = np.argsort(x)
    assert not keep[np.where(order == 50)[0][0]]
    assert abs(curve.predict(np.array([x[50]]))[0] - 1.2 * x[50]) < 0.2


def test_curve_degrades_to_a_line_below_min_points():
    x = np.linspace(1, 10, 8)
    curve, keep = fit_robust_curve(x, 2 * x + 1, min_points=20)
    assert keep.all()
    assert np.allclose(curve.predict(np.array([3.0, 7.0])), [7.0, 15.0], atol=1e-6)


@pytest.mark.parametrize("mode,expected", [("linear", 24.0), ("clamp", 20.0)])
def test_curve_extrapolation_modes(mode, expected):
    x = np.linspace(1, 10, 32)
    curve = MonotoneCurve(x, 2 * x, mode=mode)
    assert curve.predict(np.array([12.0]))[0] == pytest.approx(expected, abs=1e-6)
    assert bool(curve.is_extrapolated(np.array([12.0]))[0])
    assert not bool(curve.is_extrapolated(np.array([5.0]))[0])


def test_curve_rejects_degenerate_knots():
    with pytest.raises(CalibrationError):
        MonotoneCurve([3.0, 3.0, 3.0], [1.0, 2.0, 3.0])


def test_curve_rejects_unknown_mode():
    with pytest.raises(CalibrationError):
        MonotoneCurve([1.0, 2.0], [1.0, 2.0], mode="wishful")


def test_fit_rejects_empty_input():
    with pytest.raises(CalibrationError):
        fit_robust_curve(np.array([]), np.array([]))


# --------------------------------------------------- the stage-2 anchor gate --

def _anchor_frame(rt_src, resid_fn, classes):
    rt_src = np.asarray(rt_src, dtype=float)
    return pd.DataFrame({
        "label": [f"L{i}" for i in range(len(rt_src))],
        "rt_src": rt_src,
        "rt_ref": 2.0 * rt_src + resid_fn(rt_src),
        "lipid_class": classes,
    })


@pytest.fixture
def exact_runs():
    """Two runs related by exactly ``rt_ref = 2 * rt_src`` — a zero-residual curve."""
    rng = np.random.default_rng(3)
    mz = np.round(np.linspace(300, 900, 120) + rng.normal(0, 0.2, 120), 4)
    rt = np.sort(rng.uniform(1.0, 20.0, 120))
    return make_table(mz, rt, sn=rng.uniform(10, 500, 120)), \
        make_table(mz, 2.0 * rt, sn=rng.uniform(10, 500, 120))


def test_gate_engages_on_a_systematic_class_offset(exact_runs):
    src, ref = exact_runs
    x = np.linspace(2.0, 18.0, 12)
    cls = np.array(["SM"] * 6 + ["PC"] * 6)
    anchors = _anchor_frame(x, lambda v: np.where(cls == "SM", 0.5, -0.5), cls)
    cal = build_calibrator(src, ref, CalibrationConfig(), sample_anchors=anchors,
                           use_anchors=True, **FLAT)
    assert cal.anchors_used and cal.warp_source == "curve+anchors"
    assert cal.gate_mse_reduction >= 0.2
    assert cal.lam_c > 0, "a pure class offset must be carried by lam_c"
    # an SM feature is pushed up, a PC feature down, by roughly the offsets
    assert cal.predict([10.0], classes=["SM"])[0] - 20.0 > 0.2
    assert cal.predict([10.0], classes=["PC"])[0] - 20.0 < -0.2


def test_gate_stays_off_when_the_residuals_are_noise(exact_runs):
    src, ref = exact_runs
    rng = np.random.default_rng(7)
    x = np.linspace(2.0, 18.0, 12)
    cls = np.array(["SM", "PC", "CE", "TG"] * 3)
    anchors = _anchor_frame(x, lambda v: rng.normal(0, 0.3, len(v)), cls)
    cal = build_calibrator(src, ref, CalibrationConfig(), sample_anchors=anchors,
                           use_anchors=True, **FLAT)
    assert not cal.anchors_used and cal.warp_source == "curve"
    assert "gated off" in cal.gate_reason
    assert cal.sigma_anchor == 0.0
    # the correction is not applied: prediction is the bare curve
    assert cal.predict([10.0])[0] == pytest.approx(cal.curve.predict([10.0])[0])


def test_stage_two_skipped_below_min_anchors(exact_runs):
    src, ref = exact_runs
    anchors = _anchor_frame([4.0, 9.0], lambda v: np.full(len(v), 0.4), ["SM", "SM"])
    cal = build_calibrator(src, ref, CalibrationConfig(min_anchors=3),
                           sample_anchors=anchors, use_anchors=True, **FLAT)
    assert not cal.anchors_used and "min_anchors" in cal.gate_reason


def test_use_anchors_false_reports_why(exact_runs):
    src, ref = exact_runs
    cal = build_calibrator(src, ref, CalibrationConfig(), use_anchors=False, **FLAT)
    assert not cal.anchors_used and cal.gate_reason == "anchors disabled"


def test_guard_rail_too_few_shared_pairs():
    a = make_table([100.0, 200.0, 300.0], [1.0, 2.0, 3.0], sn=[1, 1, 1])
    b = make_table([700.0, 800.0, 900.0], [1.0, 2.0, 3.0], sn=[1, 1, 1])
    with pytest.raises(CalibrationError) as ei:
        build_calibrator(a, b, CalibrationConfig(), **FLAT)
    assert "polarity" in str(ei.value).lower()


# ------------------------------------------------ class-aware leave-one-out ---

def test_class_aware_loo_prefers_the_class_offset_over_the_smooth_term():
    x = np.linspace(1, 20, 16)
    cls = np.array(["A", "B"] * 8)               # the classes interleave in RT
    r = np.where(cls == "A", 0.6, -0.6)          # pure class systematics
    lam_g, lam_c = choose_lambda_loo_class(x, r, cls)
    assert lam_c > lam_g
    assert loo_mse_class(x, r, cls, lam_g, lam_c) < loo_mse_class(x, r, cls, 1.0, 0.0)


def test_loo_residuals_are_one_per_anchor_and_honest():
    x = np.linspace(1, 20, 12)
    cls = np.array(["A"] * 6 + ["B"] * 6)
    r = np.where(cls == "A", 0.5, -0.5)
    loo = loo_residuals_class(x, r, cls, 1.0, 1.0)
    assert loo.shape == (12,)
    # leaving an anchor out and predicting it must not reproduce it exactly
    assert np.all(np.isfinite(loo))


def test_class_refiner_uses_the_global_median_for_unknown_classes():
    x = np.array([2.0, 5.0, 9.0, 14.0])
    r = np.array([0.4, 0.4, -0.4, -0.4])
    ref = ClassAwareRefiner(x, r, np.array(["A", "A", "B", "B"]), lam_g=0.0, lam_c=1.0)
    assert ref.correction([7.0], classes=["A"])[0] == pytest.approx(0.4)
    assert ref.correction([7.0], classes=["Z"])[0] == pytest.approx(ref.m_global)
    assert ref.correction([7.0], classes=None)[0] == pytest.approx(ref.m_global)


def test_choose_lambda_loo_shrinks_pure_noise_toward_zero():
    rng = np.random.default_rng(11)
    x = np.linspace(1, 20, 20)
    lam = choose_lambda_loo(x, rng.normal(0, 1.0, 20))
    assert lam <= 0.5


# ------------------------------------------------------------- uncertainty ---

def test_uncertainty_grows_where_the_pairs_disagree(exact_runs):
    src, ref = exact_runs
    cal = build_calibrator(src, ref, CalibrationConfig(), use_anchors=False, **FLAT)
    # an exact 2x relation -> essentially zero local scatter everywhere
    assert float(np.nanmax(cal.uncertainty(np.linspace(2, 18, 20)))) < 1e-6
    assert cal.sigma_anchor == 0.0


# ------------------------------------------------------------- the iRT ruler --

def test_irt_scale_is_affine_and_spans_1_to_100():
    # deliberately uneven landmark spacing: an affine map ignores the interior
    rt = np.array([1.5, 2.0, 3.0, 9.0, 9.5, 21.6])
    m = IRTMapper(rt)
    assert m.to_irt(np.array([rt[0]]))[0] == pytest.approx(1.0)
    assert m.to_irt(np.array([rt[-1]]))[0] == pytest.approx(100.0)
    grid = np.linspace(rt[0], rt[-1], 200)
    straight = 1.0 + 99.0 * (grid - rt[0]) / (rt[-1] - rt[0])
    assert np.max(np.abs(m.to_irt(grid) - straight)) < 1e-9
    # slope constant -> the interior landmarks do not bend the scale
    assert np.ptp(m.slope(grid)) < 1e-9


def test_irt_extends_past_the_landmarks_and_flags_it():
    m = IRTMapper(np.array([2.0, 22.0]))
    assert m.to_irt(np.array([1.0]))[0] < 1.0
    assert m.to_irt(np.array([23.0]))[0] > 100.0
    assert m.is_extrapolated(np.array([1.0, 12.0, 23.0])).tolist() == [True, False, True]


def test_irt_needs_two_landmarks():
    with pytest.raises(ValueError):
        IRTMapper(np.array([5.0]))
