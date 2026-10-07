"""Stage 1b — the homologous-series term (:mod:`rt_anchor.series`).

Unit coverage of the algorithm (Kendrick phase, clustering, carbon-order
validation, leave-own-family-out, the gate) lives here on synthetic series;
the shipped-dataset behaviour is pinned against the reference implementation's
numbers in ``tests/data/expected_series_repo_datasets.json``, and the v1.2.1
reproduction is pinned by the numeric goldens in
``tests/data/expected_v121_orbitrap_golden.json`` (compared to 1e-8 — a text
hash would break on another numpy/pandas build although the engine is right).

Both goldens also record a **fingerprint** of the stage-1 fit as the reference
environment produced it (pairs matched, pairs kept, and a sha1 of the kept
pairs' (m/z, RT)): which pairs survive MAD trimming depends on how tied rows
happen to be ordered, and numpy/pandas order ties differently across versions
and CPUs. A platform whose fit legitimately differs skips the pinned-number
comparison instead of failing (the fingerprint helpers live in ``conftest.py``,
shared with ``test_lattice.py``).
"""

from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate
from rt_anchor.cli import main as cli_main
from rt_anchor.config import CalibrationConfig
from rt_anchor.crosscolumn import AnchorRefiner, ColumnCalibrator, MonotoneCurve
from rt_anchor.helpers import _sanitize_json
from rt_anchor.io.schema import RESULT_COLUMNS
from rt_anchor.pipeline import _series_columns, _warp_rows
from rt_anchor.series import CH2_MASS, KENDRICK_FACTOR, PERIOD, SeriesTerm, kendrick_phase

from conftest import DATA, _skip_unless_reference_stage1_fit

EXPECTED_PATH = DATA / "expected_series_repo_datasets.json"


# ---------------------------------------------------------------------------
# Synthetic homologous series
# ---------------------------------------------------------------------------

def _homologue_masses(phase_target, ks):
    """m/z values whose Kendrick phase is exactly ``phase_target``.

    ``k * CH2 + offset`` with the offset chosen so the phase lands on the
    target — the mass-defect signature of one homologous series.
    """
    off = CH2_MASS * phase_target / PERIOD
    return np.array([k * CH2_MASS + off for k in ks], dtype=float)


def _two_series():
    """Two homologous series with opposite, slowly-varying curve residuals.

    Series A: ks 40..50 with 42 missing (a gap to interpolate across), phase
    3.0. Series B: ks 60..69, phase 8.0. RT rises 0.8 min per CH2 on both
    columns; residuals run +side along A and -side along B, with a bump on one
    member of A so a leave-own-out prediction is visibly not self-fitting.
    """
    ks_a = [40, 41, 43, 44, 45, 46, 47, 48, 49, 50]
    ks_b = list(range(60, 70))
    mz = np.concatenate([_homologue_masses(3.0, ks_a), _homologue_masses(8.0, ks_b)])
    xa = np.concatenate([[5.0 + 0.8 * (k - 40) for k in ks_a],
                         [5.3 + 0.8 * (k - 60) for k in ks_b]])
    r = np.concatenate([0.10 + 0.02 * np.arange(10), -(0.08 + 0.015 * np.arange(10))])
    r[5] += 0.05                      # the bump
    xb = xa + r                       # identity curve + residual
    return mz, xa, xb, r, float(xa.min()), float(xa.max())


# ---------------------------------------------------------------------------
# kendrick_phase
# ---------------------------------------------------------------------------

def test_ch2_homologues_share_a_phase():
    p0 = kendrick_phase(700.0)
    for k in (1, 2, 5):
        assert kendrick_phase(700.0 + k * CH2_MASS)[()] == pytest.approx(p0, abs=1e-9)


def test_a_double_bond_moves_the_phase():
    # one more double bond: -2H; one more oxygen: a different class. Both must
    # land far from the original phase on the 14-Da circle.
    p0 = float(kendrick_phase(700.0))

    def circ(p, q):
        d = abs(p - q)
        return min(d, PERIOD - d)

    assert circ(float(kendrick_phase(700.0 - 2 * 1.007825)), p0) > 0.1
    assert circ(float(kendrick_phase(700.0 + 15.994915)), p0) > 0.1


def test_phase_lives_on_the_circle():
    phases = kendrick_phase(np.linspace(100.0, 1500.0, 500))
    assert ((phases >= 0.0) & (phases < PERIOD)).all()


def test_clustering_wraps_around_the_zero_fourteen_seam():
    # six members of one series straddling the seam: three just below 14,
    # three wrapped just above 0. Neither side alone reaches min_members (4),
    # so a series is found only if the wrap-around merge fires.
    ks = list(range(40, 46))
    base = _homologue_masses(13.995, ks)
    mz = base + np.array([-0.0035, -0.002, -0.0005, 0.006, 0.0075, 0.009])
    phases = kendrick_phase(mz)
    assert (phases > 13.98).sum() == 3 and (phases < 0.01).sum() == 3
    xa = 5.0 + np.arange(6)
    r = 0.1 + 0.01 * np.arange(6)
    term = SeriesTerm.fit(mz, xa, xa + r, r, 5.0, 10.0)
    assert term.n_series == 1
    assert term.member_mask.all()
    # and a homologue query beyond the seam-side end is corrected
    corr, n = term.raw_correction(_homologue_masses(13.995, [47]), [11.5])
    assert n[0] == 6 and corr[0] == pytest.approx(r[-1], abs=1e-12)


# ---------------------------------------------------------------------------
# fit + predict behaviour
# ---------------------------------------------------------------------------

def test_two_opposite_series_are_found_and_the_gate_engages():
    mz, xa, xb, r, x0, x1 = _two_series()
    term = SeriesTerm.fit(mz, xa, xb, r, x0, x1)
    assert term.n_series == 2
    assert term.member_mask.all()
    assert term.engaged
    assert term.gate_mse_reduction >= term.gate_threshold
    assert term.n_pairs_covered == 20
    assert "engaged" in term.gate_reason

    # a query in A's gap (k=42), halfway in RT between k=41 and k=43, is
    # interpolated from exactly those two neighbours
    mq = _homologue_masses(3.0, [42])[0]
    corr, n = term.raw_correction([mq], [6.6])
    assert n[0] == 10
    assert corr[0] == pytest.approx((r[1] + r[2]) / 2, abs=1e-12)

    # a member queried at its own coordinates is corrected by its NEIGHBOURS:
    # its prediction is the interpolation of the members on either side
    for i in range(1, 9):
        corr_i, n_i = term.raw_correction(mz[[i]], xa[[i]])
        a, b = i - 1, i + 1
        t = (xa[i] - xa[a]) / (xa[b] - xa[a])
        assert n_i[0] == 9                       # itself excluded from M
        assert corr_i[0] == pytest.approx(r[a] + t * (r[b] - r[a]), abs=1e-12)

    # the engaged correction is the raw correction; a gated-off term is zeros
    assert np.array_equal(term.correction(mz, xa)[0], term.raw_correction(mz, xa)[0])


def test_leave_own_family_out():
    mz, xa, xb, r, x0, x1 = _two_series()
    dirty = r.copy()
    dirty[5] += 5.0                            # one pair's residual goes wild
    clean = SeriesTerm.fit(mz, xa, xb, r, x0, x1)
    wild = SeriesTerm.fit(mz, xa, xb, dirty, x0, x1)
    # a query at that pair's own (m/z, RT) is unaffected by its residual ...
    c0, _ = clean.raw_correction(mz[[5]], xa[[5]])
    c1, _ = wild.raw_correction(mz[[5]], xa[[5]])
    assert c0[0] == c1[0]
    # ... but its neighbours' corrections move
    c0, _ = clean.raw_correction(mz[[4]], xa[[4]])
    c1, _ = wild.raw_correction(mz[[4]], xa[[4]])
    assert c0[0] != c1[0]


def test_coeluting_family_members_are_excluded_from_the_correction():
    # two members of one series 0.02 min apart (< exclude_rt_min): a query at
    # the first one's coordinates must not see either of them
    ks = [40, 41, 42, 43, 44, 45]
    mz = _homologue_masses(3.0, ks)
    xa = np.array([5.0, 5.02, 6.0, 7.0, 8.0, 9.0])
    r = 0.1 + 0.02 * np.arange(6)
    term = SeriesTerm.fit(mz, xa, xa + r, r, 5.0, 9.0)
    assert term.exclude_rt_min == pytest.approx(0.03)
    corr, n = term.raw_correction(mz[[0]], xa[[0]])
    assert n[0] == 4                           # k=40 (self) and k=41 excluded
    assert corr[0] == pytest.approx(r[2], abs=1e-12)   # the next member up
    # even if the co-eluting member's residual is an outlier
    dirty = r.copy()
    dirty[1] += 5.0
    term2 = SeriesTerm.fit(mz, xa, xa + dirty, dirty, 5.0, 9.0)
    corr2, _ = term2.raw_correction(mz[[0]], xa[[0]])
    assert corr2[0] == pytest.approx(r[2], abs=1e-12)


def test_carbon_order_must_hold_on_both_columns():
    ks = list(range(40, 46))
    mz = _homologue_masses(3.0, ks)
    xa = 5.0 + np.arange(6)                    # rises with m/z on the source
    r = np.full(6, 0.2)
    # falling with m/z on the reference -> no chain -> no series
    term = SeriesTerm.fit(mz, xa, xa[::-1].copy(), r, 5.0, 10.0)
    assert term.n_series == 0 and not term.member_mask.any()
    # and the mirror image: falling on the source
    term2 = SeriesTerm.fit(mz, xa[::-1].copy(), xa, r, 5.0, 10.0)
    assert term2.n_series == 0
    assert term2.gate_reason == ("only 0 matched pairs sit in a validated homologous "
                                 "series (need 20) — stage-1 curve only")


def test_end_reach_extends_two_ch2_but_not_four():
    ks = list(range(40, 46))
    mz = _homologue_masses(3.0, ks)
    xa = 5.0 + np.arange(6)
    r = 0.1 + 0.02 * np.arange(6)
    term = SeriesTerm.fit(mz, xa, xa + r, r, 5.0, 10.0)
    assert term.n_series == 1
    # one CH2 beyond the last member: the end member's residual, carried out
    corr, n = term.raw_correction(_homologue_masses(3.0, [46]), [11.5])
    assert n[0] == 6 and corr[0] == pytest.approx(r[-1], abs=1e-12)
    # two CH2 is still within reach (2 * CH2 + 1 Da of slack) ...
    corr, n = term.raw_correction(_homologue_masses(3.0, [47]), [11.5])
    assert n[0] == 6
    # ... three or four are not
    for k in (48, 49):
        corr, n = term.raw_correction(_homologue_masses(3.0, [k]), [11.5])
        assert n[0] == 0 and corr[0] == 0.0
    # the same below the series: one CH2 under the first member, later RTs off
    corr, n = term.raw_correction(_homologue_masses(3.0, [39]), [4.0])
    assert n[0] == 6 and corr[0] == pytest.approx(r[0], abs=1e-12)
    # ... but only when the query sits outside the series' elution span
    corr, n = term.raw_correction(_homologue_masses(3.0, [39]), [6.0])
    assert n[0] == 0


# ---------------------------------------------------------------------------
# the gate
# ---------------------------------------------------------------------------

def test_gate_declines_on_pure_noise():
    ks_a, ks_b = list(range(40, 50)), list(range(60, 70))
    mz = np.concatenate([_homologue_masses(3.0, ks_a), _homologue_masses(8.0, ks_b)])
    xa = np.concatenate([5.0 + 0.8 * np.arange(10), 5.3 + 0.8 * np.arange(10)])
    rng = np.random.default_rng(5)
    r = rng.normal(0.0, 0.3, 20)
    term = SeriesTerm.fit(mz, xa, xa + r, r, 5.0, 13.0)
    assert term.n_series == 2 and term.n_pairs_covered == 20
    assert not term.engaged
    assert term.gate_mse_reduction < term.gate_threshold
    assert "gated off" in term.gate_reason
    # gated off -> correction() is all zeros, raw stays inspectable
    corr, n = term.correction(mz, xa)
    assert (corr == 0.0).all() and (n == 0).all()
    assert (term.raw_correction(mz, xa)[1] > 0).all()


def test_gate_declines_when_too_few_pairs_are_covered():
    ks = list(range(40, 45))                   # one series of five
    mz = _homologue_masses(3.0, ks)
    xa = 5.0 + np.arange(5)
    r = 0.1 + 0.02 * np.arange(5)
    term = SeriesTerm.fit(mz, xa, xa + r, r, 5.0, 9.0)
    assert term.n_series == 1 and term.n_pairs_covered == 5
    assert not term.engaged and math.isnan(term.gate_mse_reduction)
    assert term.gate_reason == ("only 5 matched pairs sit in a validated homologous "
                                "series (need 20) — stage-1 curve only")


def test_self_calibration_reports_no_correction_needed():
    # r == 0 everywhere: the gate ratio would be 0/0 — say why, not a percentage
    ks_a, ks_b = list(range(40, 50)), list(range(60, 70))
    mz = np.concatenate([_homologue_masses(3.0, ks_a), _homologue_masses(8.0, ks_b)])
    xa = np.concatenate([5.0 + 0.8 * np.arange(10), 5.3 + 0.8 * np.arange(10)])
    term = SeriesTerm.fit(mz, xa, xa.copy(), np.zeros(20), 5.0, 13.0)
    assert term.n_pairs_covered == 20
    assert not term.engaged and term.gate_mse_reduction == 0.0
    assert term.gate_reason == ("no correction needed: the stage-1 curve already "
                                "reproduces the series pairs")


def test_non_finite_pairs_and_queries_are_ignored():
    mz, xa, xb, r, x0, x1 = _two_series()
    mz_bad = np.append(mz, np.nan)
    xa_bad = np.append(xa, 7.0)
    xb_bad = np.append(xb, 7.0)
    r_bad = np.append(r, 0.0)
    term = SeriesTerm.fit(mz_bad, xa_bad, xb_bad, r_bad, x0, x1)
    assert term.member_mask.shape == mz_bad.shape
    assert not term.member_mask[-1]
    assert term.n_pairs_covered == 20
    corr, n = term.raw_correction([np.nan, mz[0]], [1.0, np.nan])
    assert (corr == 0.0).all() and (n == 0).all()


# ---------------------------------------------------------------------------
# §4 item 3: predict(rt) without mz never includes the series term
# ---------------------------------------------------------------------------

def test_predict_without_mz_never_includes_the_series_term():
    mz, xa, xb, r, x0, x1 = _two_series()
    term = SeriesTerm.fit(mz, xa, xb, r, x0, x1)
    assert term.engaged
    grid = np.linspace(x0 - 0.5, x1 + 0.5, 32)
    curve = MonotoneCurve(grid, 1.0 * grid)
    cal = ColumnCalibrator(curve=curve, series=term)
    assert cal.series_used and cal.warp_source == "curve+series"
    qrt = np.linspace(x0 + 0.4, x1 - 0.4, 50)
    assert np.array_equal(cal.predict(qrt), curve.predict(qrt))
    assert np.array_equal(cal.stage1_predict(qrt), curve.predict(qrt))
    qmz = np.full_like(qrt, _homologue_masses(3.0, [44])[0])
    with_mz = cal.predict(qrt, mz=qmz)
    assert not np.allclose(with_mz, curve.predict(qrt))
    assert np.allclose(with_mz, curve.predict(qrt) + term.correction(qmz, qrt)[0])


def test_warp_rows_and_series_columns_track_what_was_applied():
    mz, xa, xb, r, x0, x1 = _two_series()
    term = SeriesTerm.fit(mz, xa, xb, r, x0, x1)
    grid = np.linspace(x0 - 0.5, x1 + 0.5, 32)
    curve = MonotoneCurve(grid, 1.0 * grid)
    cal = ColumnCalibrator(curve=curve, series=term)
    n_mem = np.array([0, 4, 0, 5])
    assert list(_warp_rows(cal, n_mem)) == ["curve", "curve+series", "curve", "curve+series"]
    cal_anchored = ColumnCalibrator(curve=curve, series=term,
                                    refiner=AnchorRefiner([6.0, 12.0], [0.1, 0.1]))
    assert cal_anchored.warp_source == "curve+series+anchors"
    assert list(_warp_rows(cal_anchored, n_mem)) == [
        "curve+anchors", "curve+series+anchors", "curve+anchors", "curve+series+anchors"]
    # _series_columns returns the gated correction, zeros when the gate is off
    corr, n = _series_columns(cal, mz, xa)
    assert (n > 0).all() and np.allclose(corr, term.raw_correction(mz, xa)[0])
    cal_off = ColumnCalibrator(curve=curve)
    corr, n = _series_columns(cal_off, mz, xa)
    assert (corr == 0.0).all() and (n == 0).all()


# ---------------------------------------------------------------------------
# §4 item 4: the shipped datasets against the reference implementation
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module", params=["Orbitrap", "QTOF", "QTOF_full"])
def dataset_result(request):
    name = request.param
    res = calibrate(str(DATA / name / "samples" / "aligned_feature_table.txt"),
                    "positive",
                    standards_table=str(DATA / name / "standards" / "aligned_feature_table.txt"),
                    panel="mix15")
    return name, res


def test_shipped_datasets_match_the_reference_implementation(dataset_result):
    name, res = dataset_result
    exp = json.loads(EXPECTED_PATH.read_text())[name]
    _skip_unless_reference_stage1_fit(res, exp["stage1_fingerprint"])
    st = res.calibrator.series
    assert st is not None
    assert res.calibrator.curve.x0 == pytest.approx(exp["x0"], abs=1e-9)
    assert res.calibrator.curve.x1 == pytest.approx(exp["x1"], abs=1e-9)
    assert st.exclude_rt_min == pytest.approx(exp["exclude_rt_min"], abs=1e-9)
    assert st.n_series == exp["n_series"]
    assert st.n_pairs_covered == exp["n_pairs_covered"]
    assert st.gate_mse_reduction == pytest.approx(exp["gate_mse_reduction"], abs=1e-9)
    assert st.engaged == exp["engaged"]
    # the model block reports the same fit
    s = res.model["series"]
    assert s["enabled"] and s["engaged"] == exp["engaged"]
    assert s["n_series"] == exp["n_series"] and s["n_pairs_covered"] == exp["n_pairs_covered"]
    assert s["gate_mse_reduction"] == pytest.approx(exp["gate_mse_reduction"], abs=1e-9)
    if name == "Orbitrap":
        # a gated-off run whose leave-own-out error is a *loss* (reduction
        # -0.345 here) must not read "only -35% below" — say it plainly
        assert "35% above curve-only" in st.gate_reason
        assert "-35" not in st.gate_reason

    mz = pd.to_numeric(res.table[res.mz_col], errors="coerce").to_numpy(dtype=float)
    rt = res.rt_minutes().to_numpy(dtype=float)
    assert len(mz) == exp["n_features"]
    corr, n = st.raw_correction(mz, rt)
    assert int((n > 0).sum()) == exp["n_features_raw_covered"]
    assert float(np.abs(corr).sum()) == pytest.approx(exp["sum_abs_raw_correction"], abs=1e-6)
    assert float(np.abs(corr).max()) == pytest.approx(exp["max_abs_raw_correction"], abs=1e-6)
    for row in exp["samples"]:
        i = row["row"]
        assert mz[i] == pytest.approx(row["mz"], abs=1e-6)
        assert rt[i] == pytest.approx(row["rt_min"], abs=1e-9)
        assert corr[i] == pytest.approx(row["raw_correction_min"], abs=1e-9)
        assert int(n[i]) == row["n_members"]

    # none of the three engages the series term -> it corrects nothing: no row
    # reads "+series" and both series columns are all zero (whether a row
    # carries "+lattice" — QTOF_full's do — is test_lattice.py's business)
    assert not res.series_used
    assert not res.table[res.col("warp_source")].astype(str).str.contains(
        "+series", regex=False).any()
    assert (res.values("series_n_members") == 0).all()
    assert (res.values("series_correction_min") == 0.0).all()


# ---------------------------------------------------------------------------
# §4 items 1-2: v1.2.1 reproduction, pinned by numeric goldens
# ---------------------------------------------------------------------------

#: The golden is ``tests/data/expected_v121_orbitrap_golden.json``, generated
#: from this tree with the series term off (the configuration the reviewer
#: verified reproduces v1.2.1 exactly on six independent datasets). Numbers are
#: compared to 1e-8 — a text hash would break on another numpy/pandas build
#: although the engine is right.
GOLDEN_V121_PATH = DATA / "expected_v121_orbitrap_golden.json"
GOLDEN_SAMPLE_STEP = 97
GOLDEN_NUMERIC_COLUMNS = ["Cal_RT_min", "Cal_RT_uncertainty_min", "iRT",
                          "iRT_uncertainty", "RI_spread", "n_contributing",
                          "series_correction_min", "series_n_members"]
GOLDEN_TEXT_COLUMNS = ["iRT_reliability", "calibration_scope", "warp_source",
                       "is_extrapolated"]


def _golden_snapshot(res):
    """The v1.2.1-comparable numbers of one run, as a JSON-safe dict."""
    n = len(res.table)
    rows = np.arange(0, n, GOLDEN_SAMPLE_STEP)
    numeric = {}
    for c in GOLDEN_NUMERIC_COLUMNS:
        v = pd.to_numeric(res.table[res.col(c)], errors="coerce").to_numpy(dtype=float)
        numeric[c] = [None if not math.isfinite(x) else float(x) for x in v[rows]]
    text = {}
    for c in GOLDEN_TEXT_COLUMNS:
        vc = res.table[res.col(c)].astype(str).value_counts()
        text[c] = {str(k): int(v) for k, v in sorted(vc.items())}
    anch_num = {}
    for c in sorted(res.anchors.select_dtypes(include=[np.number]).columns):
        v = pd.to_numeric(res.anchors[c], errors="coerce").to_numpy(dtype=float)
        anch_num[c] = [None if not math.isfinite(x) else float(x) for x in v]
    pairs = res.pairs
    pair_counts = {"total": int(len(pairs)),
                   "kept": int(pairs["kept"].astype(bool).sum()) if len(pairs) else 0,
                   "standards": int((pairs["source"] == "standards").sum()) if len(pairs) else 0,
                   "sample": int((pairs["source"] == "sample").sum()) if len(pairs) else 0}
    m = res.model
    model = {"curve": m["curve"], "irt": m["irt"],
             "anchors": {k: v for k, v in m["anchors"].items() if k != "table"}}
    return {"n_rows": n, "sample_step": GOLDEN_SAMPLE_STEP,
            "result_columns": numeric, "text_value_counts": text,
            "anchors_numeric": anch_num, "pairs": pair_counts,
            "model": _sanitize_json(model)}


def _assert_golden(actual, want, path="golden"):
    """Numbers to 1e-8 absolute; counts, strings and structure exact."""
    if isinstance(want, dict):
        assert isinstance(actual, dict) and set(actual) == set(want), path
        for k in want:
            _assert_golden(actual[k], want[k], f"{path}.{k}")
        return
    if isinstance(want, list):
        assert isinstance(actual, (list, tuple)) and len(actual) == len(want), path
        for i, (a, w) in enumerate(zip(actual, want)):
            _assert_golden(a, w, f"{path}[{i}]")
        return
    if want is None:
        assert actual is None or (isinstance(actual, float) and math.isnan(actual)), path
        return
    if isinstance(want, bool):
        assert actual is want or actual == want, path
        return
    if isinstance(want, str):
        assert actual == want, path
        return
    if isinstance(want, int):
        assert int(actual) == want, path
        return
    np.testing.assert_allclose(float(actual), float(want), rtol=0, atol=1e-8,
                               err_msg=path)


@pytest.fixture(scope="module")
def v121_golden():
    return json.loads(GOLDEN_V121_PATH.read_text())["configurations"]


def test_v121_default_configuration_matches_the_golden_numbers(
        orbitrap_samples, orbitrap_standards, v121_golden):
    """v1.2.1's default == series term off + stage-2 anchors on (§4 item 1).

    Lattice term off explicitly too: v1.2.1 has no stage 1c, so its gate must not decide this.
    """
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=CalibrationConfig.orbitrap(
                        use_sample_anchors=True, use_series_term=False,
                        use_lattice_term=False))
    want = dict(v121_golden["v121_default"])
    _skip_unless_reference_stage1_fit(res, want.pop("stage1_fingerprint"))
    _assert_golden(_golden_snapshot(res), want)


def test_v121_no_sample_anchors_configuration_matches_the_golden_numbers(
        orbitrap_samples, orbitrap_standards, v121_golden):
    """v1.2.1's ``--no-sample-anchors`` == both stages off (§4 item 2).

    Lattice term off explicitly too: v1.2.1 has no stage 1c, so its gate must not decide this.
    """
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=CalibrationConfig.orbitrap(
                        use_sample_anchors=False, use_series_term=False,
                        use_lattice_term=False))
    want = dict(v121_golden["v121_no_sample_anchors"])
    _skip_unless_reference_stage1_fit(res, want.pop("stage1_fingerprint"))
    _assert_golden(_golden_snapshot(res), want)


# ---------------------------------------------------------------------------
# §4 item 5: the bundled reference calibrated against itself
# ---------------------------------------------------------------------------

def test_self_calibration_of_the_reference_is_exact(bundled_reference):
    res = calibrate(bundled_reference.sample_path, "positive",
                    standards_table=bundled_reference.standards_path, panel="mix21")
    rt = res.rt_minutes().to_numpy()
    cal = res.values("Cal_RT_min").to_numpy()
    ok = np.isfinite(rt) & np.isfinite(cal)
    assert np.max(np.abs(cal[ok] - rt[ok])) <= 1e-9
    s = res.model["series"]
    assert s["engaged"] is False and s["gate_mse_reduction"] == 0.0
    assert s["gate_reason"] == ("no correction needed: the stage-1 curve already "
                                "reproduces the series pairs")
    assert not res.series_used


# ---------------------------------------------------------------------------
# table shape, per-sample tier, CLI
# ---------------------------------------------------------------------------

def test_table_shape_and_series_columns(orbitrap_samples, orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=CalibrationConfig.orbitrap())
    t = res.table
    # the two series columns are followed by the lattice term's two (stage 1c)
    assert RESULT_COLUMNS[-4:] == ["series_correction_min", "series_n_members",
                                   "lattice_correction_min", "lattice_n_members"]
    assert res.col("series_correction_min") in t.columns
    assert res.col("series_n_members") in t.columns
    # Orbitrap does not engage either term -> zeros, and the bare-curve warp source
    assert (res.values("series_n_members") == 0).all()
    assert (res.values("series_correction_min") == 0.0).all()
    assert set(t[res.col("warp_source")]) == {"curve"}
    # the model block is always present
    s = res.model["series"]
    for key in ("enabled", "engaged", "gate_mse_reduction", "gate_threshold",
                    "gate_reason", "n_series", "n_pairs_covered",
                    "n_features_corrected", "kmd_tol", "min_members",
                    "end_reach_ch2", "exclude_rt_min"):
        assert key in s, f"model['series'] is missing '{key}'"
    assert s["enabled"] and not s["engaged"]
    assert s["n_features_corrected"] == 0
    # the pairs frame carries the membership flag
    assert "in_series" in res.pairs.columns
    assert res.pairs["in_series"].dtype == bool
    assert int(res.pairs["in_series"].sum()) >= 4 * s["n_series"]
    assert any("series term:" in line for line in res.log)


def test_per_sample_tier_reports_the_series_columns(qtof_full_samples,
                                                    qtof_full_standards,
                                                    qtof_full_single_files):
    res = calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                    panel="mix15", single_files=list(qtof_full_single_files))
    assert res.model["calibration_scope"] == "sample"
    corr = res.values("series_correction_min").to_numpy(dtype=float)
    n_mem = res.values("series_n_members").to_numpy(dtype=int)
    assert corr.shape == n_mem.shape == (len(res.table),)
    assert (n_mem >= 0).all()
    assert (corr[n_mem == 0] == 0.0).all()
    ws = res.table[res.col("warp_source")].to_numpy()
    # the per-injection columns aggregate, so a row can carry both stage-1
    # terms ("curve+series+lattice"); "+anchors" is appended on every row when
    # stage 2 engages
    assert set(ws) <= {"curve", "curve+series", "curve+lattice", "curve+series+lattice",
                       "curve+anchors", "curve+series+anchors", "curve+lattice+anchors",
                       "curve+series+lattice+anchors"}
    # the series flag decides whether "+series" is there (stage 2 is off by
    # default; the lattice tag beside it is test_lattice.py's business)
    assert np.array_equal(np.array(["+series" in w for w in ws]), n_mem > 0)
    assert res.model["series"]["n_features_corrected"] == int((n_mem > 0).sum())


def test_cli_no_series_term(tmp_path, orbitrap_samples, orbitrap_standards):
    prefix = str(tmp_path / "nos")
    rc = cli_main(["calibrate", "--samples", orbitrap_samples,
                   "--standards", orbitrap_standards, "--polarity", "positive",
                   "--panel", "mix15", "--no-series-term", "--no-report",
                   "--out", prefix])
    assert rc == 0
    model = json.loads((tmp_path / "nos_model.json").read_text())
    s = model["series"]
    assert s["enabled"] is False and s["engaged"] is False
    assert s["gate_reason"] == "disabled (use_series_term=False)"
    assert s["n_series"] == 0 and s["n_pairs_covered"] == 0
    assert s["gate_mse_reduction"] is None
    assert model["config"]["use_series_term"] is False
    pairs = pd.read_csv(f"{prefix}_pairs.csv")
    assert not pairs["in_series"].any()
    assert ("series term: disabled (use_series_term=False) — stage-1 curve only"
            in (tmp_path / "nos_log.txt").read_text(encoding="utf-8"))


def test_cli_sample_anchors_is_opt_in(tmp_path, orbitrap_samples, orbitrap_standards):
    # default: stage 2 off, and the anchors CSV is empty
    prefix = str(tmp_path / "off")
    rc = cli_main(["calibrate", "--samples", orbitrap_samples,
                   "--standards", orbitrap_standards, "--polarity", "positive",
                   "--panel", "mix15", "--no-report", "--out", prefix])
    assert rc == 0
    model = json.loads((tmp_path / "off_model.json").read_text())
    assert model["config"]["use_sample_anchors"] is False
    assert len(pd.read_csv(f"{prefix}_anchors.csv")) == 0

    # --sample-anchors opts in
    prefix = str(tmp_path / "on")
    rc = cli_main(["calibrate", "--samples", orbitrap_samples,
                   "--standards", orbitrap_standards, "--polarity", "positive",
                   "--panel", "mix15", "--sample-anchors", "--no-report",
                   "--out", prefix])
    assert rc == 0
    model = json.loads((tmp_path / "on_model.json").read_text())
    assert model["config"]["use_sample_anchors"] is True
    assert len(pd.read_csv(f"{prefix}_anchors.csv")) > 0

    # --no-sample-anchors stays accepted and wins over --sample-anchors
    prefix = str(tmp_path / "both")
    rc = cli_main(["calibrate", "--samples", orbitrap_samples,
                   "--standards", orbitrap_standards, "--polarity", "positive",
                   "--panel", "mix15", "--sample-anchors", "--no-sample-anchors",
                   "--no-report", "--out", prefix])
    assert rc == 0
    model = json.loads((tmp_path / "both_model.json").read_text())
    assert model["config"]["use_sample_anchors"] is False
