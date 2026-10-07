"""Stage 1c — the lattice term (:mod:`rt_anchor.lattice`).

Unit coverage of the algorithm (the H2 step on the Kendrick circle, lattice
coordinates, order pruning on both columns, leave-own-family-out with a reach
limit, the gate's basis, and the precedence behind the series term) lives here
on synthetic classes whose ladders are deliberately too short for the series
term; the shipped-dataset behaviour is pinned against the reference
implementation's numbers in ``tests/data/expected_lattice_repo_datasets.json``,
and the v1.2.2 reproduction — with the lattice term off, 1.2.3 must reproduce
1.2.2 exactly — is pinned by the numeric goldens in
``tests/data/expected_v122_golden.json`` (compared to 1e-8, like the v1.2.1
goldens in ``test_series.py``).

Both goldens carry a stage-1 fingerprint (see ``conftest.py``): a platform
whose fit legitimately differs skips the pinned-number comparison instead of
failing.
"""

from __future__ import annotations

import json
import math
import warnings

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate
from rt_anchor.cli import main as cli_main
from rt_anchor.config import CalibrationConfig
from rt_anchor.crosscolumn import ColumnCalibrator, MonotoneCurve
from rt_anchor.errors import ConfigError
from rt_anchor.helpers import _sanitize_json
from rt_anchor.io.schema import RESULT_COLUMNS
from rt_anchor.lattice import H2_MASS, H2_PHASE, LatticeTerm
from rt_anchor.series import CH2_MASS, PERIOD, SeriesTerm, kendrick_phase

from conftest import DATA, _skip_unless_reference_stage1_fit
from test_series import _assert_golden

EXPECTED_PATH = DATA / "expected_lattice_repo_datasets.json"
GOLDEN_V122_PATH = DATA / "expected_v122_golden.json"


# ---------------------------------------------------------------------------
# Synthetic lattice classes
# ---------------------------------------------------------------------------

def _lattice_class(phase_target, na, nb, k0=40, rt0=5.0,
                   offset=0.15, slope_a=0.02, slope_b=-0.01):
    """A synthetic lipid class on the lattice ``mz = off + a*CH2 + b*H2``.

    ``a`` counts CH2 (carbons), ``b`` counts H2 — one more b is one double bond
    fewer. The b=0 ladder sits at Kendrick phase ``phase_target``; with
    ``na < 4`` every ladder is too short for the series term. RT rises with
    both coordinates (reversed phase: more carbons and/or fewer double bonds
    elutes later) and the residual is the class-wide offset plus a smooth tilt,
    so no member ever elutes out of order on either column.
    """
    off = CH2_MASS * phase_target / PERIOD + k0 * CH2_MASS
    aa = np.repeat(np.arange(na), nb)
    bb = np.tile(np.arange(nb), na)
    mz = off + aa * CH2_MASS + bb * H2_MASS
    xa = rt0 + 0.8 * aa + 0.35 * bb
    r = offset + slope_a * aa + slope_b * bb
    return mz, xa, xa + r, r, aa, bb


def _class21():
    """The 21-member class: 7 ladders of 3 members, one double bond apart."""
    mz, xa, xb, r, aa, bb = _lattice_class(3.0, 3, 7)
    return mz, xa, xb, r, aa, bb, float(xa.min()), float(xa.max())


# ---------------------------------------------------------------------------
# the H2 step
# ---------------------------------------------------------------------------

def test_one_double_bond_is_one_h2_step_on_the_circle():
    def circ(d):
        m = d % PERIOD
        return min(m, PERIOD - m)
    p0 = float(kendrick_phase(700.0))
    # one more double bond: -H2 -> exactly H2_PHASE away on the circle
    p1 = float(kendrick_phase(700.0 - H2_MASS))
    assert circ(p1 - p0) == pytest.approx(H2_PHASE, abs=1e-9)
    # a different class (different heteroatom mass defect — one more oxygen) is
    # not a whole number of H2 steps away, so it can never link
    p_ox = float(kendrick_phase(700.0 + 15.994915))
    assert min(circ(p_ox - p0 - s * H2_PHASE) for s in (1, 2, -1, -2)) > 0.008


def test_ladders_link_only_when_one_h2_step_apart():
    # two ladders of one class, one double bond apart: one family of 6
    mz, xa, xb, r, _, _ = _lattice_class(3.0, 3, 2, offset=0.1)
    term = LatticeTerm.fit(mz, xa, xb, r, float(xa.min()), float(xa.max()))
    assert term.n_families == 1 and term.family_sizes == [6]
    # add a third ladder at a phase that is no whole number of H2 steps away:
    # it stays a lone ladder — the series term's business, not the lattice's
    mz2, xa2, xb2, r2, _, _ = _lattice_class(4.5, 3, 1, offset=0.1)
    mz_c = np.concatenate([mz, mz2])
    xa_c = np.concatenate([xa, xa2 + 12.0])         # keep the RT spans disjoint
    r_c = np.concatenate([r, r2])
    xb_c = xa_c + r_c
    term = LatticeTerm.fit(mz_c, xa_c, xb_c, r_c, float(xa_c.min()), float(xa_c.max()))
    assert term.n_families == 1 and term.family_sizes == [6]
    assert term.member_mask.sum() == 6 and not term.member_mask[6:].any()


# ---------------------------------------------------------------------------
# fit + predict behaviour
# ---------------------------------------------------------------------------

def test_short_ladders_join_into_one_lattice_family():
    mz, xa, xb, r, aa, bb, x0, x1 = _class21()
    # every ladder has 3 members < series_min_members: the series term is silent
    assert SeriesTerm.fit(mz, xa, xb, r, x0, x1).n_series == 0
    term = LatticeTerm.fit(mz, xa, xb, r, x0, x1)
    assert term.n_families == 1 and term.family_sizes == [21]
    assert term.member_mask.all()
    # the coordinates are exactly the synthetic (a, b), and m0 the class root
    idx, fa, fb, m0 = term._families[0]
    assert sorted(zip(fa.tolist(), fb.tolist())) == [(a, b)
                                                     for a in range(3) for b in range(7)]
    root = CH2_MASS * 3.0 / PERIOD + 40 * CH2_MASS
    assert m0 == pytest.approx(root, abs=1e-9)
    # the gate engages on the class-wide offset
    assert term.engaged and term.n_pairs_covered == 21
    assert term.gate_mse_reduction >= term.gate_threshold
    assert "lattice term engaged" in term.gate_reason
    # a member's correction comes from its lattice neighbours: member (1,1) is
    # predicted from its four nearest donors, (0,1)/(2,1) at d=0.5 and
    # (1,0)/(1,2) at d=1, with w = 1/(d + 0.5)^2 — itself excluded
    i = int(np.flatnonzero((aa == 1) & (bb == 1))[0])
    corr, n = term.raw_correction(mz[[i]], xa[[i]])
    r01, r21 = 0.15 - 0.01, 0.15 + 0.04 - 0.01
    r10, r12 = 0.15 + 0.02, 0.15 + 0.02 - 0.02
    expected = (1.0 * (r01 + r21) + (4.0 / 9.0) * (r10 + r12)) / (2.0 + 8.0 / 9.0)
    assert n[0] == 12                     # every member within d <= 3, minus itself
    assert corr[0] == pytest.approx(expected, abs=1e-12)


def test_two_classes_with_opposite_offsets_stay_two_families():
    mzA, xaA, xbA, rA, _, _ = _lattice_class(3.0, 3, 7)
    mzB, xaB, xbB, rB, _, _ = _lattice_class(9.0, 3, 7, offset=-0.12,
                                             slope_a=0.015, slope_b=-0.008)
    mz = np.concatenate([mzA, mzB])
    xa = np.concatenate([xaA, xaB])
    r = np.concatenate([rA, rB])
    term = LatticeTerm.fit(mz, xa, xa + r, r, float(xa.min()), float(xa.max()))
    assert term.n_families == 2 and term.family_sizes == [21, 21]
    # the families are found in phase order: class B's seam-crossing ladder
    # sorts first, so B's family is found first
    rootA = CH2_MASS * 3.0 / PERIOD + 40 * CH2_MASS
    rootB = CH2_MASS * 9.0 / PERIOD + 40 * CH2_MASS
    assert term.family_offsets[0] == pytest.approx(rootB, abs=1e-9)
    assert term.family_offsets[1] == pytest.approx(rootA, abs=1e-9)
    assert term.engaged
    # and each class's correction carries its own offset's sign
    corrA, _ = term.raw_correction(mzA[[10]], xaA[[10]])
    corrB, _ = term.raw_correction(mzB[[10]], xaB[[10]])
    assert corrA[0] > 0.0 and corrB[0] < 0.0


def test_order_pruning_drops_the_out_of_order_member():
    mz, xa, xb, r, aa, bb, x0, x1 = _class21()
    iX = int(np.flatnonzero((aa == 2) & (bb == 3))[0])
    # member (2,3) elutes out of order on the SOURCE column only: earlier than
    # every member below it on the lattice
    xa_p = xa.copy()
    xa_p[iX] = 5.4
    clean = LatticeTerm.fit(mz, xa, xb, r, x0, x1)
    pruned = LatticeTerm.fit(mz, xa_p, xa_p + r, r, x0, x1)
    assert pruned.n_families == 1 and pruned.family_sizes == [20]
    assert not pruned.member_mask[iX] and pruned.member_mask.sum() == 20
    # and the dropped member is not a donor: a query at its neighbour (2,4)
    # loses exactly that one donor
    i24 = int(np.flatnonzero((aa == 2) & (bb == 4))[0])
    c_clean, n_clean = clean.raw_correction(mz[[i24]], xa[[i24]])
    c_pruned, n_pruned = pruned.raw_correction(mz[[i24]], xa[[i24]])
    assert n_pruned[0] == n_clean[0] - 1
    assert c_pruned[0] != c_clean[0]


def test_isomers_on_one_lattice_point_do_not_conflict():
    mz, xa, xb, r, aa, bb, x0, x1 = _class21()
    i11 = int(np.flatnonzero((aa == 1) & (bb == 1))[0])
    # a second feature on (1,1)'s lattice point, 0.2 min later — same (a, b),
    # so it can never be in conflict, and it donates like any other member
    mz_i = np.append(mz, mz[i11] + 0.001)
    xa_i = np.append(xa, xa[i11] + 0.2)
    r_i = np.append(r, r[i11])
    term = LatticeTerm.fit(mz_i, xa_i, xa_i + r_i, r_i, x0, x1)
    assert term.n_families == 1 and term.family_sizes == [22]
    corr, n = term.raw_correction(mz[[i11]], xa[[i11]])
    assert n[0] == 13                     # the 12 neighbours plus the isomer
    # the isomer is the nearest donor (d=0): the four averaged are the isomer,
    # (0,1) and (2,1) at d=0.5, and (1,0) at d=1 (member order breaks the tie
    # against (1,2))
    expected = (4.0 * 0.16 + 0.14 + 0.18 + (4.0 / 9.0) * 0.17) / (4.0 + 1.0 + 1.0 + 4.0 / 9.0)
    assert corr[0] == pytest.approx(expected, abs=1e-12)


def test_leave_own_family_out():
    mz, xa, xb, r, aa, bb, x0, x1 = _class21()
    i13 = int(np.flatnonzero((aa == 1) & (bb == 3))[0])
    dirty = r.copy()
    dirty[i13] += 5.0                            # one pair's residual goes wild
    clean = LatticeTerm.fit(mz, xa, xb, r, x0, x1)
    wild = LatticeTerm.fit(mz, xa, xa + dirty, dirty, x0, x1)
    # a query at that pair's own (m/z, RT) is unaffected by its residual ...
    c0, _ = clean.raw_correction(mz[[i13]], xa[[i13]])
    c1, _ = wild.raw_correction(mz[[i13]], xa[[i13]])
    assert c0[0] == c1[0]
    # ... but its neighbours' corrections move
    i03 = int(np.flatnonzero((aa == 0) & (bb == 3))[0])
    c0, _ = clean.raw_correction(mz[[i03]], xa[[i03]])
    c1, _ = wild.raw_correction(mz[[i03]], xa[[i03]])
    assert c0[0] != c1[0]


def test_reach_limits_the_donors():
    mz, xa, xb, r, _, _ = _lattice_class(3.0, 3, 2)      # 6 members, b in {0, 1}
    term = LatticeTerm.fit(mz, xa, xb, r, float(xa.min()), float(xa.max()))
    assert term.n_families == 1
    m0 = term.family_offsets[0]
    xq = [float(xa.max()) + 3.0]
    # (a=7, b=0): two members in reach — (2,0) at d=2.5, (1,0) at d=3.0
    corr, n = term.raw_correction([m0 + 7 * CH2_MASS], xq)
    assert n[0] == 2 and corr[0] != 0.0
    # (a=8, b=0): a single member in reach (d=3.0) — fewer than min_neighbours
    corr, n = term.raw_correction([m0 + 8 * CH2_MASS], xq)
    assert n[0] == 0 and corr[0] == 0.0
    # (a=10, b=0): nobody in reach
    corr, n = term.raw_correction([m0 + 10 * CH2_MASS], xq)
    assert n[0] == 0 and corr[0] == 0.0


@pytest.mark.parametrize("overrides,expected_n", [
    (dict(lattice_min_neighbours=0), [2, 1, 0, 0]),
    (dict(lattice_max_neighbours=0), [2, 0, 0, 0]),
    (dict(lattice_min_neighbours=0, lattice_max_neighbours=0), [2, 1, 0, 0]),
])
def test_zero_neighbour_settings_are_clamped_and_never_give_nan(overrides, expected_n):
    # A correction needs at least one donor. A bound of 0 used to let a query
    # that sits on a lattice point with nobody in reach divide 0 by 0 (every
    # weight zero), or with max 0 average no donor at all, and come out NaN;
    # both bounds are now clamped to 1 where they are stored.
    mz, xa, xb, r, _, _ = _lattice_class(3.0, 3, 2)      # 6 members, b in {0, 1}
    x0, x1 = float(xa.min()), float(xa.max())
    term = LatticeTerm.fit(mz, xa, xb, r, x0, x1, CalibrationConfig(**overrides))
    assert term.n_families == 1
    assert term.min_neighbours >= 1 and term.max_neighbours >= 1
    m0 = term.family_offsets[0]
    # lattice points (a, b) = (7, 0): two donors in reach; (8, 0): one donor;
    # (4, 4): located, but nobody within three steps; (10, 0): too far out to
    # be located at all (see test_reach_limits_the_donors)
    mq = [m0 + a * CH2_MASS + b * H2_MASS for a, b in ((7, 0), (8, 0), (4, 4), (10, 0))]
    xq = [x1 + 3.0] * len(mq)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)   # a 0/0 also warns
        corr, n = term.raw_correction(mq, xq)
        # every pair of the fit is a leave-own-out query too — the gate's basis
        fit_corr, _ = term.raw_correction(mz, xa)
    assert np.isfinite(corr).all() and np.isfinite(fit_corr).all()
    assert n.tolist() == expected_n
    assert (corr[n == 0] == 0.0).all() and (corr[n > 0] != 0.0).all()
    # the same through correction(): never NaN, whatever the gate decided
    assert np.isfinite(term.correction(mq, xq)[0]).all()
    # the defaults are what they were: two donors needed, four averaged
    default = LatticeTerm.fit(mz, xa, xb, r, x0, x1)
    assert (default.min_neighbours, default.max_neighbours) == (2, 4)
    assert default.raw_correction(mq, xq)[1].tolist() == [2, 0, 0, 0]


def test_non_finite_pairs_and_queries_are_ignored():
    mz, xa, xb, r, aa, bb, x0, x1 = _class21()
    mz_bad = np.append(mz, np.nan)
    xa_bad = np.append(xa, 7.0)
    xb_bad = np.append(xb, 7.0)
    r_bad = np.append(r, 0.0)
    term = LatticeTerm.fit(mz_bad, xa_bad, xb_bad, r_bad, x0, x1)
    assert term.member_mask.shape == mz_bad.shape
    assert not term.member_mask[-1]
    assert term.n_pairs_covered == 21
    corr, n = term.raw_correction([np.nan, mz[0]], [1.0, np.nan])
    assert (corr == 0.0).all() and (n == 0).all()


# ---------------------------------------------------------------------------
# the gate and its basis
# ---------------------------------------------------------------------------

def test_gate_declines_on_pure_noise():
    mz, xa, _, _, _, _, x0, x1 = _class21()
    rng = np.random.default_rng(5)
    r = rng.normal(0.0, 0.3, len(mz))
    term = LatticeTerm.fit(mz, xa, xa + r, r, x0, x1)
    assert term.n_families == 1 and term.n_pairs_covered == 21
    assert not term.engaged
    assert term.gate_mse_reduction < term.gate_threshold
    assert "gated off" in term.gate_reason
    # gated off -> correction() is all zeros, raw stays inspectable
    corr, n = term.correction(mz, xa)
    assert (corr == 0.0).all() and (n == 0).all()
    assert (term.raw_correction(mz, xa)[1] > 0).all()


def test_gate_declines_when_too_few_pairs_are_covered():
    mz, xa, xb, r, _, _ = _lattice_class(3.0, 3, 2)      # one family of six
    term = LatticeTerm.fit(mz, xa, xb, r, float(xa.min()), float(xa.max()))
    assert term.n_families == 1 and term.n_pairs_covered == 6
    assert not term.engaged and math.isnan(term.gate_mse_reduction)
    assert term.gate_reason == ("only 6 matched pairs the series term leaves alone "
                                "sit in a lattice family (need 20) — no lattice correction")


def test_series_covered_pairs_do_not_count_towards_the_gate():
    mz, xa, xb, r, _, _, x0, x1 = _class21()
    covered = np.zeros(len(mz), dtype=bool)
    covered[:5] = True
    term = LatticeTerm.fit(mz, xa, xb, r, x0, x1, series_covered=covered)
    assert term.n_pairs_covered == 16            # 21 lattice pairs minus the 5
    assert not term.engaged and math.isnan(term.gate_mse_reduction)
    assert term.gate_reason == ("only 16 matched pairs the series term leaves alone "
                                "sit in a lattice family (need 20) — no lattice correction")
    covered_all = np.ones(len(mz), dtype=bool)
    term = LatticeTerm.fit(mz, xa, xb, r, x0, x1, series_covered=covered_all)
    assert term.n_pairs_covered == 0
    # the masked pairs still sit in the family — the basis is the gate's, not
    # the fit's
    assert term.n_families == 1 and term.member_mask.all()


def test_self_calibration_reports_no_correction_needed():
    # r == 0 everywhere: the gate ratio would be 0/0 — say why, not a percentage
    mz, xa, _, _, _, _, x0, x1 = _class21()
    term = LatticeTerm.fit(mz, xa, xa.copy(), np.zeros(len(mz)), x0, x1)
    assert term.n_pairs_covered == 21
    assert not term.engaged and term.gate_mse_reduction == 0.0
    assert term.gate_reason == ("no correction needed: the stage-1 curve already "
                                "reproduces the lattice pairs")


# ---------------------------------------------------------------------------
# precedence (§2.4) and the bare-curve guarantee (§4 item 2)
# ---------------------------------------------------------------------------

def _series_and_lattice_class():
    """Three ladders of 8 members: long enough for the series term AND one family."""
    mz, xa, xb, r, aa, bb = _lattice_class(3.0, 8, 3, offset=0.0,
                                           slope_a=0.01, slope_b=0.1)
    return mz, xa, xb, r, aa, bb, float(xa.min()), float(xa.max())


def test_series_term_has_precedence_in_stage1_terms():
    mz, xa, xb, r, aa, bb, x0, x1 = _series_and_lattice_class()
    sterm = SeriesTerm.fit(mz, xa, xb, r, x0, x1)
    assert sterm.engaged and sterm.n_pairs_covered == 24
    # the lattice gate is judged only on what the series term leaves alone
    series_covered = sterm.correction(mz, xa)[1] > 0
    masked = LatticeTerm.fit(mz, xa, xb, r, x0, x1, series_covered=series_covered)
    assert masked.n_families == 1 and masked.n_pairs_covered == 0
    assert not masked.engaged
    # fitted on its own, the same family engages
    lattice = LatticeTerm.fit(mz, xa, xb, r, x0, x1)
    assert lattice.engaged and lattice.n_pairs_covered == 24

    grid = np.linspace(x0 - 0.5, x1 + 0.5, 32)
    curve = MonotoneCurve(grid, 1.0 * grid)
    cal = ColumnCalibrator(curve=curve, series=sterm, lattice=lattice)
    assert cal.series_used and cal.lattice_used
    assert cal.warp_source == "curve+series+lattice"

    # a member both terms cover: the series correction stands, lattice is zeroed
    i = int(np.flatnonzero((aa == 3) & (bb == 0))[0])
    assert lattice.raw_correction(mz[[i]], xa[[i]])[1][0] > 0
    s_corr, s_n, l_corr, l_n = cal.stage1_terms(xa[[i]], mz[[i]])
    assert s_n[0] > 0 and l_n[0] == 0 and l_corr[0] == 0.0
    # three CH2 beyond the series' end (past its reach of two): the series term
    # is silent, the lattice term — max_distance 3 — still corrects
    mq = np.array([CH2_MASS * 3.0 / PERIOD + 40 * CH2_MASS + 10 * CH2_MASS])
    xq = np.array([xa.max() + 2.4])
    assert sterm.raw_correction(mq, xq)[1][0] == 0
    s_corr, s_n, l_corr, l_n = cal.stage1_terms(xq, mq)
    assert s_n[0] == 0 and l_n[0] == 6 and l_corr[0] != 0.0
    # stage1_predict with mz is the curve plus both applied corrections ...
    qrt = np.linspace(x0 + 0.4, x1 - 0.4, 50)
    qmz = np.full_like(qrt, mz[i])
    s_c, _, l_c, _ = cal.stage1_terms(qrt, qmz)
    assert np.allclose(cal.stage1_predict(qrt, mz=qmz), curve.predict(qrt) + s_c + l_c)
    # ... and without mz it is the bare curve, exactly (§4 item 2)
    assert np.array_equal(cal.predict(qrt), curve.predict(qrt))
    assert np.array_equal(cal.stage1_predict(qrt), curve.predict(qrt))


# ---------------------------------------------------------------------------
# §4 item 4: the shipped datasets against the reference implementation
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lattice_results():
    """The three shipped datasets under the defaults (both stage-1 terms on)."""
    return {name: calibrate(str(DATA / name / "samples" / "aligned_feature_table.txt"),
                            "positive",
                            standards_table=str(DATA / name / "standards"
                                                / "aligned_feature_table.txt"),
                            panel="mix15")
            for name in ("Orbitrap", "QTOF", "QTOF_full")}


@pytest.fixture(scope="module")
def lattice_off_results():
    """The same three runs with ``use_lattice_term=False`` (the 1.2.2 engine)."""
    return {name: calibrate(str(DATA / name / "samples" / "aligned_feature_table.txt"),
                            "positive",
                            standards_table=str(DATA / name / "standards"
                                                / "aligned_feature_table.txt"),
                            panel="mix15",
                            config=CalibrationConfig(use_lattice_term=False))
            for name in ("Orbitrap", "QTOF", "QTOF_full")}


@pytest.mark.parametrize("name", ["Orbitrap", "QTOF", "QTOF_full"])
def test_shipped_datasets_match_the_lattice_reference(name, lattice_results):
    res = lattice_results[name]
    exp = json.loads(EXPECTED_PATH.read_text())[name]
    _skip_unless_reference_stage1_fit(res, exp["stage1_fingerprint"])
    lt = res.calibrator.lattice
    assert lt is not None
    assert res.calibrator.curve.x0 == pytest.approx(exp["x0"], abs=1e-9)
    assert res.calibrator.curve.x1 == pytest.approx(exp["x1"], abs=1e-9)
    assert lt.exclude_rt_min == pytest.approx(exp["exclude_rt_min"], abs=1e-9)
    assert lt.n_families == exp["n_families"]
    assert lt.family_sizes == exp["family_sizes"]
    assert lt.family_offsets == pytest.approx(exp["family_m0"], abs=1e-9)
    assert int(lt.member_mask.sum()) == exp["n_pairs_in_families"]
    # the gate's basis: the pairs the engaged series term leaves alone
    kept = res.pairs[res.pairs["kept"].astype(bool)]
    s_cov = res.calibrator.series.correction(kept["mz_src"].to_numpy(dtype=float),
                                             kept["rt_src"].to_numpy(dtype=float))[1] > 0
    assert int(s_cov.sum()) == exp["n_pairs_series_covered"]
    assert res.model["series"]["engaged"] == exp["series_engaged"]
    assert lt.n_pairs_covered == exp["n_pairs_covered"]
    assert lt.gate_mse_reduction == pytest.approx(exp["gate_mse_reduction"], abs=1e-9)
    assert lt.engaged == exp["engaged"]
    # the model block reports the same fit
    l = res.model["lattice"]
    assert l["enabled"] and l["engaged"] == exp["engaged"]
    assert l["n_families"] == exp["n_families"]
    assert l["n_pairs_in_families"] == exp["n_pairs_in_families"]
    assert l["n_pairs_covered"] == exp["n_pairs_covered"]
    assert l["gate_mse_reduction"] == pytest.approx(exp["gate_mse_reduction"], abs=1e-9)
    assert l["exclude_rt_min"] == pytest.approx(exp["exclude_rt_min"], abs=1e-9)

    mz = pd.to_numeric(res.table[res.mz_col], errors="coerce").to_numpy(dtype=float)
    rt = res.rt_minutes().to_numpy(dtype=float)
    assert len(mz) == exp["n_features"]
    corr, n = lt.raw_correction(mz, rt)
    assert int((n > 0).sum()) == exp["n_features_raw_covered"]
    assert float(np.abs(corr).sum()) == pytest.approx(exp["sum_abs_raw_correction"], abs=1e-6)
    assert float(np.abs(corr).max()) == pytest.approx(exp["max_abs_raw_correction"], abs=1e-6)
    for row in exp["samples"]:
        i = row["row"]
        assert mz[i] == pytest.approx(row["mz"], abs=1e-6)
        assert rt[i] == pytest.approx(row["rt_min"], abs=1e-9)
        assert corr[i] == pytest.approx(row["raw_correction_min"], abs=1e-9)
        assert int(n[i]) == row["n_members"]

    # the applied columns: QTOF_full engages (the series term does not, so the
    # lattice term corrects every feature it covers); the other two do not
    n_applied = res.values("lattice_n_members").to_numpy(dtype=int)
    assert int((n_applied > 0).sum()) == exp["n_features_corrected"]
    if name == "QTOF_full":
        assert exp["engaged"] and exp["n_features_corrected"] == exp["n_features_raw_covered"]
    else:
        assert exp["n_features_corrected"] == 0


# ---------------------------------------------------------------------------
# §4 item 1: use_lattice_term=False reproduces 1.2.2 exactly
# ---------------------------------------------------------------------------

#: The golden is ``tests/data/expected_v122_golden.json``, generated by the
#: reviewer from the untouched 1.2.2 tree. The snapshot below is exactly the
#: generating script's — the columns / keys 1.2.3 adds are simply not read.
GOLDEN_V122_STEP = 97
GOLDEN_V122_NUMERIC = ["Cal_RT_min", "Cal_RT_uncertainty_min", "iRT",
                       "iRT_uncertainty", "RI_spread", "n_contributing",
                       "series_correction_min", "series_n_members"]
GOLDEN_V122_TEXT = ["iRT_reliability", "calibration_scope", "warp_source",
                    "is_extrapolated"]


def _v122_snapshot(res):
    """The 1.2.2-comparable numbers of one run, as a JSON-safe dict."""
    n = len(res.table)
    rows = np.arange(0, n, GOLDEN_V122_STEP)
    numeric = {}
    for c in GOLDEN_V122_NUMERIC:
        v = pd.to_numeric(res.table[res.col(c)], errors="coerce").to_numpy(dtype=float)
        numeric[c] = [None if not math.isfinite(x) else float(x) for x in v[rows]]
    text = {}
    for c in GOLDEN_V122_TEXT:
        vc = res.table[res.col(c)].astype(str).value_counts()
        text[c] = {str(k): int(v) for k, v in sorted(vc.items())}
    anch = {}
    for c in sorted(res.anchors.select_dtypes(include=[np.number]).columns):
        v = pd.to_numeric(res.anchors[c], errors="coerce").to_numpy(dtype=float)
        anch[c] = [None if not math.isfinite(x) else float(x) for x in v]
    p = res.pairs
    pairs = {"total": int(len(p)), "kept": int(p["kept"].astype(bool).sum()),
             "standards": int((p["source"] == "standards").sum()),
             "sample": int((p["source"] == "sample").sum()),
             "in_series": int(p["in_series"].astype(bool).sum())}
    m = res.model
    model = {"curve": m["curve"], "series": m["series"], "irt": m["irt"],
             "anchors": {k: v for k, v in m["anchors"].items() if k != "table"}}
    return {"n_rows": n, "sample_step": GOLDEN_V122_STEP,
            "result_columns": numeric, "text_value_counts": text,
            "anchors_numeric": anch, "pairs": pairs, "model": _sanitize_json(model)}


def _run_v122(name, single=False, **cfg_overrides):
    """The golden script's ``run()``, with the lattice term switched off."""
    kw = {}
    if single:
        d = DATA / name / "samples" / "single_files"
        kw["single_files"] = sorted(str(p) for p in d.glob("*.txt"))
    return calibrate(str(DATA / name / "samples" / "aligned_feature_table.txt"),
                     "positive",
                     standards_table=str(DATA / name / "standards"
                                         / "aligned_feature_table.txt"),
                     panel="mix15",
                     config=CalibrationConfig(use_lattice_term=False, **cfg_overrides),
                     **kw)


@pytest.fixture(scope="module")
def v122_golden():
    return json.loads(GOLDEN_V122_PATH.read_text())["configurations"]


@pytest.mark.parametrize("key,args", [
    ("orbitrap_default", dict(name="Orbitrap")),
    ("orbitrap_sample_anchors", dict(name="Orbitrap", use_sample_anchors=True)),
    ("qtof_full_default", dict(name="QTOF_full")),
    ("qtof_full_per_sample", dict(name="QTOF_full", single=True)),
])
def test_v122_reproduced_with_lattice_off(key, args, v122_golden):
    res = _run_v122(**args)
    want = dict(v122_golden[key])
    _skip_unless_reference_stage1_fit(res, want.pop("stage1_fingerprint"))
    _assert_golden(_v122_snapshot(res), want)


# ---------------------------------------------------------------------------
# §4 item 3: precedence on the shipped datasets
# ---------------------------------------------------------------------------

def test_series_precedence_on_shipped_datasets(lattice_results, lattice_off_results):
    for name, res in lattice_results.items():
        off = lattice_off_results[name]
        s_n = res.values("series_n_members").to_numpy(dtype=int)
        # the lattice term never changes what the series term does ...
        assert np.array_equal(s_n, off.values("series_n_members").to_numpy(dtype=int))
        sel = s_n > 0
        assert np.allclose(res.values("series_correction_min").to_numpy(dtype=float)[sel],
                           off.values("series_correction_min").to_numpy(dtype=float)[sel])
        assert np.allclose(res.values("Cal_RT_min").to_numpy(dtype=float)[sel],
                           off.values("Cal_RT_min").to_numpy(dtype=float)[sel])
        # ... a per-project row never carries both terms ...
        l_n = res.values("lattice_n_members").to_numpy(dtype=int)
        assert (s_n[l_n > 0] == 0).all()
        # ... and the lattice term only ever adds to the 1.2.2 output
        cal_on = res.values("Cal_RT_min").to_numpy(dtype=float)
        cal_off = off.values("Cal_RT_min").to_numpy(dtype=float)
        l_corr = res.values("lattice_correction_min").to_numpy(dtype=float)
        ok = np.isfinite(cal_on)
        bare = ok & (l_n == 0)
        assert np.allclose(cal_on[bare], cal_off[bare], atol=1e-9)
        touched = ok & (l_n > 0)
        assert np.allclose(cal_on[touched], cal_off[touched] + l_corr[touched], atol=1e-9)


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
    l = res.model["lattice"]
    assert l["enabled"] and l["engaged"] is False
    assert l["gate_mse_reduction"] == 0.0
    assert l["gate_reason"] == ("no correction needed: the stage-1 curve already "
                                "reproduces the lattice pairs")
    assert not res.lattice_used


# ---------------------------------------------------------------------------
# table shape, per-sample tier, CLI
# ---------------------------------------------------------------------------

def test_table_shape_and_lattice_columns(orbitrap_samples, orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                    panel="mix15", config=CalibrationConfig.orbitrap())
    t = res.table
    assert RESULT_COLUMNS[-2:] == ["lattice_correction_min", "lattice_n_members"]
    assert res.col("lattice_correction_min") in t.columns
    assert res.col("lattice_n_members") in t.columns
    # Orbitrap does not engage -> zeros
    assert (res.values("lattice_n_members") == 0).all()
    assert (res.values("lattice_correction_min") == 0.0).all()
    assert not res.lattice_used
    # the model block is always present
    l = res.model["lattice"]
    for key in ("enabled", "engaged", "gate_mse_reduction", "gate_threshold",
                "gate_reason", "n_families", "n_pairs_in_families",
                "n_pairs_covered", "n_features_corrected", "kmd_tol",
                "min_members", "max_distance", "min_neighbours",
                "max_neighbours", "exclude_rt_min"):
        assert key in l, f"model['lattice'] is missing '{key}'"
    assert l["enabled"] and not l["engaged"]
    assert l["n_features_corrected"] == 0
    # the pairs frame carries the membership flag
    assert "in_lattice" in res.pairs.columns
    assert res.pairs["in_lattice"].dtype == bool
    assert int(res.pairs["in_lattice"].sum()) == l["n_pairs_in_families"]
    assert any(line.startswith("lattice term:") for line in res.log)


def test_qtof_full_engages_and_corrects(lattice_results):
    res = lattice_results["QTOF_full"]
    assert res.lattice_used
    n = res.values("lattice_n_members").to_numpy(dtype=int)
    corr = res.values("lattice_correction_min").to_numpy(dtype=float)
    assert (n > 0).any() and (corr[n == 0] == 0.0).all()
    # per-row warp_source carries "+lattice" exactly where a correction applied
    # (the series term is gated off on this dataset, anchors are off)
    ws = res.table[res.col("warp_source")].to_numpy()
    assert set(ws) <= {"curve", "curve+lattice"}
    assert np.array_equal(np.array(["+lattice" in w for w in ws]), n > 0)
    # Cal_RT is the bare 1.2.2 output plus the applied lattice correction
    rt = res.rt_minutes().to_numpy(dtype=float)
    ok = np.isfinite(rt)
    assert np.allclose(res.values("Cal_RT_min").to_numpy(dtype=float)[ok],
                       res.curve.predict(rt[ok]) + corr[ok], atol=1e-9)


def test_per_sample_tier_reports_the_lattice_columns(qtof_full_samples,
                                                     qtof_full_standards,
                                                     qtof_full_single_files):
    res = calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                    panel="mix15", single_files=list(qtof_full_single_files))
    assert res.model["calibration_scope"] == "sample"
    corr = res.values("lattice_correction_min").to_numpy(dtype=float)
    n_mem = res.values("lattice_n_members").to_numpy(dtype=int)
    assert corr.shape == n_mem.shape == (len(res.table),)
    assert (n_mem >= 0).all()
    assert (corr[n_mem == 0] == 0.0).all()
    # the columns aggregate over injections, so a row may show both terms; the
    # lattice flag tracks lattice_n_members exactly
    ws = res.table[res.col("warp_source")].to_numpy()
    assert np.array_equal(np.array(["+lattice" in w for w in ws]), n_mem > 0)
    assert res.model["lattice"]["n_features_corrected"] == int((n_mem > 0).sum())
    # on QTOF_full the per-injection lattice terms engage: some rows carry it
    assert (n_mem > 0).any()


def test_cli_no_lattice_term(tmp_path, orbitrap_samples, orbitrap_standards):
    prefix = str(tmp_path / "nol")
    rc = cli_main(["calibrate", "--samples", orbitrap_samples,
                   "--standards", orbitrap_standards, "--polarity", "positive",
                   "--panel", "mix15", "--no-lattice-term", "--no-report",
                   "--out", prefix])
    assert rc == 0
    model = json.loads((tmp_path / "nol_model.json").read_text())
    l = model["lattice"]
    assert l["enabled"] is False and l["engaged"] is False
    assert l["gate_reason"] == "disabled (use_lattice_term=False)"
    assert l["n_families"] == 0 and l["n_pairs_covered"] == 0
    assert l["n_pairs_in_families"] == 0 and l["n_features_corrected"] == 0
    assert l["gate_mse_reduction"] is None
    # the parameters still come from the config
    assert l["kmd_tol"] == 0.008 and l["min_members"] == 6
    assert l["max_distance"] == 3.0
    assert l["min_neighbours"] == 2 and l["max_neighbours"] == 4
    assert model["config"]["use_lattice_term"] is False
    # the series term is untouched
    assert model["series"]["enabled"] is True
    pairs = pd.read_csv(f"{prefix}_pairs.csv")
    assert not pairs["in_lattice"].any()
    assert ("lattice term: disabled (use_lattice_term=False)"
            in (tmp_path / "nol_log.txt").read_text(encoding="utf-8"))


def test_lattice_config_defaults_and_from_dict():
    c = CalibrationConfig()
    assert c.use_lattice_term is True
    assert c.lattice_kmd_tol == 0.008 and c.lattice_min_members == 6
    assert c.lattice_max_distance == 3.0
    assert c.lattice_min_neighbours == 2 and c.lattice_max_neighbours == 4
    assert c.lattice_gate_min_mse_reduction == 0.2 and c.lattice_min_covered_pairs == 20
    assert CalibrationConfig.from_dict(c.to_dict()) == c
    with pytest.raises(ConfigError):
        CalibrationConfig.from_dict({"lattice_kmd_tol": 0.01, "bogus_key": 1})
