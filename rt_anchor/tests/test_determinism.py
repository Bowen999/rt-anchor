"""Calibration must not depend on how a platform's sort orders tied keys.

Retention times are heavily tied (the Orbitrap fixture has 595 matched pairs
but only ~390 distinct source RTs), and numpy's default ``np.argsort`` is not
stable: its x86-64 SIMD paths and its arm64 paths order equal keys
differently. Before this was fixed, the same input produced a different set of
trimmed pairs, a ~3% different ``Cal_RT_uncertainty_min`` and — on the Orbitrap
fixture — a different stage-2 gate decision on macOS arm64 and on Linux /
Windows x64.

These tests stand in for "another CPU": they replace ``np.argsort`` with a
*valid* sort that deliberately orders ties differently (reversed, or at random)
and require the calibration output to be unchanged. A sort that breaks ties by
the data (``np.lexsort``) or asks for ``kind="stable"`` is passed through
untouched, so only code that still leans on an unspecified tie order can fail.
"""

from __future__ import annotations

import numpy as np
import pytest

from rt_anchor import CalibrationConfig, calibrate
from rt_anchor.crosscolumn import fit_robust_curve_sorted

_REAL_ARGSORT = np.argsort


def _adversarial_argsort(mode: str, seed: int = 0):
    """A correct ``np.argsort`` whose order for *equal* keys is chosen by ``mode``.

    ``reverse``: tied rows come out in reverse input order. ``shuffle``: tied
    rows come out in a random (seeded) order. A stable request (``kind="stable"``
    / ``"mergesort"``), a multi-dimensional input or a structured ``order`` is
    delegated to the real implementation, as a stable sort is *defined* to keep
    input order for ties.
    """
    rng = np.random.default_rng(seed)

    def argsort(a, axis=-1, kind=None, order=None, **kwargs):
        arr = np.asarray(a)
        if kind in ("stable", "mergesort") or arr.ndim != 1 or order is not None:
            return _REAL_ARGSORT(arr, axis=axis, kind=kind, order=order, **kwargs)
        n = arr.shape[0]
        tie = np.arange(n)[::-1] if mode == "reverse" else rng.permutation(n)
        return np.lexsort((tie, arr))

    return argsort


# --------------------------------------------------------------- the engine ---

def test_the_adversary_really_reorders_ties():
    """Guard the guard: the stand-in sort differs from a stable one on ties."""
    a = np.array([2.0, 1.0, 2.0, 1.0, 2.0])
    stable = np.argsort(a, kind="stable")
    mp = pytest.MonkeyPatch()
    try:
        mp.setattr(np, "argsort", _adversarial_argsort("reverse"))
        assert list(np.argsort(a)) == [3, 1, 4, 2, 0]       # valid, ties reversed
        assert list(np.argsort(a, kind="stable")) == list(stable)
    finally:
        mp.undo()
    assert list(stable) == [1, 3, 0, 2, 4]


@pytest.mark.parametrize("mode", ["reverse", "shuffle"])
def test_robust_curve_ignores_tie_order(mode, monkeypatch):
    """The stage-1 curve and its trimmed set depend on the data, not tie order."""
    rng = np.random.default_rng(3)
    x = np.round(rng.uniform(1.0, 20.0, 600), 1)           # ~190 distinct, many ties
    y = 1.6 * x ** 1.15 + rng.normal(0, 0.05, x.size)
    y[rng.choice(x.size, 40, replace=False)] += rng.normal(0, 1.5, 40)   # outliers

    curve0, keep0, order0 = fit_robust_curve_sorted(x, y)
    monkeypatch.setattr(np, "argsort", _adversarial_argsort(mode))
    curve1, keep1, order1 = fit_robust_curve_sorted(x, y)

    assert np.array_equal(order0, order1)
    assert np.array_equal(keep0, keep1)
    grid = np.linspace(x.min(), x.max(), 200)
    assert np.allclose(curve0.predict(grid), curve1.predict(grid), rtol=0, atol=1e-9)


# ------------------------------------------------------ the shipped fixtures ---

def _snapshot(res) -> dict:
    """The calibration output a user sees, in a comparable form."""
    a = res.model["anchors"]
    return {
        "cal": res.values("Cal_RT_min").to_numpy(dtype=float),
        "unc": res.values("Cal_RT_uncertainty_min").to_numpy(dtype=float),
        "kept": res.pairs["kept"].to_numpy(dtype=bool),
        "n_kept": int(res.pairs["kept"].astype(bool).sum()),
        "anchors_engaged": bool(a["engaged"]),
        "gate_reduction": a["gate_mse_reduction"],
        "series_engaged": bool(res.model["series"]["engaged"]),
        "lattice_engaged": bool(res.model["lattice"]["engaged"]),
    }


def _assert_same(got: dict, want: dict) -> None:
    assert got["n_kept"] == want["n_kept"]
    assert np.array_equal(got["kept"], want["kept"]), "the trimmed pair set changed"
    assert got["anchors_engaged"] == want["anchors_engaged"], "the anchor gate flipped"
    assert got["series_engaged"] == want["series_engaged"]
    assert got["lattice_engaged"] == want["lattice_engaged"]
    assert got["gate_reduction"] == pytest.approx(want["gate_reduction"], abs=1e-9)
    np.testing.assert_allclose(got["cal"], want["cal"], rtol=0, atol=1e-9,
                               equal_nan=True, err_msg="Cal_RT_min")
    np.testing.assert_allclose(got["unc"], want["unc"], rtol=0, atol=1e-9,
                               equal_nan=True, err_msg="Cal_RT_uncertainty_min")


@pytest.fixture(scope="module")
def orbitrap_run(orbitrap_samples, orbitrap_standards):
    """Stage 2 requested on the Orbitrap fixture — the dataset whose gate flipped."""
    def run():
        return calibrate(orbitrap_samples, "positive",
                         standards_table=orbitrap_standards, panel="mix15",
                         config=CalibrationConfig.orbitrap(use_sample_anchors=True))
    return run


@pytest.fixture(scope="module")
def qtof_run(qtof_full_samples, qtof_full_standards):
    """Defaults plus stage 2 on QTOF_full (the lattice term engages there)."""
    def run():
        return calibrate(qtof_full_samples, "positive",
                         standards_table=qtof_full_standards, panel="mix15",
                         config=CalibrationConfig(use_sample_anchors=True))
    return run


@pytest.fixture(scope="module")
def orbitrap_baseline(orbitrap_run):
    return _snapshot(orbitrap_run())


@pytest.fixture(scope="module")
def qtof_baseline(qtof_run):
    return _snapshot(qtof_run())


@pytest.mark.parametrize("mode", ["reverse", "shuffle"])
def test_orbitrap_calibration_ignores_tie_order(mode, orbitrap_run, orbitrap_baseline,
                                                monkeypatch):
    monkeypatch.setattr(np, "argsort", _adversarial_argsort(mode))
    _assert_same(_snapshot(orbitrap_run()), orbitrap_baseline)


@pytest.mark.parametrize("mode", ["reverse", "shuffle"])
def test_qtof_calibration_ignores_tie_order(mode, qtof_run, qtof_baseline, monkeypatch):
    monkeypatch.setattr(np, "argsort", _adversarial_argsort(mode))
    _assert_same(_snapshot(qtof_run()), qtof_baseline)


def test_orbitrap_stage2_gate_decision_is_the_deterministic_one(orbitrap_baseline):
    """Pin the decision that used to differ between platforms (threshold 0.2)."""
    assert orbitrap_baseline["anchors_engaged"] is True
    assert orbitrap_baseline["gate_reduction"] == pytest.approx(0.2677, abs=1e-3)
    assert orbitrap_baseline["n_kept"] == 428
