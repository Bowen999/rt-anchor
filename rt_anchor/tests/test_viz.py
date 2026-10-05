"""Visualisation / report regression tests.

``viz/`` reads the v2 ``CalibrationResult``: the model blocks (curve, series,
anchors, irt), the appended result columns and the companion frames. The report
must agree with ``model.json`` by construction — the series-term tile, radar
axis and facts row all read the ``series`` block rather than recomputing
anything, and that is what the state-matrix tests below pin down.
"""

from __future__ import annotations

import copy
import dataclasses
import re

import numpy as np
import pytest

from rt_anchor import CalibrationConfig, calibrate
from rt_anchor.viz import metrics, performance, report, repeatability, tic

TILE_LABELS = ["Features", "Matched pairs", "Curve residual", "Series term",
               "iRT landmarks", "Output"]


@pytest.fixture(scope="module")
def proj_result(orbitrap_samples, orbitrap_standards):
    """The 1.2.2 defaults: series term on (gated off on Orbitrap), anchors off."""
    return calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                     panel="mix15", config=CalibrationConfig.orbitrap())


@pytest.fixture(scope="module")
def anchors_result(orbitrap_samples, orbitrap_standards):
    """Series term on and stage-2 anchors requested (the opt-in combination)."""
    return calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                     panel="mix15",
                     config=CalibrationConfig.orbitrap(use_sample_anchors=True))


@pytest.fixture(scope="module")
def no_series_result(orbitrap_samples, orbitrap_standards):
    """The series term disabled by configuration."""
    return calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                     panel="mix15",
                     config=CalibrationConfig.orbitrap(use_series_term=False))


@pytest.fixture(scope="module")
def sample_result(qtof_full_samples, qtof_full_standards, qtof_full_single_files):
    return calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                     panel="mix15", single_files=list(qtof_full_single_files))


def _engaged_result(res):
    """A copy of a real result whose series block reads as engaged.

    The report reads the model's ``series`` block rather than recomputing, so
    flipping the block exercises the engaged rendering path without inventing
    an engine state (none of the shipped datasets engages the gate).
    """
    model = copy.deepcopy(res.model)
    model["series"].update(
        enabled=True, engaged=True, gate_mse_reduction=0.27, gate_threshold=0.2,
        gate_reason="series term engaged: leave-own-out MSE 27% below curve-only "
                    "(threshold 20%)",
        n_series=20, n_pairs_covered=106, n_features_corrected=405)
    return dataclasses.replace(res, model=model)


# ------------------------------------------------------------- series headline -

def test_series_headline_states():
    base = dict(series_enabled=True, series_engaged=False, series_reduction=np.nan,
                series_threshold=0.2, series_reason="", n_series=0,
                n_series_pairs=0, n_features_series=0, n_features=405)
    assert metrics.series_headline({**base, "series_enabled": False}) == ("off", "disabled")
    v, n = metrics.series_headline({**base, "series_engaged": True,
                                    "series_reduction": 0.27,
                                    "n_series": 20, "n_features_series": 106})
    assert v == "on · 27%" and n == "20 series · 106 of 405 features corrected"
    v, n = metrics.series_headline({**base, "series_reduction": 0.09})
    assert v == "off · 9%" and n == "below the 20% gate · stage-1 curve only"
    v, n = metrics.series_headline({**base, "series_reduction": 0.0,
                                    "series_reason": "no correction needed: the stage-1 "
                                    "curve already reproduces the series pairs"})
    assert v == "not needed"
    v, n = metrics.series_headline(
        {**base, "series_reason": "only 5 matched pairs sit in a validated homologous "
                                  "series (need 20) — stage-1 curve only"})
    assert v == "off" and "only 5 matched pairs" in n


def test_kpi_tiles(proj_result):
    tiles = metrics.kpi_tiles(proj_result)
    assert len(tiles) == 6 and all(len(t) == 3 for t in tiles)
    assert [t[0] for t in tiles] == TILE_LABELS
    # the fourth tile is the series term; on Orbitrap the gate declines
    label, value, note = tiles[3]
    assert label == "Series term"
    assert value.startswith("off") and "stage-1 curve only" in note


def test_series_tile_carries_the_anchor_state_when_requested(anchors_result):
    label, value, note = metrics.kpi_tiles(anchors_result)[3]
    assert label == "Series term"
    # stage 2 was requested and gated off on Orbitrap: the note says so
    assert "anchors off" in note


def test_series_tile_when_disabled(no_series_result):
    assert metrics.kpi_tiles(no_series_result)[3] == ("Series term", "off", "disabled")


def test_series_tile_when_engaged(proj_result):
    label, value, note = metrics.kpi_tiles(_engaged_result(proj_result))[3]
    assert (label, value) == ("Series term", "on · 27%")
    assert "20 series · 405 of" in note


# ------------------------------------------------------------------ radar ------

def test_radar_builds(proj_result):
    axes = metrics.radar_axes(proj_result)
    assert len(axes) == 5 and all(0.0 <= v <= 1.0 for _, v in axes)
    assert [a for a, _ in axes] == ["Pair support", "Curve fit", "Series term",
                                    "Landmarks", "Reliability"]
    assert len(metrics.radar_plotly(proj_result).data) > 0


def test_radar_series_axis(proj_result):
    # gated off -> 0; engaged -> the clipped MSE reduction
    axes = dict(metrics.radar_axes(proj_result))
    assert axes["Series term"] == 0.0
    axes = dict(metrics.radar_axes(_engaged_result(proj_result)))
    assert axes["Series term"] == pytest.approx(0.27)


# ------------------------------------------------------------- method facts ----

def test_method_facts_series_row(proj_result):
    facts = report.method_facts(proj_result)
    keys = [k for k, _ in facts]
    assert keys.index("Series term") < keys.index("Stage 2 anchors")
    d = dict(facts)
    assert d["Series term"] == proj_result.model["series"]["gate_reason"]
    # anchors off by configuration: the stage-2 row says "not requested"
    assert d["Stage 2 anchors"].startswith("not requested")


def test_method_facts_anchors_requested(anchors_result):
    d = dict(report.method_facts(anchors_result))
    assert not d["Stage 2 anchors"].startswith("not requested")


# ------------------------------------------------------------ the render matrix -

def test_all_mpl_figures_build(proj_result):
    import matplotlib.pyplot as plt
    for f in (metrics.figure_mpl, lambda r: tic.figure_mpl(r, "clean"),
              lambda r: tic.figure_mpl(r, "realistic"), performance.figure_mpl):
        assert f(proj_result) is not None
    plt.close("all")


def test_all_plotly_figures_build(proj_result):
    for f in (metrics.figure_plotly, lambda r: tic.figure_plotly(r, "clean"),
              performance.figure_plotly):
        fig = f(proj_result)
        assert fig is not None and len(fig.data) > 0


def test_figures_build_per_sample(sample_result):
    assert sample_result.model["calibration_scope"] == "sample"
    assert metrics.figure_mpl(sample_result) is not None
    assert tic.figure_plotly(sample_result, "clean") is not None


@pytest.mark.parametrize("state", ["gated_off", "disabled", "engaged", "anchors_on"])
def test_report_renders_in_every_series_state(state, proj_result, anchors_result,
                                              no_series_result, tmp_path):
    """HTML + PDF for: series engaged / gated off / disabled / anchors requested."""
    res = {"gated_off": proj_result, "disabled": no_series_result,
           "engaged": _engaged_result(proj_result),
           "anchors_on": anchors_result}[state]
    paths = report.write_report(res, str(tmp_path / state))
    assert {"report_html", "report_pdf"} <= set(paths)


def test_report_html_is_self_contained(proj_result, tmp_path):
    paths = report.write_report(proj_result, str(tmp_path / "r"), tic_style="both")
    assert {"report_html", "report_pdf"} <= set(paths)
    html = open(paths["report_html"], encoding="utf-8").read()
    assert "plotly" in html.lower()
    # no external <script src> / <link href> that would fetch over the network
    assert not re.search(r'<(script[^>]*src|link[^>]*href)="https?', html)
    assert open(paths["report_pdf"], "rb").read(5) == b"%PDF-"


def test_tic_styles_differ(proj_result):
    d_clean = tic.compute_tic(proj_result, "clean")
    d_real = tic.compute_tic(proj_result, "realistic")
    assert not np.allclose(d_clean["yb"], d_real["yb"])


def test_report_default_on_write_results(proj_result, tmp_path):
    from rt_anchor import write_results
    paths = write_results(proj_result, str(tmp_path / "w"))  # report defaults True
    assert "report_html" in paths and "report_pdf" in paths


def test_repeatability_per_sample_only(proj_result, sample_result):
    assert repeatability.is_applicable(sample_result)
    d = repeatability.compute_repro(sample_result)
    assert d["n"] > 0 and d["median"] >= 0
    assert repeatability.figure_mpl(sample_result) is not None
    assert len(repeatability.figure_plotly(sample_result).data) > 0
    # per-project: not applicable (RI_spread is NaN)
    assert not repeatability.is_applicable(proj_result)


def test_repeatability_in_per_sample_report(sample_result, tmp_path):
    html = report.build_html(sample_result)
    assert "fig_repeatability" in html
