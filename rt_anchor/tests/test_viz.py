"""Visualisation / report regression tests."""

import re

import numpy as np
import pytest

from rt_anchor import CalibrationConfig, calibrate
from rt_anchor.viz import metrics, performance, report, repeatability, tic


@pytest.fixture(scope="module")
def proj_result(orbitrap_samples, orbitrap_standards):
    return calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                     config=CalibrationConfig.orbitrap())


@pytest.fixture(scope="module")
def sample_result(qtof_full_samples, qtof_full_standards, qtof_full_single_files):
    return calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                     single_files=list(qtof_full_single_files))


def test_kpi_tiles(proj_result):
    tiles = metrics.kpi_tiles(proj_result)
    labels = [t[0] for t in tiles]
    assert len(tiles) == 6 and all(len(t) == 3 for t in tiles)
    # deliberately excluded per redesign
    assert "Calibration" not in labels and "High reliability" not in labels


def test_radar_builds(proj_result):
    axes = metrics.radar_axes(proj_result)
    assert len(axes) == 5 and all(0.0 <= v <= 1.0 for _, v in axes)
    assert len(metrics.radar_plotly(proj_result).data) > 0


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


def test_metrics_covers_panel_and_samples(proj_result):
    m = metrics.compute_metrics(proj_result)
    assert m["n_panel_detected"] >= 1 and m["n_sample_detected"] >= 1
    assert not np.isnan(m["rt_range_panel"][0])          # panel RT range present
    labels = [t[0] for t in metrics.kpi_tiles(proj_result)]
    assert "Standards" in labels and "RT offset" in labels


def test_report_default_on_write_results(proj_result, tmp_path):
    from rt_anchor import write_results
    paths = write_results(proj_result, str(tmp_path / "w"))  # report defaults True
    assert "report_html" in paths and "report_pdf" in paths


def test_repeatability_per_sample_only(proj_result, sample_result):
    # per-sample: applicable, figures build, RI_spread present
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
