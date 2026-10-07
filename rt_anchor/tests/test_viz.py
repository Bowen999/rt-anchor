"""Visualisation / report regression tests.

``viz/`` reads the v2 ``CalibrationResult``: the model blocks (curve, series,
lattice, anchors, irt), the appended result columns and the companion frames.
The report must agree with ``model.json`` by construction — the tile that pairs
the series and lattice terms, the radar axis and the facts rows all read the
``series`` / ``lattice`` blocks rather than recomputing anything, and that is
what the state-matrix tests below pin down. QTOF_full is the shipped dataset on
which the lattice gate really engages (the series gate engages on none of them,
so the series-engaged state is a real result with its block flipped); Orbitrap
is the gated-off run.
"""

from __future__ import annotations

import copy
import dataclasses
import re

import numpy as np
import pytest

from rt_anchor import CalibrationConfig, calibrate
from rt_anchor.viz import metrics, performance, report, repeatability, tic

TILE_LABELS = ["Features", "Matched pairs", "Curve residual", "Series + lattice",
               "iRT landmarks", "Output"]


@pytest.fixture(scope="module")
def proj_result(orbitrap_samples, orbitrap_standards):
    """The 1.2.3 defaults on the gated-off run: both stage-1 terms on and
    declining on Orbitrap, anchors off."""
    return calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                     panel="mix15", config=CalibrationConfig.orbitrap())


@pytest.fixture(scope="module")
def anchors_result(orbitrap_samples, orbitrap_standards):
    """Both terms on and stage-2 anchors requested (the opt-in combination)."""
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
def no_lattice_result(orbitrap_samples, orbitrap_standards):
    """The lattice term disabled by configuration."""
    return calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                     panel="mix15",
                     config=CalibrationConfig.orbitrap(use_lattice_term=False))


@pytest.fixture(scope="module")
def lattice_result(qtof_full_samples, qtof_full_standards):
    """The engaged run: QTOF_full under the defaults (the lattice gate engages,
    the series gate does not)."""
    return calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                     panel="mix15")


@pytest.fixture(scope="module")
def lattice_anchors_result(qtof_full_samples, qtof_full_standards):
    """Anchors requested together with an engaged lattice term: the run that
    puts three gate decisions on the curve figure."""
    return calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                     panel="mix15", config=CalibrationConfig(use_sample_anchors=True))


@pytest.fixture(scope="module")
def sample_result(qtof_full_samples, qtof_full_standards, qtof_full_single_files):
    return calibrate(qtof_full_samples, "positive", standards_table=qtof_full_standards,
                     panel="mix15", single_files=list(qtof_full_single_files))


def _with_block(res, name, **fields):
    """A copy of a real result whose model block ``name`` has ``fields`` overwritten.

    The report reads the model's blocks rather than recomputing, so flipping a
    block exercises a rendering path without inventing an engine state.
    """
    model = copy.deepcopy(res.model)
    model[name].update(fields)
    return dataclasses.replace(res, model=model)


_SERIES_ENGAGED = dict(
    enabled=True, engaged=True, gate_mse_reduction=0.27, gate_threshold=0.2,
    gate_reason="series term engaged: leave-own-out MSE 27% below curve-only "
                "(threshold 20%)",
    n_series=20, n_pairs_covered=106, n_features_corrected=405)

_LATTICE_ENGAGED = dict(
    enabled=True, engaged=True, gate_mse_reduction=0.43, gate_threshold=0.2,
    gate_reason="lattice term engaged: leave-own-out MSE 43% below curve-only "
                "(threshold 20%)",
    n_families=17, n_pairs_in_families=312, n_pairs_covered=270,
    n_features_corrected=300)


def _engaged_result(res, series=True, lattice=False):
    """A copy of a real result whose series (and optionally lattice) block reads
    as engaged.

    No shipped dataset engages the series gate, so the series-engaged state is
    always a flipped block; QTOF_full engages the lattice gate for real, and
    ``lattice=True`` is for the "both engaged" combination.
    """
    out = res
    if series:
        out = _with_block(out, "series", **_SERIES_ENGAGED)
    if lattice:
        out = _with_block(out, "lattice", **_LATTICE_ENGAGED)
    return out


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
    # a measured *loss* (the Orbitrap run: reduction -0.345) is not "off · -35%"
    v, n = metrics.series_headline({**base, "series_reduction": -0.34535720991053154})
    assert v == "off" and n == "35% worse than the curve alone · stage-1 curve only"
    v, n = metrics.series_headline({**base, "series_reduction": 0.0,
                                    "series_reason": "no correction needed: the stage-1 "
                                    "curve already reproduces the series pairs"})
    assert v == "not needed"
    v, n = metrics.series_headline(
        {**base, "series_reason": "only 5 matched pairs sit in a validated homologous "
                                  "series (need 20) — stage-1 curve only"})
    assert v == "off" and "only 5 matched pairs" in n


# ----------------------------------------------------------- lattice headline --

def test_lattice_headline_states():
    base = dict(lattice_enabled=True, lattice_engaged=False, lattice_reduction=np.nan,
                lattice_threshold=0.2, lattice_reason="", n_lattice_families=0,
                n_lattice_pairs=0, n_features_lattice=0, n_features=694)
    assert metrics.lattice_headline({**base, "lattice_enabled": False}) == ("off", "disabled")
    v, n = metrics.lattice_headline({**base, "lattice_engaged": True,
                                     "lattice_reduction": 0.43,
                                     "n_lattice_families": 17, "n_features_lattice": 300})
    assert v == "on · 43%" and n == "17 families · 300 of 694 features corrected"
    v, n = metrics.lattice_headline({**base, "lattice_engaged": True,
                                     "lattice_reduction": 0.43,
                                     "n_lattice_families": 1, "n_features_lattice": 300})
    assert n == "1 family · 300 of 694 features corrected"
    v, n = metrics.lattice_headline({**base, "lattice_reduction": 0.09})
    assert v == "off · 9%" and n == "below the 20% gate · no lattice correction"
    # a measured *loss* is not "off · -14%": it is "off", said plainly
    v, n = metrics.lattice_headline({**base, "lattice_reduction": -0.14})
    assert v == "off" and n == "14% worse than the curve alone · no lattice correction"
    v, n = metrics.lattice_headline({**base, "lattice_reduction": 0.0,
                                     "lattice_reason": "no correction needed: the stage-1 "
                                     "curve already reproduces the lattice pairs"})
    assert v == "not needed"
    reason = ("only 6 matched pairs the series term leaves alone sit in a lattice "
              "family (need 20) — no lattice correction")
    v, n = metrics.lattice_headline({**base, "lattice_reason": reason})
    assert (v, n) == ("off", reason)


def test_lattice_headline_reads_the_model_block(lattice_result, proj_result):
    # the engaged run (QTOF_full): the headline carries the block's own numbers
    m = metrics.compute_metrics(lattice_result)
    blk = lattice_result.model["lattice"]
    assert m["lattice_engaged"] and m["lattice_enabled"]
    assert m["lattice_reduction"] == pytest.approx(blk["gate_mse_reduction"])
    assert m["lattice_threshold"] == blk["gate_threshold"]
    assert m["lattice_reason"] == blk["gate_reason"]
    assert m["n_lattice_families"] == blk["n_families"]
    assert m["n_lattice_pairs"] == blk["n_pairs_covered"]
    assert m["n_features_lattice"] == blk["n_features_corrected"] > 0
    v, n = metrics.lattice_headline(m)
    assert v == f"on · {100 * blk['gate_mse_reduction']:.0f}%"
    assert n == (f"{blk['n_families']} families · {blk['n_features_corrected']} of "
                 f"{m['n_features']} features corrected")
    # the gated-off run (Orbitrap): "off", and the note says no correction was made
    m = metrics.compute_metrics(proj_result)
    assert not m["lattice_engaged"]
    v, n = metrics.lattice_headline(m)
    assert v.startswith("off") and "no lattice correction" in n


def test_run_level_warp_source_names_every_engaged_stage(lattice_result, proj_result):
    assert metrics.compute_metrics(proj_result)["warp_source"] == "curve"
    lattice_only = _with_block(lattice_result, "series", engaged=False)
    assert metrics.compute_metrics(lattice_only)["warp_source"] == "curve+lattice"
    both = _engaged_result(lattice_result, series=True, lattice=True)
    assert metrics.compute_metrics(both)["warp_source"] == "curve+series+lattice"
    both.model["anchors"]["engaged"] = True
    assert metrics.compute_metrics(both)["warp_source"] == "curve+series+lattice+anchors"


# ------------------------------------------------------------------- tiles ------

def test_kpi_tiles(proj_result):
    tiles = metrics.kpi_tiles(proj_result)
    assert len(tiles) == 6 and all(len(t) == 3 for t in tiles)
    assert [t[0] for t in tiles] == TILE_LABELS
    # the fourth tile stacks the two stage-1 terms, one headline each; on
    # Orbitrap both gates decline
    label, value, note = tiles[3]
    assert label == "Series + lattice"
    assert [k for k, _ in value] == ["series", "lattice"]
    assert all(v.startswith("off") for _, v in value)
    assert note == f"0 + 0 of {len(proj_result.table):,} features corrected"


def test_terms_tile_carries_the_anchor_state_when_requested(anchors_result):
    label, value, note = metrics.kpi_tiles(anchors_result)[3]
    assert label == "Series + lattice"
    # stage 2 was requested and gated off on Orbitrap: the note says so
    assert note.endswith(" · anchors off")


def test_terms_tile_when_a_term_is_disabled(no_series_result, no_lattice_result):
    # one term disabled: its row reads "off" (the headline's note says why)
    label, value, _ = metrics.kpi_tiles(no_series_result)[3]
    assert label == "Series + lattice"
    assert value[0] == ("series", "off") and value[1][0] == "lattice"
    assert metrics.series_headline(metrics.compute_metrics(no_series_result)) == ("off", "disabled")
    _, value, _ = metrics.kpi_tiles(no_lattice_result)[3]
    assert value[1] == ("lattice", "off") and value[0][0] == "series"
    assert metrics.lattice_headline(metrics.compute_metrics(no_lattice_result)) == ("off", "disabled")
    # both disabled
    both_off = _with_block(no_lattice_result, "series", enabled=False, engaged=False,
                           gate_mse_reduction=None,
                           gate_reason="disabled (use_series_term=False)")
    assert metrics.kpi_tiles(both_off)[3] == (
        "Series + lattice", [("series", "off"), ("lattice", "off")],
        f"0 + 0 of {len(both_off.table):,} features corrected")


def test_terms_tile_when_the_series_term_is_engaged(proj_result):
    label, value, note = metrics.kpi_tiles(_engaged_result(proj_result))[3]
    assert label == "Series + lattice"
    assert value[0] == ("series", "on · 27%") and value[1][0] == "lattice"
    assert note == f"405 + 0 of {len(proj_result.table):,} features corrected"


def test_terms_tile_when_the_lattice_term_is_engaged(lattice_result):
    # a real engaged run (QTOF_full): the lattice row says "on" with the gate's
    # own figure, the series row is whatever the series gate decided
    label, value, note = metrics.kpi_tiles(lattice_result)[3]
    ser, blk = lattice_result.model["series"], lattice_result.model["lattice"]
    n = len(lattice_result.table)
    assert label == "Series + lattice"
    assert value[1] == ("lattice", f"on · {100 * blk['gate_mse_reduction']:.0f}%")
    assert value[0] == ("series",
                        metrics.series_headline(metrics.compute_metrics(lattice_result))[0])
    # the tile's counts are the table's: rows that carry a correction from each term
    k_series = int((lattice_result.values("series_n_members") > 0).sum())
    k = int((lattice_result.values("lattice_n_members") > 0).sum())
    assert k == blk["n_features_corrected"] > 0
    assert k_series == ser["n_features_corrected"]
    assert note == f"{k_series:,} + {k:,} of {n:,} features corrected"


def test_terms_tile_when_both_terms_are_engaged_and_anchors_requested(anchors_result):
    both = _engaged_result(anchors_result, series=True, lattice=True)
    label, value, note = metrics.kpi_tiles(both)[3]
    assert value == [("series", "on · 27%"), ("lattice", "on · 43%")]
    assert note == f"405 + 300 of {len(both.table):,} features corrected · anchors off"
    both.model["anchors"].update(engaged=True, gate_mse_reduction=0.72)
    assert metrics.kpi_tiles(both)[3][2].endswith(" · anchors on 72%")


def test_the_pdf_cover_and_the_html_draw_the_stacked_tile(lattice_result):
    # both tile renderers take a list of (caption, value) rows: the two rows of
    # the fourth tile must reach each of them
    html = report.build_html(lattice_result)
    assert ">Series + lattice<" in html
    assert re.search(r'class="kpi-k">series</span><span class="kpi-n">[^<]+</span>', html)
    assert re.search(r'class="kpi-k">lattice</span><span class="kpi-n">on · \d+%</span>', html)
    import matplotlib.pyplot as plt
    fig = report._cover_page(lattice_result)
    texts = [t.get_text() for t in fig.texts]
    plt.close(fig)
    assert "SERIES + LATTICE" in texts and "series" in texts and "lattice" in texts


# ------------------------------------------------------------------ radar ------

def test_radar_builds(proj_result):
    axes = metrics.radar_axes(proj_result)
    assert len(axes) == 5 and all(0.0 <= v <= 1.0 for _, v in axes)
    assert [a for a, _ in axes] == ["Pair support", "Curve fit", "Series + lattice",
                                    "Landmarks", "Reliability"]
    assert len(metrics.radar_plotly(proj_result).data) > 0


def test_radar_terms_axis(proj_result, lattice_result):
    # both gates declined -> 0
    axes = dict(metrics.radar_axes(proj_result))
    assert axes["Series + lattice"] == 0.0
    # the engaged run (QTOF_full): the gate's own reduction (the larger one, if
    # the series gate engaged as well)
    axes = dict(metrics.radar_axes(lattice_result))
    gates = [lattice_result.model[k]["gate_mse_reduction"] for k in ("series", "lattice")
             if lattice_result.model[k]["engaged"]]
    assert lattice_result.model["lattice"]["engaged"]
    assert axes["Series + lattice"] == pytest.approx(max(gates))
    assert 0.0 < axes["Series + lattice"] < 1.0
    # series engaged alone -> the series gate's reduction
    axes = dict(metrics.radar_axes(_engaged_result(proj_result)))
    assert axes["Series + lattice"] == pytest.approx(0.27)


def test_radar_terms_axis_takes_the_larger_engaged_gate(proj_result):
    def axis(res):
        return dict(metrics.radar_axes(res))["Series + lattice"]

    both = _engaged_result(proj_result, series=True, lattice=True)      # 0.27 and 0.43
    assert axis(both) == pytest.approx(0.43)
    both = _with_block(both, "series", gate_mse_reduction=0.61)
    assert axis(both) == pytest.approx(0.61)
    # a gate that declined does not count, however large its measured figure
    declined = _with_block(both, "lattice", engaged=False, gate_mse_reduction=0.19)
    assert axis(declined) == pytest.approx(0.61)
    declined = _with_block(declined, "series", engaged=False)
    assert axis(declined) == 0.0
    # clipped to [0, 1]
    assert axis(_with_block(both, "series", gate_mse_reduction=1.3)) == 1.0
    assert axis(_with_block(both, "lattice", gate_mse_reduction=-0.2,
                            engaged=True)) == pytest.approx(0.61)


# ------------------------------------------------------------- method facts ----

def test_method_facts_term_rows(proj_result, lattice_result):
    facts = report.method_facts(proj_result)
    keys = [k for k, _ in facts]
    # the lattice row sits right after the series row, before stage 2
    assert keys.index("Lattice term") == keys.index("Series term") + 1
    assert keys.index("Lattice term") < keys.index("Stage 2 anchors")
    d = dict(facts)
    assert d["Series term"] == proj_result.model["series"]["gate_reason"]
    assert d["Lattice term"] == proj_result.model["lattice"]["gate_reason"]
    # anchors off by configuration: the stage-2 row says "not requested"
    assert d["Stage 2 anchors"].startswith("not requested")
    # the engaged run says so in the model block's own words
    d = dict(report.method_facts(lattice_result))
    assert d["Lattice term"] == lattice_result.model["lattice"]["gate_reason"]
    assert d["Lattice term"].startswith("lattice term engaged")


def test_method_facts_lattice_row_when_disabled_or_missing(no_lattice_result):
    d = dict(report.method_facts(no_lattice_result))
    assert d["Lattice term"] == "disabled (use_lattice_term=False)"
    # a model with no lattice block at all (not produced by this engine) says
    # "not fitted" rather than failing
    model = copy.deepcopy(no_lattice_result.model)
    del model["lattice"]
    bare = dataclasses.replace(no_lattice_result, model=model)
    assert dict(report.method_facts(bare))["Lattice term"] == "not fitted"


def test_method_facts_anchors_requested(anchors_result):
    d = dict(report.method_facts(anchors_result))
    assert not d["Stage 2 anchors"].startswith("not requested")


# ------------------------------------------------------ the curve figure's notes -

def _curve_note_geometry(res):
    """``(note lines, title bottom, note top, note bottom, axes top)`` of the
    curve figure, in figure fractions, as the renderer lays it out."""
    import matplotlib.pyplot as plt
    fig = performance.figure_mpl(res)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    height = fig.get_figheight() * fig.dpi
    title_text = next(t for t in fig.texts if t.get_text().startswith("Cross-column"))
    note = next(t for t in fig.texts if t is not title_text)
    title = title_text.get_window_extent(renderer)
    box = note.get_window_extent(renderer)
    axes_top = fig.axes[0].get_position().y1
    lines = note.get_text().split("\n")
    plt.close(fig)
    return lines, title.y0 / height, box.y1 / height, box.y0 / height, axes_top


def test_curve_data_carries_both_terms_decisions(lattice_result, proj_result):
    d = performance.compute_curve(lattice_result)
    blk = lattice_result.model["lattice"]
    assert d["lattice_engaged"] is True and d["lattice_reason"] == blk["gate_reason"]
    ser = lattice_result.model["series"]
    assert d["series_engaged"] is bool(ser["engaged"]) and d["series_reason"] == ser["gate_reason"]
    d = performance.compute_curve(proj_result)
    assert d["lattice_engaged"] is False
    assert d["lattice_reason"] == proj_result.model["lattice"]["gate_reason"]


def test_curve_figure_prints_both_terms_reasons_without_stuttering(lattice_result,
                                                                   no_lattice_result):
    lines, *_ = _curve_note_geometry(lattice_result)
    # series first, lattice on the next line, each named once
    assert lines[0].startswith("series term")
    assert lines[1].startswith("lattice term engaged")
    assert not any(line.startswith(("series term: series", "lattice term: lattice"))
                   for line in lines)
    # a reason that does not open with the term's name gets the name in front
    lines, *_ = _curve_note_geometry(no_lattice_result)
    assert lines[1] == "lattice term: disabled (use_lattice_term=False)"


def test_curve_figure_notes_clear_the_title_and_the_plot(lattice_result,
                                                         lattice_anchors_result,
                                                         proj_result):
    """Two or three gate decisions hang under the title; none may touch either."""
    for res, n_lines in ((lattice_result, 2), (lattice_anchors_result, 3),
                         (proj_result, 2)):
        lines, title_y0, note_y1, note_y0, axes_top = _curve_note_geometry(res)
        assert len(lines) == n_lines
        assert note_y1 < title_y0, "the notes run into the title"
        assert note_y0 > axes_top, "the notes run into the plot"


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


@pytest.mark.parametrize("state", [
    "gated_off", "series_disabled", "lattice_disabled", "series_engaged",
    "lattice_engaged", "both_engaged", "anchors_on", "lattice_and_anchors", "per_sample"])
def test_report_renders_in_every_term_state(state, proj_result, anchors_result,
                                            no_series_result, no_lattice_result,
                                            lattice_result, lattice_anchors_result,
                                            sample_result, tmp_path):
    """HTML + PDF for: both terms gated off (Orbitrap) / each term disabled /
    series engaged / lattice engaged (QTOF_full, for real) / both engaged /
    anchors requested with both terms / the per-sample tier."""
    res = {"gated_off": proj_result,
           "series_disabled": no_series_result,
           "lattice_disabled": no_lattice_result,
           "series_engaged": _engaged_result(proj_result),
           "lattice_engaged": lattice_result,
           "both_engaged": _engaged_result(lattice_result, series=True, lattice=True),
           "anchors_on": anchors_result,
           "lattice_and_anchors": lattice_anchors_result,
           "per_sample": sample_result}[state]
    paths = report.write_report(res, str(tmp_path / state))
    assert {"report_html", "report_pdf"} <= set(paths)
    html = open(paths["report_html"], encoding="utf-8").read()
    # the tile, the facts rows and the radar axis all made it into the page
    assert ">Series + lattice<" in html
    assert ">Lattice term<" in html and ">Series term<" in html
    assert open(paths["report_pdf"], "rb").read(5) == b"%PDF-"


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
