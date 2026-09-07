"""CLI + output writers."""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from rt_anchor import calibrate, write_results
from rt_anchor.cli import main as cli_main
from rt_anchor.config import CalibrationConfig

DATA_OUTPUTS = {"calibrated_csv", "model_json", "anchors_csv", "pairs_csv",
                "landmarks_csv", "log_txt"}


# ------------------------------------------------------------- describe -------

def test_cli_describe_prints_json(orbitrap_samples, capsys):
    rc = cli_main(["describe", orbitrap_samples])
    assert rc == 0
    info = json.loads(capsys.readouterr().out)   # must be valid JSON
    assert info["source_format"] == "masscube"
    assert info["n_features"] > 0


# ----------------------------------------------------------- references -------

def test_cli_references_lists_the_bundle(capsys):
    rc = cli_main(["references"])
    assert rc == 0
    refs = json.loads(capsys.readouterr().out)
    assert any(r["key"] == "col35" and r["available"] for r in refs)


# ------------------------------------------------------------- calibrate ------

def test_cli_calibrate_writes_all_outputs(tmp_path, orbitrap_samples,
                                          orbitrap_standards):
    prefix = str(tmp_path / "run")
    rc = cli_main([
        "calibrate", "--samples", orbitrap_samples,
        "--standards", orbitrap_standards,
        "--polarity", "positive", "--panel", "mix15",
        "--mz-tol-ppm", "8", "--rt-window", "0.3",
        "--no-report", "--out", prefix,
    ])
    assert rc == 0
    for suffix in ("_calibrated.csv", "_model.json", "_anchors.csv",
                   "_pairs.csv", "_landmarks.csv", "_log.txt"):
        f = tmp_path / f"run{suffix}"
        assert f.exists(), f"missing {suffix}"
    df = pd.read_csv(f"{prefix}_calibrated.csv")
    assert {"Cal_RT_min", "iRT", "iRT_reliability", "calibration_scope"} <= set(df.columns)
    model = json.loads((tmp_path / "run_model.json").read_text())
    assert model["method"] == "cross-column-v2"
    assert model["calibration_scope"] == "project"
    assert model["config"]["mz_tol_ppm"] == 8.0
    # the pairs file is the curve's raw evidence, not a summary
    pairs = pd.read_csv(f"{prefix}_pairs.csv")
    assert len(pairs) == model["curve"]["n_pairs"]


def test_cli_panel_none_and_flag_plumbing(tmp_path, orbitrap_samples,
                                          orbitrap_standards):
    prefix = str(tmp_path / "none")
    rc = cli_main([
        "calibrate", "--samples", orbitrap_samples, "--standards", orbitrap_standards,
        "--polarity", "positive", "--panel", "none",
        "--no-sample-anchors", "--no-sample-pairs", "--curve-frac", "0.2",
        "--mz-tol-da", "0.01", "--extrapolate-mode", "clamp",
        "--no-report", "--out", prefix,
    ])
    assert rc == 0
    model = json.loads((tmp_path / "none_model.json").read_text())
    cfg = model["config"]
    assert cfg["use_sample_anchors"] is False and cfg["use_sample_pairs"] is False
    assert cfg["curve_frac"] == 0.2 and cfg["match_mz_tol_da"] == 0.01
    assert cfg["extrapolate_mode"] == "clamp"
    assert model["panel"]["key"] == "none"
    assert "fallback" in model["irt"]["panel_used"]


@pytest.mark.parametrize("flag,arg,replacement", [
    ("--reference", "somefile.csv", "--reference-sample"),
    ("--no-stds-fallback", None, "required input"),
])
def test_cli_removed_flags_exit_2_naming_the_replacement(
        tmp_path, orbitrap_samples, orbitrap_standards, flag, arg, replacement, caplog):
    argv = ["calibrate", "--samples", orbitrap_samples, "--standards", orbitrap_standards,
            "--polarity", "positive", "--out", str(tmp_path / "x"), flag]
    if arg:
        argv.append(arg)
    with caplog.at_level("ERROR"):
        rc = cli_main(argv)
    assert rc == 2
    assert flag in caplog.text and replacement in caplog.text


def test_cli_bad_input_exits_nonzero(tmp_path):
    junk = tmp_path / "junk.txt"
    junk.write_text("[build-system]\nrequires=1\n")
    rc = cli_main([
        "calibrate", "--samples", str(junk), "--standards", str(junk),
        "--polarity", "positive", "--out", str(tmp_path / "bad"),
    ])
    assert rc != 0


def test_cli_bad_input_exits_nonzero_subprocess(tmp_path):
    # exercise the real console-script + exit status end-to-end
    junk = tmp_path / "junk.txt"
    junk.write_text("nope,not,a,table\n1,2,3,4\n")
    proc = subprocess.run(
        [sys.executable, "-m", "rt_anchor.cli", "calibrate",
         "--samples", str(junk), "--standards", str(junk), "--polarity", "positive",
         "--out", str(tmp_path / "bad")],
        capture_output=True, text=True,
    )
    assert proc.returncode != 0
    assert (proc.stderr + proc.stdout).strip() != ""


# --------------------------------------------------------- write_results ------

@pytest.fixture(scope="module")
def result(orbitrap_samples, orbitrap_standards):
    return calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                     panel="mix15", config=CalibrationConfig.orbitrap())


def test_write_results_roundtrip(tmp_path, result):
    prefix = str(tmp_path / "out")
    paths = write_results(result, prefix, report=False)  # data-only roundtrip
    assert set(paths) == DATA_OUTPUTS

    df = pd.read_csv(paths["calibrated_csv"])
    assert len(df) == len(result.table)
    assert set(result.table.columns) <= set(df.columns)

    with open(paths["model_json"]) as fh:
        model = json.load(fh)            # valid JSON, no NaN tokens
    assert model["n_features"] == len(result.table)

    anchors = pd.read_csv(paths["anchors_csv"])
    assert len(anchors) == len(result.anchors)
    landmarks = pd.read_csv(paths["landmarks_csv"])
    assert len(landmarks) == result.irt.n_landmarks


def test_write_results_json_has_no_raw_nan(tmp_path, result):
    prefix = str(tmp_path / "out")
    write_results(result, prefix, report=False)
    raw = (tmp_path / "out_model.json").read_text()
    # bare NaN is invalid JSON; the writer must emit null instead
    assert "NaN" not in raw
    json.loads(raw)


def test_a_failing_report_never_costs_the_data(tmp_path, result, monkeypatch):
    """The numbers are the deliverable; a broken chart is not a reason to lose them."""
    import rt_anchor.viz.report as report_mod

    def boom(*a, **k):
        raise RuntimeError("chart exploded")

    monkeypatch.setattr(report_mod, "write_report", boom)
    paths = write_results(result, str(tmp_path / "r"), report=True)
    assert DATA_OUTPUTS <= set(paths)
    assert "report_html" not in paths
    assert any("report skipped" in line for line in result.log)
