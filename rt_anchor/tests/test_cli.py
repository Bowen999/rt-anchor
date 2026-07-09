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


# ------------------------------------------------------------- describe -------

def test_cli_describe_prints_json(orbitrap_samples, capsys):
    rc = cli_main(["describe", orbitrap_samples])
    assert rc == 0
    out = capsys.readouterr().out
    info = json.loads(out)               # must be valid JSON
    assert info["source_format"] == "masscube"
    assert info["n_features"] > 0


# ------------------------------------------------------------- calibrate ------

def test_cli_calibrate_writes_all_outputs(tmp_path, orbitrap_samples,
                                          orbitrap_standards):
    prefix = str(tmp_path / "run")
    rc = cli_main([
        "calibrate", "--samples", orbitrap_samples,
        "--standards", orbitrap_standards,
        "--polarity", "positive", "--mz-tol-ppm", "8", "--rt-window", "0.3",
        "--out", prefix,
    ])
    assert rc == 0
    for suffix in ("_calibrated.csv", "_model.json", "_anchors.csv", "_log.txt"):
        f = tmp_path / f"run{suffix}"
        assert f.exists() and f.stat().st_size > 0, f"missing/empty {suffix}"

    # calibrated CSV reloads with all original columns + result columns
    df = pd.read_csv(f"{prefix}_calibrated.csv")
    assert "RI" in df.columns and "calibration_scope" in df.columns
    # model JSON parses
    model = json.loads((tmp_path / "run_model.json").read_text())
    assert model["calibration_scope"] == "project"
    assert model["n_anchors_used"] >= 6


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

def test_write_results_roundtrip(tmp_path, orbitrap_samples, orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive",
                    standards_table=orbitrap_standards,
                    config=CalibrationConfig.orbitrap())
    prefix = str(tmp_path / "out")
    paths = write_results(res, prefix, report=False)  # data-only roundtrip
    assert set(paths) == {"calibrated_csv", "model_json", "anchors_csv", "log_txt"}

    df = pd.read_csv(paths["calibrated_csv"])
    assert len(df) == len(res.table)
    assert set(res.table.columns) <= set(df.columns)

    with open(paths["model_json"]) as fh:
        model = json.load(fh)            # valid JSON, no NaN tokens
    assert model["n_features"] == len(res.table)

    anchors = pd.read_csv(paths["anchors_csv"])
    assert len(anchors) == res.model["n_anchors_used"]


def test_write_results_json_has_no_raw_nan(tmp_path, orbitrap_samples,
                                           orbitrap_standards):
    res = calibrate(orbitrap_samples, "positive",
                    standards_table=orbitrap_standards,
                    config=CalibrationConfig.orbitrap())
    prefix = str(tmp_path / "out")
    paths = write_results(res, prefix)
    raw = (tmp_path / "out_model.json").read_text()
    # bare NaN is invalid JSON; the writer must emit null instead
    assert "NaN" not in raw
    json.loads(raw)
