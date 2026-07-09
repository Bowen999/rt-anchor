"""Shared fixtures for the rt_anchor QA suite.

Provides:
* paths to the shipped example inputs and the real LipidScreener files,
* builders that synthesise small-but-valid MS-DIAL / MZmine v2 / MZmine v3
  tables seeded with real panel-standard m/z at plausible RTs (so calibration
  can actually run), and
* a helper to write a minimal MassCube ``m/z``/``RT`` table.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

from rt_anchor.config import CalibrationConfig
from rt_anchor.panel import build_panel

HERE = Path(__file__).resolve().parent
PKG_ROOT = HERE.parent                      # .../rt_anchor
EXAMPLES = PKG_ROOT / "example_input"
PROCESSED = PKG_ROOT.parent / "processed no RT"   # real LipidScreener files


# ----------------------------------------------------------------- paths ------

@pytest.fixture(scope="session")
def examples() -> Path:
    assert EXAMPLES.is_dir(), f"missing example_input at {EXAMPLES}"
    return EXAMPLES


@pytest.fixture(scope="session")
def orbitrap_samples(examples) -> str:
    return str(examples / "Orbitrap" / "samples" / "aligned_feature_table.txt")


@pytest.fixture(scope="session")
def orbitrap_standards(examples) -> str:
    return str(examples / "Orbitrap" / "standards" / "aligned_feature_table.txt")


@pytest.fixture(scope="session")
def qtof_samples(examples) -> str:
    return str(examples / "QTOF" / "samples" / "aligned_feature_table.txt")


@pytest.fixture(scope="session")
def qtof_standards(examples) -> str:
    return str(examples / "QTOF" / "standards" / "aligned_feature_table.txt")


@pytest.fixture(scope="session")
def qtof_full(examples) -> Path:
    return examples / "QTOF_full"


@pytest.fixture(scope="session")
def qtof_full_samples(qtof_full) -> str:
    return str(qtof_full / "samples" / "aligned_feature_table.txt")


@pytest.fixture(scope="session")
def qtof_full_standards(qtof_full) -> str:
    return str(qtof_full / "standards" / "aligned_feature_table.txt")


@pytest.fixture(scope="session")
def qtof_full_single_files(qtof_full):
    d = qtof_full / "samples" / "single_files"
    files = sorted(str(p) for p in d.glob("*.txt"))
    assert files, "no single_files found"
    return files


@pytest.fixture(scope="session")
def lipidscreener_pos() -> str:
    p = PROCESSED / "2_data_cleansing_complete_positive.csv"
    if not p.exists():
        pytest.skip(f"real LipidScreener file not present: {p}")
    return str(p)


@pytest.fixture(scope="session")
def lipidscreener_neg() -> str:
    p = PROCESSED / "2_data_cleansing_complete_negative.csv"
    if not p.exists():
        pytest.skip(f"real LipidScreener file not present: {p}")
    return str(p)


# --------------------------------------------------------------- panel data ---

@pytest.fixture(scope="session")
def pos_targets() -> pd.DataFrame:
    """The 15 positive-mode panel targets (name, mz, rt_ref_min, irt, adduct)."""
    return build_panel("positive", CalibrationConfig()).targets.copy()


# ------------------------------------------------------- synthetic builders ---

def _write(path, header, rows, delim="\t"):
    lines = [delim.join(map(str, header))]
    lines += [delim.join(map(str, r)) for r in rows]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


@pytest.fixture
def make_msdial(tmp_path, pos_targets):
    """Synthesise a realistic MS-DIAL alignment result.

    4 metadata rows (Class / File type / Injection order / blank) whose leading
    descriptor cell sits above the *metadata* columns, a real header with
    ``Alignment ID`` / ``Average Rt(min)`` / ``Average Mz``, then two sample
    columns. Seeded with every panel standard at its reference RT.
    """
    def _build(name="msdial.txt", n_extra=3):
        cols = ["Alignment ID", "Average Rt(min)", "Average Mz",
                "Metabolite name", "Adduct type", "S1", "S2"]
        meta = [
            ["Class", "", "", "", "", "Serum", "QC"],
            ["File type", "", "", "", "", "Sample", "QC"],
            ["Injection order", "", "", "", "", "1", "2"],
            ["", "", "", "", "", "", ""],
        ]
        rows = []
        i = 0
        for _, r in pos_targets.iterrows():
            rows.append([i, f"{r['rt_ref_min']:.3f}", f"{r['mz']:.4f}",
                         r["name"], r["adduct"], 12000 + i, 11000 + i])
            i += 1
        for k in range(n_extra):     # non-anchor filler features
            rows.append([i, f"{2.0 + k:.3f}", "300.1234", "Unknown", "[M+H]+",
                         5000, 4000])
            i += 1
        path = tmp_path / name
        lines = ["\t".join(map(str, m)) for m in meta]
        lines.append("\t".join(cols))
        lines += ["\t".join(map(str, rr)) for rr in rows]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return str(path)
    return _build


@pytest.fixture
def make_mzmine2(tmp_path, pos_targets):
    """Synthesise an MZmine v2 alignment CSV."""
    def _build(name="mzmine2.csv"):
        header = ["row ID", "row m/z", "row retention time",
                  "SampleA.mzML Peak area", "SampleB.mzML Peak area"]
        rows = []
        for i, (_, r) in enumerate(pos_targets.iterrows()):
            rows.append([i, f"{r['mz']:.4f}", f"{r['rt_ref_min']:.3f}", 8000, 7000])
        return _write(tmp_path / name, header, rows, delim=",")
    return _build


@pytest.fixture
def make_mzmine3(tmp_path, pos_targets):
    """Synthesise an MZmine v3 alignment CSV."""
    def _build(name="mzmine3.csv"):
        header = ["id", "mz", "rt",
                  "datafile:SampleA.mzML:area", "datafile:SampleB.mzML:area"]
        rows = []
        for i, (_, r) in enumerate(pos_targets.iterrows()):
            rows.append([i, f"{r['mz']:.4f}", f"{r['rt_ref_min']:.3f}", 8000, 7000])
        return _write(tmp_path / name, header, rows, delim=",")
    return _build


@pytest.fixture
def make_masscube(tmp_path, pos_targets):
    """Minimal MassCube m/z + RT table, optionally with extra columns / rows.

    ``rt_scale`` multiplies the RT (use 60 to emit a seconds table declared in
    minutes, to trigger the range-sanity warning). ``extra_cols`` is a dict of
    column-name -> constant value appended after m/z/RT.
    """
    def _build(name="masscube.txt", extra_rows=None, rt_scale=1.0,
               extra_cols=None, include_targets=True):
        header = ["m/z", "RT"]
        extra_cols = extra_cols or {}
        header += list(extra_cols.keys())
        rows = []
        if include_targets:
            for _, r in pos_targets.iterrows():
                row = [f"{r['mz']:.4f}", f"{r['rt_ref_min'] * rt_scale:.4f}"]
                row += [extra_cols[c] for c in extra_cols]
                rows.append(row)
        for er in (extra_rows or []):
            row = [f"{er[0]:.4f}", f"{er[1] * rt_scale:.4f}"]
            row += [extra_cols[c] for c in extra_cols]
            rows.append(row)
        return _write(tmp_path / name, header, rows, delim="\t")
    return _build
