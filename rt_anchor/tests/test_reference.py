"""Reference-dataset resolution and the bundled standard mixtures.

The reference pair is the second axis of every v2 calibration, so which pair was
used has to be unambiguous, recorded, and impossible to half-specify.
"""

from __future__ import annotations

import os

import pandas as pd
import pytest

from rt_anchor import mixtures
from rt_anchor.config import CalibrationConfig
from rt_anchor.errors import ConfigError
from rt_anchor.panel import build_panel
from rt_anchor.reference import (
    DEFAULT_REFERENCE_KEY,
    available_references,
    load_reference,
    resolve_reference,
)


# ------------------------------------------------------- reference resolution -

def test_default_reference_is_the_bundled_col35():
    pair = resolve_reference()
    assert pair.key == DEFAULT_REFERENCE_KEY == "col35"
    assert pair.is_default
    assert os.path.isfile(pair.sample_path) and os.path.isfile(pair.standards_path)
    assert pair.provenance, "the bundled pair must carry its provenance"


def test_reference_block_is_serialisable(bundled_reference):
    d = bundled_reference.to_dict()
    assert set(d) == {"key", "label", "sample", "standards", "is_default", "provenance"}
    assert d["is_default"] is True


def test_half_specified_reference_is_refused(bundled_reference):
    with pytest.raises(ConfigError) as ei:
        resolve_reference(sample=bundled_reference.sample_path)
    assert "together" in str(ei.value)
    with pytest.raises(ConfigError):
        resolve_reference(standards=bundled_reference.standards_path)


def test_custom_reference_pair_is_marked_not_default(bundled_reference):
    pair = resolve_reference(sample=bundled_reference.sample_path,
                             standards=bundled_reference.standards_path)
    assert pair.key == "custom" and not pair.is_default
    assert "not comparable" in pair.provenance.get("note", "")


def test_missing_custom_file_is_refused(bundled_reference, tmp_path):
    with pytest.raises(ConfigError):
        resolve_reference(sample=str(tmp_path / "nope.csv"),
                          standards=bundled_reference.standards_path)


def test_unknown_reference_key_is_refused():
    with pytest.raises(ConfigError) as ei:
        resolve_reference(key="col99")
    assert "col99" in str(ei.value)


def test_available_references_lists_the_bundle():
    refs = available_references()
    keys = {r["key"] for r in refs}
    assert "col35" in keys
    col35 = next(r for r in refs if r["key"] == "col35")
    assert col35["available"] is True and col35["label"]


def test_load_reference_parses_both_runs():
    ref = load_reference(polarity="positive")
    assert ref.sample.n_features() > 100 and ref.standards.n_features() > 100
    # the §3 loader fix: Average / Stdev are not sample columns
    assert len(ref.sample.sample_cols) == 4
    assert all("average" not in c.lower() and "stdev" not in c.lower()
               for c in ref.sample.sample_cols)


def test_calibrating_against_a_custom_reference_pair_still_runs(
        orbitrap_samples, orbitrap_standards, bundled_reference):
    """Passing the bundled pair explicitly must equal using it implicitly."""
    from rt_anchor import calibrate

    cfg = CalibrationConfig.orbitrap()
    a = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                  panel="mix15", config=cfg)
    b = calibrate(orbitrap_samples, "positive", standards_table=orbitrap_standards,
                  panel="mix15", config=cfg,
                  reference_sample=bundled_reference.sample_path,
                  reference_standards=bundled_reference.standards_path)
    assert b.model["reference"]["key"] == "custom"
    assert b.model["reference"]["is_default"] is False
    pd.testing.assert_series_equal(a.values("Cal_RT_min"), b.values("Cal_RT_min"))


# ------------------------------------------------------------------ mixtures --

def test_offered_mixtures_are_the_three_cards():
    assert mixtures.offered_keys() == ["mix15", "mix21", "none"]
    assert mixtures.DEFAULT_MIXTURE == "mix21"
    # mix21lpc stays registered for reproducibility but is not offered
    assert "mix21lpc" in mixtures.mixture_keys()
    assert "mix21lpc" not in mixtures.offered_keys()


def test_none_has_no_manifest_rather_than_an_empty_one():
    assert mixtures.get_manifest("none") is None
    assert mixtures.n_standards("none") == 0


@pytest.mark.parametrize("key,n", [("mix15", 15), ("mix21", 21)])
def test_manifests_build_a_positive_panel(key, n):
    man = mixtures.get_manifest(key)
    assert len(man) == n
    panel = build_panel("positive", CalibrationConfig(), manifest=man)
    assert 0 < len(panel) <= n


def test_build_panel_no_longer_stamps_a_v1_irt():
    """§7: Panel.targets.irt is NaN — the v2 scale comes from the landmarks."""
    targets = build_panel("positive", CalibrationConfig()).targets
    assert "irt" in targets.columns
    assert targets["irt"].isna().all()


def test_preview_payload_is_json_safe():
    payload = mixtures.preview_payload()
    assert [p["key"] for p in payload] == mixtures.offered_keys()
    for p in payload:
        assert {"key", "label", "sub", "n_standards", "rows"} <= set(p)
