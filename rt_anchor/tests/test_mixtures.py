"""Diacylglycerols are positive-mode standards only.

DGs ionise very poorly in negative ESI (the [M+HCOO]- adduct is >1000x weaker
than [M+NH4]+ in positive), so they must never be resolved as negative-mode
standards or iRT landmarks, while positive mode keeps both of them.
"""

from __future__ import annotations

import pytest

from rt_anchor import mixtures
from rt_anchor.config import CalibrationConfig
from rt_anchor.irt import landmark_panel
from rt_anchor.panel import DEFAULT_MANIFEST, build_panel

DGS = {"DG 16:0_18:0", "DG 18:0_18:0"}
PANELS = ["mix15", "mix21"]


def _manifest(key):
    return mixtures.get_manifest(key)


@pytest.mark.parametrize("key", PANELS)
def test_manifest_marks_dg_negative_as_nd(key):
    man = _manifest(key)
    dg = man[man["class"] == "DG"]
    assert set(dg["name"]) == DGS
    assert (dg["adduct_neg"] == "ND").all()
    assert (dg["adduct_pos"] == "[M+NH4]+").all()


def test_mix15_is_the_default_manifest():
    dg = DEFAULT_MANIFEST[DEFAULT_MANIFEST["class"] == "DG"]
    assert (dg["adduct_neg"] == "ND").all()


@pytest.mark.parametrize("key", PANELS)
def test_build_panel_negative_has_no_dg(key):
    panel = build_panel("negative", CalibrationConfig(), manifest=_manifest(key))
    assert len(panel) > 0
    assert not (set(panel.targets["name"]) & DGS)
    assert "DG" not in set(panel.targets["class"])


@pytest.mark.parametrize("key", PANELS)
def test_build_panel_positive_keeps_both_dgs_as_ammonium(key):
    panel = build_panel("positive", CalibrationConfig(), manifest=_manifest(key))
    t = panel.targets
    dg = t[t["name"].isin(DGS)]
    assert set(dg["name"]) == DGS
    assert (dg["adduct"] == "[M+NH4]+").all()


@pytest.mark.parametrize("key", PANELS)
def test_landmark_panel_negative_has_no_dg(key):
    panel = landmark_panel(_manifest(key), "negative")
    assert len(panel) > 0
    assert not (set(panel.targets["name"]) & DGS)


@pytest.mark.parametrize("key", PANELS)
def test_landmark_panel_positive_keeps_both_dgs(key):
    panel = landmark_panel(_manifest(key), "positive")
    t = panel.targets
    dg = t[t["name"].isin(DGS)]
    assert set(dg["name"]) == DGS
    assert (dg["adduct"] == "[M+NH4]+").all()


@pytest.mark.parametrize("key", PANELS)
def test_preview_payload_dg_has_positive_mz_only(key):
    card = next(c for c in mixtures.preview_payload([key]) if c["key"] == key)
    dg = [r for r in card["rows"] if r["class"] == "DG"]
    assert {r["name"] for r in dg} == DGS
    for r in dg:
        assert r["mz_pos"] is not None
        assert r["mz_neg"] is None
    # the other classes are untouched: PCs still list a negative m/z
    pc = [r for r in card["rows"] if r["class"] == "PC"]
    assert pc and all(r["mz_neg"] is not None for r in pc)
