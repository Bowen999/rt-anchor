"""rt_anchor — cross-column retention-time calibration for LC-MS lipidomics.

Unifies feature tables from MS-DIAL, MZmine, MassCube and LipidScreener, then
maps each feature's retention time onto a **reference column's** time axis
(``Cal_RT_min``) and onto a dimensionless 1-100 retention index (``iRT``).

The method needs two runs from your column — the sample and a standards run on
the same gradient — and a reference pair, which defaults to the bundled Column 25
reference. Stage 1 matches features between the runs by accurate m/z alone (no
identities, so an unknown or undetectable mixture costs nothing) and fits a
robust monotone curve. Stage 2 refines it with endogenous plasma lipids found in
both sample runs, and only if that provably helps.

Quick start
-----------
>>> from rt_anchor import calibrate, write_results
>>> res = calibrate("samples.txt", polarity="positive",
...                 standards_table="std.txt", panel="mix21")
>>> write_results(res, "out/run1")
"""

from __future__ import annotations

from .config import CalibrationConfig
from .crosscolumn import (
    AnchorRefiner,
    ClassAwareRefiner,
    ColumnCalibrator,
    MonotoneCurve,
    build_calibrator,
    exclude_pairs_near_mz,
    fit_robust_curve,
    match_features_by_mz,
)
from .errors import (
    AnchorIdentificationError,
    CalibrationError,
    ColumnResolutionError,
    ConfigError,
    InputFormatError,
    PanelError,
    RtAnchorError,
    RTUnitError,
)
from .helpers import describe_input, setup_logger, which_format, write_results
from .identify import build_native_template, identify_anchors
from .io import FeatureTable, detect_format, load_feature_table
from .irt import IRTMapper, build_irt, detect_landmarks
from .mixtures import (
    DEFAULT_MIXTURE,
    get_manifest,
    get_meta,
    mixture_keys,
    offered_keys,
    preview_payload,
)
from .panel import DEFAULT_MANIFEST, Panel, build_panel, load_manifest_csv, load_reference_csv
from .pipeline import CalibrationResult, calibrate
from .plasma_lipids import plasma_lipid_candidates
from .reference import (
    DEFAULT_REFERENCE_KEY,
    REFERENCE_SETS,
    available_references,
    load_reference,
    resolve_reference,
)

__version__ = "1.2.0"

__all__ = [
    # entry points
    "calibrate", "CalibrationResult", "CalibrationConfig",
    "describe_input", "which_format", "write_results", "setup_logger",
    # io
    "load_feature_table", "detect_format", "FeatureTable",
    # the cross-column engine
    "build_calibrator", "ColumnCalibrator", "MonotoneCurve", "fit_robust_curve",
    "match_features_by_mz", "exclude_pairs_near_mz",
    "AnchorRefiner", "ClassAwareRefiner", "plasma_lipid_candidates",
    # the iRT ruler
    "IRTMapper", "build_irt", "detect_landmarks",
    # reference datasets + panels
    "load_reference", "resolve_reference", "available_references",
    "REFERENCE_SETS", "DEFAULT_REFERENCE_KEY",
    "mixture_keys", "offered_keys", "get_manifest", "get_meta",
    "preview_payload", "DEFAULT_MIXTURE",
    # panel / identification machinery
    "build_panel", "Panel", "DEFAULT_MANIFEST", "load_manifest_csv", "load_reference_csv",
    "build_native_template", "identify_anchors",
    # errors
    "RtAnchorError", "InputFormatError", "ColumnResolutionError", "RTUnitError",
    "PanelError", "AnchorIdentificationError", "CalibrationError", "ConfigError",
]
