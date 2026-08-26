"""rt_anchor — standard-panel retention-index (iRT) calibration for LC-MS lipidomics.

Unifies feature tables from MS-DIAL, MZmine, MassCube and LipidScreener, then
warps retention time onto a dimensionless iRT scale anchored by a spiked lipid
standard panel. Outputs are CSV / JSON / log (no visualisation).

Quick start
-----------
>>> from rt_anchor import calibrate, write_results
>>> res = calibrate("samples.txt", polarity="positive", standards_table="std.txt")
>>> write_results(res, "out/run1")
"""

from __future__ import annotations

from .calibrate import MonotoneWarp, apply_warp, fit_warp
from .config import CalibrationConfig
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
from .panel import DEFAULT_MANIFEST, Panel, build_panel, load_manifest_csv, load_reference_csv
from .pipeline import CalibrationResult, calibrate

__version__ = "1.1.0"

__all__ = [
    "calibrate", "CalibrationResult", "CalibrationConfig",
    "load_feature_table", "detect_format", "FeatureTable",
    "build_panel", "Panel", "DEFAULT_MANIFEST", "load_manifest_csv", "load_reference_csv",
    "build_native_template", "identify_anchors",
    "fit_warp", "apply_warp", "MonotoneWarp",
    "describe_input", "which_format", "write_results", "setup_logger",
    "RtAnchorError", "InputFormatError", "ColumnResolutionError", "RTUnitError",
    "PanelError", "AnchorIdentificationError", "CalibrationError", "ConfigError",
]
