"""Typed exceptions for rt_anchor.

Every error carries an actionable message: what went wrong, which input it
concerns, and (where possible) how to fix it. Callers can catch the base
``RtAnchorError`` to handle any package-level failure.
"""

from __future__ import annotations


class RtAnchorError(Exception):
    """Base class for all rt_anchor errors."""


class ConfigError(RtAnchorError):
    """The configuration is invalid (unknown keys, bad values)."""


class InputFormatError(RtAnchorError):
    """The input file could not be parsed or its format could not be recognised."""


class ColumnResolutionError(InputFormatError):
    """A required logical column (m/z or RT) could not be located in the table."""


class RTUnitError(RtAnchorError):
    """Retention-time units are ambiguous or inconsistent with the declared unit."""


class PanelError(RtAnchorError):
    """The standard manifest / reference baseline is malformed or incomplete."""


class AnchorIdentificationError(RtAnchorError):
    """Too few standards could be identified to build a calibration."""


class CalibrationError(RtAnchorError):
    """The monotone warp could not be fitted (e.g. degenerate anchors)."""
