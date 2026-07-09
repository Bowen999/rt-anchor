"""User-facing helpers: format inspection, result writers, logging setup.

None of these produce visualisations — outputs are CSV / JSON / plain-text log.
"""

from __future__ import annotations

import json
import logging
import os
from typing import Dict, List, Optional

import pandas as pd

from .io.adapters import detect_format
from .io.loader import load_feature_table
from .pipeline import CalibrationResult


def describe_input(path: str, source_format: Optional[str] = None,
                   polarity: Optional[str] = None) -> Dict:
    """Detect + summarise an input table without calibrating (dry-run helper)."""
    ft = load_feature_table(path, source_format=source_format, polarity=polarity)
    info = ft.summary()
    info["path"] = path
    if ft.meta:
        info["meta_keys"] = list(ft.meta.keys())
    return info


def which_format(path: str) -> str:
    fmt, _, _ = detect_format(path)
    return fmt


def setup_logger(logfile: Optional[str] = None, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger("rt_anchor")
    logger.setLevel(level)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def write_results(result: CalibrationResult, out_prefix: str,
                  write_log: bool = True, report: bool = True,
                  tic_style: str = "clean", report_formats=("html", "pdf")) -> Dict[str, str]:
    """Write the calibrated table (CSV), model (JSON), log (txt), and — by
    default — a visual Report (interactive HTML + static PDF).

    Set ``report=False`` for data-only output. Returns {kind: path}.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)
    paths: Dict[str, str] = {}

    csv_path = f"{out_prefix}_calibrated.csv"
    result.table.to_csv(csv_path, index=False)
    paths["calibrated_csv"] = csv_path

    model_path = f"{out_prefix}_model.json"
    with open(model_path, "w") as fh:
        # allow_nan=False + a NaN->null sanitiser so the file is *valid* JSON
        # (bare ``NaN`` is rejected by strict parsers, e.g. JS ``JSON.parse``).
        json.dump(_sanitize_json(result.model), fh, indent=2,
                  default=_json_default, allow_nan=False)
    paths["model_json"] = model_path

    anchors_path = f"{out_prefix}_anchors.csv"
    result.anchors.to_csv(anchors_path, index=False)
    paths["anchors_csv"] = anchors_path

    if write_log:
        log_path = f"{out_prefix}_log.txt"
        with open(log_path, "w") as fh:
            fh.write("\n".join(result.log) + "\n")
        paths["log_txt"] = log_path

    if report:
        try:
            from .viz.report import write_report
            paths.update(write_report(result, out_prefix, tic_style=tic_style,
                                      formats=report_formats))
        except ImportError as e:
            result.log.append(f"report skipped (missing viz deps: {e}); "
                              f"install rt_anchor[report] for matplotlib+plotly")
    return paths


def _sanitize_json(o):
    """Recursively replace non-finite floats (NaN/inf) with ``None``.

    Needed because numpy floats are ``float`` subclasses, so the json C-encoder
    serialises them natively (emitting a bare, invalid ``NaN``) and never calls
    ``default``. We walk the structure and null-out non-finite values so the
    written file is valid JSON.
    """
    import math

    if isinstance(o, dict):
        return {k: _sanitize_json(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_sanitize_json(v) for v in o]
    if isinstance(o, float):  # includes np.float64 (a float subclass)
        return None if not math.isfinite(o) else float(o)
    return o


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and (o != o):  # NaN
        return None
    return str(o)
