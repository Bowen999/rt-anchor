"""Command-line interface.

    rt-anchor describe   <table>
    rt-anchor references
    rt-anchor calibrate  --samples S --standards T --polarity positive
                         [--panel {mix15,mix21,none}] [--reference-sample R
                          --reference-standards RT] [...] --out PREFIX

``calibrate`` writes every companion output of the spec's §5 (calibrated CSV,
model JSON, anchors / pairs / landmarks CSVs, log) plus the visual report unless
``--no-report`` is given.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from typing import List, Optional

from .config import CalibrationConfig
from .errors import RtAnchorError
from .helpers import describe_input, setup_logger, write_results
from .mixtures import DEFAULT_MIXTURE, offered_keys
from .panel import load_manifest_csv
from .pipeline import calibrate
from .reference import available_references

#: Flags the v2 method retired, and what replaced them. Kept parseable so the
#: user gets a sentence instead of argparse's "unrecognized arguments".
REMOVED_FLAGS = {
    "reference": ("--reference", "the v1 reference-baseline CSV is gone — the v2 "
                  "method calibrates against a reference *dataset*: use "
                  "--reference-sample / --reference-standards, or omit both for "
                  "the bundled Column 25 reference"),
    "no_stds_fallback": ("--no-stds-fallback", "the standards run is now a required "
                         "input of every calibration, not a fallback; there is "
                         "nothing to switch off"),
}


def _collect_single_files(spec: Optional[str]) -> Optional[List[str]]:
    if not spec:
        return None
    if os.path.isdir(spec):
        files = sorted(glob.glob(os.path.join(spec, "*.txt")) +
                       glob.glob(os.path.join(spec, "*.csv")))
    else:
        files = sorted(glob.glob(spec))
    return files or None


def _build_config(args) -> CalibrationConfig:
    if args.config:
        with open(args.config) as fh:
            cfg = CalibrationConfig.from_dict(json.load(fh))
    else:
        cfg = CalibrationConfig()
    if args.mz_tol_ppm is not None:
        cfg.mz_tol_ppm = args.mz_tol_ppm
    if args.match_mz_tol_ppm is not None:
        cfg.match_mz_tol_ppm = args.match_mz_tol_ppm
    if args.mz_tol_da is not None:
        cfg.match_mz_tol_da = args.mz_tol_da
    if args.rt_window is not None:
        cfg.rt_window_min = args.rt_window
    if args.curve_frac is not None:
        cfg.curve_frac = args.curve_frac
    if args.min_anchors is not None:
        cfg.min_anchors = args.min_anchors
    if args.no_sample_anchors:
        cfg.use_sample_anchors = False
    if args.no_sample_pairs:
        cfg.use_sample_pairs = False
    if args.no_extrapolate:
        cfg.extrapolate = False
    if args.extrapolate_mode is not None:
        cfg.extrapolate_mode = args.extrapolate_mode
    return cfg


def _check_removed(args) -> Optional[str]:
    for dest, (flag, why) in REMOVED_FLAGS.items():
        val = getattr(args, dest, None)
        if val:
            return f"{flag} was removed: {why}."
    return None


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="rt-anchor",
        description="Cross-column retention-time calibration for LC-MS lipidomics.")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("describe", help="Detect format + summarise a table (no calibration).")
    d.add_argument("table")
    d.add_argument("--format", dest="source_format", default=None)
    d.add_argument("--polarity", default=None)

    sub.add_parser("references", help="List the bundled reference datasets.")

    c = sub.add_parser("calibrate", help="Calibrate a feature table onto the reference column.")
    c.add_argument("--samples", required=True, help="Feature table to calibrate.")
    c.add_argument("--standards", required=True,
                   help="Standards run from the SAME column and gradient as --samples.")
    c.add_argument("--polarity", required=True, choices=["positive", "negative", "pos", "neg"])
    c.add_argument("--panel", default=DEFAULT_MIXTURE, choices=offered_keys() + ["mix21lpc"],
                   help=f"Standard mixture used (default {DEFAULT_MIXTURE}); sets the iRT "
                        f"landmark panel and the detection QC. 'none' still calibrates — "
                        f"the iRT ruler falls back to irt_landmark_panel.")
    c.add_argument("--reference-sample", dest="reference_sample", default=None,
                   help="Custom reference sample run (needs --reference-standards too).")
    c.add_argument("--reference-standards", dest="reference_standards", default=None,
                   help="Custom reference standards run (needs --reference-sample too).")
    c.add_argument("--reference-key", dest="reference_key", default="col35",
                   help="Bundled reference dataset key (default col35).")
    c.add_argument("--single-files", default=None,
                   help="Dir or glob of per-injection tables -> per-sample tier.")
    c.add_argument("--format", dest="source_format", default=None,
                   help="Force input format (masscube/lipidscreener/msdial/mzmine).")
    c.add_argument("--rt-unit", dest="rt_unit", default=None, choices=["min", "sec"])
    c.add_argument("--manifest", default=None,
                   help="Custom standard manifest CSV; overrides --panel.")
    c.add_argument("--config", default=None,
                   help="JSON config file (base; the flags below override individual fields).")
    # --- matching / robustness -------------------------------------------
    c.add_argument("--match-mz-tol-ppm", dest="match_mz_tol_ppm", type=float,
                   default=None,
                   help="ppm window for ANONYMOUS feature matching "
                        "(match_mz_tol_ppm, default 15). This is the window the "
                        "cross-column curve is built with.")
    c.add_argument("--mz-tol-da", dest="mz_tol_da", type=float, default=None,
                   help="Absolute floor beneath the anonymous matching ppm "
                        "window, Da (match_mz_tol_da, default 0.008).")
    c.add_argument("--mz-tol-ppm", dest="mz_tol_ppm", type=float, default=None,
                   help="m/z tolerance in ppm for TARGETED panel-standard "
                        "identification: landmarks and detection QC (default 15).")
    c.add_argument("--rt-window", dest="rt_window", type=float, default=None,
                   help="RT search window around the native standard RT, minutes (default 0.5).")
    c.add_argument("--curve-frac", dest="curve_frac", type=float, default=None,
                   help="LOESS fraction of the stage-1 curve (default 0.1, CV-chosen).")
    c.add_argument("--min-anchors", dest="min_anchors", type=int, default=None,
                   help="Minimum stage-2 plasma-lipid anchors (default 3); below this "
                        "the stage-1 curve is used.")
    c.add_argument("--no-sample-anchors", dest="no_sample_anchors", action="store_true",
                   help="Skip stage 2 entirely (stage-1 curve only).")
    c.add_argument("--no-sample-pairs", dest="no_sample_pairs", action="store_true",
                   help="Do not merge the sample-run m/z pairs into the stage-1 curve.")
    c.add_argument("--no-extrapolate", dest="no_extrapolate", action="store_true",
                   help="Do not assign Cal_RT beyond the matched-pair span "
                        "(default: assign it and flag is_extrapolated).")
    c.add_argument("--extrapolate-mode", dest="extrapolate_mode", default=None,
                   choices=["linear", "clamp"],
                   help="Beyond-span behaviour: 'linear' extends at the terminal slope "
                        "(default), 'clamp' holds the edge value.")
    # --- output ------------------------------------------------------------
    c.add_argument("--no-report", action="store_true", help="Skip the HTML/PDF report (data only).")
    c.add_argument("--tic-style", default="clean", choices=["clean", "realistic", "both"],
                   help="Reconstructed-TIC style in the report (default clean).")
    c.add_argument("--report-format", default="html,pdf",
                   help="Comma list of report formats: html,pdf.")
    c.add_argument("--out", required=True, help="Output path prefix.")
    # --- retired flags: parsed only so they can be refused by name ---------
    c.add_argument("--reference", default=None, help=argparse.SUPPRESS)
    c.add_argument("--no-stds-fallback", dest="no_stds_fallback",
                   action="store_true", help=argparse.SUPPRESS)

    args = p.parse_args(argv)

    if args.cmd == "describe":
        try:
            print(json.dumps(describe_input(args.table, source_format=args.source_format,
                                            polarity=args.polarity), indent=2, default=str))
            return 0
        except RtAnchorError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
        except Exception as e:  # e.g. FileNotFoundError -> clean message, no traceback
            print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
            return 3

    if args.cmd == "references":
        print(json.dumps(available_references(), indent=2, default=str))
        return 0

    # calibrate
    log = setup_logger()
    removed = _check_removed(args)
    if removed:
        log.error(removed)
        return 2
    try:
        cfg = _build_config(args)
        manifest = load_manifest_csv(args.manifest) if args.manifest else None
        single = _collect_single_files(args.single_files)
        if single:
            log.info(f"per-sample tier: {len(single)} injection tables")
        result = calibrate(
            sample_table=args.samples, polarity=args.polarity,
            standards_table=args.standards, panel=args.panel,
            reference_sample=args.reference_sample,
            reference_standards=args.reference_standards,
            reference_key=args.reference_key,
            single_files=single, config=cfg, manifest=manifest,
            source_format=args.source_format, rt_unit=args.rt_unit,
        )
        for line in result.log:
            log.info(line)
        fmts = tuple(f.strip() for f in args.report_format.split(",") if f.strip())
        paths = write_results(result, args.out, report=not args.no_report,
                              tic_style=args.tic_style, report_formats=fmts)
        log.info(f"wrote: {json.dumps(paths)}")
        return 0
    except RtAnchorError as e:
        log.error(str(e))
        return 2
    except Exception as e:  # unexpected -> surface with type
        log.error(f"{type(e).__name__}: {e}")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
