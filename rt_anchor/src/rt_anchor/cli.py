"""Command-line interface.

    rt-anchor describe  <table>
    rt-anchor calibrate --samples <table> --standards <run> --polarity positive [options] --out <prefix>

Outputs are CSV / JSON / log only (no plots).
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
from .panel import load_manifest_csv, load_reference_csv
from .pipeline import calibrate


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
        cfg = CalibrationConfig()        # sensible defaults; tune with the flags below
    if args.mz_tol_ppm is not None:
        cfg.mz_tol_ppm = args.mz_tol_ppm
    if args.rt_window is not None:
        cfg.rt_window_min = args.rt_window
    if args.min_anchors is not None:
        cfg.min_anchors = args.min_anchors
    if args.extrapolate:
        cfg.extrapolate = True
    return cfg


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(prog="rt-anchor", description="Standard-panel iRT calibration for LC-MS lipidomics.")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("describe", help="Detect format + summarise a table (no calibration).")
    d.add_argument("table")
    d.add_argument("--format", dest="source_format", default=None)
    d.add_argument("--polarity", default=None)

    c = sub.add_parser("calibrate", help="Calibrate a feature table to iRT.")
    c.add_argument("--samples", required=True, help="Feature table to calibrate.")
    c.add_argument("--standards", required=True,
                   help="Standard-panel run (required; builds the native anchor template).")
    c.add_argument("--polarity", required=True, choices=["positive", "negative", "pos", "neg"])
    c.add_argument("--single-files", default=None,
                   help="Dir or glob of per-injection tables -> per-sample tier.")
    c.add_argument("--format", dest="source_format", default=None,
                   help="Force input format (masscube/lipidscreener/msdial/mzmine).")
    c.add_argument("--rt-unit", dest="rt_unit", default=None, choices=["min", "sec"])
    c.add_argument("--manifest", default=None, help="Custom standard manifest CSV (Mode B).")
    c.add_argument("--reference", default=None, help="Custom reference-baseline CSV.")
    c.add_argument("--config", default=None,
                   help="JSON config file (base; the --mz-tol-ppm/--rt-window/--min-anchors flags override individual fields).")
    # matching / robustness parameters (replace any notion of an 'instrument type')
    c.add_argument("--mz-tol-ppm", dest="mz_tol_ppm", type=float, default=None,
                   help="m/z match tolerance in ppm (default 15).")
    c.add_argument("--rt-window", dest="rt_window", type=float, default=None,
                   help="RT search window around the native anchor RT, minutes (default 0.5).")
    c.add_argument("--min-anchors", dest="min_anchors", type=int, default=None,
                   help="Minimum anchors for a confident warp (default 6).")
    c.add_argument("--extrapolate", action="store_true",
                   help="Assign RI beyond the anchor span (capped); default off -> NaN beyond span.")
    c.add_argument("--no-report", action="store_true", help="Skip the HTML/PDF report (data only).")
    c.add_argument("--tic-style", default="clean", choices=["clean", "realistic", "both"],
                   help="Reconstructed-TIC style in the report (default clean; realistic adds simulated noise).")
    c.add_argument("--report-format", default="html,pdf",
                   help="Comma list of report formats: html,pdf.")
    c.add_argument("--out", required=True, help="Output path prefix.")

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

    # calibrate
    log = setup_logger()
    try:
        cfg = _build_config(args)
        manifest = load_manifest_csv(args.manifest) if args.manifest else None
        reference = load_reference_csv(args.reference) if args.reference else None
        single = _collect_single_files(args.single_files)
        if single:
            log.info(f"per-sample tier: {len(single)} injection tables")
        result = calibrate(
            sample_table=args.samples, polarity=args.polarity,
            standards_table=args.standards, single_files=single,
            config=cfg, manifest=manifest, reference=reference,
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
