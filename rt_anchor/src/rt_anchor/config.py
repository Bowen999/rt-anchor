"""Configuration for the v2 cross-column calibration.

Defaults follow the **validated** method (see ``RI_CALIBRATION_SPEC_V2.md``
§7 and §9) with two deliberate changes since:

* the default m/z window is **15 ppm** for every match (details below);
* from **1.2.2** the homologous-series term (stage 1b) is on by default and
  the plasma-lipid anchor stage (stage 2) is **opt-in**
  (``use_sample_anchors=False``) — pass ``use_sample_anchors=True`` or the
  CLI's ``--sample-anchors`` to enable it.

Two values deserve to be read before they are changed:

* ``curve_frac = 0.1`` was chosen by 5-fold CV on the panel-masked training
  pairs of all five validation columns. Do not change it casually.
* ``match_mz_tol_ppm = 15.0`` is the m/z window of every *anonymous*
  operation. The validation runs used a flat 0.008 Da window instead;
  ``match_mz_tol_ppm=0, match_mz_tol_da=0.008`` reproduces them.

**One default, two knobs.** Both m/z tolerances default to 15 ppm with no
absolute floor, but they do different jobs and are tuned separately:

``match_mz_tol_ppm`` / ``match_mz_tol_da``
    The *anonymous* matching window, ``max(ppm window, absolute floor)`` —
    reciprocal feature matching between two runs, the stage-2 panel-mask
    exclusion, and the plasma-lipid m/z windows, all the places where nothing
    is claimed about a feature's identity. The floor is 0 (off) by default.

``mz_tol_ppm`` / ``mz_tol_min_da``
    The same window for *targeted* identification of named panel standards —
    the iRT landmarks and the detection QC, both of which go through
    :mod:`rt_anchor.identify`. The ``orbitrap()`` preset narrows it to 8 ppm.

Everything is override-able (constructor kwargs, the :meth:`qtof` /
:meth:`orbitrap` presets, or :meth:`from_dict`). Nothing is instrument-locked.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional, Tuple

#: v1 fields the v2 method retired, mapped to what replaced them. ``from_dict``
#: refuses a config carrying any of these rather than silently ignoring it — an
#: old config file that still says ``stds_fallback: false`` is asking for a
#: behaviour that no longer exists, and quietly running something else is worse
#: than failing.
REMOVED_FIELDS = {
    "scale_rt_lo_min":
        "the v2 iRT scale comes from landmarks detected on the reference "
        "standards run, not from scale constants — see irt_landmark_panel",
    "scale_rt_hi_min":
        "the v2 iRT scale comes from landmarks detected on the reference "
        "standards run, not from scale constants — see irt_landmark_panel",
    "coverage_gap_min":
        "the v1 anchor-gap guard is gone; curve quality is reported as "
        "curve.residual_min in model.json and as Cal_RT_uncertainty_min per row",
    "min_anchors_hard":
        "use min_anchors — under v2 too few anchors skips stage 2, it is never fatal",
    "sigma_gap_factor":
        "use sigma_window_pairs — uncertainty is now the local scatter of the "
        "matched pairs about the curve (spec §5.2)",
    "stds_fallback":
        "removed — the standards run is a required input of the v2 method, not a "
        "fallback for a sample whose own anchors ran out",
    "max_extrapolation_min":
        "use extrapolate_mode ('linear' or 'clamp') — v2 always assigns a value "
        "beyond the pair span and flags the row with is_extrapolated",
}


@dataclass
class CalibrationConfig:
    # ---- anonymous m/z matching (spec §2.1, §7) ----
    match_mz_tol_ppm: float = 15.0      # the matching window, ppm
    match_mz_tol_da: float = 0.0        # optional absolute floor, Da (0 = off)

    # ---- targeted identification: landmarks + detection QC (v1 meaning) ----
    mz_tol_ppm: float = 15.0            # QTOF default (Orbitrap ~8)
    mz_tol_min_da: float = 0.0          # optional absolute floor on the ppm window
    rt_window_min: float = 0.5          # window around the native-template RT
    rt_window_seed_min: float = 2.0     # window when no native template exists

    # ---- stage 1: the monotone curve ----
    curve_frac: float = 0.1             # LOESS fraction — CV-chosen, do not change
    curve_iter: int = 3                 # MAD outlier-removal rounds
    curve_mad_k: float = 3.0            # |resid - median| > k * 1.4826 * MAD -> dropped
    curve_min_points: int = 20          # below this the curve degrades to a line
    use_sample_pairs: bool = True       # merge the sample-run pairs into the fit

    # ---- stage 1b: homologous-series term ----
    use_series_term: bool = True        # stage 1b on/off (see rt_anchor.series)
    series_kmd_tol: float = 0.008       # tolerance on the Kendrick phase
    series_min_members: int = 4         # smallest cluster/chain that counts as a series
    series_end_reach_ch2: int = 2       # how many CH2 past a series end a query may sit
    series_exclude_rt_frac: float = 0.0015   # co-elution exclusion, fraction of knot span
    series_exclude_rt_floor_min: float = 0.03    # ... with this absolute floor, minutes
    series_gate_min_mse_reduction: float = 0.2   # the "do no harm" gate
    series_min_covered_pairs: int = 20  # below this the gate declines unread

    # ---- stage 2: class-aware sample anchors ----
    use_sample_anchors: bool = False    # stage 2 on/off (opt-in since 1.2.2)
    class_aware: bool = True            # class-offset residual decomposition
    min_anchors: int = 3                # below this stage 2 is skipped (not fatal)
    anchor_gate_min_mse_reduction: float = 0.2   # the "do no harm" gate
    anchor_tie_tol_min: float = 0.15    # isomer re-pick tie window (minutes)
    isomer_sn_frac: float = 0.10        # isomer census: S/N fraction of the pick

    # ---- uncertainty + confidence (spec §5.2) ----
    sigma_window_pairs: int = 50        # pairs in the local-scatter window
    # Banded against the method's own validated accuracy rather than against
    # v1's numbers. Held-out cross-column error is median 0.09-0.19 min and P90
    # 0.23-0.44 min; on the bundled reference's ruler (~4.9 iRT per minute)
    # that is a median of ~0.4-0.9 iRT and a P90 of ~1.1-2.1. v1's thresholds
    # (0.6 / 3.0) were set for the far smaller within-run warp residual, and
    # under v2 they label essentially every real feature "medium" — a column
    # that never varies tells the reader nothing.
    conf_high_irt: float = 1.0          # iRT_uncertainty < this -> high
    conf_low_irt: float = 2.0           # iRT_uncertainty > this -> low

    # ---- the iRT ruler ----
    irt_landmark_panel: str = "mix21"   # landmarks used when panel="none"

    # ---- extrapolation ----
    # v2 always assigns a value outside the matched-pair span and flags the row.
    extrapolate: bool = True
    extrapolate_mode: str = "linear"    # "linear" = terminal slope, "clamp" = hold

    # ---- monotonicity guard, read by identify._drop_nonmonotone ----
    tie_epsilon_min: float = 0.005      # equal-RT tolerance
    max_drop_anchors: int = 3           # cap on anchors removed by the guard

    # ---- optional peak-quality gate (applied only if the columns exist) ----
    min_gaussian_similarity: float = 0.0            # 0 disables
    asymmetry_range: Tuple[float, float] = (0.0, 1e9)   # widen -> disabled

    # ---- RT unit handling ----
    rt_unit: Optional[str] = None       # None = auto (per-format default + range check)

    # ---- misc ----
    verbose: bool = True

    # ---- serialisation ----
    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "CalibrationConfig":
        from .errors import ConfigError

        known = {f for f in cls.__dataclass_fields__}  # type: ignore[attr-defined]
        retired = sorted(set(d) & set(REMOVED_FIELDS))
        if retired:
            detail = "; ".join(f"'{k}' -> {REMOVED_FIELDS[k]}" for k in retired)
            raise ConfigError(
                f"Config keys removed by the v2 cross-column method: {retired}. {detail}."
            )
        unknown = set(d) - known
        if unknown:
            raise ConfigError(f"Unknown config keys: {sorted(unknown)}. Known: {sorted(known)}")
        cfg = cls(**{k: v for k, v in d.items() if k in known})
        if cfg.asymmetry_range is not None:
            cfg.asymmetry_range = tuple(cfg.asymmetry_range)  # json gives lists
        return cfg

    # ---- instrument presets ----
    @classmethod
    def qtof(cls, **overrides) -> "CalibrationConfig":
        return cls(**overrides)

    @classmethod
    def orbitrap(cls, **overrides) -> "CalibrationConfig":
        base = dict(mz_tol_ppm=8.0, rt_window_min=0.3)
        base.update(overrides)
        return cls(**base)
