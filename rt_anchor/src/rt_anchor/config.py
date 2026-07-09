"""Configuration for the calibration.

Defaults are tuned for **QTOF** data (the common case for this project):
wider m/z tolerance and RT windows than an Orbitrap would need. Every value is
override-able (constructor kwargs, ``CalibrationConfig.qtof()`` /
``.orbitrap()`` presets, or ``from_dict``). Nothing here is instrument-locked.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Optional, Tuple


@dataclass
class CalibrationConfig:
    # ---- m/z matching ----
    mz_tol_ppm: float = 15.0            # QTOF default (Orbitrap ~8)
    mz_tol_min_da: float = 0.0          # optional absolute floor on the ppm window

    # ---- RT search window around the native-template anchor RT (minutes) ----
    rt_window_min: float = 0.5          # QTOF default (Orbitrap ~0.3)
    # window used when NO native template exists and we must seed from reference RT
    rt_window_seed_min: float = 2.0

    # ---- iRT scale (fixed method constants; dimensionless, minute-derived) ----
    scale_rt_lo_min: float = 1.0
    scale_rt_hi_min: float = 18.4

    # ---- anchor acceptance / warp robustness ----
    min_anchors: int = 6                # below this -> low-confidence / uncalibrated
    min_anchors_hard: int = 3           # below this -> refuse to calibrate
    coverage_gap_min: float = 3.0       # max in-domain inter-anchor RT gap (minutes)
    tie_epsilon_min: float = 0.005      # equal-RT tolerance for monotonicity guard
    max_drop_anchors: int = 3           # cap on anchors removed by the guard

    # ---- extrapolation ----
    extrapolate: bool = False           # if False: RI = NaN beyond anchor span
    max_extrapolation_min: float = 1.0  # only used when extrapolate=True

    # ---- confidence thresholds (iRT units) ----
    conf_high_irt: float = 0.6          # RI_uncertainty < this  -> high
    conf_low_irt: float = 3.0           # RI_uncertainty > this  -> low
    sigma_gap_factor: float = 0.1       # sigma_local ~= factor * nearest-anchor gap (min) * slope

    # ---- optional peak-quality gate (applied only if the columns exist) ----
    min_gaussian_similarity: float = 0.0    # 0 disables; spec suggests ~0.7
    asymmetry_range: Tuple[float, float] = (0.0, 1e9)  # widen -> disabled

    # ---- RT unit handling ----
    rt_unit: Optional[str] = None       # None = auto (per-format default + range check)

    # ---- misc ----
    verbose: bool = True

    def slope_irt_per_min(self) -> float:
        return 100.0 / (self.scale_rt_hi_min - self.scale_rt_lo_min)

    def irt_from_rt(self, rt_min: float) -> float:
        return 100.0 * (rt_min - self.scale_rt_lo_min) / (self.scale_rt_hi_min - self.scale_rt_lo_min)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "CalibrationConfig":
        from .errors import ConfigError
        known = {f for f in cls.__dataclass_fields__}  # type: ignore[attr-defined]
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
