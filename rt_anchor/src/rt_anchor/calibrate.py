"""The calibration itself: a monotone RT->iRT warp plus numeric uncertainty.

* :class:`MonotoneWarp` — shape-preserving monotone interpolation (PCHIP) with
  ``extrapolate=False`` so features beyond the anchor span become NaN (never a
  fabricated value). Leave-one-anchor-out residuals give the honest in-domain
  interpolation error.
* :func:`apply_warp` — turn a warp + a config into per-feature ``RI``,
  ``RI_uncertainty``, ``RI_reliability`` and ``is_extrapolated``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from .config import CalibrationConfig
from .errors import CalibrationError


@dataclass
class MonotoneWarp:
    rt: np.ndarray          # anchor observed RT (minutes), strictly increasing
    irt: np.ndarray         # anchor iRT values
    _pchip: PchipInterpolator
    loo_resid: np.ndarray   # per-anchor leave-one-out residual (iRT units)

    @property
    def rt_min(self) -> float:
        return float(self.rt[0])

    @property
    def rt_max(self) -> float:
        return float(self.rt[-1])

    def predict(self, rt_min: np.ndarray, extrapolate: bool = False,
                max_extrap_min: float = 1.0) -> np.ndarray:
        rt_min = np.asarray(rt_min, dtype=float)
        out = self._pchip(rt_min)  # NaN outside [rt_min, rt_max] (extrapolate=False)
        if extrapolate:
            below = rt_min < self.rt_min
            above = rt_min > self.rt_max
            # linear terminal-slope extrapolation, capped
            if below.any():
                s = (self.irt[1] - self.irt[0]) / (self.rt[1] - self.rt[0])
                d = np.clip(self.rt_min - rt_min[below], 0, max_extrap_min)
                out[below] = self.irt[0] - s * d
            if above.any():
                s = (self.irt[-1] - self.irt[-2]) / (self.rt[-1] - self.rt[-2])
                d = np.clip(rt_min[above] - self.rt_max, 0, max_extrap_min)
                out[above] = self.irt[-1] + s * d
        return out


def fit_warp(anchor_rt_min: np.ndarray, anchor_irt: np.ndarray,
             config: CalibrationConfig) -> MonotoneWarp:
    rt = np.asarray(anchor_rt_min, dtype=float)
    irt = np.asarray(anchor_irt, dtype=float)
    order = np.argsort(rt)
    rt, irt = rt[order], irt[order]

    # strict-increase guard (PCHIP requires strictly increasing x)
    keep = np.concatenate([[True], np.diff(rt) > config.tie_epsilon_min])
    rt, irt = rt[keep], irt[keep]
    if rt.size < max(2, config.min_anchors_hard):
        raise CalibrationError(
            f"Only {rt.size} usable anchors after the monotonicity guard "
            f"(need >= {max(2, config.min_anchors_hard)}). Loosen mz_tol_ppm / RT "
            f"window, provide a standard-panel run, or check the polarity."
        )
    if not np.all(np.diff(irt) > 0):
        raise CalibrationError("Anchor iRT values are not strictly increasing — reference baseline is inconsistent.")

    pchip = PchipInterpolator(rt, irt, extrapolate=False)
    loo = _loo_residuals(rt, irt)
    return MonotoneWarp(rt=rt, irt=irt, _pchip=pchip, loo_resid=loo)


def _loo_residuals(rt: np.ndarray, irt: np.ndarray) -> np.ndarray:
    """Leave-one-anchor-out residual (iRT) at each anchor; NaN for the two ends
    (holding out an endpoint is extrapolation, not interpolation)."""
    n = rt.size
    res = np.full(n, np.nan)
    for i in range(1, n - 1):  # interior only
        m = np.ones(n, bool); m[i] = False
        f = PchipInterpolator(rt[m], irt[m], extrapolate=True)
        res[i] = float(f(rt[i]) - irt[i])
    return res


def apply_warp(rt_min: np.ndarray, warp: MonotoneWarp, config: CalibrationConfig,
               scope: str, warp_source: str = "self") -> pd.DataFrame:
    """Compute RI + uncertainty for an array of feature RTs (minutes)."""
    rt_min = np.asarray(rt_min, dtype=float)
    ri = warp.predict(rt_min, extrapolate=config.extrapolate,
                      max_extrap_min=config.max_extrapolation_min)
    in_span = (rt_min >= warp.rt_min) & (rt_min <= warp.rt_max)
    is_extrap = ~in_span & np.isfinite(rt_min)

    sigma = _sigma_ri(rt_min, warp, config)
    conf = _reliability(sigma, is_extrap, config)

    # keep sigma / confidence consistent with RI: an uncalibrated feature (RI NaN,
    # e.g. beyond the anchor span with extrapolate=False) gets sigma NaN + 'none'.
    uncal = ~np.isfinite(ri)
    sigma[uncal] = np.nan
    conf[uncal] = "none"

    return pd.DataFrame({
        "RI": ri,
        "RI_uncertainty": sigma,
        "RI_reliability": conf,
        "is_extrapolated": is_extrap,
        "calibration_scope": scope,
        "warp_source": warp_source,
    })


def _sigma_ri(rt_min: np.ndarray, warp: MonotoneWarp, config: CalibrationConfig) -> np.ndarray:
    slope = config.slope_irt_per_min()
    # interior LOO residual as a function of RT (|resid|, iRT units)
    interior = np.isfinite(warp.loo_resid)
    if interior.any():
        loo_interp = np.interp(rt_min, warp.rt[interior], np.abs(warp.loo_resid[interior]),
                               left=np.abs(warp.loo_resid[interior][0]),
                               right=np.abs(warp.loo_resid[interior][-1]))
        sigma_scale = float(np.nanmedian(np.abs(warp.loo_resid[interior])))
    else:
        loo_interp = np.zeros_like(rt_min)
        sigma_scale = 0.0
    # nearest-anchor gap term (minutes -> iRT)
    nearest_gap = np.min(np.abs(rt_min[:, None] - warp.rt[None, :]), axis=1)
    gap_term = slope * config.sigma_gap_factor * nearest_gap
    sigma_local = np.maximum(loo_interp, gap_term)
    sigma = np.sqrt(sigma_local ** 2 + sigma_scale ** 2)
    sigma[~np.isfinite(rt_min)] = np.nan
    return sigma


def _reliability(sigma: np.ndarray, is_extrap: np.ndarray, config: CalibrationConfig) -> np.ndarray:
    conf = np.full(sigma.shape, "medium", dtype=object)
    conf[sigma < config.conf_high_irt] = "high"
    conf[sigma > config.conf_low_irt] = "low"
    conf[is_extrap] = "low"
    conf[~np.isfinite(sigma)] = "none"
    return conf
