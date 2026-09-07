"""Cross-column RT calibration: the two-stage engine.

This is the port of the validated v2 method (``rt_calibrate.py``). It maps a
feature's retention time from *the user's column* onto *the reference column's*
time axis. Two stages, and the order matters:

**Stage 1 — the anchor-free monotone curve.** Every feature of the user's
standards run is matched against every feature of the reference standards run
by accurate m/z, reciprocal best match only, so the pairing is 1-to-1. No
standard identities are used at any point, which is the whole reason the method
survives a mixture that is unknown, mis-spiked or simply undetectable. Pairs
from the two *sample* runs are merged in as well, because serum features
densely cover the very early and very late elution regions where a standards
mixture is sparse, and that is exactly where an unconstrained curve goes wrong.
A LOESS smooth, an isotonic projection to guarantee monotonicity, and a PCHIP
through the isotonic knots give the curve; up to ``curve_iter`` rounds of
MAD-based trimming remove pairs the curve cannot explain (removal only, never
resurrection). Beyond the knots the curve extends at its terminal slope, and
every such RT is flagged extrapolated.

``curve_frac = 0.1`` is not a guess. It was chosen by 5-fold CV on the
panel-masked training pairs of all five validation columns. Do not change it.

**Stage 2 — class-aware anchor refinement, gated.** Anchor residuals on long
gradients are class-systematic: SM and CE sit high while PC and TG sit low, and
the classes interleave along the RT axis. A single piecewise-linear correction
therefore both misses the sign flips and gets shrunk toward nothing by
cross-class contamination in the shrinkage selection. So the residual is
decomposed instead::

    resid_i            = m[class_i] + g(x_i) + eps
    correction(x, cls) = lam_c * m[cls] + lam_g * g(x)
    Cal_RT             = f(RT) + correction(RT, class)

with ``m_c`` the median residual within class ``c`` and ``g`` a piecewise-linear
fit through the class-detrended residuals. ``(lam_g, lam_c)`` come from a 2-D
grid scored by leave-one-anchor-out MSE. Then the gate: the correction is
engaged only if its LOO MSE is at least ``anchor_gate_min_mse_reduction`` below
the zero-correction MSE. This is the "do no harm" rule, and it is what
correctly switches column 90 off — there the anchor residuals are
class-incoherent single-lipid deviations that do not generalise, and applying
them would make the calibration worse while looking like it was doing work.

Everything reads its inputs through :class:`~rt_anchor.io.schema.FeatureTable`,
never through literal MS-DIAL column names, so MZmine / MassCube /
LipidScreener exports go through the same path. Where a table has no adduct or
no Fill%, the corresponding filter is skipped rather than faked.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from sklearn.isotonic import IsotonicRegression
from statsmodels.nonparametric.smoothers_lowess import lowess

from .errors import CalibrationError
from .io.schema import FeatureTable

# Defaults for the spec §7 fields. They are read off the config when it carries
# them (which it does from v2 onwards) and used from here otherwise, so this
# module stays importable against an older CalibrationConfig.
DEFAULTS = {
    "mz_tol_ppm": 15.0,
    "mz_tol_min_da": 0.008,
    "curve_frac": 0.1,
    "curve_iter": 3,
    "curve_mad_k": 3.0,
    "curve_min_points": 20,
    "class_aware": True,
    "min_anchors": 3,
    "anchor_gate_min_mse_reduction": 0.2,
    "sigma_window_pairs": 50,
    "extrapolate_mode": "linear",
}

#: Shrinkage grid searched by the leave-one-anchor-out lambda selection.
LAMBDA_GRID = np.round(np.arange(0, 1.05, 0.1), 2)

PAIR_COLUMNS = ["idx_a", "idx_b", "mz_a", "rt_a", "rt_b"]

# Anchor-residual MSE below which there is nothing for stage 2 to correct
# (1e-12 min^2 = 1 microsecond RMS — far below any real chromatographic scale,
# so this only ever catches the exactly-zero case, e.g. self-calibration).
_NEGLIGIBLE_RESID_MSE = 1e-12


def cfg(config: Any, name: str, override: Any = None) -> Any:
    """One config value: explicit override, else the config, else the default."""
    if override is not None:
        return override
    if config is not None:
        val = getattr(config, name, None)
        if val is not None:
            return val
    return DEFAULTS[name]


def mz_window(mz, config: Any = None, *, mz_tol_ppm: Optional[float] = None,
              mz_tol_min_da: Optional[float] = None) -> np.ndarray:
    """Half-width of the m/z match window, per feature, in Da.

    ``max(ppm window, absolute floor)``. Anonymous matching uses
    ``match_mz_tol_ppm`` / ``match_mz_tol_da`` (defaults 15 ppm, 0.008 Da), so
    the floor dominates below m/z 533 and the ppm term above it; targeted
    identification uses ``mz_tol_ppm`` / ``mz_tol_min_da``.
    """
    ppm = float(cfg(config, "mz_tol_ppm", mz_tol_ppm))
    floor = float(cfg(config, "mz_tol_min_da", mz_tol_min_da))
    return np.maximum(np.asarray(mz, dtype=float) * ppm / 1e6, floor)


# ---------------------------------------------------------------------------
# Feature matching
# ---------------------------------------------------------------------------

def match_features_by_mz(ft_a: FeatureTable,
                         ft_b: FeatureTable,
                         config: Any = None,
                         *,
                         mz_tol_ppm: Optional[float] = None,
                         mz_tol_min_da: Optional[float] = None,
                         reciprocal: bool = True) -> pd.DataFrame:
    """Match features between two runs by accurate m/z alone.

    For each feature of A the strongest in-window candidate of B wins, and vice
    versa; with ``reciprocal`` only mutual best matches survive, which makes the
    pairing 1-to-1 and drops the ambiguous windows where a wrong pair would bend
    the curve. "Strongest" is ``FeatureTable.sn()`` — the S/N column where the
    format has one, summed abundance where it does not (see ``_match_arrays``).

    ``idx_a`` / ``idx_b`` are **row positions** in the respective tables (not
    labels), so a caller can go back to the original rows. Rows with a missing
    m/z or RT cannot be matched and are excluded.

    Returns a frame with columns ``idx_a, idx_b, mz_a, rt_a, rt_b``.
    """
    mz_a, rt_a, rank_a, keep_a = _match_arrays(ft_a)
    mz_b, rt_b, rank_b, keep_b = _match_arrays(ft_b)
    if mz_a.size == 0 or mz_b.size == 0:
        return pd.DataFrame(columns=PAIR_COLUMNS)

    tol_a = mz_window(mz_a, config, mz_tol_ppm=mz_tol_ppm, mz_tol_min_da=mz_tol_min_da)
    tol_b = mz_window(mz_b, config, mz_tol_ppm=mz_tol_ppm, mz_tol_min_da=mz_tol_min_da)

    a_to_b = _best_match(mz_a, tol_a, mz_b, rank_b)
    b_to_a = _best_match(mz_b, tol_b, mz_a, rank_a)

    pairs = []
    for i, j in enumerate(a_to_b):
        if j < 0:
            continue
        if reciprocal and b_to_a[j] != i:
            continue
        pairs.append((i, j))
    if not pairs:
        return pd.DataFrame(columns=PAIR_COLUMNS)
    ia, ib = np.array(pairs).T
    return pd.DataFrame({
        "idx_a": keep_a[ia],
        "idx_b": keep_b[ib],
        "mz_a": mz_a[ia],
        "rt_a": rt_a[ia],
        "rt_b": rt_b[ib],
    })


def exclude_pairs_near_mz(pairs: pd.DataFrame,
                          mz_values: Iterable[float],
                          config: Any = None,
                          *,
                          mz_tol_ppm: Optional[float] = None,
                          mz_tol_min_da: Optional[float] = None) -> pd.DataFrame:
    """Drop pairs whose source m/z is within tolerance of any given m/z.

    Used to keep the stage-2 anchor candidates' own feature pairs out of the
    stage-1 curve. Without it the anchors would partly be fitting themselves,
    and their residuals — the quantity the gate is judging — would be
    optimistically small.
    """
    mzs = np.asarray(list(mz_values), dtype=float)
    if pairs.empty or mzs.size == 0:
        return pairs.reset_index(drop=True)
    src = pairs["mz_a"].to_numpy(dtype=float)
    tol = mz_window(src, config, mz_tol_ppm=mz_tol_ppm, mz_tol_min_da=mz_tol_min_da)
    hit = np.any(np.abs(src[:, None] - mzs[None, :]) < tol[:, None], axis=1)
    return pairs.loc[~hit].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Stage 1: robust monotone curve
# ---------------------------------------------------------------------------

class MonotoneCurve:
    """Monotone RT mapping through knots, with defined out-of-range behaviour.

    ``mode="linear"`` extends at the terminal slope of the end segment (the v2
    behaviour: a chromatographic gradient does not stop being a gradient one
    knot past the last matched pair). ``mode="clamp"`` holds the edge value
    instead, for callers who would rather a wrong answer be obviously flat.
    Either way :meth:`is_extrapolated` flags every RT outside the knot range.
    """

    def __init__(self, x_knots: np.ndarray, y_knots: np.ndarray, mode: str = "linear"):
        x_knots = np.asarray(x_knots, dtype=float)
        y_knots = np.asarray(y_knots, dtype=float)
        ux, idx = np.unique(x_knots, return_index=True)
        uy = y_knots[idx]
        if ux.size < 2:
            raise CalibrationError(
                "Cannot build a retention-time curve from fewer than two distinct "
                "source RTs. This usually means the two runs share almost no "
                "features — check the polarity and the m/z tolerance."
            )
        if mode not in ("linear", "clamp"):
            raise CalibrationError(
                f"Unknown extrapolate_mode '{mode}' (expected 'linear' or 'clamp')."
            )
        self.mode = mode
        self.x_knots, self.y_knots = ux, uy
        self.x0, self.x1 = float(ux[0]), float(ux[-1])
        self.y0, self.y1 = float(uy[0]), float(uy[-1])
        self._pchip = PchipInterpolator(ux, uy, extrapolate=False)
        self._m_lo = (uy[1] - uy[0]) / max(ux[1] - ux[0], 1e-9)
        self._m_hi = (uy[-1] - uy[-2]) / max(ux[-1] - ux[-2], 1e-9)

    def predict(self, x) -> np.ndarray:
        x = np.atleast_1d(np.asarray(x, dtype=float))
        y = self._pchip(np.clip(x, self.x0, self.x1))
        lo = x < self.x0
        hi = x > self.x1
        if self.mode == "linear":
            y[lo] = self.y0 + (x[lo] - self.x0) * self._m_lo
            y[hi] = self.y1 + (x[hi] - self.x1) * self._m_hi
        else:
            y[lo] = self.y0
            y[hi] = self.y1
        return y

    def is_extrapolated(self, x) -> np.ndarray:
        x = np.atleast_1d(np.asarray(x, dtype=float))
        return (x < self.x0) | (x > self.x1)

    def slope(self, x) -> np.ndarray:
        """d(reference RT)/d(source RT); the terminal slope outside the knots."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        d = self._pchip.derivative()(np.clip(x, self.x0, self.x1))
        d = np.asarray(d, dtype=float)
        if self.mode == "linear":
            d[x < self.x0] = self._m_lo
            d[x > self.x1] = self._m_hi
        else:
            d[(x < self.x0) | (x > self.x1)] = 0.0
        return d


def fit_robust_curve(rt_x: np.ndarray,
                     rt_y: np.ndarray,
                     frac: float = DEFAULTS["curve_frac"],
                     n_iter: int = DEFAULTS["curve_iter"],
                     mad_k: float = DEFAULTS["curve_mad_k"],
                     min_points: int = DEFAULTS["curve_min_points"],
                     mode: str = DEFAULTS["extrapolate_mode"],
                     ) -> Tuple[MonotoneCurve, np.ndarray]:
    """Robust monotone curve ``rt_x -> rt_y``.

    LOESS with robustifying iterations, then isotonic regression to make
    monotonicity a guarantee rather than a hope, then up to ``n_iter`` rounds of
    MAD-based outlier removal. Below ``min_points`` pairs there is not enough
    signal to justify a smoother and it degrades to a straight line.

    Returns ``(curve, keep_mask)``, where ``keep_mask`` is indexed in the order
    of the **RT-sorted** pairs — see :func:`fit_robust_curve_sorted` if you need
    to map the mask back onto the input order.
    """
    curve, keep, _ = fit_robust_curve_sorted(rt_x, rt_y, frac=frac, n_iter=n_iter,
                                             mad_k=mad_k, min_points=min_points,
                                             mode=mode)
    return curve, keep


def fit_robust_curve_sorted(rt_x: np.ndarray,
                            rt_y: np.ndarray,
                            frac: float = DEFAULTS["curve_frac"],
                            n_iter: int = DEFAULTS["curve_iter"],
                            mad_k: float = DEFAULTS["curve_mad_k"],
                            min_points: int = DEFAULTS["curve_min_points"],
                            mode: str = DEFAULTS["extrapolate_mode"],
                            ) -> Tuple[MonotoneCurve, np.ndarray, np.ndarray]:
    """:func:`fit_robust_curve` plus the sort order that the mask is indexed by.

    ``order`` satisfies ``keep_in_input_order[order] == keep``, which is what a
    caller needs to write the per-pair ``kept`` flag into ``_pairs.csv``.
    """
    rt_x = np.asarray(rt_x, dtype=float)
    rt_y = np.asarray(rt_y, dtype=float)
    if rt_x.size == 0:
        raise CalibrationError(
            "No matched feature pairs — a retention-time curve cannot be fitted. "
            "The usual causes are a polarity mismatch between the two runs, an "
            "m/z tolerance that is far too tight, or the wrong reference dataset."
        )
    order = np.argsort(rt_x)
    x, y = rt_x[order], rt_y[order]
    keep = np.ones(len(x), dtype=bool)

    if len(x) < min_points:
        coef = np.polyfit(x, y, 1)
        grid = np.linspace(x.min(), x.max(), 32)
        return MonotoneCurve(grid, np.polyval(coef, grid), mode=mode), keep, order

    curve = None
    for _ in range(n_iter):
        xs, ys = x[keep], y[keep]
        est = lowess(ys, xs, frac=frac, it=3, return_sorted=True)
        iso = IsotonicRegression(increasing=True, out_of_bounds="clip")
        y_iso = iso.fit_transform(est[:, 0], est[:, 1])
        curve = MonotoneCurve(est[:, 0], y_iso, mode=mode)
        resid = y - curve.predict(x)
        mad = np.median(np.abs(resid[keep] - np.median(resid[keep])))
        thresh = mad_k * 1.4826 * max(mad, 1e-6)
        new_keep = np.abs(resid - np.median(resid[keep])) <= thresh
        new_keep &= keep  # removal only, no resurrection
        if new_keep.sum() < min_points or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return curve, keep, order


# ---------------------------------------------------------------------------
# Stage 2: anchor refinement
# ---------------------------------------------------------------------------

class AnchorRefiner:
    """Piecewise-linear residual correction from anchor pairs, shrunk by ``lam``.

    ``lam`` runs from 0 (no correction) to 1 (full interpolation) and is chosen
    by leave-one-out over the anchors. The shrinkage is what stops anchor
    residual *noise* from being injected into every calibrated RT when the noise
    rivals the true local bias.

    Used when the anchors carry no class labels; otherwise see
    :class:`ClassAwareRefiner`.
    """

    def __init__(self, anchor_rt_src: np.ndarray, residuals: np.ndarray, lam: float = 1.0):
        order = np.argsort(anchor_rt_src)
        self._x = np.asarray(anchor_rt_src, dtype=float)[order]
        self._r = np.asarray(residuals, dtype=float)[order]
        self.lam = float(lam)

    def correction(self, x, classes=None) -> np.ndarray:
        # np.interp holds the end residuals constant beyond the anchor range
        return self.lam * np.interp(np.asarray(x, dtype=float), self._x, self._r)


class ClassAwareRefiner:
    """Residual correction that separates a class offset from a smooth drift.

    ``m_c`` is the median residual within class ``c`` — the systematic shift of
    a whole lipid class relative to the anchor-free curve, e.g. SM and CE
    against PC and TG on a long gradient. ``g`` is piecewise-linear through the
    class-detrended residuals, i.e. what is left once those offsets are taken
    out::

        correction(x, cls) = lam_c * m[cls] + lam_g * g(x)

    Features whose class is unknown get the global median offset instead of a
    class offset, which is the honest thing to do: they are not evidence for any
    particular class's shift.
    """

    def __init__(self, x, resid, classes, lam_g: float = 1.0, lam_c: float = 1.0):
        x = np.asarray(x, dtype=float)
        r = np.asarray(resid, dtype=float)
        cls = np.asarray(classes)
        order = np.argsort(x)
        self._x = x[order]
        self.m_global = float(np.median(r))
        self.m_class = {c: float(np.median(r[cls == c])) for c in np.unique(cls)}
        detrended = r - np.array([self.m_class[c] for c in cls])
        self._g = detrended[order]
        self.lam_g, self.lam_c = float(lam_g), float(lam_c)

    def correction(self, x, classes=None) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        base = self.lam_g * np.interp(x, self._x, self._g)
        if classes is None:
            return base + self.lam_c * self.m_global
        off = np.array([self.m_class.get(c, self.m_global) for c in classes])
        return base + self.lam_c * off


def loo_residuals_global(x, r, lam: float) -> np.ndarray:
    """Leave-one-anchor-out prediction errors of the plain shrunk correction."""
    x = np.asarray(x, dtype=float)
    r = np.asarray(r, dtype=float)
    order = np.argsort(x)
    x, r = x[order], r[order]
    n = len(x)
    if n < 3:
        return r.copy()
    errs = np.empty(n, dtype=float)
    for i in range(n):
        keep = np.arange(n) != i
        errs[i] = lam * np.interp(x[i], x[keep], r[keep]) - r[i]
    return errs


def loo_mse_global(x, r, lam: float) -> float:
    """LOO MSE of the plain shrunk piecewise-linear correction at ``lam``."""
    return float(np.mean(loo_residuals_global(x, r, lam) ** 2))


def choose_lambda_loo(anchor_rt_src: np.ndarray,
                      residuals: np.ndarray,
                      grid: Optional[np.ndarray] = None) -> float:
    """Pick the shrinkage factor lambda by leave-one-anchor-out MSE."""
    if grid is None:
        grid = LAMBDA_GRID
    x = np.asarray(anchor_rt_src, dtype=float)
    r = np.asarray(residuals, dtype=float)
    if len(x) < 3:
        return 1.0
    best_lam, best_mse = 1.0, np.inf
    for lam in grid:
        mse = loo_mse_global(x, r, lam)
        if mse < best_mse - 1e-12:
            best_mse, best_lam = mse, float(lam)
    return best_lam


def loo_residuals_class(x, r, cls, lam_g: float, lam_c: float) -> np.ndarray:
    """Leave-one-anchor-out prediction errors of the class-aware correction.

    Each fold refits both the class offsets and the detrended piecewise-linear
    component without the held-out anchor, so a held-out singleton class falls
    back to the global median offset — the same situation a real unknown-class
    feature is in.
    """
    x = np.asarray(x, dtype=float)
    r = np.asarray(r, dtype=float)
    cls = np.asarray(cls)
    order = np.argsort(x)
    x, r, cls = x[order], r[order], cls[order]
    n = len(x)
    if n < 3:
        return r.copy()
    errs = np.empty(n, dtype=float)
    for i in range(n):
        keep = np.arange(n) != i
        rk, ck = r[keep], cls[keep]
        m_glob = float(np.median(rk))
        m_cls = {c: float(np.median(rk[ck == c])) for c in np.unique(ck)}
        det = rk - np.array([m_cls[c] for c in ck])
        g = np.interp(x[i], x[keep], det)
        errs[i] = lam_g * g + lam_c * m_cls.get(cls[i], m_glob) - r[i]
    return errs


def loo_mse_class(x, r, cls, lam_g: float, lam_c: float) -> float:
    """LOO MSE of the class-aware correction at fixed ``(lam_g, lam_c)``."""
    return float(np.mean(loo_residuals_class(x, r, cls, lam_g, lam_c) ** 2))


def choose_lambda_loo_class(anchor_rt_src: np.ndarray,
                            residuals: np.ndarray,
                            classes: np.ndarray,
                            grid: Optional[np.ndarray] = None) -> Tuple[float, float]:
    """Pick ``(lam_g, lam_c)`` for the :class:`ClassAwareRefiner` by LOO MSE."""
    if grid is None:
        grid = LAMBDA_GRID
    x = np.asarray(anchor_rt_src, dtype=float)
    r = np.asarray(residuals, dtype=float)
    cls = np.asarray(classes)
    if len(x) < 3:
        return 1.0, 1.0
    best, best_mse = (1.0, 1.0), np.inf
    for lg in grid:
        for lc in grid:
            mse = loo_mse_class(x, r, cls, lg, lc)
            if mse < best_mse - 1e-12:
                best_mse, best = mse, (float(lg), float(lc))
    return best


@dataclass
class ColumnCalibrator:
    """A complete ``RT_source -> RT_reference`` calibrator for one run pair."""

    curve: MonotoneCurve
    source_label: str = ""
    refiner: Optional[object] = None
    anchors: pd.DataFrame = field(default_factory=pd.DataFrame)
    anchors_dropped: pd.DataFrame = field(default_factory=pd.DataFrame)  # sanity filter
    pairs: pd.DataFrame = field(default_factory=pd.DataFrame)   # incl. `source`, `kept`
    n_pairs: int = 0
    n_pairs_kept: int = 0
    n_pairs_standards: int = 0
    n_pairs_sample: int = 0
    lam_g: float = 1.0                  # shrinkage on the smooth component
    lam_c: float = 0.0                  # shrinkage on the class-offset component
    gate_mse_reduction: float = np.nan  # LOO MSE reduction at the chosen lambdas
    gate_threshold: float = DEFAULTS["anchor_gate_min_mse_reduction"]
    gate_reason: str = ""
    n_anchors_validated: int = 0
    loo_residuals: np.ndarray = field(default_factory=lambda: np.empty(0))
    sigma_anchor: float = 0.0
    sigma_window_pairs: int = DEFAULTS["sigma_window_pairs"]

    @property
    def anchors_used(self) -> bool:
        return self.refiner is not None

    @property
    def class_aware(self) -> bool:
        return isinstance(self.refiner, ClassAwareRefiner)

    @property
    def lam(self) -> float:
        """Single shrinkage factor (the smooth component), for v1-shaped callers."""
        return self.lam_g

    @property
    def warp_source(self) -> str:
        return "curve+anchors" if self.anchors_used else "curve"

    def predict(self, rt, classes: Optional[Sequence] = None) -> np.ndarray:
        rt = np.asarray(rt, dtype=float)
        y = self.curve.predict(rt)
        if self.refiner is not None:
            y = y + self.refiner.correction(rt, classes)
        return y

    def is_extrapolated(self, rt) -> np.ndarray:
        return self.curve.is_extrapolated(rt)

    def residuals(self) -> np.ndarray:
        """Curve residuals of every matched pair, ``rt_ref - f(rt_src)``, in minutes."""
        if self.pairs.empty:
            return np.empty(0)
        return (self.pairs["rt_b"].to_numpy(dtype=float)
                - self.curve.predict(self.pairs["rt_a"].to_numpy(dtype=float)))

    def sigma_curve(self, rt) -> np.ndarray:
        """Local 1-sigma scatter of the matched pairs about the curve, in minutes.

        Robust (1.4826 x MAD) over the ``sigma_window_pairs`` pairs nearest in
        RT, so the uncertainty grows where the pairs thin out or disagree rather
        than being one global number that flatters the sparse ends.
        """
        rt = np.atleast_1d(np.asarray(rt, dtype=float))
        if self.pairs.empty:
            return np.full(rt.shape, np.nan)
        kept = self.pairs
        if "kept" in kept.columns:
            sub = kept[kept["kept"].astype(bool)]
            if len(sub) >= 3:
                kept = sub
        x = kept["rt_a"].to_numpy(dtype=float)
        r = kept["rt_b"].to_numpy(dtype=float) - self.curve.predict(x)
        order = np.argsort(x)
        x, r = x[order], r[order]
        n = len(x)
        w = int(min(max(self.sigma_window_pairs, 3), n))
        # position of each query in the sorted pair RTs -> a window of w pairs
        pos = np.searchsorted(x, rt)
        lo = np.clip(pos - w // 2, 0, max(n - w, 0))
        out = np.empty(rt.shape, dtype=float)
        for k, s in enumerate(lo):
            seg = r[int(s):int(s) + w]
            out[k] = 1.4826 * float(np.median(np.abs(seg - np.median(seg))))
        return out

    def uncertainty(self, rt) -> np.ndarray:
        """Total 1-sigma on ``Cal_RT_min``: curve scatter and anchor LOO in quadrature."""
        sc = self.sigma_curve(rt)
        return np.sqrt(np.nan_to_num(sc, nan=0.0) ** 2 + self.sigma_anchor ** 2)

    def anchor_table(self) -> pd.DataFrame:
        """The stage-2 anchors with their curve residuals and LOO errors."""
        if self.anchors is None or self.anchors.empty:
            return pd.DataFrame()
        out = self.anchors.copy()
        rt_src = out["rt_src"].to_numpy(dtype=float)
        out["residual_min"] = out["rt_ref"].to_numpy(dtype=float) - self.curve.predict(rt_src)
        if self.loo_residuals.size == len(out):
            # loo_residuals are computed in RT-sorted order
            order = np.argsort(rt_src)
            loo = np.empty(len(out), dtype=float)
            loo[order] = self.loo_residuals
            out["loo_residual_min"] = loo
        return out

    def to_model(self) -> dict:
        """The ``curve`` + ``anchors`` blocks of ``model.json`` (spec §5.1)."""
        resid = np.abs(self.residuals())
        rt_src = self.pairs["rt_a"].to_numpy(dtype=float) if not self.pairs.empty else np.empty(0)
        rt_ref = self.pairs["rt_b"].to_numpy(dtype=float) if not self.pairs.empty else np.empty(0)
        loo = np.abs(self.loo_residuals)
        return {
            "curve": {
                "n_pairs": int(self.n_pairs),
                "n_pairs_kept": int(self.n_pairs_kept),
                "n_pairs_standards": int(self.n_pairs_standards),
                "n_pairs_sample": int(self.n_pairs_sample),
                "rt_span_src_min": [float(rt_src.min()), float(rt_src.max())] if rt_src.size else [0.0, 0.0],
                "rt_span_ref_min": [float(rt_ref.min()), float(rt_ref.max())] if rt_ref.size else [0.0, 0.0],
                "residual_min": {
                    "median_abs": float(np.median(resid)) if resid.size else 0.0,
                    "p90_abs": float(np.percentile(resid, 90)) if resid.size else 0.0,
                },
            },
            "anchors": {
                "engaged": bool(self.anchors_used),
                "class_aware": bool(self.class_aware),
                "n_validated": int(self.n_anchors_validated),
                "n_used": int(len(self.anchors)) if self.anchors is not None else 0,
                "lam_g": float(self.lam_g),
                "lam_c": float(self.lam_c),
                "gate_mse_reduction": (None if not np.isfinite(self.gate_mse_reduction)
                                       else float(self.gate_mse_reduction)),
                "gate_threshold": float(self.gate_threshold),
                "gate_reason": self.gate_reason,
                "sigma_anchor_min": float(self.sigma_anchor),
                "loo_residual_min": {
                    "median": float(np.median(loo)) if loo.size else 0.0,
                    "p90": float(np.percentile(loo, 90)) if loo.size else 0.0,
                    "max": float(loo.max()) if loo.size else 0.0,
                },
            },
        }


def build_calibrator(std_src: FeatureTable,
                     std_ref: FeatureTable,
                     config: Any = None,
                     *,
                     source_label: str = "",
                     use_anchors: bool = True,
                     sample_anchors: Optional[pd.DataFrame] = None,
                     exclude_anchor: Optional[str] = None,
                     extra_pairs: Optional[pd.DataFrame] = None,
                     class_aware: Optional[bool] = None,
                     min_anchors: Optional[int] = None,
                     min_mse_reduction: Optional[float] = None,
                     frac: Optional[float] = None,
                     mz_tol_ppm: Optional[float] = None,
                     mz_tol_min_da: Optional[float] = None) -> ColumnCalibrator:
    """Build a source-to-reference calibrator from two standards runs.

    Stage 1 always runs: the two standards runs are matched by m/z and a robust
    monotone curve is fitted. ``extra_pairs`` (a frame with ``rt_a``/``rt_b``,
    typically the m/z-matched *sample* features) is merged into that fit, which
    is what extends coverage into the elution regions the mixture does not
    reach.

    Stage 2 is optional. ``sample_anchors`` is a frame of ``label``, ``rt_src``,
    ``rt_ref`` and optionally ``lipid_class`` — endogenous lipids validated in
    both sample runs. With class labels a :class:`ClassAwareRefiner` is fitted,
    otherwise a plain :class:`AnchorRefiner`; ``exclude_anchor`` drops one named
    anchor, which is how honest leave-one-out evaluation is done. The
    correction is then gated: engaged only if its LOO MSE is at least
    ``min_mse_reduction`` below the zero-correction MSE, so a column whose
    anchor residuals do not generalise keeps the anchor-free curve.
    """
    frac = float(cfg(config, "curve_frac", frac))
    n_iter = int(cfg(config, "curve_iter"))
    mad_k = float(cfg(config, "curve_mad_k"))
    min_points = int(cfg(config, "curve_min_points"))
    mode = str(cfg(config, "extrapolate_mode"))
    min_anchors = int(cfg(config, "min_anchors", min_anchors))
    gate = float(cfg(config, "anchor_gate_min_mse_reduction", min_mse_reduction))
    want_class = bool(cfg(config, "class_aware", class_aware))

    std_pairs = match_features_by_mz(std_src, std_ref, config,
                                     mz_tol_ppm=mz_tol_ppm, mz_tol_min_da=mz_tol_min_da)
    if len(std_pairs) < 5:
        raise CalibrationError(
            f"Only {len(std_pairs)} feature pairs matched between the standards run "
            f"and the reference standards run. That is almost always a polarity "
            f"mismatch, or a reference dataset from a different experiment — check "
            f"both before widening mz_tol_ppm."
        )
    std_pairs = std_pairs.assign(source="standards")

    if extra_pairs is not None and len(extra_pairs) > 0:
        extra = extra_pairs.loc[:, PAIR_COLUMNS].assign(source="sample")
        pairs = pd.concat([std_pairs, extra], ignore_index=True)
    else:
        pairs = std_pairs.reset_index(drop=True)

    if len(pairs) < min_points and (extra_pairs is None or len(extra_pairs) == 0):
        raise CalibrationError(
            f"Only {len(pairs)} matched pairs (need {min_points}) and no sample-run "
            f"pairs to fall back on. Check that both runs are the same polarity, "
            f"and that mz_tol_ppm / mz_tol_min_da are not unreasonably tight."
        )

    curve, keep, order = fit_robust_curve_sorted(
        pairs["rt_a"].to_numpy(), pairs["rt_b"].to_numpy(),
        frac=frac, n_iter=n_iter, mad_k=mad_k, min_points=min_points, mode=mode)
    kept_flag = np.zeros(len(pairs), dtype=bool)
    kept_flag[order] = keep
    pairs = pairs.assign(kept=kept_flag)

    cal = ColumnCalibrator(
        curve=curve,
        source_label=source_label,
        pairs=pairs,
        n_pairs=len(pairs),
        n_pairs_kept=int(keep.sum()),
        n_pairs_standards=int((pairs["source"] == "standards").sum()),
        n_pairs_sample=int((pairs["source"] == "sample").sum()),
        gate_threshold=gate,
        sigma_window_pairs=int(cfg(config, "sigma_window_pairs")),
    )

    if not use_anchors or sample_anchors is None or len(sample_anchors) == 0:
        cal.lam_g, cal.lam_c = 1.0, 0.0
        cal.gate_reason = ("anchors disabled" if not use_anchors
                           else "no validated sample anchors")
        return cal

    anchors = sample_anchors.reset_index(drop=True)
    cal.n_anchors_validated = len(anchors)
    if exclude_anchor is not None:
        anchors = anchors[anchors["label"] != exclude_anchor].reset_index(drop=True)
    if len(anchors) < min_anchors:
        cal.anchors = anchors
        cal.lam_g, cal.lam_c = 1.0, 0.0
        cal.gate_reason = (f"only {len(anchors)} anchors validated "
                           f"(min_anchors={min_anchors}) — stage-1 curve only")
        return cal

    resid = anchors["rt_ref"].to_numpy(dtype=float) - curve.predict(
        anchors["rt_src"].to_numpy(dtype=float))
    has_class = want_class and "lipid_class" in anchors.columns

    # Robust anchor sanity filter: drop anchors the curve cannot explain at all.
    # Class-systematic shifts are detrended first, so a genuine class offset is
    # not mistaken for an outlier and thrown away.
    anchors = anchors.assign(dropped_by_sanity_filter=False)
    if len(anchors) >= 5:
        if has_class:
            cls = anchors["lipid_class"].to_numpy()
            detr = resid - np.array([np.median(resid[cls == c]) for c in cls])
        else:
            detr = resid - np.median(resid)
        mad = 1.4826 * np.median(np.abs(detr - np.median(detr))) + 1e-9
        ok = np.abs(detr - np.median(detr)) / mad <= 3.5
        if ok.sum() >= min_anchors and not ok.all():
            dropped = anchors[~ok].reset_index(drop=True)
            cal.anchors_dropped = dropped.assign(dropped_by_sanity_filter=True)
            anchors = anchors[ok].reset_index(drop=True)
            resid = resid[ok]

    x_a = anchors["rt_src"].to_numpy(dtype=float)
    mse0 = float(np.mean(resid ** 2))  # the zero-correction error the gate judges

    # Degenerate case: the curve already lands on every anchor (the clearest
    # instance being a run calibrated against itself, where every residual is
    # exactly 0). The gate ratio 1 - mse1/mse0 is then 0/0, and the guard
    # `max(mse0, 1e-12)` would report a triumphant "100% below no-correction"
    # for a correction that is doing nothing at all. There is genuinely nothing
    # to correct, so say that instead.
    if mse0 < _NEGLIGIBLE_RESID_MSE:
        cal.anchors = anchors
        cal.lam_g, cal.lam_c = 1.0, 0.0
        cal.gate_mse_reduction = 0.0
        cal.gate_reason = (
            f"no correction needed: the stage-1 curve already reproduces all "
            f"{len(anchors)} anchors to within "
            f"{np.sqrt(mse0):.2e} min RMS")
        return cal

    if has_class:
        classes = anchors["lipid_class"].to_numpy()
        lam_g, lam_c = choose_lambda_loo_class(x_a, resid, classes)
        loo = loo_residuals_class(x_a, resid, classes, lam_g, lam_c)
        mse1 = float(np.mean(loo ** 2))
        gate_red = 1.0 - mse1 / max(mse0, 1e-12)
        if gate_red >= gate:
            cal.refiner = ClassAwareRefiner(x_a, resid, classes, lam_g, lam_c)
            cal.lam_g, cal.lam_c = lam_g, lam_c
            cal.gate_reason = (f"class-aware correction engaged: LOO MSE "
                               f"{100 * gate_red:.0f}% below no-correction "
                               f"(threshold {100 * gate:.0f}%)")
        else:
            lam_g, lam_c = 1.0, 0.0  # gated off: pure stage-1 curve
            cal.lam_g, cal.lam_c = lam_g, lam_c
            cal.gate_reason = (f"gated off: LOO MSE only {100 * gate_red:.0f}% below "
                               f"no-correction (threshold {100 * gate:.0f}%) — the "
                               f"anchor residuals do not generalise")
    else:
        lam_g = choose_lambda_loo(x_a, resid)
        loo = loo_residuals_global(x_a, resid, lam_g)
        mse1 = float(np.mean(loo ** 2))
        gate_red = 1.0 - mse1 / max(mse0, 1e-12)
        if gate_red >= gate:
            cal.refiner = AnchorRefiner(x_a, resid, lam=lam_g)
            cal.lam_g, cal.lam_c = lam_g, 0.0
            cal.gate_reason = (f"correction engaged (no class labels): LOO MSE "
                               f"{100 * gate_red:.0f}% below no-correction")
        else:
            cal.lam_g, cal.lam_c = 1.0, 0.0
            cal.gate_reason = (f"gated off: LOO MSE only {100 * gate_red:.0f}% below "
                               f"no-correction (threshold {100 * gate:.0f}%)")

    cal.anchors = anchors
    cal.gate_mse_reduction = gate_red
    cal.loo_residuals = loo
    cal.sigma_anchor = float(np.sqrt(np.mean(loo ** 2))) if cal.anchors_used else 0.0
    return cal


# ---------------------------------------------------------------- internals ---

def _match_arrays(ft: FeatureTable):
    """``(mz, rt_min, rank, row_positions)`` for the matchable rows of a table.

    Candidates are ranked by :meth:`~rt_anchor.io.schema.FeatureTable.sn`, which
    is ``S/N average`` where the format has it and summed abundance where it does
    not. S/N is the quantity the validated implementation ranked by, and it is
    the better one to rank by: it is a property of the peak, whereas a summed
    intensity also carries how much of the sample was injected, so a noisy
    shoulder in a concentrated run can out-rank the real peak.
    """
    mz = ft.mz().to_numpy(dtype=float)
    rt = ft.rt_minutes().to_numpy(dtype=float)
    rank = ft.sn().to_numpy(dtype=float)
    good = np.isfinite(mz) & np.isfinite(rt)
    idx = np.where(good)[0]
    rank = np.where(np.isfinite(rank[idx]), rank[idx], -np.inf)
    return mz[idx], rt[idx], rank, idx


def _best_match(mz_from: np.ndarray, tol_from: np.ndarray,
                mz_to: np.ndarray, rank_to: np.ndarray) -> np.ndarray:
    """For each element of ``from``, the highest-ranked in-window index of ``to``."""
    order = np.argsort(mz_to)
    mz_sorted = mz_to[order]
    out = np.full(len(mz_from), -1, dtype=int)
    for i, (m, tol) in enumerate(zip(mz_from, tol_from)):
        lo = np.searchsorted(mz_sorted, m - tol)
        hi = np.searchsorted(mz_sorted, m + tol, side="right")
        if hi > lo:
            cand = order[lo:hi]
            out[i] = cand[np.argmax(rank_to[cand])]
    return out
