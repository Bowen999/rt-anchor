"""The iRT scale: a dimensionless 1-100 retention index on the reference column.

Landmarks are the chosen panel's standards, located in the **reference**
standards run. The scale is then::

    iRT(rt_ref) = 1 + 99 * (rt_ref - landmark_first) / (landmark_last - landmark_first)

Be clear about what that is. The reference implementation spaces iRT *linearly
in RT* across the landmarks, so the landmark points are collinear by
construction and the interpolation between them is a straight line: **the v2 iRT
mapping is affine**. The interior landmarks do not bend the scale. They only
certify that the panel was found — that the standards eluted where a panel of
standards should elute. Only the earliest and the latest detected landmark set
the endpoints and therefore the numbers.

This is the deliberate v2 definition, chosen so iRT values stay comparable with
every number produced so far. The alternative — rank-spacing the landmarks
(landmark *k* of *n* maps to ``1 + 99k/(n-1)``) — would make the interior
landmarks load-bearing and the scale robust to gradient shape, at the cost of
breaking that comparability. See §10.1 of the port specification.

Two further consequences worth keeping in mind:

* iRT is a property of the reference column plus the landmark panel, and of
  nothing else. It does **not** require the user's own column to contain the
  panel: the engine computes ``iRT(Cal_RT_min)``, after the cross-column curve
  has already moved the feature onto the reference axis.
* Which panel is chosen moves the scale, because mix15 and mix21 have different
  earliest and latest detected landmarks on the reference run. iRT values are
  comparable within a panel choice, not across one.

Landmark detection reuses :func:`rt_anchor.identify.build_native_template`,
which already does m/z window, intensity pick, and longest-increasing-
subsequence monotonicity enforcement. If fewer than two landmarks survive there
is no scale: iRT is NaN for every feature and reliability is ``"none"``.
``Cal_RT_min`` is unaffected — the two outputs are independent, and a missing
ruler is not a reason to withhold a calibrated retention time.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .identify import build_native_template
from .io.schema import FeatureTable
from .panel import Panel, adduct_mz

IRT_MIN = 1.0
IRT_MAX = 100.0

LANDMARK_COLUMNS = ["name", "class", "mz", "rt_ref_run_min", "iRT"]


class IRTMapper:
    """Landmark iRT scale on the reference column (affine — see the module docstring).

    Constructed from the landmark RTs observed in the reference standards run.
    Outside the landmark span the scale continues at the terminal slope, which
    for an affine scale is simply the same slope: an RT before the first
    landmark gets an iRT below 1, and after the last, above 100. Callers flag
    those rows via the curve's own ``is_extrapolated``.
    """

    def __init__(self, landmark_rt, irt_min: float = IRT_MIN, irt_max: float = IRT_MAX):
        rt = np.sort(np.unique(np.asarray(landmark_rt, dtype=float)))
        rt = rt[np.isfinite(rt)]
        if rt.size < 2:
            raise ValueError(
                f"An iRT scale needs at least 2 distinct landmark RTs, got {rt.size}."
            )
        self.rt = rt
        self.irt = irt_min + (irt_max - irt_min) * (rt - rt[0]) / (rt[-1] - rt[0])
        self._m_lo = (self.irt[1] - self.irt[0]) / (rt[1] - rt[0])
        self._m_hi = (self.irt[-1] - self.irt[-2]) / (rt[-1] - rt[-2])

    # ---- construction ----
    @classmethod
    def from_landmarks(cls, landmarks: pd.DataFrame, **kw) -> "IRTMapper":
        return cls(landmarks["rt_ref_run_min"].to_numpy(dtype=float), **kw)

    # ---- use ----
    def to_irt(self, rt) -> np.ndarray:
        rt = np.atleast_1d(np.asarray(rt, dtype=float))
        v = np.interp(rt, self.rt, self.irt)
        lo = rt < self.rt[0]
        hi = rt > self.rt[-1]
        v[lo] = self.irt[0] + (rt[lo] - self.rt[0]) * self._m_lo
        v[hi] = self.irt[-1] + (rt[hi] - self.rt[-1]) * self._m_hi
        return v

    def slope(self, rt) -> np.ndarray:
        """``d(iRT)/d(rt_ref)`` — the factor that carries an RT sigma into iRT units."""
        rt = np.atleast_1d(np.asarray(rt, dtype=float))
        out = np.empty(rt.shape, dtype=float)
        idx = np.clip(np.searchsorted(self.rt, rt) - 1, 0, len(self.rt) - 2)
        seg = (self.irt[idx + 1] - self.irt[idx]) / (self.rt[idx + 1] - self.rt[idx])
        out[:] = seg
        out[rt < self.rt[0]] = self._m_lo
        out[rt > self.rt[-1]] = self._m_hi
        return out

    def is_extrapolated(self, rt) -> np.ndarray:
        rt = np.atleast_1d(np.asarray(rt, dtype=float))
        return (rt < self.rt[0]) | (rt > self.rt[-1])

    # ---- reporting ----
    @property
    def n_landmarks(self) -> int:
        return int(self.rt.size)

    @property
    def rt_span(self) -> Tuple[float, float]:
        return float(self.rt[0]), float(self.rt[-1])

    def to_model(self, panel_used: str = "") -> dict:
        """The ``irt`` block of ``model.json`` (spec §5.1)."""
        return {
            "n_landmarks": self.n_landmarks,
            "landmark_rt_span_min": [round(self.rt_span[0], 4), round(self.rt_span[1], 4)],
            "irt_range": [float(self.irt[0]), float(self.irt[-1])],
            "definition": "affine on the reference-column RT axis",
            "panel_used": panel_used,
        }


def landmark_panel(manifest: pd.DataFrame,
                   polarity: str,
                   config: Any = None) -> Panel:
    """Resolve a manifest into per-polarity ion targets for landmark detection.

    Deliberately does **not** go through :func:`rt_anchor.panel.build_panel`:
    that function stamps a v1 fixed-affine ``irt`` onto every target from
    ``config.irt_from_rt``, and under v2 the iRT scale is derived from the
    landmarks themselves, not from a pair of scale constants. Landmark detection
    only ever reads ``name`` / ``class`` / ``mz`` / ``rt_ref_min``, so the ``irt``
    column is carried as NaN rather than filled with a number that would mean
    nothing.
    """
    pol = str(polarity).strip().lower()
    col = "adduct_pos" if pol in ("positive", "pos", "+", "p") else "adduct_neg"
    rows = []
    for _, r in manifest.iterrows():
        adduct = r[col]
        mz = adduct_mz(r["exact_mass"], adduct)
        if mz is None:
            continue  # not ionisable in this polarity
        rows.append({
            "name": r["name"], "class": r.get("class", ""), "adduct": adduct,
            "mz": float(mz), "rt_ref_min": float(r["rt_ref_min"]), "irt": np.nan,
            "endogenous": bool(r.get("endogenous", False)),
            "void": bool(r.get("void", False)),
        })
    targets = pd.DataFrame(rows).sort_values("rt_ref_min").reset_index(drop=True)
    return Panel(polarity=("positive" if col == "adduct_pos" else "negative"),
                 targets=targets)


def detect_landmarks(ref_std_table: FeatureTable,
                     panel: Union[Panel, pd.DataFrame],
                     config: Any,
                     polarity: Optional[str] = None) -> pd.DataFrame:
    """Locate the landmark panel in the reference standards run.

    ``panel`` may be a resolved :class:`~rt_anchor.panel.Panel` or a raw manifest
    frame (in which case ``polarity`` — or the table's own — resolves the
    adducts). Detection is :func:`rt_anchor.identify.build_native_template`: an
    m/z window, the most intense candidate across the whole gradient, and a
    longest-increasing-subsequence pass that discards any landmark eluting out
    of panel order.

    Returns ``name, class, mz, rt_ref_run_min``, sorted by RT. ``iRT`` is added
    by :func:`build_irt` once the scale exists — a landmark's own iRT is a
    consequence of the landmark set, not an input to it.
    """
    if not isinstance(panel, Panel):
        pol = polarity or ref_std_table.polarity or "positive"
        panel = landmark_panel(panel, pol, config)
    template = build_native_template(ref_std_table, panel, config)
    meta = panel.targets.set_index("name")
    rows = [{
        "name": name,
        "class": meta.loc[name, "class"] if name in meta.index else "",
        "mz": float(meta.loc[name, "mz"]) if name in meta.index else np.nan,
        "rt_ref_run_min": float(rt),
    } for name, rt in template.items()]
    if not rows:
        return pd.DataFrame(columns=LANDMARK_COLUMNS[:-1])
    return pd.DataFrame(rows).sort_values("rt_ref_run_min").reset_index(drop=True)


def build_irt(ref_std_table: FeatureTable,
              panel: Union[Panel, pd.DataFrame, None],
              config: Any,
              polarity: Optional[str] = None,
              ) -> Tuple[Optional[IRTMapper], pd.DataFrame, str]:
    """Detect landmarks and build the scale in one step.

    Returns ``(mapper, landmarks, reason)``. ``mapper`` is ``None`` when the
    panel is absent or fewer than two landmarks were found; ``reason`` says
    which, in words fit for the run log and the model JSON. ``landmarks`` always
    comes back — an empty or single-row landmark table is itself the diagnosis.
    """
    if panel is None:
        return None, pd.DataFrame(columns=LANDMARK_COLUMNS), \
            "no landmark panel supplied — iRT not computed"
    landmarks = detect_landmarks(ref_std_table, panel, config, polarity=polarity)
    if len(landmarks) < 2:
        return None, landmarks, (
            f"only {len(landmarks)} landmark(s) detected in the reference "
            f"standards run (2 needed) — iRT is NaN; Cal_RT is unaffected")
    mapper = IRTMapper.from_landmarks(landmarks)
    landmarks = landmarks.assign(
        iRT=mapper.to_irt(landmarks["rt_ref_run_min"].to_numpy()))
    return mapper, landmarks, (
        f"{len(landmarks)} landmarks, RT {mapper.rt_span[0]:.2f}-"
        f"{mapper.rt_span[1]:.2f} min on the reference standards run")
