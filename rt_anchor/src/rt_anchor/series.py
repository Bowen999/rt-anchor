"""Stage 1b — the homologous-series term.

Why a single monotone curve is not enough. Stage 1 fits one mapping
``RT_source -> RT_reference`` from anonymous m/z-matched pairs, and that is the
right backbone — but two methods can differ in how strongly they retain a
double bond relative to a CH2 group. Two lipids that co-elute on the source
column can then sit more than a minute apart on the reference column, and no
single curve can place both: the error it leaves behind is systematic within a
**homologous series** (same class, same number of double bonds, different chain
length), shared by every member of the series and varying smoothly along RT.

Members of a homologous series differ by CH2, so they can be recognised from
m/z alone — no identities — through the Kendrick mass: rescaling the mass axis
so that CH2 is exactly 14.0 makes every member of a series share one mass
defect, i.e. one position on a 14-Da circle (the *phase*). The series term
corrects a feature with the curve residual of the *other* matched pairs in its
own homologous series, interpolated along the series' elution order.

Two disciplines keep it honest:

* **Carbon-order validation.** A cluster of same-phase pairs only becomes a
  series if its longest chain in which RT rises with carbon number on *both*
  columns is long enough — a cluster whose elution order disagrees between the
  columns is not a series the reference column recognises, and is ignored.
* **Leave-own-family-out.** A feature is never corrected with pairs that
  co-elute with it (within ``exclude_rt_min``), so a pair's own residual — or
  a mis-assigned co-eluter's — cannot leak into its own correction. The end
  members extrapolate at most two CH2 beyond the series (``end_reach_ch2``).

And finally the gate, the same "do no harm" rule as stage 2: every matched
pair is predicted as a leave-own-out query, and the term is applied only if
that cuts the pairs' MSE by at least ``series_gate_min_mse_reduction``. The
raw correction is always available regardless of the gate, so a reviewer can
inspect exactly what the gate was judging.
"""

from __future__ import annotations

from typing import Any, List, Optional, Tuple

import numpy as np

#: Kendrick rescaling: CH2 becomes exactly ``PERIOD`` Da, so the members of one
#: homologous series share a mass defect and ``kendrick_phase`` puts them at one
#: position of the 14-Da circle.
KENDRICK_FACTOR = 14.0 / 14.01565006
CH2_MASS = 14.01565006
PERIOD = 14.0

#: Curve-residual MSE below which there is nothing for the series term to
#: correct (mirrors ``crosscolumn._NEGLIGIBLE_RESID_MSE``: it only ever catches
#: the exactly-zero case, e.g. self-calibration).
_NEGLIGIBLE_RESID_MSE = 1e-12


def kendrick_phase(mz) -> np.ndarray:
    """Position on the Kendrick 14-Da circle, in [0, 14).

    Two CH2 homologues map to the same phase; a different double-bond count or
    a different class (different heteroatom mass defect) maps elsewhere.
    """
    return np.mod(np.asarray(mz, dtype=float) * KENDRICK_FACTOR, PERIOD)


def _circ(p, q):
    """Circular distance between two phases on the 14-Da circle."""
    d = np.abs(p - q)
    return np.minimum(d, PERIOD - d)


def _resolve(config: Any, name: str, override: Any = None) -> Any:
    """One parameter: explicit override, else the config, else the default.

    Defers to :func:`rt_anchor.crosscolumn.cfg` (imported lazily — crosscolumn
    imports this module) so the series parameters behave exactly like the rest
    of the engine's settings.
    """
    from .crosscolumn import cfg
    return cfg(config, name, override)


class SeriesTerm:
    """The fitted homologous-series correction for one run pair.

    Built by :meth:`fit` from the stage-1 pairs that survived MAD trimming.
    ``engaged`` is the gate's decision; :meth:`raw_correction` answers what the
    series would do whether or not the gate engaged, and :meth:`correction` is
    the gated version the engine actually applies.
    """

    def __init__(self, mz: np.ndarray, xa: np.ndarray, xb: np.ndarray,
                 r: np.ndarray, *, kmd_tol: float, min_members: int,
                 end_reach_ch2: int, exclude_rt_min: float,
                 gate_min_mse_reduction: float, min_covered_pairs: int):
        self.kmd_tol = float(kmd_tol)
        self.min_members = int(min_members)
        self.end_reach_ch2 = int(end_reach_ch2)
        self.exclude_rt_min = float(exclude_rt_min)
        self.gate_threshold = float(gate_min_mse_reduction)
        self.min_covered_pairs = int(min_covered_pairs)
        # a query may hang at most `end_reach_ch2` CH2 off a series end
        # (+1 Da of slack for the mass-defect spread within a cluster)
        self._reach = self.end_reach_ch2 * CH2_MASS + 1.0

        self._mz = mz                      # source m/z of the kept pairs
        self._xa = xa                      # source RT
        self._r = r                        # curve residual, rt_ref - f(rt_src)
        n = len(mz)
        # one entry per validated series: (member pair indices in ascending
        # m/z, the series' centre phase)
        self._series: List[Tuple[np.ndarray, float]] = []
        self._member = np.zeros(n, dtype=bool)
        self._find_series(xb)
        self.member_mask = self._member

        self.n_pairs_covered = 0
        self.gate_mse_reduction = float("nan")
        self.engaged = False
        self.gate_reason = ""
        self._gate()

    # ---------------------------------------------------------------- fit ---

    @classmethod
    def fit(cls, mz, rt_src, rt_ref, residual, x0, x1, config=None,
            **overrides) -> "SeriesTerm":
        """Fit the series term on the stage-1 pairs that survived MAD trimming.

        ``mz`` / ``rt_src`` / ``rt_ref`` are the pairs' source m/z, source RT
        and reference RT; ``residual`` is the curve residual
        ``rt_ref - curve.predict(rt_src)``; ``x0`` / ``x1`` are the curve's
        source-RT knot span (they set the co-elution exclusion window).
        Parameters come from an explicit ``series_*`` override, else the
        config, else :data:`rt_anchor.crosscolumn.DEFAULTS`.
        """
        def p(name):
            return _resolve(config, name, overrides.get(name))

        exclude_rt_min = max(float(p("series_exclude_rt_floor_min")),
                             float(p("series_exclude_rt_frac")) * (float(x1) - float(x0)))
        mz = np.asarray(mz, dtype=float)
        xa = np.asarray(rt_src, dtype=float)
        xb = np.asarray(rt_ref, dtype=float)
        r = np.asarray(residual, dtype=float)
        finite = np.isfinite(mz) & np.isfinite(xa) & np.isfinite(xb) & np.isfinite(r)
        term = cls(mz[finite], xa[finite], xb[finite], r[finite],
                   kmd_tol=float(p("series_kmd_tol")),
                   min_members=int(p("series_min_members")),
                   end_reach_ch2=int(p("series_end_reach_ch2")),
                   exclude_rt_min=exclude_rt_min,
                   gate_min_mse_reduction=float(p("series_gate_min_mse_reduction")),
                   min_covered_pairs=int(p("series_min_covered_pairs")))
        # align the membership mask with the arrays as passed in (non-finite
        # pairs were dropped above and are never members)
        member = np.zeros(len(mz), dtype=bool)
        member[np.flatnonzero(finite)] = term._member
        term.member_mask = member
        return term

    def _find_series(self, xb: np.ndarray) -> None:
        """Phase-cluster the pairs, then keep the clusters with a validated chain."""
        mz, xa = self._mz, self._xa
        n = len(mz)
        if n == 0:
            return
        phase = kendrick_phase(mz)
        order = np.lexsort((np.arange(n), phase))       # by (phase, index)
        sp = phase[order]
        new_cluster = np.ones(n, dtype=bool)
        new_cluster[1:] = np.diff(sp) > self.kmd_tol
        starts = np.flatnonzero(new_cluster)
        clusters = [order[s:e] for s, e in zip(starts, np.append(starts[1:], n))]
        # wrap-around: a cluster straddling the 0/14 seam is split into the
        # first and last cluster of the sorted order; rejoin them
        if len(clusters) >= 2 and sp[0] + PERIOD - sp[-1] <= self.kmd_tol:
            clusters = [np.concatenate([clusters[-1], clusters[0]])] + clusters[1:-1]
        for cl in clusters:
            if len(cl) < self.min_members:
                continue
            o = cl[np.lexsort((cl, mz[cl]))]            # by (mz, index)
            chain = _longest_carbon_chain(o, mz, xa, xb)
            if len(chain) < self.min_members:
                continue
            self._series.append((chain, _series_centre(phase, chain)))
            self._member[chain] = True

    # ------------------------------------------------------------- the gate ---

    def _gate(self) -> None:
        """Judge the term by its leave-own-out error on the fitted pairs.

        Every pair is predicted as a query at its own (m/z, RT) — which the
        co-elution exclusion turns into a leave-own-family-out prediction — and
        the term engages only if that beats the bare curve by the configured
        MSE reduction.
        """
        n = len(self._mz)
        if n:
            pred, n_mem = self.raw_correction(self._mz, self._xa)
        else:
            pred = np.empty(0)
            n_mem = np.empty(0, dtype=int)
        cov = n_mem > 0
        self.n_pairs_covered = int(cov.sum())
        if self.n_pairs_covered < self.min_covered_pairs:
            self.gate_mse_reduction = float("nan")
            self.gate_reason = (f"only {self.n_pairs_covered} matched pairs sit in a "
                                f"validated homologous series (need {self.min_covered_pairs})"
                                f" — stage-1 curve only")
            return
        mse0 = float(np.mean(self._r[cov] ** 2))
        mse1 = float(np.mean((self._r[cov] - pred[cov]) ** 2))
        if mse0 < _NEGLIGIBLE_RESID_MSE:
            # the self-calibration case: the gate ratio would be 0/0, and there
            # is genuinely nothing to correct — say that, not a percentage
            self.gate_mse_reduction = 0.0
            self.gate_reason = ("no correction needed: the stage-1 curve already "
                                "reproduces the series pairs")
            return
        red = 1.0 - mse1 / mse0
        self.gate_mse_reduction = red
        self.engaged = bool(red >= self.gate_threshold)
        if self.engaged:
            self.gate_reason = (f"series term engaged: leave-own-out MSE {100 * red:.0f}% "
                                f"below curve-only (threshold {100 * self.gate_threshold:.0f}%)")
        else:
            self.gate_reason = (f"series term gated off: leave-own-out MSE only "
                                f"{100 * red:.0f}% below curve-only "
                                f"(threshold {100 * self.gate_threshold:.0f}%)")

    # --------------------------------------------------------- prediction ---

    @property
    def n_series(self) -> int:
        return len(self._series)

    def raw_correction(self, mz, rt) -> Tuple[np.ndarray, np.ndarray]:
        """The ungated series correction: ``(correction_min, n_members)`` per query.

        "No correction" is ``(0.0, 0)``. The queries are grouped by their
        nearest series centre, so correcting a whole feature table is a handful
        of vectorised passes, not a Python loop over features against series.
        """
        mq = np.atleast_1d(np.asarray(mz, dtype=float))
        xq = np.atleast_1d(np.asarray(rt, dtype=float))
        corr = np.zeros(mq.shape, dtype=float)
        n_mem = np.zeros(mq.shape, dtype=int)
        if not self._series:
            return corr, n_mem
        ok = np.isfinite(mq) & np.isfinite(xq)
        if not ok.any():
            return corr, n_mem
        pq = kendrick_phase(mq[ok])
        centres = np.array([c for _, c in self._series])
        d = _circ(pq[:, None], centres[None, :])
        nearest = np.argmin(d, axis=1)                  # first on ties
        close = d[np.arange(len(pq)), nearest] <= self.kmd_tol
        qidx_all = np.flatnonzero(ok)
        for s in range(len(self._series)):
            sel = close & (nearest == s)
            if not sel.any():
                continue
            qidx = qidx_all[sel]
            c, m = self._correct_with_series(self._series[s][0], mq[qidx], xq[qidx])
            corr[qidx] = c
            n_mem[qidx] = m
        return corr, n_mem

    def _correct_with_series(self, idx: np.ndarray, mq: np.ndarray,
                             xq: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """The §2.2 rule for the queries assigned to one series.

        The members are the series' chain in ascending m/z; ``M`` is those that
        do not co-elute with the query (leave-own-family-out). A query between
        two members is linearly interpolated in RT; past an end member it takes
        that member's residual, but only within ``end_reach_ch2`` CH2.
        """
        mmz = self._mz[idx]
        mxa = self._xa[idx]
        mr = self._r[idx]
        nm = len(idx)
        keep = np.abs(mxa[None, :] - xq[:, None]) > self.exclude_rt_min   # (nq, nm)
        n_m = keep.sum(axis=1)
        enough = n_m >= self.min_members - 1
        lo = keep & (mmz[None, :] < mq[:, None] - 7.0)
        hi = keep & (mmz[None, :] > mq[:, None] + 7.0)
        has_lo = lo.any(axis=1)
        has_hi = hi.any(axis=1)
        a_idx = np.where(has_lo, nm - 1 - np.argmax(lo[:, ::-1], axis=1), 0)  # last lo
        b_idx = np.where(has_hi, np.argmax(hi, axis=1), 0)                    # first hi

        corr = np.zeros(len(mq), dtype=float)
        n_mem = np.zeros(len(mq), dtype=int)

        both = has_lo & has_hi & enough
        if both.any():
            sel = np.flatnonzero(both)
            xa_a, xa_b = mxa[a_idx[both]], mxa[b_idx[both]]
            xqb = xq[both]
            sits = (xa_a < xqb) & (xqb < xa_b)      # in the series' elution order
            t = (xqb[sits] - xa_a[sits]) / (xa_b[sits] - xa_a[sits])
            val = mr[a_idx[both]][sits] + t * (mr[b_idx[both]][sits] - mr[a_idx[both]][sits])
            corr[sel[sits]] = val
            n_mem[sel[sits]] = n_m[both][sits]

        only_lo = has_lo & ~has_hi & enough
        if only_lo.any():
            sel = np.flatnonzero(only_lo)
            a = a_idx[only_lo]
            ok = (mq[only_lo] - mmz[a] <= self._reach) & (xq[only_lo] > mxa[a])
            corr[sel[ok]] = mr[a[ok]]
            n_mem[sel[ok]] = n_m[only_lo][ok]

        only_hi = ~has_lo & has_hi & enough
        if only_hi.any():
            sel = np.flatnonzero(only_hi)
            b = b_idx[only_hi]
            ok = (mmz[b] - mq[only_hi] <= self._reach) & (xq[only_hi] < mxa[b])
            corr[sel[ok]] = mr[b[ok]]
            n_mem[sel[ok]] = n_m[only_hi][ok]

        return corr, n_mem

    def correction(self, mz, rt) -> Tuple[np.ndarray, np.ndarray]:
        """The gated correction the engine applies: raw when engaged, zeros otherwise."""
        mq = np.atleast_1d(np.asarray(mz, dtype=float))
        if not self.engaged:
            return np.zeros(mq.shape, dtype=float), np.zeros(mq.shape, dtype=int)
        return self.raw_correction(mz, rt)

    def to_model(self) -> dict:
        """The fit-level part of ``model["series"]`` (the pipeline adds the rest)."""
        return {
            "enabled": True,
            "engaged": bool(self.engaged),
            "gate_mse_reduction": (None if not np.isfinite(self.gate_mse_reduction)
                                   else float(self.gate_mse_reduction)),
            "gate_threshold": float(self.gate_threshold),
            "gate_reason": self.gate_reason,
            "n_series": int(self.n_series),
            "n_pairs_covered": int(self.n_pairs_covered),
            "kmd_tol": float(self.kmd_tol),
            "min_members": int(self.min_members),
            "end_reach_ch2": int(self.end_reach_ch2),
            "exclude_rt_min": float(self.exclude_rt_min),
        }


# ---------------------------------------------------------------- internals ---

def _longest_carbon_chain(o: np.ndarray, mz: np.ndarray, xa: np.ndarray,
                          xb: np.ndarray) -> np.ndarray:
    """Longest chain rising with carbon number on BOTH columns, in ascending m/z.

    Dynamic programming over the cluster's pairs ordered by ``(mz, index)``:
    ``best[i]`` is the longest chain ending at ``o[i]``, extended from the
    first ``j`` that strictly improves it — a homologous series must gain
    carbon (> 7 Da per step), and elute later with every CH2 on the source and
    the reference column alike, or it is not a series the mapping can trust.
    """
    m = len(o)
    omz, oxa, oxb = mz[o], xa[o], xb[o]
    best = np.ones(m, dtype=int)
    prev = np.full(m, -1, dtype=int)
    for i in range(m):
        for j in range(i):
            if (omz[i] - omz[j] > 7.0 and oxa[i] > oxa[j] and oxb[i] > oxb[j]
                    and best[j] + 1 > best[i]):
                best[i] = best[j] + 1
                prev[i] = j
    k = int(np.argmax(best))            # first index attaining the max
    chain = []
    while k != -1:
        chain.append(k)
        k = int(prev[k])
    return o[chain[::-1]]


def _series_centre(phase: np.ndarray, chain: np.ndarray) -> float:
    """The series' centre phase: first member's phase plus the median wrapped offset.

    Offsets are measured from the chain's first member and wrapped into
    (-7, 7], so a series straddling the 0/14 seam centres correctly instead of
    averaging to the far side of the circle.
    """
    p0 = phase[chain[0]]
    d = phase[chain] - p0
    d = np.where(d > 7.0, d - PERIOD, d)
    d = np.where(d <= -7.0, d + PERIOD, d)
    return float(np.mod(p0 + np.median(d), PERIOD))
