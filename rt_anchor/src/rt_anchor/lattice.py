"""Stage 1c — the lattice term.

The homologous-series term (stage 1b, :mod:`rt_anchor.series`) corrects a
feature with the curve residuals of the other matched pairs in its own CH2
ladder — same class, same number of double bonds. That needs at least
``series_min_members`` members per ladder, and some lipid classes never offer
it: serum cholesteryl esters, for instance, have at most three members per
double-bond count (CE 16:1 / 17:1 / 18:1), so the series term never touches
them and their class-wide offset against the stage-1 curve stays in the
calibrated RT. Yet the class as a whole has a dozen members — spread over
*several* ladders that differ by H2 (one double bond).

The lattice term uses that. CH2 ladders whose Kendrick phases differ by one H2
step are joined into a **lattice family**, and every member gets integer
coordinates ``(a, b)`` — carbons and H2 count relative to the family's
lightest member — read straight off the mass: ``mz = m0 + a * CH2 + b * H2``.
A feature's correction is then a distance-weighted mean of the curve residuals
of the family's members near its own lattice point, with one double bond
counting like two carbons: ``d = |Δa| / 2 + |Δb|``.

It is a **fallback behind the series term, never a replacement**: a feature
the series term corrects keeps that correction and is not touched, and the
lattice gate is judged only on the pairs the series term leaves alone. The two
terms answer the same question at different resolutions — the series term
interpolates along one ladder, which is the better evidence wherever a ladder
is long enough to trust; the lattice term borrows across ladders where it is
not.

The two disciplines are the series term's, one dimension up:

* **Order validation on both columns.** On a reversed-phase column more
  carbons and/or fewer double bonds elutes later; a candidate member whose
  elution order disagrees with its lattice coordinates — on either column —
  is dropped before it can donate a residual.
* **Leave-own-family-out with a reach limit.** A feature is never corrected
  with pairs that co-elute with it (the series term's exclusion window — one
  notion of "co-eluting" for both terms), it borrows only from members within
  ``lattice_max_distance`` lattice steps, and only when at least
  ``lattice_min_neighbours`` of them are in reach.

And finally the same "do no harm" gate: every matched pair the series term
leaves alone is predicted as a leave-own-out query, and the term is applied
only if that cuts those pairs' MSE by at least
``lattice_gate_min_mse_reduction``. The raw correction is always available
regardless of the gate, so a reviewer can inspect exactly what the gate was
judging.
"""

from __future__ import annotations

from typing import Any, List, Optional, Tuple

import numpy as np

from .series import CH2_MASS, KENDRICK_FACTOR, PERIOD, _resolve, kendrick_phase

#: One H2 step (one double bond) in Da, and its length on the Kendrick 14-Da
#: circle: ladders whose phases differ by exactly this are neighbouring rows of
#: one lipid class's lattice.
H2_MASS = CH2_MASS - 12.0
H2_PHASE = H2_MASS * KENDRICK_FACTOR

#: Curve-residual MSE below which there is nothing for the lattice term to
#: correct (the series term's ``_NEGLIGIBLE_RESID_MSE``: it only ever catches
#: the exactly-zero case, e.g. self-calibration).
_NEGLIGIBLE_RESID_MSE = 1e-12


class LatticeTerm:
    """The fitted lattice correction for one run pair.

    Built by :meth:`fit` from the stage-1 pairs that survived MAD trimming.
    ``engaged`` is the gate's decision; :meth:`raw_correction` answers what the
    lattice would do whether or not the gate engaged, and :meth:`correction` is
    the gated version the engine actually applies.
    """

    def __init__(self, mz: np.ndarray, xa: np.ndarray, xb: np.ndarray,
                 r: np.ndarray, series_covered: np.ndarray, *, kmd_tol: float,
                 min_members: int, max_distance: float, min_neighbours: int,
                 max_neighbours: int, exclude_rt_min: float,
                 gate_min_mse_reduction: float, min_covered_pairs: int):
        self.kmd_tol = float(kmd_tol)
        self.min_members = int(min_members)
        self.max_distance = float(max_distance)
        # A correction needs at least one donor. With min_neighbours at 0 a
        # query with nobody in reach would pass the "enough donors" test, and
        # with max_neighbours at 0 no donor would be averaged: either way
        # raw_correction divides 0 by 0 (every weight zero). Both are clamped
        # where they are stored; the defaults (2 and 4) are untouched, and
        # to_model() reports what actually ran.
        self.min_neighbours = max(1, int(min_neighbours))
        self.max_neighbours = max(1, int(max_neighbours))
        self.exclude_rt_min = float(exclude_rt_min)
        self.gate_threshold = float(gate_min_mse_reduction)
        self.min_covered_pairs = int(min_covered_pairs)

        self._mz = mz                      # source m/z of the kept pairs
        self._xa = xa                      # source RT
        self._r = r                        # curve residual, rt_ref - f(rt_src)
        self._series_covered = series_covered   # the gate never judges these
        n = len(mz)
        # one entry per accepted family: (member pair indices in ascending
        # (m/z, index), their a coordinates, their b coordinates, the family's
        # m0 intercept mass)
        self._families: List[Tuple[np.ndarray, np.ndarray, np.ndarray, float]] = []
        self._member = np.zeros(n, dtype=bool)
        self._find_families(xb)
        self.member_mask = self._member

        self.n_pairs_covered = 0
        self.gate_mse_reduction = float("nan")
        self.engaged = False
        self.gate_reason = ""
        self._gate()

    # ---------------------------------------------------------------- fit ---

    @classmethod
    def fit(cls, mz, rt_src, rt_ref, residual, x0, x1, config=None,
            series_covered=None, **overrides) -> "LatticeTerm":
        """Fit the lattice term on the stage-1 pairs that survived MAD trimming.

        ``mz`` / ``rt_src`` / ``rt_ref`` are the pairs' source m/z, source RT
        and reference RT; ``residual`` is the curve residual
        ``rt_ref - curve.predict(rt_src)``; ``x0`` / ``x1`` are the curve's
        source-RT knot span (they set the co-elution exclusion window, which is
        the series term's setting — one notion of "co-eluting" for both terms).
        ``series_covered`` marks the pairs the *engaged* series term already
        corrects (``None`` = none); those pairs may sit in a family and donate,
        but the gate is judged only on the pairs left alone. Parameters come
        from an explicit ``lattice_*`` override, else the config, else
        :data:`rt_anchor.crosscolumn.DEFAULTS`.
        """
        def p(name):
            return _resolve(config, name, overrides.get(name))

        exclude_rt_min = max(float(p("series_exclude_rt_floor_min")),
                             float(p("series_exclude_rt_frac")) * (float(x1) - float(x0)))
        mz = np.asarray(mz, dtype=float)
        xa = np.asarray(rt_src, dtype=float)
        xb = np.asarray(rt_ref, dtype=float)
        r = np.asarray(residual, dtype=float)
        n_in = len(mz)
        if series_covered is None:
            series_covered = np.zeros(n_in, dtype=bool)
        series_covered = np.asarray(series_covered, dtype=bool)
        finite = np.isfinite(mz) & np.isfinite(xa) & np.isfinite(xb) & np.isfinite(r)
        term = cls(mz[finite], xa[finite], xb[finite], r[finite],
                   series_covered[finite],
                   kmd_tol=float(p("lattice_kmd_tol")),
                   min_members=int(p("lattice_min_members")),
                   max_distance=float(p("lattice_max_distance")),
                   min_neighbours=int(p("lattice_min_neighbours")),
                   max_neighbours=int(p("lattice_max_neighbours")),
                   exclude_rt_min=exclude_rt_min,
                   gate_min_mse_reduction=float(p("lattice_gate_min_mse_reduction")),
                   min_covered_pairs=int(p("lattice_min_covered_pairs")))
        # align the membership mask with the arrays as passed in (non-finite
        # pairs were dropped above and are never members)
        member = np.zeros(n_in, dtype=bool)
        member[np.flatnonzero(finite)] = term._member
        term.member_mask = member
        return term

    def _find_families(self, xb: np.ndarray) -> None:
        """Cluster by phase, link clusters one H2 step apart, prune by elution order."""
        mz, xa = self._mz, self._xa
        n = len(mz)
        if n == 0:
            return
        # -- clusters: exactly the series term's step, but every cluster is
        # kept — a ladder too short to be a series can still be one b-row of a
        # lattice family
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
        centres = [_cluster_centre(phase, cl) for cl in clusters]
        cluster_of = np.zeros(n, dtype=int)
        for k, cl in enumerate(clusters):
            cluster_of[cl] = k

        # -- components of the H2-link graph, carrying the H2 index along ----
        seen = np.zeros(len(clusters), dtype=bool)
        b_cluster = np.zeros(len(clusters), dtype=int)
        for s0 in range(len(clusters)):
            if seen[s0]:
                continue
            seen[s0] = True
            queue = [s0]
            comp = [s0]
            while queue:
                u = queue.pop(0)
                for v in range(len(clusters)):
                    if seen[v]:
                        continue
                    s = _h2_step(centres[u], centres[v], self.kmd_tol)
                    if s is None:
                        continue
                    seen[v] = True
                    b_cluster[v] = b_cluster[u] + s
                    queue.append(v)
                    comp.append(v)
            if len(comp) < 2:
                continue                    # a lone ladder is the series term's business
            self._add_family(comp, cluster_of, b_cluster, xb)

    def _add_family(self, comp: List[int], cluster_of: np.ndarray,
                    b_cluster: np.ndarray, xb: np.ndarray) -> None:
        """Read the component's lattice coordinates off the masses and keep the
        family if enough members survive with their elution order intact."""
        mz, xa, r = self._mz, self._xa, self._r
        mem = np.flatnonzero(np.isin(cluster_of, comp))
        o = mem[np.lexsort((mem, mz[mem]))]             # by (mz, index)
        b = b_cluster[cluster_of[o]]
        b = b - b[0]                    # relative to the lightest member
        root = mz[o[0]]
        a = np.rint((mz[o] - root - b * H2_MASS) / CH2_MASS).astype(int)
        err = mz[o] - root - a * CH2_MASS - b * H2_MASS
        ok = np.abs(err - np.median(err)) <= self.kmd_tol
        o, a, b = o[ok], a[ok], b[ok]
        alive = _prune_disorder(xa[o], xb[o], r[o], a, b)
        o, a, b = o[alive], a[alive], b[alive]
        if len(o) < self.min_members or len(np.unique(b)) < 2:
            return
        m0 = float(np.median(mz[o] - a * CH2_MASS - b * H2_MASS))
        self._families.append((o, a, b, m0))
        self._member[o] = True

    # ------------------------------------------------------------- the gate ---

    def _gate(self) -> None:
        """Judge the term by its leave-own-out error on the pairs it would correct.

        Every fitted pair the series term leaves alone is predicted as a query
        at its own (m/z, RT) — which the co-elution exclusion turns into a
        leave-own-family-out prediction — and the term engages only if that
        beats the bare curve by the configured MSE reduction.
        """
        n = len(self._mz)
        if n:
            pred, n_mem = self.raw_correction(self._mz, self._xa)
        else:
            pred = np.empty(0)
            n_mem = np.empty(0, dtype=int)
        cov = (n_mem > 0) & ~self._series_covered
        self.n_pairs_covered = int(cov.sum())
        if self.n_pairs_covered < self.min_covered_pairs:
            self.gate_mse_reduction = float("nan")
            self.gate_reason = (f"only {self.n_pairs_covered} matched pairs the series "
                                f"term leaves alone sit in a lattice family (need "
                                f"{self.min_covered_pairs}) — no lattice correction")
            return
        mse0 = float(np.mean(self._r[cov] ** 2))
        mse1 = float(np.mean((self._r[cov] - pred[cov]) ** 2))
        if mse0 < _NEGLIGIBLE_RESID_MSE:
            # the self-calibration case: the gate ratio would be 0/0, and there
            # is genuinely nothing to correct — say that, not a percentage
            self.gate_mse_reduction = 0.0
            self.gate_reason = ("no correction needed: the stage-1 curve already "
                                "reproduces the lattice pairs")
            return
        red = 1.0 - mse1 / mse0
        self.gate_mse_reduction = red
        self.engaged = bool(red >= self.gate_threshold)
        if self.engaged:
            self.gate_reason = (f"lattice term engaged: leave-own-out MSE {100 * red:.0f}% "
                                f"below curve-only (threshold {100 * self.gate_threshold:.0f}%)")
        else:
            p = round(100 * red)
            if p >= 0:
                self.gate_reason = (f"lattice term gated off: leave-own-out MSE only "
                                    f"{p}% below curve-only "
                                    f"(threshold {100 * self.gate_threshold:.0f}%)")
            else:
                # a measured loss is not "only -35% below" — say it plainly
                self.gate_reason = (f"lattice term gated off: leave-own-out MSE "
                                    f"{-p}% above curve-only "
                                    f"(threshold: {100 * self.gate_threshold:.0f}% below)")

    # --------------------------------------------------------- prediction ---

    @property
    def n_families(self) -> int:
        return len(self._families)

    @property
    def family_sizes(self) -> List[int]:
        """Members per family, in the order found."""
        return [int(len(f[0])) for f in self._families]

    @property
    def family_offsets(self) -> List[float]:
        """The ``m0`` intercept mass per family, in the order found."""
        return [float(f[3]) for f in self._families]

    def raw_correction(self, mz, rt) -> Tuple[np.ndarray, np.ndarray]:
        """The ungated lattice correction: ``(correction_min, n_members)`` per query.

        "No correction" is ``(0.0, 0)``. The location step is vectorised over
        the queries per family and H2 row, and the donor step over each
        family's queries, so correcting a whole feature table never loops in
        Python over features.
        """
        mq = np.atleast_1d(np.asarray(mz, dtype=float))
        xq = np.atleast_1d(np.asarray(rt, dtype=float))
        corr = np.zeros(mq.shape, dtype=float)
        n_mem = np.zeros(mq.shape, dtype=int)
        if not self._families:
            return corr, n_mem
        ok = np.isfinite(mq) & np.isfinite(xq)
        if not ok.any():
            return corr, n_mem
        qidx = np.flatnonzero(ok)
        mqf = mq[ok]
        D = self.max_distance
        reach_b = int(np.floor(D))              # bq reach beyond the family's rows
        reach_a = int(np.floor(2.0 * D))        # aq reach beyond its columns

        # -- locate the lattice point: the smallest |m/z error| over all
        # families (the first found wins ties: family order, then ascending bq)
        best_e = np.full(len(mqf), np.inf)
        best_fam = np.full(len(mqf), -1, dtype=int)
        best_aq = np.zeros(len(mqf), dtype=int)
        best_bq = np.zeros(len(mqf), dtype=int)
        for k, (_, fa, fb, m0) in enumerate(self._families):
            a_lo = int(fa.min()) - reach_a
            a_hi = int(fa.max()) + reach_a
            for bq in range(int(fb.min()) - reach_b, int(fb.max()) + reach_b + 1):
                aq = np.rint((mqf - m0 - bq * H2_MASS) / CH2_MASS).astype(int)
                e = np.abs(mqf - m0 - aq * CH2_MASS - bq * H2_MASS)
                cand = (e <= self.kmd_tol) & (aq >= a_lo) & (aq <= a_hi) & (e < best_e)
                if cand.any():
                    best_e[cand] = e[cand]
                    best_fam[cand] = k
                    best_aq[cand] = aq[cand]
                    best_bq[cand] = bq

        # -- donors: the family's members in reach that do not co-elute ------
        for k, (idx, fa, fb, _) in enumerate(self._families):
            sel = best_fam == k
            if not sel.any():
                continue
            qi = qidx[sel]
            d = (np.abs(fa[None, :] - best_aq[sel][:, None]) / 2.0
                 + np.abs(fb[None, :] - best_bq[sel][:, None]))
            donor = ((np.abs(self._xa[idx][None, :] - xq[qi][:, None]) > self.exclude_rt_min)
                     & (d <= D))
            n_donors = donor.sum(axis=1)
            enough = n_donors >= self.min_neighbours
            if not enough.any():
                continue
            qi = qi[enough]
            d = np.where(donor[enough], d[enough], np.inf)
            nearest = np.argsort(d, axis=1, kind="stable")[:, :self.max_neighbours]
            d_near = np.take_along_axis(d, nearest, axis=1)     # inf past the donors
            w = np.where(np.isfinite(d_near), 1.0 / (d_near + 0.5) ** 2, 0.0)
            corr[qi] = (w * self._r[idx][nearest]).sum(axis=1) / w.sum(axis=1)
            # n_members counts every donor, not only the ones averaged
            n_mem[qi] = n_donors[enough]
        return corr, n_mem

    def correction(self, mz, rt) -> Tuple[np.ndarray, np.ndarray]:
        """The gated correction the engine applies: raw when engaged, zeros otherwise."""
        mq = np.atleast_1d(np.asarray(mz, dtype=float))
        if not self.engaged:
            return np.zeros(mq.shape, dtype=float), np.zeros(mq.shape, dtype=int)
        return self.raw_correction(mz, rt)

    def to_model(self) -> dict:
        """The fit-level part of ``model["lattice"]`` (the pipeline adds the rest)."""
        return {
            "enabled": True,
            "engaged": bool(self.engaged),
            "gate_mse_reduction": (None if not np.isfinite(self.gate_mse_reduction)
                                   else float(self.gate_mse_reduction)),
            "gate_threshold": float(self.gate_threshold),
            "gate_reason": self.gate_reason,
            "n_families": int(self.n_families),
            "n_pairs_covered": int(self.n_pairs_covered),
            "kmd_tol": float(self.kmd_tol),
            "min_members": int(self.min_members),
            "max_distance": float(self.max_distance),
            "min_neighbours": int(self.min_neighbours),
            "max_neighbours": int(self.max_neighbours),
            "exclude_rt_min": float(self.exclude_rt_min),
        }


# ---------------------------------------------------------------- internals ---

def _cluster_centre(phase: np.ndarray, cl: np.ndarray) -> float:
    """A cluster's centre phase: first member's phase plus the median wrapped offset.

    The same construction as the series term's centre, but over every pair of
    the cluster — at this point there is no validated chain to restrict to.
    Offsets are measured from the first member and wrapped into (-7, 7], so a
    cluster straddling the 0/14 seam centres correctly instead of averaging to
    the far side of the circle.
    """
    p0 = phase[cl[0]]
    d = phase[cl] - p0
    d = np.where(d > 7.0, d - PERIOD, d)
    d = np.where(d <= -7.0, d + PERIOD, d)
    return float(np.mod(p0 + np.median(d), PERIOD))


def _circ(d: float) -> float:
    """Circular distance on the 14-Da Kendrick circle (scalar)."""
    m = d % PERIOD
    return min(m, PERIOD - m)


def _h2_step(cu: float, cv: float, tol: float) -> Optional[int]:
    """The H2-index step from cluster ``u`` to cluster ``v``, or ``None``.

    The first ``s`` in (+1, +2, -1, -2) that puts ``v`` exactly ``s`` H2 steps
    from ``u`` on the Kendrick circle: ±1 links neighbouring ladders of one
    class, ±2 catches a ladder whose middle row is missing from the data.
    """
    for s in (1, 2, -1, -2):
        if _circ(cv - cu - s * H2_PHASE) <= tol:
            return s
    return None


def _prune_disorder(xa: np.ndarray, xb: np.ndarray, r: np.ndarray,
                    a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Which candidate members to keep, by elution order on BOTH columns.

    On a reversed-phase column more carbons and/or fewer double bonds elutes
    later, so member ``i`` must elute after member ``j`` on both columns
    whenever ``(a_i, b_i) >= (a_j, b_j)`` with at least one strict. The loop
    removes the most-conflicted member — the one whose residual sits furthest
    from the family median on ties — until no conflict remains. Members on the
    same lattice point never conflict with each other: isomers are allowed.
    """
    alive = np.ones(len(xa), dtype=bool)
    if len(xa) == 0:
        return alive
    above = ((a[:, None] >= a[None, :]) & (b[:, None] >= b[None, :])
             & ((a[:, None] > a[None, :]) | (b[:, None] > b[None, :])))
    elutes_later = (xa[:, None] > xa[None, :]) & (xb[:, None] > xb[None, :])
    conflict = above & ~elutes_later
    conflict = conflict | conflict.T          # the relation is symmetric
    while True:
        n_conf = conflict[np.ix_(alive, alive)].sum(axis=1)
        if n_conf.max() == 0:
            return alive
        cand = np.flatnonzero(alive)[n_conf == n_conf.max()]    # member order
        med = np.median(r[alive])
        drop = cand[np.argmax(np.abs(r[cand] - med))]           # first on ties
        alive[drop] = False
