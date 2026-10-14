"""Match the manually vetted lipid list against the six positive-mode serum tables.

For every row of ``lipid_vetted_list.xlsx`` the positive-mode precursor m/z is
computed from the formula, and the row counts as detected on a column when any
feature of ``../column_*/samples.csv`` lies within the m/z tolerance (15 ppm by
default). RT plays no part in detection. The script also prints the checks the
README quotes: parsing and mass self-checks, tolerance sensitivity, a decoy
baseline and a coarse RT-consistency check.

    python lipid_detection.py [--tol 15] [--out lipid_detection.csv]

Needs numpy, pandas and openpyxl.
"""
import argparse
import collections
import re
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

HERE = Path(__file__).resolve().parent
COLUMNS = ['column_12', 'column_16', 'column_20', 'column_22', 'column_41', 'column_85']
MAIN, SEPARATE_GROUP, FA_GROUP = 'main', 'DMPE/MMPE (separate)', 'FA (not searchable in +)'

# monoisotopic atomic masses and the electron
AM = {'C': 12.0, 'H': 1.00782503207, 'N': 14.0030740048, 'O': 15.99491461956,
      'P': 30.97376163, 'S': 31.97207100, 'Na': 22.9897692809}
E = 0.00054857990946
PROTON = AM['H'] - E
ADDUCT = {  # m/z = M + shift
    '[M+H]+': PROTON,
    '[M+NH4]+': AM['N'] + 4 * AM['H'] - E,
    '[M+Na]+': AM['Na'] - E,
    '[M+H-H2O]+': PROTON - (2 * AM['H'] + AM['O']),
    '[M-H]-': -PROTON,
    '[M+HCOO]-': AM['H'] + AM['C'] + 2 * AM['O'] + E,
}

# positive-mode precursor for each class of the negative-mode sheet
POS_ION = {c: '[M+H]+' for c in (
    'PC', 'EtherPC', 'LPC', 'EtherLPC', 'SM', 'PE', 'EtherPE', 'LPE', 'EtherLPE', 'PS',
    'Cer_NS', 'Cer_NDS', 'Cer_NP', 'HexCer_NS', 'HexCer_HS', 'HexCer_HDS', 'Hex2Cer',
    'BileAcid', 'DMPE', 'MMPE')}
POS_ION.update({c: '[M+NH4]+' for c in ('PA', 'LPA', 'PG', 'LPG', 'PI', 'LPI', 'SSulfate')})
# most co-elute with the same-chain PC (PC [M-CH3]- in-source fragments?), so they
# are reported on their own
SEPARATE = {'DMPE', 'MMPE'}

# neutral formula from sum composition (n carbons, m double bonds in the chains)
SPECIES_FORMULA = {
    'FA': lambda n, m: f'C{n}H{2*n - 2*m}O2',
    'CAR': lambda n, m: f'C{n + 7}H{2*n - 2*m + 13}NO4',
    'MG': lambda n, m: f'C{n + 3}H{2*n - 2*m + 6}O4',
    'DG': lambda n, m: f'C{n + 3}H{2*n - 2*m + 4}O5',
    'PE P-': lambda n, m: f'C{n + 5}H{2*n - 2*m + 10}NO7P',
    'CE': lambda n, m: f'C{n + 27}H{2*n - 2*m + 44}O2',
    'TG': lambda n, m: f'C{n + 3}H{2*n - 2*m + 2}O6',
}


def composition(formula):
    counts = collections.Counter()
    for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        counts[el] += int(n) if n else 1
    return counts


def mass(formula):
    return sum(AM[el] * n for el, n in composition(formula).items())


def ppm(a, b):
    return (a - b) / b * 1e6


def read_negative_sheet(ws):
    """The 492-row MS-DIAL sheet (negative ESI): one positive ion per row, by class."""
    rows, checks = [], []
    for i, r in enumerate(ws.iter_rows(values_only=True), start=1):
        if not (isinstance(r[0], (int, float)) and r[2] is not None):
            continue   # header, banner, blank and m/z-range rows
        rt, mz, name, adduct, ref_mz, formula, cls = r[1], float(r[2]), r[3], r[4], float(r[9]), r[10], r[11]
        M = mass(formula)
        calc = M + ADDUCT[adduct]
        checks.append((ppm(calc, ref_mz), ppm(mz, calc)))
        ion = POS_ION[cls]
        rows.append(dict(sheet=ws.title, xlsx_row=i, name=name, cls=cls,
                         group=SEPARATE_GROUP if cls in SEPARATE else MAIN, formula=formula,
                         neg_adduct=adduct, rt_vetted=rt, ions={ion: M + ADDUCT[ion]}))
    return rows, np.array(checks)


def read_holpie_sheet(ws):
    """The 182-row manuscript sheet: 14 FA in negative ESI, then 168 in positive ESI.

    Every precursor ion the sheet lists is a target (TG and MG list two); fragment
    and neutral-loss columns are skipped.
    """
    table = [list(r) + [None] * 20 for r in ws.iter_rows(values_only=True)]
    rows, checks = [], []
    precursors, species, polarity = [], None, 'neg'
    for i, r in enumerate(table, start=1):
        if r[0] == 'POSITIVE ESI':
            polarity = 'pos'
        elif r[0] == 'Lipid species level ':
            sub = table[i]   # the row below: m/z measured / m/z theoretical / accuracy
            precursors = []
            for c, h in enumerate(r):
                m = re.match(r'Ion \d+: (\[M[^\]]*\][+-])', str(h)) if h else None
                if m:
                    theo = next(k for k in range(c, len(sub)) if sub[k] == 'm/z theoretical')
                    meas = c if sub[c] == 'm/z measured' else None
                    precursors.append((m.group(1), meas, theo))
        elif isinstance(r[3], (int, float)) and isinstance(r[5], (int, float)) and r[1]:
            species = (r[0] or species).strip()   # isomer rows leave the species blank
            cls, n, m = re.match(r'(FA|CAR|MG|DG|PE P-|CE|TG) ?(\d+):(\d+)', species).groups()
            formula = SPECIES_FORMULA[cls](int(n), int(m))
            M = mass(formula)
            ions = {}
            for adduct, meas, theo in precursors:
                calc = M + ADDUCT[adduct]
                checks.append((species, adduct, ppm(float(r[theo]), calc),
                               ppm(float(r[meas]), calc) if meas is not None else np.nan))
                ions[adduct] = calc
            row = dict(sheet=ws.title, xlsx_row=i, name=f'{species} | {r[1]}', cls=cls, group=MAIN,
                       formula=formula, neg_adduct=None, rt_vetted=r[3], ions=ions)
            if polarity == 'neg':   # FA: [M-H]- only, no positive-mode ion to look for
                row.update(group=FA_GROUP, neg_adduct='[M-H]-', ions={})
            rows.append(row)
    return rows, pd.DataFrame(checks, columns=['species', 'adduct', 'theo_vs_calc', 'meas_vs_calc'])


def load_columns(root):
    feats = {}
    for c in COLUMNS:
        t = pd.read_csv(root / c / 'samples.csv', skiprows=4, low_memory=False)
        feats[c] = (t['Average Mz'].to_numpy(float), t['Average Rt(min)'].to_numpy(float))
    return feats


def candidates(mz, ions, tol, shift_ppm=0.0):
    """Indices of the features within tol ppm of any target ion (shifted for decoys)."""
    idx = set()
    for target in ions.values():
        t = target * (1 + shift_ppm * 1e-6)
        idx.update(np.flatnonzero(np.abs(mz - t) / t * 1e6 <= tol).tolist())
    return sorted(idx)


def n_columns(lip, feats, tol, shift_ppm=0.0):
    return np.array([sum(bool(candidates(feats[c][0], ions, tol, shift_ppm)) for c in COLUMNS)
                     for ions in lip.ions])


def brute_force_n_columns(lip, feats, tol):
    """The same count with plain loops, as an independent check of n_columns."""
    out = []
    for ions in lip.ions:
        out.append(sum(any(abs(f - t) / t * 1e6 <= tol for t in ions.values() for f in feats[c][0].tolist())
                       for c in COLUMNS))
    return np.array(out)


def tally(n):
    return {k: int((np.asarray(n) >= k).sum()) for k in (3, 4, 5, 6)}


def coelutes_with_pc(lip, max_drt=0.1):
    """DMPE/MMPE rows whose same-chain PC (one / two more CH2) is listed within max_drt min."""
    pcs = lip[lip.cls.isin(['PC', 'EtherPC'])]
    pc_comp = [composition(f) for f in pcs.formula]
    sep = lip[lip.group == SEPARATE_GROUP]
    n = 0
    for _, r in sep.iterrows():
        k = 1 if r.cls == 'DMPE' else 2
        pc = composition(r.formula) + collections.Counter({'C': k, 'H': 2 * k})
        rts = pcs.rt_vetted.to_numpy(float)[[c == pc for c in pc_comp]]
        n += bool(len(rts)) and np.abs(rts - float(r.rt_vetted)).min() <= max_drt
    return n, len(sep)


def lowess(x, y, xq, frac=0.3, iters=4):
    """Robust local-linear smoother (tricube kernel, bisquare robustness weights)."""
    k = max(int(np.ceil(frac * len(x))), 5)
    w_rob = np.ones(len(x))

    def fit_at(x0, wr):
        d = np.abs(x - x0)
        h = np.sort(d)[k - 1] + 1e-9
        w = np.clip(1 - (d / h) ** 3, 0, None) ** 3 * wr
        mx, my = (w * x).sum() / w.sum(), (w * y).sum() / w.sum()
        sxx = (w * (x - mx) ** 2).sum()
        b = (w * (x - mx) * (y - my)).sum() / sxx if sxx > 0 else 0
        return my + b * (x0 - mx)

    for _ in range(iters):
        r = y - np.array([fit_at(x0, w_rob) for x0 in x])
        u = r / (6 * (np.median(np.abs(r)) + 1e-9))
        w_rob = np.clip(1 - u ** 2, 0, None) ** 2
    return np.array([fit_at(x0, w_rob) for x0 in xq]), r, w_rob


def rt_check(rows, feats, tol, seed=0):
    """Coarse RT plausibility of the m/z hits, per column.

    The vetted RT is mapped onto the column with a robust LOWESS fitted on the
    unambiguous matches (exactly one feature within 5 ppm). A hit is RT-consistent
    when one of its candidates lies within +-max(3 sigma, 0.3 min) of the predicted
    RT. The chance rate asks the same of a randomly drawn other row's predicted RT.
    """
    rng = np.random.default_rng(seed)
    x_all = rows.rt_vetted.to_numpy(float)
    ok = np.zeros((len(rows), len(COLUMNS)), bool)
    report = []
    for j, c in enumerate(COLUMNS):
        mz, rt = feats[c]
        cand = [candidates(mz, ions, tol) for ions in rows.ions]
        c5 = [candidates(mz, ions, 5.0) for ions in rows.ions]
        anchors = [(i, f[0]) for i, f in enumerate(c5) if len(f) == 1]
        pred, r, w = lowess(x_all[[a for a, _ in anchors]], rt[[f for _, f in anchors]], x_all)
        sigma = 1.4826 * np.median(np.abs(r[w > 0.5]))
        win = max(3 * sigma, 0.3)
        hits = chance = 0
        for i, f in enumerate(cand):
            if not f:
                continue
            hits += 1
            ok[i, j] = np.any(np.abs(rt[f] - pred[i]) <= win)
            chance += np.any(np.abs(rt[f] - pred[rng.integers(len(rows))]) <= win)
        report.append(dict(column=c, anchors=len(anchors), sigma_min=sigma, window_min=win,
                           mz_hits=hits, rt_consistent=int(ok[:, j].sum()), chance=chance / hits))
    return ok, pd.DataFrame(report)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--tol', type=float, default=15.0, help='m/z tolerance in ppm (default 15)')
    ap.add_argument('--out', type=Path, default=HERE / 'lipid_detection.csv')
    args = ap.parse_args()

    wb = openpyxl.load_workbook(HERE / 'lipid_vetted_list.xlsx', read_only=True, data_only=True)
    neg_rows, neg_chk = read_negative_sheet(wb.worksheets[0])
    pos_rows, pos_chk = read_holpie_sheet(wb.worksheets[1])
    lip = pd.DataFrame(neg_rows + pos_rows)
    feats = load_columns(HERE.parent)

    print(f'rows: {len(neg_rows)} + {len(pos_rows)} = {len(lip)};', lip.group.value_counts().to_dict())
    print('rows per class, 182 sheet:', lip[lip.sheet == wb.worksheets[1].title].cls.value_counts().to_dict())
    print(f'check, 492 sheet: computed vs MS-DIAL reference m/z max |d| = {np.abs(neg_chk[:, 0]).max():.2f} ppm; '
          f'measured vs computed: 95% within {np.percentile(np.abs(neg_chk[:, 1]), 95):.2f} ppm')
    off = pos_chk[pos_chk.theo_vs_calc.abs() > 1]
    print(f'check, 182 sheet: sheet theoretical vs computed max |d| = {pos_chk.theo_vs_calc.abs().max():.2f} ppm; '
          f'> 1 ppm: {off[["species", "adduct"]].values.tolist()}')

    lip['n_columns_mz'] = n_columns(lip, feats, args.tol)
    same = (brute_force_n_columns(lip, feats, args.tol) == lip.n_columns_mz).all()
    print(f'check, brute-force recount: {"identical" if same else "DIFFERENT"}')
    print('check, DMPE/MMPE co-eluting (<= 0.1 min) with the same-chain PC: %d of %d' % coelutes_with_pc(lip))
    for c in COLUMNS:
        lip[f'{c}_mz'] = [int(bool(candidates(feats[c][0], ions, args.tol))) for ions in lip.ions]
    lip['n_candidate_features_per_column'] = [
        ';'.join(str(len(candidates(feats[c][0], ions, args.tol))) for c in COLUMNS) for ions in lip.ions]
    lip['target_key'] = [tuple(sorted(round(v, 4) for v in ions.values())) for ions in lip.ions]

    print(f'\ndetected on >= k of 6 columns, {args.tol:g} ppm, m/z only')
    for g, sub in lip.groupby('group', sort=False):
        print(f'  {g:26s} rows={len(sub):4d}', tally(sub.n_columns_mz))
    main_rows = lip[lip.group == MAIN].reset_index(drop=True)
    uniq = main_rows.drop_duplicates('target_key')
    print(f'  {"main, distinct targets":26s} n={len(uniq):4d}', tally(uniq.n_columns_mz))

    print('\ntolerance sensitivity and decoy baseline (targets shifted by +-50/75/100 ppm), % of rows')
    for g in (MAIN, SEPARATE_GROUP):
        sub = lip[lip.group == g]
        for tol in (5.0, 10.0, 15.0):
            real = tally(n_columns(sub, feats, tol))
            decoy = collections.Counter()
            for s in (-100, -75, -50, 50, 75, 100):
                decoy.update(tally(n_columns(sub, feats, tol, s)))
            print(f'  {g:22s} {tol:4.0f} ppm  real ' + ' '.join(f'>={k}:{v / len(sub) * 100:5.1f}' for k, v in real.items())
                  + '   decoy ' + ' '.join(f'>={k}:{decoy[k] / 6 / len(sub) * 100:5.1f}' for k in (3, 4, 5, 6)))

    nf = np.array([[int(v) for v in s.split(';')] for s in main_rows.n_candidate_features_per_column])
    print(f'\nmain detections (row x column) with >= 2 candidate features: {(nf >= 2).sum()} of {(nf > 0).sum()}')

    ok, report = rt_check(main_rows, feats, args.tol)
    print('\nRT check (main rows)')
    print(report.round(2).to_string(index=False))
    print('  m/z + RT,', tally(ok.sum(1)))

    main_rows['n_columns_mz_and_RT'] = ok.sum(1)
    for j, c in enumerate(COLUMNS):
        main_rows[f'{c}_rt'] = ok[:, j].astype(int)
    rt_cols = [f'{c}_rt' for c in COLUMNS] + ['n_columns_mz_and_RT']
    out = lip.merge(main_rows[['sheet', 'xlsx_row'] + rt_cols], on=['sheet', 'xlsx_row'], how='left')
    out['positive_target_ions'] = ['; '.join(f'{k} {v:.4f}' for k, v in d.items()) for d in out.ions]
    out = out.rename(columns={'cls': 'class', 'neg_adduct': 'vetted_neg_adduct', 'rt_vetted': 'vetted_RT_min'})
    keep = (['sheet', 'xlsx_row', 'name', 'class', 'group', 'formula', 'vetted_neg_adduct', 'vetted_RT_min',
             'positive_target_ions'] + [f'{c}_mz' for c in COLUMNS]
            + ['n_columns_mz', 'n_candidate_features_per_column'] + rt_cols)
    for c in rt_cols:
        out[c] = out[c].astype('Int64')

    print('\nper class: rows, >= k columns (m/z only), m/z + RT on >= 3')
    for (g, cls), s in out.groupby(['group', 'class'], sort=False):
        rt = int((s.n_columns_mz_and_RT >= 3).sum()) if g == MAIN else '-'
        print(f'  {cls:10s} {len(s):3d}', tally(s.n_columns_mz), rt)
    out[keep].to_csv(args.out, index=False)
    print(f'\nwrote {args.out.name}: {len(out)} rows')


if __name__ == '__main__':
    main()
