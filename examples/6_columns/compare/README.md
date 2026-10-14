# compare — the vetted lipid list vs. the six columns

How many lipids of the manually vetted list can be found, by precursor m/z, in the
positive-mode serum tables of the six columns (`../column_*/samples.csv`).

| file | content |
|---|---|
| `lipid_vetted_list.xlsx` | the manually vetted (MS/MS) lipid list for this serum. Sheet *492 lipids MSDial format*: negative ESI, MS-DIAL export of the 20-min method. Sheet *182 Holpie manuscript SI format*: 14 FA in negative ESI + 168 lipids in positive ESI |
| `lipid_detection.py` | the matching; prints every number on this page |
| `lipid_detection.csv` | one row per list row (674): target ions, hit per column (`*_mz`), candidate features per column, RT check (`*_rt`) |

## Targets

The positive-mode precursor m/z is recomputed from the formula (theoretical, monoisotopic).

- **Negative-mode sheet**, one ion per class: `[M+H]+` for PC, LPC, SM, PE, LPE, PS,
  ether PC/LPC/PE/LPE, Cer, HexCer, Hex2Cer, bile acids and DMPE/MMPE; `[M+NH4]+` for
  PA, LPA, PG, LPG, PI, LPI and sterol sulfates.
- **Positive-mode sheet**, the precursor ions it lists: CAR and PE P- `[M+H]+`, DG
  `[M+Na]+`, CE `[M+NH4]+`; TG (`[M+NH4]+` or `[M+Na]+`) and MG (`[M+Na]+` or
  `[M+H-H2O]+`), where either ion counts.
- The 14 FA have no positive-mode ion and are not searched. DMPE/MMPE are reported
  separately (see below).

A row counts as detected on a column when any feature lies within ±15 ppm of a target
ion. RT plays no part.

## Result (15 ppm)

| | rows | ≥ 3 columns | ≥ 4 | ≥ 5 | all 6 |
|---|---|---|---|---|---|
| main | 611 | **542** | **524** | **513** | **490** |
| DMPE/MMPE | 49 | 49 | 48 | 48 | 45 |
| FA | 14 | 0 | 0 | 0 | 0 |
| main, distinct target m/z | 356 | 296 | 283 | 273 | 256 |
| main, m/z + RT-consistent (coarse, below) | 611 | 509 | 487 | 455 | 393 |

Rows are list rows: isomers and repeated entries share a target, so the 611 main rows
are 356 distinct target m/z.

## Checks

- Parsed rows per class match the sheets' own `COUNT` rows.
- Recomputed m/z vs. MS-DIAL's reference m/z: ≤ 0.13 ppm. Vs. the theoretical m/z of the
  positive-mode sheet: ≤ 0.1 ppm, except three MG `[M+Na]+` values that the sheet lists
  2.5–3 ppm low.
- An independent brute-force recount gives the same count for every row.

## What an m/z hit does and does not show

- **Mostly not chance.** Every target shifted by ±50/75/100 ppm (decoys), share of rows
  found on ≥ 3 / all 6 columns:

  | tolerance | real ≥ 3 | real all 6 | decoy ≥ 3 | decoy all 6 |
  |---|---|---|---|---|
  | 5 ppm | 86.7% | 78.1% | 15.3% | 7.3% |
  | 10 ppm | 88.4% | 79.2% | 21.3% | 11.9% |
  | 15 ppm | 88.7% | 80.2% | 26.9% | 14.2% |

  At 5–10 ppm almost no real hit is lost and up to half of the chance hits go.
- **Not an identification.** 73% of the detections (2337 of 3191, row × column) have two
  or more candidate features at the target m/z, at different RT. Isomers (sn-position,
  1,2- vs. 1,3-DG, TG chain combinations, PE O- vs. P-) cannot be told apart.
- **Exact isobars from the adduct rules.** PA/LPA `[M+NH4]+` is the same ion as the
  `[M+H]+` of the PE/LPE with two fewer carbons and one fewer double bond (LPA 18:1 and
  LPE 16:0 are both 454.2928), so a PA/LPA hit is no evidence of PA/LPA. DMPE n:m
  `[M+H]+` is the same ion as PE (n+2):m `[M+H]+`.
- **DMPE/MMPE.** 34 of the 49 co-elute (≤ 0.1 min) with the same-chain PC in the
  negative-mode run: more likely PC `[M−CH3]−` in-source fragments than DMPE.
- **The RT check is coarse.** Per column, a robust LOWESS maps the vetted RT onto the
  column, fitted on the unambiguous matches (exactly one feature within 15 ppm); a hit is
  RT-consistent when a candidate lies within ±max(3σ, 0.3 min) of the predicted RT. The
  window is ±0.45–0.62 min on columns 12–20 and ±0.9 min on 22, but ±2.7 and ±4.6 min on
  41 and 85, where a random RT passes 31% of the time. One curve per column also
  ignores class-specific selectivity — what rt-anchor's cross-column calibration models.
- **Not detected ≠ absent.** PA, LPA, PI, LPI, PG and sterol sulfates ionize weakly in
  positive mode; short-chain CAR and most LPE are low-abundance in serum.

<details>
<summary>Per class (15 ppm)</summary>

| class | rows | ≥ 3 | ≥ 4 | ≥ 5 | all 6 | m/z + RT ≥ 3 |
|---|---|---|---|---|---|---|
| TG | 94 | 94 | 94 | 94 | 94 | 94 |
| PC | 88 | 88 | 87 | 86 | 85 | 87 |
| EtherPE | 81 | 66 | 61 | 61 | 56 | 50 |
| EtherPC | 66 | 66 | 66 | 64 | 62 | 65 |
| SM | 39 | 39 | 38 | 38 | 38 | 38 |
| LPC | 32 | 32 | 32 | 32 | 31 | 32 |
| PE P- | 27 | 27 | 27 | 27 | 26 | 25 |
| PE | 24 | 23 | 23 | 21 | 19 | 19 |
| PI | 20 | 15 | 14 | 14 | 12 | 13 |
| DG | 20 | 19 | 19 | 19 | 19 | 19 |
| Cer_NS | 16 | 14 | 12 | 12 | 10 | 14 |
| LPE | 15 | 8 | 6 | 5 | 3 | 6 |
| CE | 13 | 13 | 13 | 13 | 12 | 13 |
| CAR | 10 | 3 | 3 | 3 | 3 | 3 |
| LPA | 10 | 2 | 1 | 1 | 1 | 2 |
| HexCer_NS | 9 | 8 | 8 | 6 | 6 | 8 |
| EtherLPC | 7 | 7 | 7 | 7 | 7 | 7 |
| LPI | 6 | 1 | 1 | 1 | 0 | 1 |
| HexCer_HS | 6 | 4 | 4 | 4 | 4 | 4 |
| EtherLPE | 5 | 4 | 2 | 2 | 0 | 3 |
| MG | 4 | 3 | 2 | 0 | 0 | 3 |
| SSulfate | 4 | 0 | 0 | 0 | 0 | 0 |
| Cer_NDS | 3 | 2 | 2 | 1 | 0 | 2 |
| BileAcid | 2 | 1 | 0 | 0 | 0 | 0 |
| Cer_NP | 2 | 0 | 0 | 0 | 0 | 0 |
| Hex2Cer | 2 | 2 | 2 | 2 | 2 | 1 |
| PA | 2 | 0 | 0 | 0 | 0 | 0 |
| HexCer_HDS | 1 | 0 | 0 | 0 | 0 | 0 |
| LPG | 1 | 0 | 0 | 0 | 0 | 0 |
| PG | 1 | 0 | 0 | 0 | 0 | 0 |
| PS | 1 | 1 | 0 | 0 | 0 | 0 |
| DMPE (separate) | 48 | 48 | 47 | 47 | 44 | — |
| MMPE (separate) | 1 | 1 | 1 | 1 | 1 | — |
| FA (not searched) | 14 | 0 | 0 | 0 | 0 | — |

</details>

## Reproduce

```bash
python examples/6_columns/compare/lipid_detection.py
```

`--tol` sets another tolerance; pass `--out` with it, or the 15 ppm CSV is overwritten.

Needs numpy, pandas and openpyxl (`pip install "rt-anchor[excel]"`).
