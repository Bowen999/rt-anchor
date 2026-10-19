# compare — the vetted lipid list vs. the six columns

How many lipids of the manually vetted list can be found, by precursor m/z, in the
positive-mode serum tables of the six columns (`../column_*/samples.csv`).

| file | content |
|---|---|
| `lipid_vetted_list.xlsx` | the manually vetted (MS/MS) lipid list for this serum. Sheet *492 lipids MSDial format*: negative ESI, MS-DIAL export of the 20-min method. Sheet *182 Holpie manuscript SI format*: 14 FA in negative ESI + 168 lipids in positive ESI |
| `lipid_detection.py` | the matching; prints every number on this page |
| `lipid_detection_high_confidence.csv` | the 601 list rows whose m/z hit can point to the listed lipid: target ions, hit per column (`*_mz`), candidate features per column |
| `lipid_detection_low_confidence.csv` | the other 73 rows, same columns plus `reason` |

## Targets

The positive-mode precursor m/z is recomputed from the formula (theoretical, monoisotopic).

- **Negative-mode sheet**, one ion per class: `[M+H]+` for PC, LPC, SM, PE, LPE, PS,
  ether PC/LPC/PE/LPE, Cer, HexCer, Hex2Cer, bile acids and DMPE/MMPE; `[M+NH4]+` for
  PA, LPA, PG, LPG, PI, LPI and sterol sulfates.
- **Positive-mode sheet**, the precursor ions it lists: CAR and PE P- `[M+H]+`, DG
  `[M+Na]+`, CE `[M+NH4]+`; TG (`[M+NH4]+` or `[M+Na]+`) and MG (`[M+Na]+` or
  `[M+H-H2O]+`), where either ion counts.
- The 14 FA have no positive-mode ion and are not searched.

A row counts as detected on a column when any feature lies within ±15 ppm of a target
ion. Only the precursor m/z is used; RT plays no part.

## Result (15 ppm)

| | rows | ≥ 3 columns | ≥ 4 | ≥ 5 | all 6 |
|---|---|---|---|---|---|
| **all list rows** | 674 | **591** | **572** | **561** | **535** |
| high confidence | 601 | 540 | 523 | 512 | 489 |
| low confidence | 73 | 51 | 49 | 49 | 46 |
| high confidence, distinct target m/z | 352 | 296 | 283 | 273 | 256 |

Rows are list rows: isomers and repeated entries share a target, so the 601
high-confidence rows are 352 distinct target m/z.

## Low confidence

A row is low confidence when an m/z hit cannot point to the listed lipid; `reason` in
the CSV says why, row by row.

| | rows | ≥ 3 | ≥ 4 | ≥ 5 | all 6 |
|---|---|---|---|---|---|
| DMPE/MMPE | 49 | 49 | 48 | 48 | 45 |
| PA/LPA | 10 | 2 | 1 | 1 | 1 |
| FA | 14 | 0 | 0 | 0 | 0 |

- **DMPE/MMPE.** The target is exactly the `[M+H]+` of a PE: DMPE n:m of PE (n+2):m,
  MMPE n:m of PE (n+1):m. 34 of the 49 also lie within 0.1 min, in the negative-mode
  run, of the PC whose `[M−CH3]−` is the same ion — more likely PC in-source fragments
  than DMPE.
- **PA/LPA.** The `[M+NH4]+` target is exactly the `[M+H]+` of the PE/LPE with two fewer
  carbons and one fewer double bond (LPA 18:1 and LPE 16:0 are both 454.2928), so a hit
  is no evidence of PA/LPA. Saturated LPA 16:0 and 18:0 have no such PE/LPE and are in
  the high-confidence file.
- **FA.** No positive-mode ion; not searched.

## Checks

- Parsed rows per class match the sheets' own `COUNT` rows.
- Recomputed m/z vs. MS-DIAL's reference m/z: ≤ 0.13 ppm. Vs. the theoretical m/z of the
  positive-mode sheet: ≤ 0.1 ppm, except three MG `[M+Na]+` values that the sheet lists
  2.5–3 ppm low.
- An independent brute-force recount gives the same count for every row.
- Every PE/LPE named in `reason` has exactly the target's m/z (asserted by the script).

## What an m/z hit does and does not show

- **Mostly not chance.** Every target shifted by ±50/75/100 ppm (decoys), share of the
  high-confidence rows found on ≥ 3 / all 6 columns:

  | tolerance | real ≥ 3 | real all 6 | decoy ≥ 3 | decoy all 6 |
  |---|---|---|---|---|
  | 5 ppm | 87.9% | 79.2% | 15.4% | 7.4% |
  | 10 ppm | 89.5% | 80.4% | 21.5% | 12.0% |
  | 15 ppm | 89.9% | 81.4% | 27.2% | 14.3% |

  At 5–10 ppm almost no real hit is lost and up to half of the chance hits go.
- **Not an identification.** 73% of the high-confidence detections (2330 of 3178,
  row × column) have two or more features within 15 ppm of the target, and m/z alone
  cannot say which one is the lipid. Isomers (sn-position, 1,2- vs. 1,3-DG, TG chain
  combinations, PE O- vs. P-) cannot be told apart.
- **Not detected ≠ absent.** PA, LPA, PI, LPI, PG and sterol sulfates ionize weakly in
  positive mode; short-chain CAR and most LPE are low-abundance in serum.

<details>
<summary>Per class (15 ppm)</summary>

| class | confidence | rows | ≥ 3 | ≥ 4 | ≥ 5 | all 6 |
|---|---|---|---|---|---|---|
| TG | high | 94 | 94 | 94 | 94 | 94 |
| PC | high | 88 | 88 | 87 | 86 | 85 |
| EtherPE | high | 81 | 66 | 61 | 61 | 56 |
| EtherPC | high | 66 | 66 | 66 | 64 | 62 |
| SM | high | 39 | 39 | 38 | 38 | 38 |
| LPC | high | 32 | 32 | 32 | 32 | 31 |
| PE P- | high | 27 | 27 | 27 | 27 | 26 |
| PE | high | 24 | 23 | 23 | 21 | 19 |
| DG | high | 20 | 19 | 19 | 19 | 19 |
| PI | high | 20 | 15 | 14 | 14 | 12 |
| Cer_NS | high | 16 | 14 | 12 | 12 | 10 |
| LPE | high | 15 | 8 | 6 | 5 | 3 |
| CE | high | 13 | 13 | 13 | 13 | 12 |
| CAR | high | 10 | 3 | 3 | 3 | 3 |
| HexCer_NS | high | 9 | 8 | 8 | 6 | 6 |
| EtherLPC | high | 7 | 7 | 7 | 7 | 7 |
| HexCer_HS | high | 6 | 4 | 4 | 4 | 4 |
| LPI | high | 6 | 1 | 1 | 1 | 0 |
| EtherLPE | high | 5 | 4 | 2 | 2 | 0 |
| MG | high | 4 | 3 | 2 | 0 | 0 |
| SSulfate | high | 4 | 0 | 0 | 0 | 0 |
| Cer_NDS | high | 3 | 2 | 2 | 1 | 0 |
| BileAcid | high | 2 | 1 | 0 | 0 | 0 |
| Cer_NP | high | 2 | 0 | 0 | 0 | 0 |
| Hex2Cer | high | 2 | 2 | 2 | 2 | 2 |
| LPA | high | 2 | 0 | 0 | 0 | 0 |
| HexCer_HDS | high | 1 | 0 | 0 | 0 | 0 |
| LPG | high | 1 | 0 | 0 | 0 | 0 |
| PG | high | 1 | 0 | 0 | 0 | 0 |
| PS | high | 1 | 1 | 0 | 0 | 0 |
| DMPE | low | 48 | 48 | 47 | 47 | 44 |
| FA | low | 14 | 0 | 0 | 0 | 0 |
| LPA | low | 8 | 2 | 1 | 1 | 1 |
| PA | low | 2 | 0 | 0 | 0 | 0 |
| MMPE | low | 1 | 1 | 1 | 1 | 1 |

</details>

## Reproduce

```bash
python examples/6_columns/compare/lipid_detection.py
```

`--tol` sets another tolerance; pass `--out-dir` with it, or the 15 ppm CSVs are
overwritten.

Needs numpy, pandas and openpyxl (`pip install "rt-anchor[excel]"`).
