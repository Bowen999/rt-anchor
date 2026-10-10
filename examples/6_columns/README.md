# 6 columns

Serum and standards-mixture runs of six LC methods on one QTOF (positive mode),
exported from MS-DIAL (area tables). The number in each folder name is roughly the
last-eluting retention time of that method, in minutes.

Each folder holds `samples.csv` (serum) and `standards.csv` (the Mix 4.4 standards
run on the same column; use `panel="mix21"`). `column_22` is the bundled reference
column.

| folder | serum features | standards features |
|---|---|---|
| `column_12` | 405 | 480 |
| `column_16` | 580 | 405 |
| `column_20` | 483 | 316 |
| `column_22` | 610 | 383 |
| `column_41` | 694 | 405 |
| `column_85` | 549 | 291 |

Lipid names in the tables are m/z-only MS-DIAL guesses; none is MS/MS-confirmed.
