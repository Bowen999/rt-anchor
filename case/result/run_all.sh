#!/usr/bin/env bash
# Run RT Anchor 1.2.3 cross-column calibration for all 6 columns with the default settings
# (stage 1 curve + homologous-series term + lattice term; the anchor stage is off by default),
# on the reprocessed serum tables of 2026-10-12 (data/<N min>/Serum/*Area_2*.csv; byte-identical to
# examples/6_columns/column_*/samples.csv on GitHub main) and the unchanged standards runs,
# using column_22 (new serum table + standards) as the reference pair.
# Folder map: column_12 = "12 min", column_16 = "17 min", column_20 = "20 min", column_22 = "22 min",
#             column_41 = "45 min", column_85 = "90 min".
# Outputs go to result_v1.2.3_reprocessed/<column>/.
set -u
PY="/c/Users/li_la/anaconda3/python.exe"
BASE="/c/Users/li_la/OneDrive/Desktop/RT_Calibration"
SRC="${RT_ANCHOR_SRC:-$BASE/software/rt-anchor/rt_anchor/src}"
OUT="$BASE/result_v1.2.3_reprocessed"
D="$BASE/data"
REF_S="$D/22 min/Serum/22 min Serum Area_2_2026_10_12_13_24_24.csv"
REF_T="$D/22 min/STDs/Res/22 min stds Area_1_2026_08_16_23_14_28.csv"

declare -A SAMPLES=(
  [column_12]="$D/12 min/Serum/12 min Serum Area_2_2026_10_12_13_18_31.csv"
  [column_16]="$D/17 min/Serum/17 Serum Area_2_2026_10_12_13_33_45.csv"
  [column_20]="$D/20 min/Serum/20 min Serum Area_2_2026_10_12_12_12_37.csv"
  [column_22]="$D/22 min/Serum/22 min Serum Area_2_2026_10_12_13_24_24.csv"
  [column_41]="$D/45 min/Serum/45 Serum Area_2_2026_10_12_12_48_51.csv"
  [column_85]="$D/90 min/Serum/90 Serum Area_2_2026_10_12_13_11_13.csv"
)
declare -A STANDARDS=(
  [column_12]="$D/12 min/stds/12 min STDs Area_1_2026_08_28_02_03_07.csv"
  [column_16]="$D/17 min/STDs/17 STDs Area_1_2026_08_29_19_34_08.csv"
  [column_20]="$D/20 min/STDs/Results/20 min STDs Area_1_2026_08_17_00_20_03.csv"
  [column_22]="$D/22 min/STDs/Res/22 min stds Area_1_2026_08_16_23_14_28.csv"
  [column_41]="$D/45 min/STDS/45 stds Area_1_2026_08_29_17_12_33.csv"
  [column_85]="$D/90 min/STDs/90 STDS Area_1_2026_08_27_23_47_07.csv"
)

echo "rt_anchor version: $(PYTHONPATH="$SRC" "$PY" -c 'import rt_anchor; print(rt_anchor.__version__, rt_anchor.__file__)')"
overall=0
for col in column_12 column_16 column_20 column_22 column_41 column_85; do
  out="$OUT/$col/rt_anchor_run"
  mkdir -p "$OUT/$col"
  echo "================ $col ================"
  PYTHONPATH="$SRC" "$PY" -m rt_anchor.cli calibrate \
    --samples "${SAMPLES[$col]}" \
    --standards "${STANDARDS[$col]}" \
    --polarity positive \
    --panel mix21 \
    --reference-sample "$REF_S" \
    --reference-standards "$REF_T" \
    --report-format html \
    --out "$out" 2>&1
  rc=$?
  echo "---- $col exit code: $rc ----"
  [ $rc -ne 0 ] && overall=1
done
exit $overall
