"""Rebuild everything in vetted_panel/ from the calibrated tables of the six columns.

Runs, in this order, the scripts that make the deliverables:

  build_vetted_panel.py      targets, picks, panel rules, vetted-RT mapping     -> work/*.csv, logs/
  make_final_table.py        the 222-row table                                  -> vetted_panel_table.csv
  build_figures.py           overview heatmap, ECDF, violins of the 222 rows    -> figures/
  make_metrics_table.py      accuracy before / after calibration                -> vetted_panel_metrics.csv
  build_rt_check_figure.py   vetted RT vs column_20 RT                          -> figures/
  build_ring_heatmaps.py     one circular heatmap per lipid category            -> figures/category_ring_heatmaps/

and then deletes work/, the intermediate tables, so that the folder keeps the two tables. Reads
../../column_*/rt_anchor_run_calibrated.csv (run_all.sh of the parent folder) and ../../../data/.

Not part of this run, since they need many calibration runs (minutes): per_ion_panel.py, holdout_cv.py and
summary_stats.py, the held-out (10-fold) statistics. They read and write work/ as well.

Run:  python run_all.py
"""
import os
import shutil
import subprocess
import sys

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(SCRIPTS)
STEPS = ["build_vetted_panel.py", "make_final_table.py", "build_figures.py", "make_metrics_table.py",
         "build_rt_check_figure.py", "build_ring_heatmaps.py"]

for step in STEPS:
    print(f"\n######## {step}", flush=True)
    subprocess.run([sys.executable, os.path.join(SCRIPTS, step)], check=True)
shutil.rmtree(os.path.join(HERE, "work"), ignore_errors=True)
print("\nrun_all: done; work/ removed")
