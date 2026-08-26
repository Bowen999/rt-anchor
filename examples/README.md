# Examples

End-to-end examples of `rt-anchor` on bundled Orbitrap datasets — one per standard
mixture panel (pick the same mixture in the desktop app's **Mixture panel**).

- **[`quickstart.ipynb`](quickstart.ipynb)** — Jupyter notebook: load → calibrate → inspect → write outputs (outputs are pre-run so it renders without executing). Runs the **Mix 15** dataset.
- **[`input_15_mixture/`](input_15_mixture/)** — inputs for the **Mix 15** panel (the Caley lipid RT-calibration mix):
  - `samples.txt` — the feature table to calibrate.
  - `standards.txt` — the standards mixture run (builds the anchor template).
- **[`input_21_mixture/input/`](input_21_mixture/input/)** — inputs for the **Mix 21** panel (the extended Mix 4.4 standards):
  - `sample.txt` — the feature table to calibrate.
  - `standards.txt` — the standards mixture run (builds the anchor template).
- **[`output/`](output/)** — what `write_results` produces for the Mix 15 run:
  - `run_calibrated.csv` — original table + appended `RI…` columns.
  - `run_model.json` · `run_anchors.csv` · `run_log.txt` — fitted model, anchors, log.
  - `run_report.html` · `run_report.pdf` — interactive + static report.

Reproduce:

```bash
pip install "rt-anchor[report]"
cd examples
jupyter lab quickstart.ipynb      # or: jupyter nbconvert --to notebook --execute quickstart.ipynb
```

The Mix 15 panel is the engine's built-in default thus the notebook needs no
extra argument. For the Mix 21 dataset, select the panel in the **desktop app**
(Mixture panel → Mix 21); with the engine API, pass the corresponding
`manifest` DataFrame (see `RT Anchor Desktop/rtad/mixtures.py`).
