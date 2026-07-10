# Examples

End-to-end example of `rt-anchor` on a bundled Orbitrap dataset.

- **[`quickstart.ipynb`](quickstart.ipynb)** — Jupyter notebook: load → calibrate → inspect → write outputs (outputs are pre-run so it renders without executing).
- **[`input/`](input/)** — the inputs:
  - `samples.txt` — the feature table to calibrate.
  - `standards.txt` — the standard-panel run (builds the anchor template).
- **[`output/`](output/)** — what `write_results` produces:
  - `run_calibrated.csv` — original table + appended `RI…` columns.
  - `run_model.json` · `run_anchors.csv` · `run_log.txt` — fitted model, anchors, log.
  - `run_report.html` · `run_report.pdf` — interactive + static report.

Reproduce:

```bash
pip install "rt-anchor[report]"
cd examples
jupyter lab quickstart.ipynb      # or: jupyter nbconvert --to notebook --execute quickstart.ipynb
```
