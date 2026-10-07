# Examples

End-to-end examples of `rt-anchor` (v2, cross-column calibration). Each input folder
holds the two runs the method needs — the sample feature table and a run of the
standards mixture from the same column — named after the panel to pass as
`panel=` (the same choice as the desktop app's **Standards panel** cards).

| folder | contents | panel |
|---|---|---|
| [`mix15/`](mix15/) | Orbitrap, positive mode: `samples.txt` (22,314 features) + `standards.txt` | `mix15` — the 15-standard Caley lipid RT-calibration mix |
| [`mix21/`](mix21/) | QTOF, positive mode: `samples.txt` (610 features) + `standards.txt` — the bundled reference column's own runs | `mix21` — Mix 4.4, the extended 21-standard panel |
| [`output/`](output/) | what `write_results` writes for the Mix 15 run | |

Under 1.2.3 the Mix 15 example shows **both gates declining**. Stage 1 matches
595 pairs (173 standards + 422 sample) and keeps 433 after trimming. The series
term finds 27 validated series holding 182 of the 433 kept pairs, and its gate,
read on the 148 of them that get a leave-own-out prediction, declines — the
leave-own-out error would be 35% *higher* than the curve alone. The lattice term
finds 12 lattice families holding 243 pairs, and its gate, read on 270 pairs
(the series term corrects none here, so every pair the lattice covers counts),
declines too — the leave-own-out error would be 2% *higher* than the curve
alone. So 0 of the 22,314 features are corrected by either term and the output
is the plain stage-1 curve; stage 2 (opt-in since 1.2.2) is not requested. The
iRT ruler is set by 14 landmarks on the reference standards run, and detection
QC locates 15 of the 15 panel standards in the example's own standards run
(reporting only).

That is the designed behaviour on data where neither the homologous series nor
the lattice families carry a consistent selectivity difference — it is not a
demonstration of either term at work.

**[`quickstart.ipynb`](quickstart.ipynb)** — load → calibrate → inspect → write
outputs, on the Mix 15 example (outputs are pre-run, so it renders without
executing). It ends with the Mix 21 example as a sanity check: those are the
bundled reference column's own sample and standards runs, so calibrating them
against the default reference returns `Cal_RT_min == RT` (to ~1e-11 min), and
both the series and the lattice term report "no correction needed". That makes
Mix 21 a quick installation check rather than a demonstration of a cross-column
shift.

## Output files

| file | content |
|---|---|
| `run_calibrated.csv` | the input table + `Cal_RT_min`, `Cal_RT_uncertainty_min`, `iRT`, `iRT_uncertainty`, `iRT_reliability`, `is_extrapolated`, `calibration_scope`, `warp_source`, `RI_spread`, `n_contributing`, `series_correction_min`, `series_n_members`, `lattice_correction_min`, `lattice_n_members` |
| `run_model.json` | the fitted model: curve, the series and lattice terms' fit and gate decisions, stage-2 anchors and gate decision, iRT definition, reference provenance, detection QC, config |
| `run_pairs.csv` | the stage-1 matched pairs, whether each survived outlier trimming, and whether it sits in a validated homologous series (`in_series`) or in a lattice family (`in_lattice`) |
| `run_landmarks.csv` | the panel standards located on the reference standards run (the iRT ruler) |
| `run_anchors.csv` | the stage-2 plasma-lipid anchors and their residuals |
| `run_log.txt` | the run log |
| `run_report.html` · `run_report.pdf` | interactive + static report |

## Reproduce

```bash
pip install "rt-anchor[report]"
cd examples
jupyter nbconvert --to notebook --execute --inplace quickstart.ipynb
```

or, with the CLI:

```bash
rt-anchor calibrate --samples mix15/samples.txt --standards mix15/standards.txt \
                    --polarity positive --panel mix15 --out output/run
```
