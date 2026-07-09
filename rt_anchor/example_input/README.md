# rt_anchor — example inputs

Three example inputs for the RI‑calibration method (all **positive mode**). The only value **required per
feature is `m/z` + `RT`**; everything else is either a built‑in default or optional.

| folder | instrument | tier | contents |
|---|---|---|---|
| `Orbitrap/` | Orbitrap | **minimal** | necessary inputs only → runs **per‑project** calibration |
| `QTOF/` | QTOF | **minimal** | necessary inputs only → runs **per‑project** calibration |
| `QTOF_full/` | QTOF | **full** | minimal + every optional input, incl. `single_files` → **per‑sample** calibration |

## Two calibration tiers

- **Per‑project (default):** needs only the aligned tables. Uses the MassCube‑aligned consensus `RT` → one
  warp → one `RI` per feature. Removes cross‑batch / cross‑instrument offset (the main goal).
- **Per‑sample (auto‑upgrade):** if `single_files/` (or `aligned_features.pkl`) is present, each injection is
  self‑anchored in raw RT space, then aggregated → adds per‑injection drift correction + an inter‑sample
  `RI_spread` QC.

## Minimal (`Orbitrap/`, `QTOF/`) — necessary inputs only

```
<folder>/
  meta.json                        # instrument + polarity (required)
  samples/aligned_feature_table.txt    # columns: m/z, RT   (the calibration input)
  standards/aligned_feature_table.txt  # columns: m/z, RT   (standard-panel run → anchors + native template)
```

The **standard manifest** and the **reference baseline** (frozen iRT ruler) are **built‑in defaults** — not
shipped here. Numeric config uses the locked defaults. No `single_files` → per‑project tier.

## Full (`QTOF_full/`) — minimal + all optional inputs

```
QTOF_full/
  meta.json                        # + scale_constants (RT_lo/RT_hi)
  manifest.csv                     # OPTIONAL (Mode B): override the default standard manifest
  reference_baseline.csv           # OPTIONAL: override the default iRT ruler (RT_ref + per-anchor σ + iRT)
  sample_table_with_time.csv       # OPTIONAL: injection timestamps (injection-order drift QC)
  samples/
    aligned_feature_table.txt      # + isotope flags, MS2, per-sample intensities
    single_files/<sample>.txt      # OPTIONAL → per-sample tier. cols: m/z, RT, + peak_height,
                                   #   Gaussian_similarity, asymmetry_factor, is_isotope,
                                   #   is_in_source_fragment, charge, MS2, MS2_scan_id
    aligned_features.pkl           # OPTIONAL: primary (index-based) join for the per-sample tier
  standards/
    aligned_feature_table.txt
    single_files/P_RT-A1.txt       # panel injection 1
    single_files/P_RT-A2.txt       # OPTIONAL replicate 2 → per-anchor σ, reference stability
```

Note: the one optional input **not** present is a **t0 / void marker** — this panel has no unretained
standard, so there is no real data to include for it.
