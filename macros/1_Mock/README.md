# Stage 1 — Mock background pipeline (overview)

Fiji/ImageJ macro that measures, on **Mock (uninfected)** immunofluorescence images, the
**Top-X % brightest cytosol pixel statistics** per marker channel. Those per-image values are
aggregated into the per-channel **background** that Stage 2 subtracts from infected images.

> **Script (current):** `01_Mock_pipeline.ijm` — `MACRO_VERSION 0.8.1` (edited in place).
> **Detailed guide:** [`USAGE.md`](USAGE.md) (every mode + setting, step by step).
> **History:** [`CHANGELOG.md`](CHANGELOG.md). **Project context:** [`../HANDOFF.md`](../HANDOFF.md).
> No extra plugins required.

---

## Where it sits

```
raw .tif (Mock + MOI mixed)
   │  ── only files containing "mock" are processed ──
   ▼
THIS SCRIPT  → per-image CSVs (one row per image × marker)
   │
   │  manual aggregation: median of the p99_9995 column
   │  per (cell line × timepoint × combo × marker)
   ▼
Stage 2 (../2_SubBg_Coloc/) — values typed into its background dialog
```

## What it does

For each Mock `.tif`: split channels → build **cell** (autofluorescence), **nucleus** (DAPI)
and **cytosol = cell − nucleus** masks → remove bright artefacts + tiny fragments → measure the
**histogram of each marker channel inside the cytosol** → write the Top-X % stats and a
percentile ladder (`p95 … p99_9999`) to a per-(timepoint × combo × marker) CSV. Also saves
binary mask TIFs, a QC PNG (marker1 + cytosol outline) and an RGB composite JPG for visual
checking.

## Modes (startup)

- **Mode** `filename` / `dialog` — metadata from the file name (with a token-mapping dialog) or
  asked per image.
- **Threshold mode** `automatic` / `manual` — auto-threshold silently, or pause on each mask
  with the Threshold slider.
- **Cell line** `Huh7` / `VeroE6` / `Other` — recorded in the CSV + used as the output prefix.

Full details, the artefact-calibration workflow, all config knobs and a QC checklist are in
[`USAGE.md`](USAGE.md).

## Output

`<input>/measure_mock/<RUN_ID>/` — the CSVs, `masks/`, `qc/`, and a copy of the IJ Log. Then
take the **median `p99_9995`** per group → Stage 2.

## Key configuration (top of the `.ijm`)

| Variable | Default | Meaning |
|---|---|---|
| `TOP_PCT` | `0.01` | Top-X % pool size for the `*_top` diagnostics |
| `CELL_THR_METHOD` / `CELL_THR_FACTOR` | `Li` / `0.5` | cell-mask threshold / permissiveness (<1 = more permissive) |
| `NUC_THR_METHOD` / `NUC_CLOSE_ITER` | `Otsu` / `2` | nucleus-mask threshold / closing |
| `ARTIFACT_UPPER_BOUND` | `4000` | artefact cut (calibrate from the logged `cyto_max_raw`) |
| `MIN_PARTICLE_SIZE` | `200` | drop cytosol fragments smaller than this (px) |
| `BLUR_SIGMA_CELL` / `BLUR_SIGMA_NUC` | `1` / `1` | Gaussian pre-blur (µm) |

After any change, bump `MACRO_VERSION` and add a `CHANGELOG.md` entry (the CSV `macro_version`
column then distinguishes runs).

---

**Contact:** Kolja Hildenbrand — Kolja.Hildenbrand@gmail.com
