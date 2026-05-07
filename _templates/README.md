# colocalisation_pipeline

Fiji/ImageJ macros for **cytosolic top‑X% intensity quantification** and downstream
**colocalisation analysis** of immunofluorescence images of Mock vs. infected cells
(SARS‑CoV‑2 / HCV systems, Bartenschlager lab).

> **Status:** Work in progress · `v0.1.0` · Production pipeline (single‑threshold) only.
> Tuning mode and colocalisation step are planned for `v0.2`+.

---

## What this does

For each Mock immunofluorescence image, the pipeline:

1. **Splits** the multi‑channel TIF into marker channels using a positional filename schema.
2. **Builds three masks** from the channels:
   - **Cell mask** from a cytoplasmic marker (HA preferred over NS4B/dsRNA in Mock).
   - **Nucleus mask** from DAPI (always channel 3).
   - **Cytosol mask** = `Cell ∧ ¬Nucleus`.
3. **Measures the top X % brightest pixels** inside the cytosol of each marker channel
   via a histogram‑based percentile estimator (default top 1 %).
4. **Writes one CSV row per (image, marker)** containing the top‑pool statistics
   (`mean_top`, `median_top`, `std_top`, `threshold_value`, `n_pixels`, percentiles)
   plus the configuration that produced them (`stat_method`, `top_pct`, threshold methods,
   `macro_version`, `run_id`).

These per‑image values are aggregated downstream into channel‑specific background /
threshold values that drive colocalisation of the infected (MOI1 / MOI5) samples.

## Quick start

1. Install [Fiji](https://imagej.net/software/fiji/downloads) (≥ 2.16.0 / ImageJ 1.54p).
2. Clone this repo:
   ```bash
   git clone https://github.com/<your-user>/colocalisation_pipeline.git
   ```
3. Place your microscopy images in a folder, named according to the
   [filename schema](docs/filename_schema.md):
   `timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif`
4. In Fiji: **Plugins → Macros → Run...** → select `macros/01_mock_top_pct.ijm`.
5. When prompted, choose your image folder. Results land in `<your_folder>/measure_mock/`.

For step‑by‑step instructions and parameter tuning, see **[USAGE.md](USAGE.md)**.

## Repository layout

```
macros/        ← Fiji IJM scripts (the pipeline itself)
configs/       ← per-cell-line parameter files
docs/          ← schemas, diagrams, examples
tests/         ← synthetic test images + expected outputs
analysis/      ← (optional) Python/R for CSV post-processing
```

Raw images and pipeline outputs are **not** stored in this repository — see
[USAGE.md](USAGE.md) for the recommended local folder layout.

## Requirements

| Tool | Version | Notes |
|---|---|---|
| Fiji | ≥ 2.16.0 | bundled ImageJ ≥ 1.54p |
| (optional) StarDist 2D | latest | for improved nucleus segmentation, see [USAGE.md](USAGE.md#optional-stardist) |
| (optional) Coloc 2 | bundled | used by the colocalisation step (planned `v0.3`) |

No additional Update Sites required for `v0.1`.

## Filename schema

Strictly positional, separated by underscores:

```
<timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<coverslip>_<imgIdx>.tif
```

Example: `12h_Huh7_Mock_HA568_dsRNA488_CS1_1.tif`

- `marker1` → channel 1, `marker2` → channel 2, **DAPI is always channel 3**.
- See [docs/filename_schema.md](docs/filename_schema.md) for full rules.

## Citation

If you use this pipeline in a publication, please cite:

> Hildenbrand, K. (2026). *colocalisation_pipeline* (v1.0.0) [Software].
> Zenodo. https://doi.org/<DOI-once-released>

A machine‑readable version is in [`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE) — free to use, modify, and redistribute with attribution.

## Contact

Kolja Hildenbrand — Bartenschlager Lab, Heidelberg University Hospital
Issues / questions: please use the [GitHub issue tracker](https://github.com/<your-user>/colocalisation_pipeline/issues).
