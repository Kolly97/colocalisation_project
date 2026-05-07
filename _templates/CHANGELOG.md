# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
versions follow [Semantic Versioning](https://semver.org/).

Unreleased changes go under `[Unreleased]`. On release: rename to the version,
add a date, and start a new `[Unreleased]` block above.

## [Unreleased]

### Planned
- Tuning mode: grid over `(cell_thr_method × nuc_thr_method × top_pct)`, results in per‑combination subfolders.
- Per‑cell‑line config files (`configs/Huh7_default.txt`, `configs/VeroE6_default.txt`).
- `02_aggregate_mock_medians.ijm` — successor of `Mock_median_value.ijm`.
- `03_coloc_infected.ijm` — colocalisation of MOI1/MOI5 using Mock thresholds.

## [0.1.0] — 2026-05-07

### Added
- Initial production pipeline `macros/01_mock_top_pct.ijm`.
- Histogram‑based top‑X% statistic (default 1 %) within a cytosol mask.
- Cell mask from cytoplasmic marker (priority: HA568 > HA488 > NS4B568 > dsRNA488).
- Nucleus mask from DAPI (Otsu, Fill Holes).
- Cytosol mask via `Cell − Nucleus` in 8‑bit (clamped subtraction).
- ROI‑based measurement (no image multiplication needed).
- Per‑(timepoint, combo, marker) CSV outputs with full provenance columns.
- Optional binary mask TIFs and QC PNG overlay output.
- Filename schema enforced positionally; only files containing `mock` are processed.
