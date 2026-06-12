# USAGE — colocalisation pipeline macros

Detailed step-by-step usage lives **next to each macro** (every dialog, mode and setting):

- **Stage 1 — Mock background:** [`1_Mock/USAGE.md`](1_Mock/USAGE.md)
  Run this first. Produces per-image cytosol-percentile CSVs; you aggregate the
  **median `p99_9995`** per (cell line × timepoint × combo × marker) → the background values.

- **Stage 2 — Subtract bg + Colocalisation:** [`2_SubBg_Coloc/USAGE.md`](2_SubBg_Coloc/USAGE.md)
  Run this second. Type the Stage-1 medians into its background dialog; it builds per-cell
  cytosol ROIs and runs the Colocalisation Threshold analysis.

## The one thing that links the two stages

Stage 1 measures background on Mock images → you take the **median of the `p99_9995`** column
per group → type those numbers into Stage 2's *Background values* dialog. Everything else is
self-contained per stage.

## Preparation that applies to both

- Work on a **small test folder first** (≈2 images per condition). For Stage 1, deliberately
  include images **with bright artefacts** so you can calibrate `ARTIFACT_UPPER_BOUND` from the
  logged `cyto_max_raw`. For Stage 2, include a **dense field with touching cells** to confirm
  the watershed splits them.
- Images must be **calibrated 16-bit `.tif`**, channel order **C1=marker1, C2=marker2,
  C3=DAPI**.
- Stage 2 automatic ROI mode needs **MorphoLibJ** (IJPB-plugins); both stages otherwise need
  only stock Fiji (+ the bundled *Colocalization Threshold* plugin for Stage 2's coloc step).

See the project [`README.md`](../README.md) for the big picture and [`HANDOFF.md`](HANDOFF.md)
for design rationale + IJM pitfalls.
