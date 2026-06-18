# CHANGELOG for 01_Mock_Background_Pipeline.ijm

---

| Version | Short description                                            |
| ------- | :----------------------------------------------------------- |
| v0.1.0  | Initial Mock Top-1% cytosol intensity pipeline with batch processing, cytosol masks, CSV export, and QC outputs. |
| v0.2.0  | Added flexible metadata handling with filename/dialog mode and fixed 16-bit Top-1% histogram calculation. |
| v0.3.0  | Added manual threshold mode for cell and nucleus masks and run-specific CSV export. |
| v0.4.0  | Added artifact exclusion, cytosol particle-size filtering, and cell-line-specific defaults. |
| v0.4.1  | Strengthened artifact handling, added extra high-percentile outputs, and saved full FIJI log files. |
| v0.5.0  | Added composite RGB JPG QC export and improved binary mask stability with repeated BlackBackground pinning. |
| v0.6.0  | Domain-list dialog (markers/combos/timepoints) is now always shown at startup, pre-filled with the workflow standards. |
| V0.6.1  | Removed cell line specific settings for `MIN_PARTICLE_SIZE`and `ARTIFACT_UPPER_BOUND` to make macro more versatile |
| v0.7.0  | Adaptive per-channel artifact bound replacing the fixed global ARTIFACT_UPPER_BOUND; per-image/per-channel threshold logged; CSV `artifact_upper_bound` → `artifact_bound`. (First tried a Q3+K·IQR fence — see v0.7.1.) |
| v0.7.1  | Top-anchored, separation-based artifact rule (`MULT·p99.9`) replacing the IQR fence; cut only when a separated bright population exists. |
| v0.7.2  | Robust trimmed `mean + k·std` artifact fence with a separation gate; CSV second artifact column → `artifact_k`. |
| v0.8.0  | Reverted to a single fixed high `ARTIFACT_UPPER_BOUND` for all channels; added `logChannelMaxInCytosol()` + CSV `cyto_max_raw` to calibrate it. |
| v0.8.1  | Flexible filename tokens: `askTokenMapping()` picks which token holds marker1/marker2/timepoint; cell line is the dialog choice, not parsed. |
| v1.0.0 (current)  | Third cytosol mode **manual draw** (freehand region(s) like the coloc macro), measured **per ROI** and **numbered** in the QC; nucleus include/exclude asked once in the mode dialog (all modes). CSV gains `roi_index`, `exclude_nuc`. |

---

## v1.0.0 (current)
**Major changes**
- **Third cytosol mode: `manual draw`** — draw the cytosol region(s) by hand, the same option the
  coloc macro offers. Useful when threshold segmentation fails (dim/uneven autofluorescence,
  debris-heavy fields) and for matching the supervisor's hand-drawn workflow.
- **Per-ROI measurement** in manual-draw mode — each drawn region is measured separately and
  **numbered** in the QC image (number == `roi_index`), so a CSV row can be matched to its region.
- **Nucleus include/exclude** is now asked **once in the startup mode dialog** (like the coloc
  macro) and applies to **all** modes.

### Added
- The startup **"Threshold mode"** radio becomes a 3-way **"Cytosol definition":**
  `automatic` | `manual threshold` | `manual draw`. The first two are the previous behaviour,
  renamed; `manual draw` is new.
- ***drawCytosolRois()*** — draw one or more FREEHAND regions on an RGB composite of all channels;
  the regions are kept in the ROI Manager for per-region measurement + numbered QC.
- ***processManualDrawRois()*** + ***measureOneDrawnRoi()*** — measure **each** drawn region
  separately: one CSV row per marker per region, with the new **`roi_index`** column (1..n).
- ***saveNumberedRoiQc()*** — one QC composite per image. The outline is the **carved cytosol**
  (union of the cleaned regions, `Cyto_All`), so excluded **artifacts and the nucleus show as
  holes**; each region is numbered (number == `roi_index`).
- **Masks saved in manual-draw mode too** — the cleaned cytosol union, the artifact mask, and the
  nucleus mask are written via `saveMasksTif` (no cell mask in draw mode). `processManualDrawRois`
  accumulates `Cyto_All` (measured-region union) + `Artifact_All` for this.
- ***makeDrawComposite()*** — builds the RGB drawing composite (duplicates the channels so the
  split windows survive the later measurement). Ported from the coloc macro.
- Global **`EXCLUDE_NUC`** (default **ON**) + a checkbox in the startup mode dialog
  (`askModeAndConfig`). Applies to all modes: `automatic` / `manual threshold` subtract the nucleus
  inside `makeCytosolMask()`; `manual draw` subtracts it from each drawn region.
- ***closeIfOpen()*** utility (also from the coloc macro).
- CSV columns **`roi_index`** (after `channel`) and **`exclude_nuc`** (after `threshold_mode`).

### Changed
- `processOneImage()` now dispatches to **`processSingleCytosol()`** (automatic / manual threshold —
  one pooled cytosol per image, `roi_index` = 1, the original v0.8 flow) or
  **`processManualDrawRois()`** (manual draw — per region). **Design:** each drawn region is turned
  into a `Cytosol_Mask` in the *same binary state* as `makeCytosolMask()`, so the shared downstream
  code (max log, artifact + particle cleanup, empty-mask guard, measurement) is **reused unchanged**.
- `makeCytosolMask()` honours `EXCLUDE_NUC` (subtract the nucleus, or keep the whole cell mask).
- CSV `threshold_mode` column also takes the value **`manual_draw`**.

### Notes
- Artifact + particle cleanup **still runs** on hand-drawn regions on purpose: it removes any
  bright dust you included and keeps the background estimate clean and consistent across all modes.
  Flip this only if you want the drawn region measured verbatim.
- A region that becomes empty after nucleus/artifact removal is skipped (no CSV row) but still
  appears numbered in the QC — so a QC number with no matching CSV row = that region was empty.
- **CSV widened by 2 columns** (`roi_index`, `exclude_nuc`). Downstream aggregation that reads by
  column *name* is unaffected; anything reading by fixed position must be updated.

### Removed
- No functionality removed (the previous `automatic` / `manual` threshold modes are unchanged) just renamed.

---

## v0.8.0
**Major change**
- **Reverted** the artifact bound to a single **fixed `ARTIFACT_UPPER_BOUND`** (set high) for all
  channels. The adaptive per-channel attempts (v0.7.0 IQR, v0.7.1 p99.9, v0.7.2 trimmed mean+std)
  were all hard to calibrate on this data.

### Rationale
A high fixed threshold (5000) removes only **extreme** dust/debris and leaves real
autofluorescence untouched, even in brighter images. We accept that a few moderate artifacts stay
in the cytosol, because:
- the **particle-size filter** (`MIN_PARTICLE_SIZE`) removes most of them, and
- the background value used downstream is the **median across ~30 images**, which is robust to a
  few contaminated pixels in any single image.
  Bonus: dropping the per-channel histogram/percentile pass makes the run **faster**.

### Added
- **`logChannelMaxInCytosol()`** — returns AND logs the brightest RAW pixel of each marker channel
  **within the cytosol** per image (computed before artifact removal, so a present artifact still
  registers). Lets you see where artifacts sit (e.g. `max=21340` vs `max=3800`) and pick
  `ARTIFACT_UPPER_BOUND` accordingly. Persisted in `macro_log_<RUN_ID>.txt`.
- **CSV column `cyto_max_raw`** — that per-channel pre-cleanup cytosol max written to every row
  (threaded `logChannelMaxInCytosol → processOneImage → measureAndWrite → appendCsvRow`), so the
  values can be sorted/plotted across the dataset to calibrate `ARTIFACT_UPPER_BOUND`.

### Removed
- `computeArtifactBound()` and the `ARTIFACT_K` / `ARTIFACT_SEP` / `ARTIFACT_TRIM_PCT` knobs.

### Changed
- `makeArtifactMaskFor` / `buildCombinedArtifactMask` / `cleanCytosolMask` no longer take a
  per-channel bound parameter; `measureAndWrite` / `appendCsvRow` now take a `cytoMax` parameter
  (for the `cyto_max_raw` column).
- CSV second artifact column back to **`artifact_upper_bound`** (= the fixed value).

### Notes
- The adaptive code is preserved in the v0.7.0–0.7.2 history above, should per-channel
  artifact handling be revisited (e.g. via the log-scale / gap-detector routes noted under v0.7.2).

---
## v0.8.1
**Major change**
- **Flexible filename tokens** — the macro no longer assumes a fixed token order.

### Added
- ***askTokenMapping()*** — at startup (filename mode) a dialog shows the first filename and lets
  you pick (via radio buttons) which underscore token holds **marker1 / marker2 / timepoint**
  (globals `TOK_M1` / `TOK_M2` / `TOK_TP`). Works with any naming layout.

### Changed
- The **cell line is no longer parsed** from the filename — it is the startup dialog choice
  (`tryParseFilename()` returns the dialog `CELL_LINE`).
- `CELL_LINE` is now an array of options (more flexible cell-line selection).

### Removed
- No functionality removed.

---

## v0.7.2
**Major change**
- Artifact bound switched to a robust **"centre + k·spread"** fence per channel (the user's
  mean+k·std idea, made robust), replacing the v0.7.1 `MULT·p99.9` anchor.

### Why
A `mean + k·std` fence is intuitive and (unlike the v0.7.0 IQR fence) its `std` "feels" the bright
tail, so it sits at a sensible height. But plain mean/std are **not robust**: the artifacts inflate
the very `std` meant to detect them (the *masking* effect), and the fence assumes ~normal data while
fluorescence is heavily skewed. v0.7.2 keeps the idea but fixes both with trimming + a gate.

### Method
```
bulk  = cytosol pixels > 0, EXCLUDING the top ARTIFACT_TRIM_PCT %   (default 1 %)
thr   = mean_bulk + ARTIFACT_K * std_bulk                            (K default 5)
cut only if  max > ARTIFACT_SEP * thr   (separated dust)  else no cut (SEP default 2)
```
- Trimming the top 1 % stops artifacts from corrupting `std_bulk` → stable image-to-image.
- The **separation gate** preserves the "artifact-free image → cut nothing" behaviour.
- **Honest limitation:** trimming makes `std_bulk` underestimate the autofluorescence spread, so
  `thr` can land inside the real tail; the gate is what protects clean images, and `ARTIFACT_K`
  lifts `thr` clear of the autofluorescence. A center+spread value sets *where* "bright" begins,
  never *whether* dust exists — the gate answers the latter.

### Added
- `ARTIFACT_K` (knob, std multiplier), `ARTIFACT_SEP` (separation guard), `ARTIFACT_TRIM_PCT`
  (bulk trim). `computeArtifactBound()` logs `mean_b / std_b / thr / max / decision` per channel.

### Changed
- CSV second artifact column `artifact_mult` → **`artifact_k`** (`artifact_bound` = `thr` or no-cut).

### Removed
- `ARTIFACT_MULT` / `ARTIFACT_ANCHOR_PCTL` (the p99.9-anchor knobs).

### Possible next steps (not implemented)
- Compute the fence on a **log scale** (fluorescence ≈ log-normal → `k` maps to a real tail
  probability).
- If center+spread still over/under-cuts, a true **gap/valley detector** between the
  autofluorescence tail and the dust population.

---

## v0.7.0

**Major change**
- **Adaptive, per-channel artifact bound** replacing the fixed global `ARTIFACT_UPPER_BOUND`.

### Problem (in the fixed-bound approach)
A single `ARTIFACT_UPPER_BOUND` was applied to every channel. But the channels live in very
different intensity regimes: **dsRNA** Mock autofluorescence is genuinely high and its real top
pixels exceeded the bound, so they were wrongly excluded as "artifact" — **deflating the dsRNA
background**. Raising the global bound to spare dsRNA then let real dust survive in the dimmer
**HA / NS4B** channels — **inflating their background**. A fixed absolute number cannot serve all
channels and does not generalise to new samples/markers.

### Fix
An artifact is defined **relative to each channel's own autofluorescence**, not as a fixed
intensity. The bound is computed per image+channel from that channel's **cytosol** pixels.

> **Formula evolution:** v0.7.0 first tried a Tukey/IQR fence `Q3 + K·(Q3−Q1)`. On real data
> this failed: the cytosol is mostly dim with a sparse bright tail, so `Q1/Q3` sat in the dark
> bulk (e.g. `Q1=10 Q3=44`), the `IQR` was tiny, and even `K=10` produced bounds of a few hundred
> — *far below* the real signal — so it cut autofluorescence and even trimmed artifact-free
> images. **v0.7.1 replaced it** with a top-anchored, separation-based rule.

**v0.7.1 rule** — anchor on the autofluorescence top (its 99.9th percentile, robust because dust
is < 0.1 % of cytosol pixels) and only cut when a **separated** bright population exists:

```
cut = ARTIFACT_MULT * p99.9
if (max_intensity > cut)  ->  bound = cut      // artifacts present, cut above
else                      ->  no cut           // bright tail is just autofluorescence
```

This gives both desired behaviours: a sensible high threshold (not buried in the dim bulk), and
**nothing is cut on artifact-free images**.

### Added
- **`ARTIFACT_MULT`** — the single knob (default `3.0`): "a pixel is an artifact if it is more than
  MULT× brighter than the channel's 99.9th percentile." Larger = more conservative.
  (`ARTIFACT_ANCHOR_PCTL`, default 99.9, is the secondary anchor.)
- **`computeArtifactBound()`** — reuses the existing `pctIndex()` histogram percentile code and
  **logs the per-image, per-channel decision**, e.g.
  `artifact bound  C1 dsRNA : p50=120 p99.9=2400 max=61000 | cut@ MULT(3)xp99.9 = 7200 -> CUT above (bound=7200)`.
- Degenerate guards: empty cytosol / `p99.9 == 0` → no cut.

### Changed
- `makeArtifactMaskFor` / `buildCombinedArtifactMask` / `cleanCytosolMask` / `measureAndWrite` /
  `appendCsvRow` now take the per-channel bound as a parameter; `makeArtifactMaskFor` builds an
  empty mask when the bound says "no cut".
- CSV: column `artifact_upper_bound` → **`artifact_bound`** (the per-channel value actually used)
  plus new column **`artifact_mult`**.

### Removed
- The fixed global `ARTIFACT_UPPER_BOUND` variable, and the abandoned `ARTIFACT_K`/IQR fence.

### Notes
- Mock-only. Script 2 (MOI) must NOT upper-cut bright pixels (there the bright pixels are the real
  specific signal), so this fence is deliberately not propagated to the coloc pipeline.

---

## v0.6.1

**Major changes**

- Removed cell line specific settings for `MIN_PARTICLE_SIZE`and `ARTIFACT_UPPER_BOUND`

### Changed

- Changed ARTIFACT_UPPER_BOUND to 2700

### Fixed

- No bug fixes in this version.

### Removed

- Cell line specific settings

---

## v0.6.0
**Major changes**
- Domain-list dialog always shown at startup

### Changed
- The setup dialog asking for **markers**, **combos**, and **timepoints** is now shown on
  every run, pre-filled with the workflow standards (a plain OK keeps the defaults).
  Previously this dialog only appeared when the startup **MODE** was `dialog`.
- **MODE** (`filename` / `dialog`) now controls ONLY the per-image metadata source,
  no longer whether the domain lists are editable.

### Fixed
- No bug fixes in this version.

### Removed
- Removed the `if (MODE == "dialog")` gate around the list-setup dialog.

---

## v0.5.0 
**Major futures**
- Composite JPG QC  
- Binary Mask Stability Fix (invert cell mask)
- Fixed mask bug (mask creation stopped after first image)

### Added
- Added **SAVE_COMPOSITE_JPG** option.
- Added ***saveCompositeJpg()** to export a merged RGB composite QC image.
- Composite JPGs include:
  - marker 1 in red
  - marker 2 in green
  - DAPI in blue
  - cytosol ROI outline
- Added stronger QC support to visually inspect whether artifacts, nuclei, and cytosol masks are handled correctly.

### Changed
- Composite JPG export is now run before mask saving to preserve access to the original channel windows.
- VeroE6 defaults updated:
  - **CELL_THR_FACTOR set to 0.5**
  - **ARTIFACT_UPPER_BOUND** set to 1000
- Huh7 defaults now explicitly set **ARTIFACT_UPPER_BOUND** to 2000.

### Fixed
- Added repeated ***setOption("BlackBackground", true)** calls to prevent Fiji from silently inverting binary masks.
- Improved stability after Merge Channels / Stack to RGB operations.
- Composite JPG generation now duplicates source channels before merging, preventing LUT/display changes from affecting original channel windows.

### Removed
- No major functionality removed in this version.

---

## v0.4.1.

** Major changes**
- Stronger Artifact Handling 
- Extended Percentiles
-  Fiji Log Export to see what happened

### Added
- Added **ARTIFACT_DILATE_ITER** to expand artifact masks around very bright artifact pixels.
- Added stronger artifact exclusion defaults:
  - **ARTIFACT_UPPER_BOUND** lowered from 4000 to 2000.
  - **ARTIFACT_DILATE_ITER** increased from 2 to 20.
- Added additional high-percentile output columns:
  - p99_9995
  - p99_9999
- Added ***saveLogToFile()*** to export the full FIJI Log window as a text file after the run.
- Added artifact settings to startup logging.

### Changed
- Updated VeroE6 defaults:
  - **CELL_THR_FACTOR** changed from 0.3 to 0.4.
  - **BLUR_SIGMA_CELL** changed from 1.5 to 1.
- Updated CSV header and row export to include the new high-percentile values.
- Updated macro version from 0.4.0 to 0.4.1.

### Fixed
- Improved reproducibility by saving the full macro log alongside the CSV results.
- Improved tracking of artifact exclusion settings in the log output.

### Removed
- No major functionality removed in this version.

---

## v0.4.0

**Major Features**

- Artifact Exclusion  via Threshold to exclude signal over **ARTIFACT_UPPER_BOUND**
- Cytosol Mask Cleanup by ecxcluding particles <= **MIN_PARTICLE_SIZE**

### Added
- Added artifact exclusion using an upper-bound intensity threshold.
- Added **ARTIFACT_UPPER_BOUND** to remove pixels with very high intensity in either marker channel from the cytosol mask.
- Added combined artifact-mask generation via ***buildCombinedArtifactMask()***
- Added per-channel artifact-mask creation via ***makeArtifactMaskFor()***
- Added particle-size filtering for the cytosol mask via ***cleanCytosolMask()***
- Added **MIN_PARTICLE_SIZE** to remove small disconnected cytosol fragments.
- Added cell-line-specific defaults via ***applyCellLineDefaults()***
- Added **CELL_THR_FACTOR** to scale the automatic cell-threshold lower bound.
- Added artifact mask saving for visual QC.
- Added CSV columns:
  - cell_thr_factor
  - artifact_upper_bound
  - min_particle_size
  - n_artifacts_excluded

### Changed
- Updated cytosol-mask generation to exclude high-intensity artifacts before measurement.
- Updated cell-mask thresholding to support cell-line-specific threshold scaling.
- Updated QC output to include artifact masks.
- Improved support for smaller VeroE6 cells through more permissive default settings.

### Fixed
- Fixed ***isValidCombo()*** by using an intermediate variable to avoid an IJM macro-language issue.
- Fixed **CELL_LINE** declaration as a single string instead of an array.

### Removed
- No major functionality removed in this version.

---

## v0.3.0

**Major changes** 

- Manual Threshold Mode
- Run-Specific CSV Export (with RUN_ID)

### Added
- Added threshold mode selector at startup:
  - automatic: fully automatic cell and nucleus thresholding.
  - manual: user can inspect and adjust cell and nucleus thresholds before mask conversion.
- Added ***makeCellMaskManual() ***for interactive cell-mask thresholding.
- Added ***makeNucleusMaskManual()*** for interactive nucleus-mask thresholding.
- Added **THR_MODE** to track whether automatic or manual thresholding was used.
- Added threshold_mode column to the CSV output.
- Added run-specific CSV output paths using **MEASURE_DIR_RUN_ID**

### Changed
- Updated ***processOneImage()*** to choose between automatic and manual mask generation depending on **THR_MODE**
- Updated CSV filenames to include the selected cell line.
- Updated startup logging to report the selected threshold mode.

### Fixed
- Improved manual QC workflow by allowing threshold adjustment before converting masks.

### Removed
- No major functionality removed in this version.

---

## v0.2.0

**Major changes**

- Dialog Mode
- Filename Fallback
- Histogram Fix

### Added

- **Added startup mode selector with two modes:**

  - **filename:** automated metadata parsing from the image filename.

  - **dialog:** manual metadata entry for more flexible analysis

- **Added startup configuration dialog if other markers, combi, cell types or timepoints are used:**

  - cell line

  - available markers

  - marker combinations

  - timepoints

- **Added per-image metadata dialog to manually define:**

  - timepoint

  - cell line

  - marker on C1

  - marker on C2

- Added fallback dialog when filename parsing fails in filename mode.
- Added ***tryParseFilename()*** to separate filename parsing from the main image-processing workflow.
- Added ***askImageMetadata()*** for per-image manual metadata input.
- Added ***askParseFailureAction()*** to choose whether to skip an image or enter metadata manually after parsing failure.
- Added utility functions:
  - ***parseCsvString()***
  - ***arrToStr()***
  - ***inArray()***
- Added **NUC_CLOSE_ITER** parameter for configurable morphological closing of the nucleus mask.

### Changed
- Updated ***processOneImage()*** to support both filename-based and dialog-based metadata handling.
- Improved cell mask generation by adding contrast enhancement and scaled Gaussian blur.
- Improved nucleus mask generation by adding contrast enhancement, scaled Gaussian blur, and morphological closing before hole filling.
- Updated top-1% histogram analysis to use the histogram bin index directly as the pixel intensity value.
- Replaced the old percentile calculation based on values[] with index-based percentile calculation.
- Updated logging to report selected mode, marker list, marker combinations, and timepoints at startup.

### Fixed
- Fixed unreliable intensity values from ***getHistogram()*** for full-range 16-bit histograms.
- Fixed top-1% statistics and percentile calculations by avoiding dependence on the values[] array.
- Improved robustness when image filenames do not match the expected naming scheme.

### Removed
- Removed the old strict filename-only parsing logic from ***processOneImage()***.
- Removed the old ***pct()*** helper function and replaced with new ***pctIndex()*** function.

---

## v0.1.0

### Description

First functional version of the **Mock image-analysis pipeline for FIJI/ImageJ**

This version automatically processes all Mock .tif images in a selected input folder. Image metadata is extracted from the filename using the expected naming scheme:

> [!Important]
>
> timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif

1. The script splits the image into individual channels and assigns them based on the marker information encoded in the filename. Channel 1 corresponds to marker 1, channel 2 to marker 2, and channel 3 to DAPI.

2. For each image, the macro generates a cell mask, a nucleus mask, and a cytosol mask. The cell mask is selected using a prioritized marker logic: HA568 -> HA488 -> NS4B568 -> dsRNA488. The cytosol mask is then used as the measurement region for both marker channels.

3. Within the cytosol ROI, the script measures the top 1% brightest pixels using a full 16-bit histogram with 65,536 bins. Exported values include the top-1% threshold, number of top pixels, total cytosol pixel count, mean, median, standard deviation, and several high-percentile intensity values (5% -> 0.01%)

4. Results are saved as separate CSV files grouped by timepoint, marker combination, and measured channel. In addition, the script can save binary masks and QC overlay images to allow visual inspection of the segmentation quality.

### Main features of this version:

- Batch processing of all Mock images in the selected folder
- Filename-based metadata extraction
- Automatic channel splitting and marker-based channel naming
- Cell, nucleus, and cytosol mask generation
- Cytosol-specific intensity measurement
- Top-1% intensity analysis using full 16-bit histogram binning
- Structured CSV export including run ID and macro version
- Optional saving of binary masks and QC overlays
- Reproducible folder structure for results, masks, and QC images
