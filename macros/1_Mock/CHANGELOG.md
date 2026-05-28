# CHANGELOG for 01_Mock_Background_Pipeline.ijm

---

| Version | Short description                                            |
| ------- | ------------------------------------------------------------ |
| v0.1.0  | Initial Mock Top-1% cytosol intensity pipeline with batch processing, cytosol masks, CSV export, and QC outputs. |
| v0.2.0  | Added flexible metadata handling with filename/dialog mode and fixed 16-bit Top-1% histogram calculation. |
| v0.3.0  | Added manual threshold mode for cell and nucleus masks and run-specific CSV export. |
| v0.4.0  | Added artifact exclusion, cytosol particle-size filtering, and cell-line-specific defaults. |
| v0.4.1  | Strengthened artifact handling, added extra high-percentile outputs, and saved full FIJI log files. |
| v0.5.0  | Added composite RGB JPG QC export and improved binary mask stability with repeated BlackBackground pinning. |

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