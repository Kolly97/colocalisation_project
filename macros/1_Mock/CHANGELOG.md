# Initial Mock Top-1% Pipeline

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

**Biggest changes**

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