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

## v0.2.0



### New function

- 

### Erased functions

- 