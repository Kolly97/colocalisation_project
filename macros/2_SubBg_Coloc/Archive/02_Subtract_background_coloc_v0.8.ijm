// ==========================================================
// Background Subtraction + Coloc Prep Pipeline  --  V0.8
// Author 	: Kolja Hildenbrand
// Date   	: 2026-06-01
// Status	: OLD (superseded by v0.14.3)
//
// Changes vs V0.7:
//	CHANGE:	Defined cell line as Array to make cell type selection more adaptive
//  AUTOMATIC cytosol redesign (fixes the two v0.7 field problems):
//   - A1: "nearest-to-centre" used to pick UNINFECTED / mitotic cells.
//         Now each candidate cell must pass an INFECTION GATE: the 99th
//         percentile of BOTH bg-subtracted marker channels inside the
//         cell must be >= INFECT_MIN_P99 (default 100). p99 (not max) is
//         robust to single hot pixels; measuring on bg-subtracted data
//         means 100 = "real puncta vs ~0 background".
//   - A2: the DAPI-dilation circle did not match the real (elongated,
//         asymmetric) cytoplasm. Replaced by a pure-IJM segmentation:
//           cell mask (autofluorescence, 568 channel)
//           ∩  nucleus-seeded VORONOI split  (cuts touching cells)
//           − nucleus
//         → a per-cell cytosol that follows the real cell shape.
//
//  MANUAL ROI overhaul (M1/M2/M3):
//   - M2: draw on a MERGED RGB composite (all channels visible), not
//         just channel 1.
//   - M3: draw MULTIPLE ROIs per image; coloc then runs on each ROI in
//         turn (one CSV row per ROI).
//   - M1: optional automatic nucleus subtraction from each drawn ROI
//         (one yes/no dialog up front).
//
//  COLOC CSV (G1):
//   - The old IJ-Log parser never worked (the "Colocalization Threshold"
//     plugin writes to its own Results TABLE, not the Log) and spammed
//     "WARN: no coloc labels parsed". Deleted. We now read the Results
//     table via getResult() after the plugin runs (auto AND manual). A
//     one-time dump of Table.headings lets you confirm the column names
//     on the first real run.
//   - CSV gains a roi_index column; nuc_dilate_um dropped; roi_mode kept.
//
// NOTE ON BIAS (publication): choosing WHICH cell to analyse by its
// marker signal is legitimate (an uninfected cell has no biology). The
// ROI BOUNDARY is still defined by morphology (cell mask + nucleus), not
// by the coloc intensity, so the measured coloc value stays unbiased.
//
// Filename schema (binding):
//   timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif
//   C1 = marker1, C2 = marker2, C3 = DAPI (always)
//
// Output (inside INPUT_DIR):
//   <RUN_ID>_bgsub/
//     bgsub_<original>.tif                   multi-ch 16-bit
//     bgsub_<original>.jpg                   8-bit RGB + scale bar (opt.)
//     qc_<original>_roi.jpg                  ROI QC composite (opt.)
//     background_values_used.md              bg values entered
//     coloc_results_<RUN_ID>.csv             one row per ROI
//     macro_log_<RUN_ID>.txt                 full IJ Log
// ==========================================================


// ============== 1. CONFIG ==================================

// Runtime state
var INPUT_DIR;
var OUTPUT_DIR;
var RUN_ID;
var MODE;        // "automatic" or "manual" (per-image metadata source)
var CELL_LINE = newArray("Huh7", "VeroE6");	// Pre-selection for cell lines gets overwritten by dialog

// Domain lists — always editable at startup (pre-filled with these)
var MARKERS       = newArray("HA568", "HA488", "dsRNA488", "NS4B568");
var ANALYSE_COMBI = newArray("HA568_dsRNA488", "NS4B568_dsRNA488", "NS4B568_HA488");
var TIMEPOINTS    = newArray("12h", "24h");

// Background lookup table (filled by askBackgroundValues)
var BG_KEYS   = newArray();
var BG_VALUES = newArray();

// --- ROI strategy for the colocalisation step ---
// "auto_central" : cell mask + Voronoi split + infection gate → one
//                  central, infected cell's cytosol.
// "manual_draw"  : user draws one or more ROIs on a composite.
var ROI_MODE = "auto_central";

// Cell mask (autofluorescence) — same family as the Mock cell mask.
var CELL_THR_METHOD = "Li";       // Li | Otsu | Triangle | Yen
var CELL_THR_FACTOR = 0.5;        // <1 = more permissive (more cytoplasm kept)
var BLUR_SIGMA_CELL = 1;          // microns (scaled Gaussian)
var CELL_MIN_SIZE   = 3000;       // px, min cell footprint in Analyze Particles

// Nucleus mask
var NUC_THR_METHOD  = "Otsu";     // Otsu | Triangle
var NUC_CLOSE_ITER  = 2;
var BLUR_SIGMA_NUC  = 1;          // microns (scaled Gaussian)

// Infection gate: a cell is "infected" only if the 99th percentile of
// BOTH bg-subtracted marker channels inside it is >= this value.
var INFECT_MIN_P99  = 100;

// A cell whose bbox comes within this many px of an image edge is treated
// as partially imaged and de-prioritised.
var BORDER_MARGIN_PX = 0;

// Coloc step. "none" | "manual" | "auto".
var COLOC_MODE = "manual";

// QC / figure exports
var SAVE_COLOC_QC   = true;       // composite + ROI outline(s) per image (JPG)
var SAVE_JPG        = false;      // bg-subtracted composite + scale bar (JPG)
var JPG_SCALEBAR_UM = 20;

// Output naming
var OUT_PREFIX = "bgsub_";

// Reproducibility
var MACRO_VERSION = "0.8.0";

// One-time flag: dump the Results-table headings on the first coloc read.
var HEADINGS_DUMPED = false;


// ============== MAIN =======================================

setOption("BlackBackground", true);

chooseInputDir();
RUN_ID = makeRunId();
askModeAndConfig();    // also calls applyCellLineDefaults()
askBackgroundValues();
buildOutputDir();
if (COLOC_MODE != "none") initColocCsv();
imageFiles = listMoiFiles();

logHeader(imageFiles.length);

for (f = 0; f < imageFiles.length; f++) {
    print("[" + (f+1) + "/" + imageFiles.length + "] " + imageFiles[f]);
    processOneImage(imageFiles[f]);
    cleanupBetweenImages();
}

print("=== DONE ===");
saveLogToFile();


// ============== 2. SETUP FUNCTIONS =========================

function chooseInputDir() {
    INPUT_DIR = getDirectory("Choose folder with MOI .tif images");
    if (INPUT_DIR == "") exit("No folder selected.");
}

function makeRunId() {
    getDateAndTime(y, m, dw, d, h, mn, s, ms);
    return "" + y + IJ.pad(m+1, 2) + IJ.pad(d, 2)
        + "_" + IJ.pad(h, 2) + IJ.pad(mn, 2);
}

function askModeAndConfig() {
    Dialog.create("Pipeline mode");

    Dialog.addMessage("How should the macro know which marker is on which channel?");
    Dialog.addRadioButtonGroup("Marker source:",
        newArray("automatic", "manual"), 2, 1, "automatic");

    Dialog.addMessage("---");
    Dialog.addMessage("What to do after subtracting the background?");
    Dialog.addRadioButtonGroup("Pipeline:",
        newArray("subtract only",
                 "subtract + manual coloc",
                 "subtract + auto coloc"),
        3, 1, "subtract + manual coloc");

    Dialog.addMessage("---");
    Dialog.addMessage("How should the colocalisation ROI be built?");
    Dialog.addRadioButtonGroup("ROI strategy:",
        newArray("automatic central cell", "manual draw"), 2, 1, "automatic central cell");

    Dialog.addMessage("---");
    Dialog.addMessage("Which cell line are you analysing?");
	Dialog.addRadioButtonGroup("Cell line:", CELL_LINE, CELL_LINE.length, 1, CELL_LINE[0]);
    Dialog.addMessage("---");
    Dialog.addCheckbox("Also save 8-bit RGB JPG with " + JPG_SCALEBAR_UM + " um scale bar", SAVE_JPG);

    Dialog.show();

    MODE          = Dialog.getRadioButton();
    pipelineScope = Dialog.getRadioButton();
    roiChoice     = Dialog.getRadioButton();
    CELL_LINE     = Dialog.getRadioButton();
    SAVE_JPG      = Dialog.getCheckbox();

    if (pipelineScope == "subtract only")              COLOC_MODE = "none";
    else if (pipelineScope == "subtract + auto coloc") COLOC_MODE = "auto";
    else                                               COLOC_MODE = "manual";

    if (roiChoice == "manual draw") ROI_MODE = "manual_draw";
    else                            ROI_MODE = "auto_central";

    applyCellLineDefaults(CELL_LINE);

    // ALWAYS show the domain-list dialog (Feature A).
    Dialog.create("Pipeline setup");
    Dialog.addMessage("Edit the lists used by the macro (CSV files, dialogs, validation).\n"
                    + "Pre-filled with the workflow standards — just click OK to keep them.\n"
                    + "Comma-separated, whitespace trimmed.");
    Dialog.addString("Markers:",    arrToStr(MARKERS),       40);
    Dialog.addString("Combos:",     arrToStr(ANALYSE_COMBI), 40);
    Dialog.addString("Timepoints:", arrToStr(TIMEPOINTS),    20);
    Dialog.show();
    MARKERS       = parseCsvString(Dialog.getString());
    ANALYSE_COMBI = parseCsvString(Dialog.getString());
    TIMEPOINTS    = parseCsvString(Dialog.getString());
}

// Cell-line tuning (strictness knobs only).
function applyCellLineDefaults(cellLine) {
    if (cellLine == "VeroE6") {
        BLUR_SIGMA_CELL = 1;
        CELL_THR_FACTOR = 0.5;
        CELL_MIN_SIZE   = 1000;   // Vero ~4x smaller
        BLUR_SIGMA_NUC  = 1;
    } else {  // Huh7
        BLUR_SIGMA_CELL = 1;
        CELL_THR_FACTOR = 0.5;
        CELL_MIN_SIZE   = 3000;
        BLUR_SIGMA_NUC  = 1;
    }
    print("Cell-line tuning for " + cellLine + ":");
    print("  BLUR_SIGMA_CELL = " + BLUR_SIGMA_CELL);
    print("  CELL_THR_FACTOR = " + CELL_THR_FACTOR);
    print("  CELL_MIN_SIZE   = " + CELL_MIN_SIZE);
    print("  INFECT_MIN_P99  = " + INFECT_MIN_P99);
}

function askBackgroundValues() {
    Dialog.create("Background values");
    Dialog.addMessage("Enter the background value for each (marker x combo x timepoint).\n"
                    + "These will be subtracted from the matching channel.", 13);

    for (t = 0; t < TIMEPOINTS.length; t++) {
        tp = TIMEPOINTS[t];
        Dialog.addMessage("--- " + tp + " samples ---", 12);
        for (c = 0; c < ANALYSE_COMBI.length; c++) {
            combo = ANALYSE_COMBI[c];
            parts = split(combo, "_");
            m1 = parts[0]; m2 = parts[1];
            Dialog.addMessage(combo);
            Dialog.addNumber(m1 + ":", 0);
            Dialog.addNumber(m2 + ":", 0);
        }
    }
    Dialog.show();

    for (t = 0; t < TIMEPOINTS.length; t++) {
        tp = TIMEPOINTS[t];
        for (c = 0; c < ANALYSE_COMBI.length; c++) {
            combo = ANALYSE_COMBI[c];
            parts = split(combo, "_");
            m1 = parts[0]; m2 = parts[1];
            setBgValue(m1, combo, tp, Dialog.getNumber());
            setBgValue(m2, combo, tp, Dialog.getNumber());
        }
    }
}

function buildOutputDir() {
    OUTPUT_DIR = INPUT_DIR + RUN_ID + "_bgsub" + File.separator;
    ensureDir(OUTPUT_DIR);
    writeBgMarkdown(OUTPUT_DIR + "background_values_used.md");
}

function initColocCsv() {
    csvPath = OUTPUT_DIR + "coloc_results_" + RUN_ID + ".csv";
    if (File.exists(csvPath)) return;
    header = "image,cell_line,timepoint,combo,channel_1,channel_2,roi_index,"
           + "Rtotal,m,b,Ch1_thresh,Ch2_thresh,"
           + "Rcoloc,R_below_thresh,M1,M2,tM1,tM2,"
           + "Ncoloc,perc_volume,perc_ch1_vol,perc_ch2_vol,"
           + "roi_mode,macro_version,run_id,status\n";
    File.append(header, csvPath);
}

function listMoiFiles() {
    files = getFileList(INPUT_DIR);
    out = newArray();
    for (i = 0; i < files.length; i++) {
        n = files[i]; ln = toLowerCase(n);
        if (endsWith(ln, ".tif") && indexOf(ln, "moi") >= 0)
            out = Array.concat(out, n);
    }
    out = Array.sort(out);
    return out;
}

function logHeader(nFiles) {
    print("=== Bg-sub + Coloc-prep V" + MACRO_VERSION + " ===");
    print("run_id            : " + RUN_ID);
    print("mode              : " + MODE);
    colocLabel = "disabled";
    if (COLOC_MODE == "manual") colocLabel = "MANUAL (pause per ROI)";
    if (COLOC_MODE == "auto")   colocLabel = "AUTO (plugin runs per ROI, Results table read)";
    print("coloc step        : " + colocLabel);
    print("ROI strategy      : " + ROI_MODE);
    if (ROI_MODE == "auto_central")
        print("infection gate    : p99 >= " + INFECT_MIN_P99 + " in BOTH channels (bg-subtracted)");
    qcLabel = "no";
    if (SAVE_COLOC_QC) qcLabel = "yes (composite + ROI outline)";
    print("coloc QC jpg      : " + qcLabel);
    jpgLabel = "no";
    if (SAVE_JPG) jpgLabel = "yes (" + JPG_SCALEBAR_UM + " um scale bar)";
    print("save jpg          : " + jpgLabel);
    print("input dir         : " + INPUT_DIR);
    print("output dir        : " + OUTPUT_DIR);
    print("MOI files         : " + nFiles);
    print("---");
}


// ============== 3. PER-IMAGE PIPELINE ======================
// ROI-first order (unbiased boundary, raw signal visible for drawing):
//   split → build ROI(s) → coloc QC → subtract bg → coloc per ROI
//   → merge → save.
// ===========================================================

function processOneImage(fname) {
    setOption("BlackBackground", true);

    open(INPUT_DIR + fname);
    title = getTitle();
    imgName = substring(title, 0, lastIndexOf(title, "."));

    cellLine = "";
    meta = resolveMetadata(imgName);
    if (meta.length == 0) { print("  -> skipped"); return; }
    tp = meta[0]; m1 = meta[1]; m2 = meta[2];
    tokens = split(imgName, "_");
    if (tokens.length >= 2) cellLine = tokens[1]; else cellLine = CELL_LINE;

    comboKey = m1 + "_" + m2;
    if (!isValidCombo(comboKey)) { print("  SKIP: combo not in ANALYSE_COMBI -> " + comboKey); return; }
    if (!inArray(tp, TIMEPOINTS)) { print("  SKIP: timepoint not in TIMEPOINTS -> " + tp); return; }
    print("  combo=" + comboKey + "  tp=" + tp + "  cellLine=" + cellLine);

    bg1 = getBgValue(m1, comboKey, tp);
    bg2 = getBgValue(m2, comboKey, tp);
    print("  bg(" + m1 + ") = " + bg1 + ", bg(" + m2 + ") = " + bg2);

    splitAndRenameChannels(title, m1, m2);

    // ---- build ROI(s) on RAW channels (BEFORE subtraction) ----
    nRois = 0;
    if (COLOC_MODE != "none") {
        nRois = buildColocRoi(m1, m2, bg1, bg2);
        if (SAVE_COLOC_QC) saveColocQc(imgName, m1, m2);
    }

    // ---- NOW subtract background -------------------------
    subtractFromChannel(m1 + "_channel", bg1);
    subtractFromChannel(m2 + "_channel", bg2);

    // ---- colocalisation step, per ROI --------------------
    if (COLOC_MODE != "none") {
        if (nRois <= 0) {
            doColocalisationStep(imgName, cellLine, tp, comboKey, m1, m2, false, -1);
        } else {
            for (r = 0; r < nRois; r++) {
                doColocalisationStep(imgName, cellLine, tp, comboKey, m1, m2, true, r);
            }
        }
    }

    // ---- merge back + save -------------------------------
    mergeChannelsBack(m1, m2);
    Stack.setDisplayMode("color");
    Stack.setChannel(1);

    outPath = OUTPUT_DIR + OUT_PREFIX + imgName + ".tif";
    saveAs("Tiff", outPath);
    print("  saved tif -> " + outPath);

    if (SAVE_JPG) {
        jpgPath = OUTPUT_DIR + OUT_PREFIX + imgName + ".jpg";
        saveJpgWithScaleBar(jpgPath);
    }

    setOption("BlackBackground", true);
}

function resolveMetadata(imgName) {
    parsed = tryParseFilename(imgName);
    result = newArray();
    if (MODE == "automatic") {
        if (parsed.length > 0) {
            result = parsed;
        } else {
            print("  Parse failed for: " + imgName);
            action = askParseFailureAction(imgName);
            if (action == "skip") { result = newArray(); }
            else                  { result = askImageMetadata("(parse failed)", "", "", ""); }
        }
    } else {
        defTp = ""; defM1 = ""; defM2 = "";
        if (parsed.length > 0) { defTp = parsed[0]; defM1 = parsed[1]; defM2 = parsed[2]; }
        result = askImageMetadata(imgName, defTp, defM1, defM2);
    }
    return result;
}

function tryParseFilename(imgName) {
    tokens = split(imgName, "_");
    if (tokens.length < 7) return newArray();
    tp_ = tokens[0]; m1_ = tokens[3]; m2_ = tokens[4];
    if (!inArray(m1_, MARKERS) || !inArray(m2_, MARKERS)) return newArray();
    if (!inArray(tp_, TIMEPOINTS))                        return newArray();
    return newArray(tp_, m1_, m2_);
}

function askImageMetadata(imgLabel, defTp, defM1, defM2) {
    if (defTp == "" || !inArray(defTp, TIMEPOINTS)) defTp = TIMEPOINTS[0];
    if (defM1 == "" || !inArray(defM1, MARKERS))    defM1 = MARKERS[0];
    if (defM2 == "" || !inArray(defM2, MARKERS))    defM2 = MARKERS[1 % MARKERS.length];
    Dialog.create("Image metadata");
    Dialog.addMessage("Image: " + imgLabel + "\nC3 is always DAPI.");
    Dialog.addRadioButtonGroup("Timepoint:",      TIMEPOINTS, 1, TIMEPOINTS.length, defTp);
    Dialog.addRadioButtonGroup("Channel 1 (C1):", MARKERS,    1, MARKERS.length,    defM1);
    Dialog.addRadioButtonGroup("Channel 2 (C2):", MARKERS,    1, MARKERS.length,    defM2);
    Dialog.show();
    tp_ = Dialog.getRadioButton();
    m1_ = Dialog.getRadioButton();
    m2_ = Dialog.getRadioButton();
    if (m1_ == m2_) { showMessage("Error", "C1 and C2 must differ."); exit("Aborted."); }
    return newArray(tp_, m1_, m2_);
}

function askParseFailureAction(imgName) {
    Dialog.create("Filename parse failed");
    Dialog.addMessage("Cannot parse: " + imgName + "\n\nWhat do you want to do?");
    Dialog.addRadioButtonGroup("Action:", newArray("skip", "manual"), 2, 1, "skip");
    Dialog.show();
    return Dialog.getRadioButton();
}


// ============== 4. CHANNEL HANDLING ========================

function isValidCombo(comboKey) {
    result = inArray(comboKey, ANALYSE_COMBI);
    return result;
}

function splitAndRenameChannels(title, m1, m2) {
    selectWindow(title);
    run("Split Channels");
    selectWindow("C1-" + title); rename(m1 + "_channel");
    selectWindow("C2-" + title); rename(m2 + "_channel");
    selectWindow("C3-" + title); rename("DAPI_channel");
}

function subtractFromChannel(chTitle, bgValue) {
    if (!isOpen(chTitle)) { print("  WARN: missing " + chTitle); return; }
    selectWindow(chTitle);
    run("Select None");
    run("Subtract...", "value=" + bgValue);
}

function mergeChannelsBack(m1, m2) {
    run("Merge Channels...",
        "c1=" + m1 + "_channel "
      + "c2=" + m2 + "_channel "
      + "c3=DAPI_channel create");
}

// Priority for the autofluorescence cell mask: 568 markers show diffuse
// cytoplasm (good outline); dsRNA is pure puncta (bad outline) → last.
function pickCellMaskSource(m1, m2) {
    priority = newArray("HA568", "NS4B568", "HA488", "dsRNA488");
    for (i = 0; i < priority.length; i++) {
        cand = priority[i] + "_channel";
        if (isOpen(cand)) return cand;
    }
    return "";
}


// ============== 5. MASKS ===================================

// Cell mask from autofluorescence (ported from Mock makeCellMask).
function makeCellMask(srcTitle) {
    selectWindow(srcTitle);
    run("Duplicate...", "title=Cell_Mask");
    run("Enhance Contrast", "saturated=0.35");
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_CELL + " scaled");
    setAutoThreshold(CELL_THR_METHOD + " dark");
    getThreshold(lo, hi);
    loNew = lo * CELL_THR_FACTOR;
    if (loNew < 0) loNew = 0;
    setThreshold(loNew, hi);
    print("  cell thr (" + CELL_THR_METHOD + "): auto=" + lo + " scaled -> " + loNew);
    run("Convert to Mask");
    run("Fill Holes");
}

function makeNucleusMask(dapiTitle) {
    selectWindow(dapiTitle);
    run("Duplicate...", "title=Nucleus_Mask");
    run("Enhance Contrast", "saturated=0.35");
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_NUC + " scaled");
    setAutoThreshold(NUC_THR_METHOD + " dark");
    run("Convert to Mask");
    if (NUC_CLOSE_ITER > 0) {
        run("Options...", "iterations=" + NUC_CLOSE_ITER + " count=1 black do=Nothing");
        run("Dilate");
        run("Fill Holes");
        run("Erode");
    } else {
        run("Fill Holes");
    }
}

// Split touching cells: Voronoi of the nuclei gives dividing ridges
// (equidistant between nuclei). Subtract those ridges from the cell mask
// → each connected component belongs to one nucleus. Needs Cell_Mask +
// Nucleus_Mask open; produces "Cells_Split".
function splitTouchingCellsVoronoi() {
    selectWindow("Nucleus_Mask");
    run("Duplicate...", "title=Voronoi_tmp");
    run("Voronoi");                      // ridges get grey value, bg = 0
    setThreshold(1, 65535);              // covers 8- and 16-bit Voronoi output
    run("Convert to Mask");              // binary ridge lines
    run("Options...", "iterations=1 count=1 black do=Nothing");
    run("Dilate");                       // thicken ridges so cuts fully separate
    imageCalculator("Subtract create", "Cell_Mask", "Voronoi_tmp");
    rename("Cells_Split");
    setThreshold(1, 255);
    run("Convert to Mask");
    waitForUser("can I close everything?");
    closeIfOpen("Voronoi_tmp");
}

// 16-bit percentile inside the current ROI of chTitle (reuses the Mock
// histogram-percentile pattern: bin index == intensity for 16-bit).
function roiPercentile(chTitle, roiIdx, p) {
    selectWindow(chTitle);
    roiManager("Select", roiIdx);
    getHistogram(values, counts, 65536, 0, 65535);
    nTotal = 0;
    for (k = 0; k < 65536; k++) nTotal += counts[k];
    if (nTotal == 0) return 0;
    r = pctIndex(counts, 65536, nTotal, p);
    return r;
}

function pctIndex(counts, nBins, nTotal, p) {
    target = nTotal * p / 100.0;
    a = 0;
    for (i = 0; i < nBins; i++) {
        a += counts[i];
        if (a >= target) return i;
    }
    return nBins - 1;
}

// Fill ROI roiIdx onto a fresh 8-bit mask (dimensions from refTitle).
function roiIndexToMask(refTitle, roiIdx, outName) {
    selectWindow(refTitle);
    getDimensions(W, H, ch, sl, fr);
    newImage(outName, "8-bit black", W, H, 1);
    roiManager("Select", roiIdx);
    setColor(255);
    fill();
    run("Select None");
}

// Bg-subtracted duplicate of a channel (for infection gating only — the
// real channels are subtracted later).
function makeBgSubDuplicate(srcTitle, outName, bg) {
    selectWindow(srcTitle);
    run("Select None");
    run("Duplicate...", "title=" + outName);
    run("Subtract...", "value=" + bg);
}


// ============== 6. ROI CONSTRUCTION ========================
// Every path leaves the coloc ROI(s) in the ROI Manager (index 0..n-1)
// and returns the count, so the per-ROI coloc loop is strategy-agnostic.
// ===========================================================

function buildColocRoi(m1, m2, bg1, bg2) {
    if (ROI_MODE == "manual_draw") {
        c = buildManualRois(m1, m2);
        return c;
    }
    c = buildAutoCytosolRoi(m1, m2, bg1, bg2);
    return c;
}

// ------------------------------
// AUTOMATIC: cell mask (autofluorescence) ∩ Voronoi split, choose the
// central INFECTED cell, cytosol = cell − nucleus. Returns 1 or 0.
// ------------------------------
function buildAutoCytosolRoi(m1, m2, bg1, bg2) {
    cellSrc = pickCellMaskSource(m1, m2);
    if (cellSrc == "") { print("  WARN: no cell-mask source — full image"); return 0; }
    print("  cell-mask source: " + cellSrc);

    makeCellMask(cellSrc);                 // Cell_Mask
    makeNucleusMask("DAPI_channel");       // Nucleus_Mask
    splitTouchingCellsVoronoi();           // Cells_Split

    // bg-subtracted duplicates for the infection gate
    makeBgSubDuplicate(m1 + "_channel", "m1_sub", bg1);
    makeBgSubDuplicate(m2 + "_channel", "m2_sub", bg2);

    // one ROI per cell
    roiManager("reset");
    selectWindow("Cells_Split");
    run("Analyze Particles...", "size=" + CELL_MIN_SIZE + "-Infinity pixel add");
    nCells = roiManager("count");
    if (nCells == 0) {
        print("  WARN: no cells found — full image");
        closeAutoMaskWindows();
        return 0;
    }

    selectWindow("Cells_Split");
    imgW = getWidth(); imgH = getHeight();
    cx = imgW / 2.0;   cy = imgH / 2.0;

    bestInfComplete = -1; dInfComplete = 1e18;
    bestInf         = -1; dInf         = 1e18;
    bestComplete    = -1; dComplete    = 1e18;
    bestAny         = -1; dAny         = 1e18;

    for (i = 0; i < nCells; i = i + 1) {
        selectWindow("Cells_Split");
        roiManager("Select", i);
        getSelectionBounds(x, y, w, h);
        ccx = x + w/2.0; ccy = y + h/2.0;
        dx = ccx - cx; dy = ccy - cy;
        dist = sqrt(dx*dx + dy*dy);
        touches = (x <= BORDER_MARGIN_PX) || (y <= BORDER_MARGIN_PX)
               || (x + w >= imgW - BORDER_MARGIN_PX) || (y + h >= imgH - BORDER_MARGIN_PX);

        p1 = roiPercentile("m1_sub", i, 99);
        p2 = roiPercentile("m2_sub", i, 99);
        infected = (p1 >= INFECT_MIN_P99) && (p2 >= INFECT_MIN_P99);

        if (dist < dAny) { dAny = dist; bestAny = i; }
        if (!touches && dist < dComplete) { dComplete = dist; bestComplete = i; }
        if (infected) {
            if (dist < dInf) { dInf = dist; bestInf = i; }
            if (!touches && dist < dInfComplete) { dInfComplete = dist; bestInfComplete = i; }
        }
    }

    chosen = bestInfComplete;
    why = "infected + fully imaged";
    if (chosen < 0) { chosen = bestInf;      why = "infected (touches border)"; }
    if (chosen < 0) { chosen = bestComplete; why = "WARN: no infected cell — uninfected, fully imaged"; }
    if (chosen < 0) { chosen = bestAny;      why = "WARN: fallback nearest-centre"; }
    if (chosen < 0) { closeAutoMaskWindows(); return 0; }
    print("  chosen cell #" + chosen + " (" + why + ")");

    // cytosol = chosen cell − nucleus
    roiIndexToMask("Cells_Split", chosen, "Chosen_Cell");
    imageCalculator("Subtract create", "Chosen_Cell", "Nucleus_Mask");
    rename("Cytosol_Mask");
    setThreshold(1, 255);
    run("Convert to Mask");
    getStatistics(area, meanCyto);
    if (meanCyto == 0) {
        print("  WARN: cytosol empty after nucleus removal — full image");
        closeAutoMaskWindows();
        closeIfOpen("Chosen_Cell");
        closeIfOpen("Cytosol_Mask");
        return 0;
    }
    run("Create Selection");
    roiManager("reset");
    roiManager("Add");
    roiManager("Select", 0);
    roiManager("Rename", "Cytosol_auto");

    closeAutoMaskWindows();
    closeIfOpen("Chosen_Cell");
    closeIfOpen("Cytosol_Mask");
    return 1;
}

function closeAutoMaskWindows() {
    closeIfOpen("Cell_Mask");
    closeIfOpen("Nucleus_Mask");
    closeIfOpen("Cells_Split");
    closeIfOpen("m1_sub");
    closeIfOpen("m2_sub");
}

// ------------------------------
// MANUAL: draw one or more ROIs on an RGB composite of all channels.
// Optional automatic nucleus subtraction. Returns the ROI count.
// ------------------------------
function buildManualRois(m1, m2) {
    Dialog.create("Manual ROI options");
    Dialog.addMessage("Draw on a merged composite of all channels.\n"
                    + "You can add several ROIs; coloc runs on each in turn.");
    Dialog.addCheckbox("Exclude nucleus from drawn ROIs (auto, via DAPI mask)", true);
    Dialog.show();
    excludeNuc = Dialog.getCheckbox();

    if (excludeNuc) makeNucleusMask("DAPI_channel");   // build once

    makeDrawComposite(m1, m2);                          // "Draw_RGB"

    roiManager("reset");
    setTool("freehand");
    more = true;
    count = 0;
    while (more) {
        selectWindow("Draw_RGB");
        run("Select None");
        waitForUser("Draw ROI #" + (count+1),
            "Draw a FREEHAND ROI around a cell on the composite, then click OK.\n"
          + "(Cancel aborts the WHOLE run.)");
        if (selectionType() != -1) {
            roiManager("Add");
            count = count + 1;
        } else {
            print("  (empty selection ignored)");
        }
        Dialog.create("More ROIs?");
        Dialog.addRadioButtonGroup("Action:", newArray("draw another", "done"), 2, 1, "done");
        Dialog.show();
        if (Dialog.getRadioButton() == "done") more = false;
    }
    closeIfOpen("Draw_RGB");

    if (count == 0) { print("  WARN: no ROI drawn — full image"); closeIfOpen("Nucleus_Mask"); return 0; }

    if (excludeNuc) {
        excludeNucleusFromAllRois(m1, count);
        closeIfOpen("Nucleus_Mask");
    }

    n = roiManager("count");
    for (i = 0; i < n; i = i + 1) { roiManager("Select", i); roiManager("Rename", "ROI_" + (i+1)); }
    print("  manual ROIs: " + n);
    return n;
}

// Build an RGB composite of the 3 RAW channels for drawing (sources
// duplicated so the split channels survive).
function makeDrawComposite(m1, m2) {
    selectWindow(m1 + "_channel");  run("Select None"); run("Duplicate...", "title=draw_m1");
    selectWindow(m2 + "_channel");  run("Select None"); run("Duplicate...", "title=draw_m2");
    selectWindow("DAPI_channel");   run("Select None"); run("Duplicate...", "title=draw_dapi");
    run("Merge Channels...", "c1=draw_m1 c2=draw_m2 c3=draw_dapi create");
    Stack.setDisplayMode("composite");
    Stack.setChannel(1); run("Enhance Contrast", "saturated=0.35");
    Stack.setChannel(2); run("Enhance Contrast", "saturated=0.35");
    Stack.setChannel(3); run("Enhance Contrast", "saturated=0.35");
    run("Stack to RGB");
    rename("Draw_RGB");
    closeIfOpen("Composite");
    setOption("BlackBackground", true);
}

// Replace each of the first n ROIs with (ROI − nucleus). Appends the
// cytosol ROIs to the end, then deletes the originals.
function excludeNucleusFromAllRois(m1, n) {
    for (i = 0; i < n; i = i + 1) {
        roiIndexToMask("Nucleus_Mask", i, "Roi_Tmp");
        imageCalculator("Subtract", "Roi_Tmp", "Nucleus_Mask");
        selectWindow("Roi_Tmp");
        setThreshold(1, 255);
        run("Convert to Mask");
        getStatistics(area, mean);
        if (mean > 0) {
            run("Create Selection");
            roiManager("Add");
        } else {
            print("  WARN: ROI " + (i+1) + " empty after nucleus removal — kept original");
            roiManager("Select", i);
            roiManager("Add");
        }
        closeIfOpen("Roi_Tmp");
    }
    del = newArray(n);
    for (i = 0; i < n; i = i + 1) del[i] = i;
    roiManager("Select", del);
    roiManager("Delete");
}

// QC: composite of the 3 RAW channels + ALL ROIs as outlines.
function saveColocQc(imgName, m1, m2) {
    if (!isOpen(m1 + "_channel") || !isOpen(m2 + "_channel") || !isOpen("DAPI_channel")) {
        print("  WARN: coloc QC skipped (missing channel window)");
        return;
    }
    selectWindow(m1 + "_channel"); run("Select None"); run("Duplicate...", "title=qc_m1");
    selectWindow(m2 + "_channel"); run("Select None"); run("Duplicate...", "title=qc_m2");
    selectWindow("DAPI_channel");  run("Select None"); run("Duplicate...", "title=qc_dapi");

    run("Merge Channels...", "c1=qc_m1 c2=qc_m2 c3=qc_dapi create");
    Stack.setDisplayMode("composite");
    Stack.setChannel(1); run("Enhance Contrast", "saturated=0.35");
    Stack.setChannel(2); run("Enhance Contrast", "saturated=0.35");
    Stack.setChannel(3); run("Enhance Contrast", "saturated=0.35");

    run("Stack to RGB");
    nR = roiManager("count");
    for (i = 0; i < nR; i = i + 1) { roiManager("Select", i); run("Add Selection..."); }
    run("Flatten");

    run("Options...", "jpeg=90");
    qcPath = OUTPUT_DIR + "qc_" + imgName + "_roi.jpg";
    saveAs("Jpeg", qcPath);
    print("  saved coloc QC -> " + qcPath);

    close();
    closeIfOpen("Composite (RGB)");
    closeIfOpen("Composite");
    setOption("BlackBackground", true);
}


// ============== 7. COLOCALISATION ==========================

function doColocalisationStep(imgName, cellLine, tp, comboKey, m1, m2, haveRoi, roiIdx) {
    if (COLOC_MODE == "auto") {
        doColocAuto(imgName, cellLine, tp, comboKey, m1, m2, haveRoi, roiIdx);
    } else {
        doColocManual(imgName, cellLine, tp, comboKey, m1, m2, haveRoi, roiIdx);
    }
}

// Apply ROI roiIdx to m1_channel as the active selection (the plugin's
// "Use ROI: Channel 1" reads that selection). roiIdx<0 = no ROI.
function applyRoiToChannel(m1, haveRoi, roiIdx) {
    if (!isOpen(m1 + "_channel")) return;
    selectWindow(m1 + "_channel");
    run("Select None");
    if (haveRoi && roiIdx >= 0 && roiManager("count") > roiIdx) roiManager("Select", roiIdx);
}

function doColocManual(imgName, cellLine, tp, comboKey, m1, m2, haveRoi, roiIdx) {
    applyRoiToChannel(m1, haveRoi, roiIdx);
    useRoi = "<None>";
    if (haveRoi) useRoi = "Channel 1";
    label = "full image";
    if (haveRoi) label = "ROI #" + (roiIdx+1);

    nBefore = nResults;
    waitForUser("Coloc " + label + "  [" + imgName + "]",
        label + " is applied to " + m1 + "_channel (bg subtracted).\n\n"
      + "Run: Analyze > Colocalisation > Colocalisation Threshold\n"
      + "  Channel 1 : " + m1 + "_channel\n"
      + "  Channel 2 : " + m2 + "_channel\n"
      + "  Use ROI   : " + useRoi + "\n"
      + "  Channel   : Red : Green     Include zero-zero: [x]\n\n"
      + "Click OK when done (or if you skipped this ROI).\n(Cancel aborts the WHOLE run.)");
    appendColocRowFromTable(imgName, cellLine, tp, comboKey, m1, m2, roiIdx, nBefore, "manual");
}

function doColocAuto(imgName, cellLine, tp, comboKey, m1, m2, haveRoi, roiIdx) {
    applyRoiToChannel(m1, haveRoi, roiIdx);

    args = "channel_1=" + m1 + "_channel "
         + "channel_2=" + m2 + "_channel "
         + "channel=[Red : Green] "
         + "include";
    if (haveRoi) args = args + " use=[Channel 1]";

    print("  auto coloc (roi " + (roiIdx+1) + "): run(\"Colocalization Threshold\", \"" + args + "\")");
    nBefore = nResults;
    run("Colocalization Threshold", args);
    appendColocRowFromTable(imgName, cellLine, tp, comboKey, m1, m2, roiIdx, nBefore, "auto_ok");
}

// Read the new row of the standard Results table written by the plugin.
// Header strings come from the plugin's own table; the first call dumps
// Table.headings to the log so the exact names can be confirmed.
function appendColocRowFromTable(imgName, cellLine, tp, comboKey, m1, m2, roiIdx, nBefore, status) {
    rtotal = ""; slope = ""; intercept = ""; ch1 = ""; ch2 = "";
    rcoloc = ""; rbelow = ""; mm1 = ""; mm2 = ""; tm1 = ""; tm2 = "";
    ncoloc = ""; pctvol = ""; pctch1 = ""; pctch2 = "";
    st = status;

    nAfter = nResults;
    if (nAfter > nBefore) {
        row = nAfter - 1;
        if (!HEADINGS_DUMPED) {
            print("  Results table headings: " + Table.headings);
            HEADINGS_DUMPED = true;
        }
        rtotal    = safeResult("Rtotal", row);
        slope     = safeResult("m", row);
        intercept = safeResult("b", row);
        ch1       = safeResult("Ch1 thresh", row);
        ch2       = safeResult("Ch2 thresh", row);
        rcoloc    = safeResult("Rcoloc", row);
        rbelow    = safeResult("R-threshold", row);
        mm1       = safeResult("M1", row);
        mm2       = safeResult("M2", row);
        tm1       = safeResult("tM1", row);
        tm2       = safeResult("tM2", row);
        ncoloc    = safeResult("Ncoloc", row);
        pctvol    = safeResult("%Volume", row);
        pctch1    = safeResult("%Ch1 Vol", row);
        pctch2    = safeResult("%Ch2 Vol", row);
    } else {
        st = status + "_noresult";
        print("  WARN: Results table did not grow for ROI " + (roiIdx+1)
            + " (skipped, or the plugin uses a custom table window).");
    }

    csvPath = OUTPUT_DIR + "coloc_results_" + RUN_ID + ".csv";
    line = imgName + "," + cellLine + "," + tp + "," + comboKey + ","
         + m1 + "," + m2 + "," + (roiIdx+1) + ","
         + rtotal + "," + slope + "," + intercept + "," + ch1 + "," + ch2 + ","
         + rcoloc + "," + rbelow + "," + mm1 + "," + mm2 + "," + tm1 + "," + tm2 + ","
         + ncoloc + "," + pctvol + "," + pctch1 + "," + pctch2 + ","
         + ROI_MODE + "," + MACRO_VERSION + "," + RUN_ID + "," + st;
    File.append(line, csvPath);
    print("  coloc row (roi " + (roiIdx+1) + ", " + st + "): Rtotal=" + rtotal
        + " M1=" + mm1 + " M2=" + mm2);
}

function safeResult(col, row) {
    v = getResult(col, row);
    if (isNaN(v)) return "";
    r = "" + v;
    return r;
}


// ============== 8. BACKGROUND LOOKUP TABLE =================

function setBgValue(marker, combo, tp, value) {
    key = marker + "_in_" + combo + "_" + tp;
    BG_KEYS   = Array.concat(BG_KEYS,   key);
    BG_VALUES = Array.concat(BG_VALUES, parseFloat(value));
}

function getBgValue(marker, combo, tp) {
    key = marker + "_in_" + combo + "_" + tp;
    for (i = 0; i < BG_KEYS.length; i++) {
        if (BG_KEYS[i] == key) {
            v = parseFloat(BG_VALUES[i]);
            return v;
        }
    }
    exit("Background value missing for key: " + key
       + "\nDid you forget to enter it in the startup dialog?");
}


// ============== 9. OUTPUT / DOCUMENTATION ==================

function writeBgMarkdown(mdPath) {
    if (File.exists(mdPath)) File.delete(mdPath);
    File.append("# Background values used\n\n", mdPath);
    File.append("Generated by `Subtract_background_coloc_v" + MACRO_VERSION + ".ijm`\n\n", mdPath);
    File.append("- Run ID        : " + RUN_ID + "\n", mdPath);
    File.append("- ROI strategy  : " + ROI_MODE + "\n", mdPath);
    if (ROI_MODE == "auto_central")
        File.append("- Infection gate: p99 >= " + INFECT_MIN_P99 + " in both channels (bg-subtracted)\n", mdPath);
    File.append("- Input folder  : " + INPUT_DIR + "\n", mdPath);
    File.append("- Output folder : " + OUTPUT_DIR + "\n\n", mdPath);

    for (t = 0; t < TIMEPOINTS.length; t++) {
        tp = TIMEPOINTS[t];
        File.append("## " + tp + " samples\n\n", mdPath);
        for (c = 0; c < ANALYSE_COMBI.length; c++) {
            combo = ANALYSE_COMBI[c];
            parts = split(combo, "_");
            m1 = parts[0]; m2 = parts[1];
            File.append("### " + combo + "\n", mdPath);
            File.append("- " + m1 + ": " + getBgValue(m1, combo, tp) + "\n", mdPath);
            File.append("- " + m2 + ": " + getBgValue(m2, combo, tp) + "\n\n", mdPath);
        }
    }
    print("Wrote background documentation: " + mdPath);
}

function saveJpgWithScaleBar(jpgPath) {
    Stack.setDisplayMode("composite");
    run("Stack to RGB");
    run("Scale Bar...",
        "width=" + JPG_SCALEBAR_UM + " height=10 thickness=20 font=100 "
      + "color=White background=None location=[Lower Right] overlay");
    run("Flatten");
    run("Options...", "jpeg=90");
    saveAs("Jpeg", jpgPath);
    print("  saved jpg -> " + jpgPath);
}

function saveLogToFile() {
    if (!isOpen("Log")) return;
    logPath = OUTPUT_DIR + "macro_log_" + RUN_ID + ".txt";
    selectWindow("Log");
    saveAs("Text", logPath);
    print("Saved log: " + logPath);
}


// ============== 10. CLEANUP & UTILS ========================

function cleanupBetweenImages() {
    while (nImages > 0) { selectImage(nImages); close(); }
    if (isOpen("ROI Manager")) roiManager("reset");
    if (isOpen("Results"))     run("Clear Results");
}

function ensureDir(path) {
    if (!File.exists(path)) File.makeDirectory(path);
}

function closeIfOpen(title) {
    if (isOpen(title)) { selectWindow(title); close(); }
}

function inArray(val, arr) {
    for (i = 0; i < arr.length; i++) if (arr[i] == val) return true;
    return false;
}

function parseCsvString(s) {
    raw = split(s, ",");
    out = newArray(raw.length);
    for (i = 0; i < raw.length; i++) out[i] = trim(raw[i]);
    return out;
}

function arrToStr(a) {
    s = "";
    for (i = 0; i < a.length; i++) {
        if (i > 0) s = s + ",";
        s = s + a[i];
    }
    return s;
}
