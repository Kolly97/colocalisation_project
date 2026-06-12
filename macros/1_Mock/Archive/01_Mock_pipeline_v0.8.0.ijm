// ==========================================================
// Mock Top-X% Pipeline  --  V0.8.0
// Author 	: Kolja Hildenbrand
// Date   	: 2026-06-02
// Status	: OLD
//
// Changes vs V0.7.2:
//  - REVERT: Back to a single FIXED ARTIFACT_UPPER_BOUND for all channels
//            (the adaptive per-channel bounds of V0.7.0–0.7.2 were hard to
//            calibrate on this data). Set HIGH (5000) so it removes only the
//            most extreme dust and leaves real autofluorescence intact even
//            in brighter images. Most artifacts are caught by the particle-
//            size filter; the cross-image MEDIAN absorbs the rest. Also
//            faster: no per-channel histogram/percentile pass.
//  - REMOVE: computeArtifactBound() and the ARTIFACT_K / ARTIFACT_SEP /
//            ARTIFACT_TRIM_PCT knobs.
//  - NEW:    logs the max pixel intensity per marker channel (within the
//            cytosol, before artifact removal) for every image — so you can
//            see where artifacts sit (e.g. 20000 vs 4000) and tune
//            ARTIFACT_UPPER_BOUND.
//	- NEW:	Max pixel value in csv file as cyto_max_raw
//	- CHANGED	CELL_LINE as newArray to make cell line selection more flexible
//  - CSV:    second artifact column back to `artifact_upper_bound`.
//
// 	- HISTORY: 
//		V0.7.0 tried a Q3+K*IQR fence; 
//		V0.7.1 MULT*p99.9; 
//		V0.7.2 a trimmed mean+k*std fence with a separation gate — all reverted here.
//
// ==========================================================


// ============== 1. CONFIG (globals via `var`) ==============
// Reproducibility
var MACRO_VERSION = "0.8.0";

// Runtime state
var INPUT_DIR;
var MEASURE_DIR;
var MEASURE_DIR_RUN_ID;
var RUN_ID;
var MODE;             // "filename" or "dialog"
var THR_MODE;         // "manual" or "automatic"
var CELL_LINE = newArray("Huh7", "VeroE6", "Other cell line");  // type in cell lines in workflow, will be later overwritten by dialog

// Lists used by CSV setup and per-image dialog
var MARKERS       = newArray("HA568", "HA488", "dsRNA488", "NS4B568");
var ANALYSE_COMBI = newArray("HA568_dsRNA488", "NS4B568_dsRNA488", "NS4B568_HA488");
var TIMEPOINTS    = newArray("12h", "24h");

// Thresholding
var CELL_THR_METHOD  = "Li";  // Li | Otsu | Triangle | Yen
var NUC_THR_METHOD   = "Otsu";      // Otsu | Triangle
var BLUR_SIGMA_CELL  = 1;           // micrometers (Gaussian "scaled")
var BLUR_SIGMA_NUC   = 1;
var NUC_CLOSE_ITER   = 2;

// CELL_THR_FACTOR (NEW V0.4): the auto-threshold's lower bound is MULTIPLIED
// by this before binarisation. 1.0 = use as-is. < 1.0 = MORE permissive
// (lower threshold = more pixels classified as cell). Critical for Vero
// where the auto method is too conservative on dim periphery.
var CELL_THR_FACTOR  = 0.5;

// ARTIFACT EXCLUSION — FIXED upper-bound (V0.7.3, reverted)
// Adaptive per-channel bounds (V0.7.0–0.7.2) proved hard to calibrate on this
// data, so we are back to a single fixed threshold for ALL channels: a pixel
// brighter than ARTIFACT_UPPER_BOUND is treated as an artifact and excluded
// from the cytosol. It is set deliberately HIGH (5000) so it removes only the
// most extreme dust/debris and leaves real autofluorescence intact — even in
// brighter images. The bulk of artifacts is caught by the particle-size filter
// below, and the downstream cross-image MEDIAN of the background is robust
// enough to absorb the few artifact pixels that remain in any single image.
// Bonus: no per-channel histogram/percentile pass -> faster.
var ARTIFACT_UPPER_BOUND = 4000;
var ARTIFACT_DILATE_ITER = 20;       // expand artifact slightly to catch the rim

// PARTICLE-SIZE FILTER (NEW V0.4)
// Cytosol regions smaller than this (in PIXELS) are dropped after the
// cytosol mask is built. Catches isolated bright dots on the coverslip
// outside any real cell that passed thresholding by accident.
var MIN_PARTICLE_SIZE = 200;

// Top-percentile measurement
var TOP_PCT     = 0.01;
var STAT_METHOD = "median_top_hist";

// Output flags
var SAVE_QC            = true;   // single-channel marker1 + cytosol outline (PNG)
var SAVE_MASKS         = true;   // binary mask TIFs (cell, nuc, cyto, artifact)
var SAVE_COMPOSITE_JPG = true;   // 3-channel composite + cytosol outline (JPG)

// CSV header — columns added in V0.4 marked with /* V0.4 */
var CSV_HEADER = "image,cell_line,timepoint,combo,channel,"
               + "stat_method,top_pct,"
               + "macro_mode,threshold_mode,"
               + "threshold_value,n_top_pixels,n_cyto_pixels,n_artifacts_excluded," /* V0.4 */
               + "mean_top,median_top,std_top,"
               + "p95,p99,p99_25,p99_5,p99_9,p99_95,p99_99,p99_995,p99_999,p99_9995,p99_9999,cyto_max_raw,"
               + "cell_thr_method,nuc_thr_method,cell_thr_factor," /* V0.4 */
               + "blur_sigma_cell,blur_sigma_nuc,"
               + "artifact_upper_bound,min_particle_size," /* V0.4 */
               + "macro_version,run_id\n";


// ============== MAIN =======================================

// CRITICAL: pin the "Black Background" binary option ON.
// Without this, Convert to Mask can silently invert (foreground=0
// instead of 255) — depending on Fiji's prior state and on side
// effects of Merge Channels / Stack to RGB. We set it once here
// AND once per iteration (see processOneImage) to be bulletproof.
setOption("BlackBackground", true);

chooseInputDir();
askModeAndConfig();   // also calls applyCellLineDefaults()
RUN_ID = makeRunId();
buildOutputDir();
mockFiles = listMockFiles();
initCsvFiles();

print("=== Mock pipeline V" + MACRO_VERSION + " | run_id=" + RUN_ID + " | mode=" + MODE + " ===");
print("Input  : " + INPUT_DIR);
print("Output : " + MEASURE_DIR_RUN_ID);
print("Cell line: " + CELL_LINE);
print("Markers   : " + arrToStr(MARKERS));
print("Combos    : " + arrToStr(ANALYSE_COMBI));
print("Timepoints: " + arrToStr(TIMEPOINTS));
print("Threshold mode: " + THR_MODE);
print("Artifact upper bound: " + ARTIFACT_UPPER_BOUND);
print("Min particle size  : " + MIN_PARTICLE_SIZE);
print("Found " + mockFiles.length + " mock .tif files.");

for (f = 0; f < mockFiles.length; f++) {
    print("[" + (f+1) + "/" + mockFiles.length + "] " + mockFiles[f]);
    processOneImage(mockFiles[f]);
    cleanupBetweenImages();
}
print("=== DONE ===");
saveLogToFile();   // persist full IJ Log next to the CSVs/masks


// ============== 2. SETUP FUNCTIONS =========================

function chooseInputDir() {
    INPUT_DIR = getDirectory("Choose folder with Mock + MOI .tif images");
    if (INPUT_DIR == "") exit("No folder selected.");
}

// Two-step startup dialog.
function askModeAndConfig() {
    Dialog.create("Pipeline mode");
    Dialog.addMessage("How should the macro know which marker is on which channel\n"
                    + "and how much influence do you want to have?");
    Dialog.addRadioButtonGroup("Mode:",
        newArray("filename", "dialog"), 2, 1, "filename");
    Dialog.addRadioButtonGroup("Threshold mode:",
        newArray("manual", "automatic"), 2, 1, "automatic");
    Dialog.addMessage("What cell line are you analysing?");
    Dialog.addRadioButtonGroup("Cell line:", CELL_LINE, CELL_LINE.length, 1, CELL_LINE[0]);
    Dialog.show();
    MODE      = Dialog.getRadioButton();
    THR_MODE  = Dialog.getRadioButton();
    CELL_LINE = Dialog.getRadioButton();

    // Apply cell-line-specific defaults BEFORE the list dialog,
    // so the user could (in theory) override them below.
    applyCellLineDefaults(CELL_LINE);

    // ALWAYS show the domain-list dialog (Feature A): markers / combos /
    // timepoints are editable at every startup, pre-filled with the
    // workflow standards so a simple OK = use the defaults. This is
    // independent of MODE (which only governs per-image metadata source).
    Dialog.create("Pipeline setup");
    Dialog.addMessage("These values define CSV files and the per-image dialog options.\n"
                    + "Pre-filled with the workflow standards — just click OK to keep them.\n"
                    + "Lists are comma-separated. Whitespace is trimmed.");
    Dialog.addString("Markers:",    arrToStr(MARKERS),       40);
    Dialog.addString("Combos:",     arrToStr(ANALYSE_COMBI), 40);
    Dialog.addString("Timepoints:", arrToStr(TIMEPOINTS),    20);
    Dialog.show();
    MARKERS       = parseCsvString(Dialog.getString());
    ANALYSE_COMBI = parseCsvString(Dialog.getString());
    TIMEPOINTS    = parseCsvString(Dialog.getString());
}

// Cell-line-specific mask parameters. Called after CELL_LINE is known.
// Add new cell lines as additional `else if` branches.
function applyCellLineDefaults(cellLine) {
    print("Cell-line defaults applied for " + cellLine + ":");
    print("  CELL_THR_METHOD   = " + CELL_THR_METHOD);
    print("  CELL_THR_FACTOR   = " + CELL_THR_FACTOR);
    print("  BLUR_SIGMA_CELL   = " + BLUR_SIGMA_CELL);
    print("  MIN_PARTICLE_SIZE = " + MIN_PARTICLE_SIZE);
    print("  ARTIFACT_UPPER_BOUND = " + ARTIFACT_UPPER_BOUND + " (fixed, all channels)");
    print("  ARTIFACT_DILATE_ITER = " + ARTIFACT_DILATE_ITER);
}

function makeRunId() {
    getDateAndTime(y, m, dw, d, h, mn, s, ms);
    return "" + y + IJ.pad(m+1, 2) + IJ.pad(d, 2)
        + "_" + IJ.pad(h, 2) + IJ.pad(mn, 2);
}

function buildOutputDir() {
	
    MEASURE_DIR = INPUT_DIR + "measure_mock" + File.separator;
    if (!File.exists(MEASURE_DIR)) File.makeDirectory(MEASURE_DIR);
    MEASURE_DIR_RUN_ID = MEASURE_DIR + RUN_ID + File.separator;
    if (!File.exists(MEASURE_DIR_RUN_ID))
        File.makeDirectory(MEASURE_DIR_RUN_ID);
    if (SAVE_MASKS && !File.exists(MEASURE_DIR_RUN_ID + "masks"))
        File.makeDirectory(MEASURE_DIR_RUN_ID + "masks");
    if (SAVE_QC && !File.exists(MEASURE_DIR_RUN_ID + "qc"))
        File.makeDirectory(MEASURE_DIR_RUN_ID + "qc");
}

function initCsvFiles() {
    for (i = 0; i < ANALYSE_COMBI.length; i++) {
        combo = ANALYSE_COMBI[i];
        parts = split(combo, "_");
        m1 = parts[0]; m2 = parts[1];
        for (t = 0; t < TIMEPOINTS.length; t++) {
            tp = TIMEPOINTS[t];
            csv1 = MEASURE_DIR_RUN_ID + CELL_LINE + "_mock_" + tp + "_" + m1 + "_in_" + combo + ".csv";
            csv2 = MEASURE_DIR_RUN_ID + CELL_LINE + "_mock_" + tp + "_" + m2 + "_in_" + combo + ".csv";
            if (!File.exists(csv1)) File.append(CSV_HEADER, csv1);
            if (!File.exists(csv2)) File.append(CSV_HEADER, csv2);
        }
    }
}

function listMockFiles() {
    files = getFileList(INPUT_DIR);
    out = newArray();
    for (i = 0; i < files.length; i++) {
        n = files[i]; ln = toLowerCase(n);
        if (endsWith(ln, ".tif") && indexOf(ln, "mock") >= 0)
            out = Array.concat(out, n);
    }
    out = Array.sort(out);
    return out;
}


// ============== 3. PER-IMAGE PIPELINE ======================

function processOneImage(fname) {
    // Re-pin Black Background per iteration — paranoid defense against
    // any side effect from save/cleanup operations that might have
    // toggled the global binary state.
    setOption("BlackBackground", true);

    open(INPUT_DIR + fname);
    title = getTitle();
    imgName = substring(title, 0, lastIndexOf(title, "."));

    parsed = tryParseFilename(imgName);

    // ---- decide metadata source -----------------------------
    tp = ""; cellLine = ""; m1 = ""; m2 = "";
    if (MODE == "filename") {
        if (parsed.length == 0) {
            print("  Parse failed for: " + imgName);
            action = askParseFailureAction(imgName);
            if (action == "skip") { print("  -> skipped"); return; }
            meta = askImageMetadata("(parse failed)", "", "", "", "");
            tp = meta[0]; cellLine = meta[1]; m1 = meta[2]; m2 = meta[3];
        } else {
            tp = parsed[0]; cellLine = parsed[1]; m1 = parsed[2]; m2 = parsed[3];
        }
    } else {
        defTp = ""; defM1 = ""; defM2 = "";
        if (parsed.length > 0) { defTp = parsed[0]; defM1 = parsed[2]; defM2 = parsed[3]; }
        meta = askImageMetadata(imgName, defTp, CELL_LINE, defM1, defM2);
        tp = meta[0]; cellLine = meta[1]; m1 = meta[2]; m2 = meta[3];
    }

    // ---- validate -----------------------------------------
    comboKey = m1 + "_" + m2;
    if (!isValidCombo(comboKey)) {
        print("  SKIP: combo not in ANALYSE_COMBI -> " + comboKey);
        return;
    }
    if (!inArray(tp, TIMEPOINTS)) {
        print("  SKIP: timepoint not in TIMEPOINTS -> " + tp);
        return;
    }
    print("  combo=" + comboKey + " tp=" + tp + " cellLine=" + cellLine);

    // ---- channels -----------------------------------------
    splitAndRenameChannels(title, m1, m2);
    cellSrc = pickCellMaskSource(m1, m2);
    if (cellSrc == "") {
        print("  SKIP: no marker channel suitable as Cell-Mask source");
        return;
    }
    print("  Cell-Mask source: " + cellSrc);

    // ---- masks --------------------------------------------
    if (THR_MODE == "automatic") {
        makeCellMask(cellSrc);
        makeNucleusMask("DAPI_channel");
    } else {
        makeCellMaskManual(cellSrc);
        makeNucleusMaskManual("DAPI_channel");
    }
    makeCytosolMask();

    // Brightest pixel per marker channel WITHIN the cytosol, BEFORE artifact
    // removal — so any artifact still shows up as the max. Logged AND written to
    // the CSV (cyto_max_raw) so you can see where artifacts sit (e.g. 20000 vs
    // 4000) and tune ARTIFACT_UPPER_BOUND.
    cytoMaxM1 = logChannelMaxInCytosol(m1 + "_channel", m1, "C1");
    cytoMaxM2 = logChannelMaxInCytosol(m2 + "_channel", m2, "C2");

    // ---- artifact + particle cleanup ---------------------
    // Fixed high upper-bound (V0.7.3): removes only extreme dust; the particle
    // filter + cross-image median handle the rest. cleanCytosolMask returns how
    // many pixels were removed by artifact exclusion (written to the CSV).
    nArtifactsExcluded = cleanCytosolMask(m1, m2);

    // Guard: empty / too-small cytosol = nothing meaningful to measure
    selectWindow("Cytosol_Mask");
    run("Invert");
    getStatistics(area, meanCyto);
    if (meanCyto == 0) {
        print("  SKIP: cytosol mask empty after cleanup");
        return;
    }

    // Make cytosol foreground a re-applicable ROI
    roiManager("reset");
    selectWindow("Cytosol_Mask");
    run("Create Selection");
    roiManager("Add");
    cytoRoiId = 0;
    nCyto = getRawStatisticsCount();
    if (nCyto < 1000) {
        print("  WARN: cytosol very small after cleanup (" + nCyto + " px) — measurements may be unreliable");
    }

    // ---- measure each marker ------------------------------
    measureAndWrite(m1, m1 + "_channel", cytoRoiId, nCyto, nArtifactsExcluded, cytoMaxM1, imgName, cellLine, tp, comboKey);
    measureAndWrite(m2, m2 + "_channel", cytoRoiId, nCyto, nArtifactsExcluded, cytoMaxM2, imgName, cellLine, tp, comboKey);

    // ---- save artefacts -----------------------------------
    // Order matters: composite JPG must come BEFORE saveMasksTif,
    // because saveMasksTif uses saveAs() which renames the mask
    // windows — by then the m1/m2/DAPI channel windows are still
    // intact, but if any later step closes them the composite
    // wouldn't have the data to merge. Doing JPG first is safe.
    if (SAVE_QC)            saveQcOverlay(imgName, m1);
    if (SAVE_COMPOSITE_JPG) saveCompositeJpg(imgName, m1, m2);
    if (SAVE_MASKS)         saveMasksTif(imgName);
}

function tryParseFilename(imgName) {
    tokens = split(imgName, "_");
    if (tokens.length < 7) return newArray();
    tp_       = tokens[0];
    cellLine_ = tokens[1];
    m1_       = tokens[3];
    m2_       = tokens[4];
    if (!inArray(m1_, MARKERS) || !inArray(m2_, MARKERS)) return newArray();
    if (!inArray(tp_, TIMEPOINTS))                       return newArray();
    return newArray(tp_, cellLine_, m1_, m2_);
}

function askImageMetadata(imgLabel, defTp, defCellLine, defM1, defM2) {
    if (defTp == "" || !inArray(defTp, TIMEPOINTS)) defTp = TIMEPOINTS[0];
    if (defM1 == "" || !inArray(defM1, MARKERS))    defM1 = MARKERS[0];
    if (defM2 == "" || !inArray(defM2, MARKERS))    defM2 = MARKERS[1 % MARKERS.length];
    if (defCellLine == "") defCellLine = "Huh7";

    Dialog.create("Image metadata");
    Dialog.addMessage("Image: " + imgLabel + "\nC3 is always DAPI.");
    Dialog.addString("Cell line:", defCellLine, 12);
    Dialog.addRadioButtonGroup("Timepoint:", TIMEPOINTS, 1, TIMEPOINTS.length, defTp);
    Dialog.addRadioButtonGroup("Channel 1 (C1):", MARKERS, 1, MARKERS.length, defM1);
    Dialog.addRadioButtonGroup("Channel 2 (C2):", MARKERS, 1, MARKERS.length, defM2);
    Dialog.show();

    cellLine_ = Dialog.getString();
    tp_       = Dialog.getRadioButton();
    m1_       = Dialog.getRadioButton();
    m2_       = Dialog.getRadioButton();

    if (m1_ == m2_) {
        showMessage("Error", "C1 and C2 must be different markers.");
        exit("Aborted: identical markers for C1 and C2.");
    }
    return newArray(tp_, cellLine_, m1_, m2_);
}

function askParseFailureAction(imgName) {
    Dialog.create("Filename parse failed");
    Dialog.addMessage("Cannot parse: " + imgName + "\n\nWhat do you want to do?");
    Dialog.addRadioButtonGroup("Action:",
        newArray("skip", "dialog"), 2, 1, "skip");
    Dialog.show();
    return Dialog.getRadioButton();
}


// ============== 4. CHANNEL HANDLING ========================

// FIX V0.4: intermediate variable — avoids the IJM
// "return userFunction(...)" type-inference bug.
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

function pickCellMaskSource(m1, m2) {
    priority = newArray("HA568", "HA488", "NS4B568", "dsRNA488");
    for (i = 0; i < priority.length; i++) {
        cand = priority[i] + "_channel";
        if (isOpen(cand)) return cand;
    }
    return "";
}


// ============== 5. MASKS ===================================

// Cell mask with auto-threshold + CELL_THR_FACTOR scaling.
// The factor lets us be more permissive without changing the method:
// the method adapts per image, factor sets global strictness.
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
    print("  cell thr (" + CELL_THR_METHOD + "): auto=" + lo
        + ", scaled by " + CELL_THR_FACTOR + " -> " + loNew);

    run("Convert to Mask");
    run("Fill Holes");
}

function makeCellMaskManual(srcTitle) {
    selectWindow(srcTitle);
    run("Duplicate...", "title=Cell_Mask");
    run("Enhance Contrast", "saturated=0.35");
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_CELL + " scaled");
    run("Threshold...");
    setAutoThreshold(CELL_THR_METHOD + " dark");
    waitForUser("Check and adjust the threshold if necessary, then click OK.");
    run("Convert to Mask");
    run("Fill Holes");
    run("Close");
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

function makeNucleusMaskManual(dapiTitle) {
    selectWindow(dapiTitle);
    run("Duplicate...", "title=Nucleus_Mask");
    run("Enhance Contrast", "saturated=0.35");
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_NUC + " scaled");
    run("Threshold...");
    setAutoThreshold(NUC_THR_METHOD + " dark");
    waitForUser("Check and adjust the threshold if necessary, then click OK.");
    run("Convert to Mask");
    run("Close");
    if (NUC_CLOSE_ITER > 0) {
        run("Options...", "iterations=" + NUC_CLOSE_ITER + " count=1 black do=Nothing");
        run("Dilate");
        run("Fill Holes");
        run("Erode");
    } else {
        run("Fill Holes");
    }
}

// Cytosol = Cell − Nucleus. 8-bit subtraction clamps at 0, giving a clean
// binary: cell-not-nucleus = 255, else 0.
function makeCytosolMask() {
    imageCalculator("Subtract create", "Cell_Mask", "Nucleus_Mask");
    rename("Cytosol_Mask");
    setThreshold(1, 255);
    run("Convert to Mask");
}


// ============== 5b. ARTIFACT + PARTICLE CLEANUP (V0.4) =====
//
// Strategy:
//   1) Build per-channel artifact masks (pixel > the channel's adaptive bound).
//   2) OR them into one combined Artifact_Mask. Pixels too bright in
//      EITHER channel are excluded — important for downstream coloc,
//      where an artifact in one channel must not be measured against
//      the other.
//   3) Subtract Artifact_Mask from Cytosol_Mask.
//   4) Run Analyze Particles to drop disconnected regions smaller
//      than MIN_PARTICLE_SIZE — handles tiny bright dots on the
//      coverslip that survived (1).
//
// Why both? Different artifact populations:
//   - Inside cells: bright dust/debris → caught by (1).
//   - Outside cells: small bright fragments that the cell threshold
//     captured by accident → caught by (4).
// They complement each other; we don't have to catch 100 % to fix
// the Mock median (we only need to drop the worst outliers).
// ===========================================================

// V0.7.3: simple FIXED upper-bound artifact exclusion. A pixel brighter than
// ARTIFACT_UPPER_BOUND (same value for every channel) is treated as an artifact.
// Set deliberately HIGH (5000) so it removes only extreme dust/debris and leaves
// real autofluorescence — even in brighter images — untouched. The bulk of
// artifacts is handled by the particle-size filter, and the downstream
// cross-image MEDIAN of the background absorbs the few that remain. Also faster:
// no per-channel histogram/percentile pass per image.
// Log the brightest pixel of each marker channel inside the cytosol. Helps
// calibrate ARTIFACT_UPPER_BOUND later (you see, per image, whether artifacts
// sit at e.g. 20000 or only 4000). Computed on the RAW channel before any
// artifact removal so a present artifact still registers as the max.
// Returns the brightest RAW pixel value of `chTitle` inside the cytosol (and
// logs it). Returns -1 if the channel/cytosol is unavailable. Called BEFORE
// artifact removal, so a present artifact still shows up as the max — this is
// the value written to the CSV as cyto_max_raw, to help tune ARTIFACT_UPPER_BOUND.
function logChannelMaxInCytosol(chTitle, marker, chLabel) {
    if (!isOpen(chTitle)) return -1;
    if (countMaskForeground("Cytosol_Mask") == 0) {
        print("  max intensity   " + chLabel + " " + marker + " : (empty cytosol)");
        return -1;
    }
    selectWindow("Cytosol_Mask");
    run("Create Selection");
    selectWindow(chTitle);
    run("Restore Selection");
    getRawStatistics(nPix, meanV, minV, maxV);   // raw (uncalibrated) intensities
    run("Select None");
    selectWindow("Cytosol_Mask"); run("Select None");
    print("  max intensity   " + chLabel + " " + marker + " : max=" + maxV
        + "  (cytosol; ARTIFACT_UPPER_BOUND=" + ARTIFACT_UPPER_BOUND + ")");
    return maxV;
}

function makeArtifactMaskFor(srcTitle, maskName) {
    selectWindow(srcTitle);
    run("Select None");
    run("Duplicate...", "title=" + maskName);
    setThreshold(ARTIFACT_UPPER_BOUND, 65535);
    run("Convert to Mask");
    if (ARTIFACT_DILATE_ITER > 0) {
        run("Options...", "iterations=" + ARTIFACT_DILATE_ITER + " count=1 black do=Nothing");
        run("Dilate");
    }
}

function buildCombinedArtifactMask(m1, m2) {
    makeArtifactMaskFor(m1 + "_channel", "Artifact_M1");
    makeArtifactMaskFor(m2 + "_channel", "Artifact_M2");
    imageCalculator("OR create", "Artifact_M1", "Artifact_M2");
    rename("Artifact_Mask");
    if (isOpen("Artifact_M1")) { selectWindow("Artifact_M1"); close(); }
    if (isOpen("Artifact_M2")) { selectWindow("Artifact_M2"); close(); }
}

// Apply artifact exclusion + particle-size filter to Cytosol_Mask.
// Returns the number of pixels removed by artifact exclusion (for CSV log).
function cleanCytosolMask(m1, m2) {
    // Count cytosol pixels BEFORE cleanup
    nBefore = countMaskForeground("Cytosol_Mask");

    // (1+2+3) Artifact subtraction (fixed upper-bound, both channels)
    buildCombinedArtifactMask(m1, m2);
    imageCalculator("Subtract create", "Cytosol_Mask", "Artifact_Mask");
    rename("Cytosol_Clean");
    if (isOpen("Cytosol_Mask")) { selectWindow("Cytosol_Mask"); close(); }
    selectWindow("Cytosol_Clean"); rename("Cytosol_Mask");
    setThreshold(1, 255);
    run("Convert to Mask");

    nAfterArtifact = countMaskForeground("Cytosol_Mask");
    nArtifactsExcluded = nBefore - nAfterArtifact;
    print("  artifact excl: " + nArtifactsExcluded + " px removed (threshold > " + ARTIFACT_UPPER_BOUND + ")");

    // (4) Particle-size filter
    selectWindow("Cytosol_Mask");
    run("Analyze Particles...", "size=" + MIN_PARTICLE_SIZE + "-Infinity show=Masks");
    // Output window: "Mask of Cytosol_Mask"
    if (isOpen("Cytosol_Mask")) { selectWindow("Cytosol_Mask"); close(); }
    selectWindow("Mask of Cytosol_Mask"); rename("Cytosol_Mask");
    nAfterParticle = countMaskForeground("Cytosol_Mask");
    print("  particle filter (>= " + MIN_PARTICLE_SIZE + " px): " + (nAfterArtifact - nAfterParticle) + " px removed");

    return nArtifactsExcluded;
}

// Count foreground pixels (value == 255) of a binary mask via its histogram.
// IJM quirk: getHistogram(values, counts, nBins, histMin, histMax) ONLY works
// for 16- or 32-bit images. On 8-bit images (our binary masks) you MUST call
// the 3-arg form — IJM then uses the fixed 0..255 range automatically.
function countMaskForeground(maskTitle) {
    selectWindow(maskTitle);
    getHistogram(values, counts, 256);   // no min/max for 8-bit!
    return counts[255];
}


// ============== 6. MEASUREMENT (TOP-X% via histogram) ======
// IJM quirk: getHistogram with nBins=65536, range 0..65535 does NOT
// populate values[]. We use bin index `i` directly as the pixel value.
// ===========================================================

function measureAndWrite(marker, chTitle, cytoRoiId, nCyto, nArtifactsExcluded, cytoMax,
                         imgName, cellLine, tp, comboKey) {
    if (!isOpen(chTitle)) {
        print("  WARN: channel window missing: " + chTitle);
        return;
    }
    selectWindow(chTitle);
    roiManager("Select", cytoRoiId);

    nBins = 65536;
    getHistogram(values, counts, nBins, 0, 65535);

    nTotal = 0;
    for (i = 0; i < nBins; i++) nTotal += counts[i];
    if (nTotal == 0) {
        print("  WARN: 0 pixels in cytosol selection for " + chTitle);
        return;
    }

    stats = computeTopStats(counts, nBins, nTotal, TOP_PCT);
    csvPath = MEASURE_DIR_RUN_ID + CELL_LINE + "_mock_" + tp + "_" + marker + "_in_" + comboKey + ".csv";
    appendCsvRow(csvPath, imgName, cellLine, tp, comboKey, marker, stats, nCyto, nArtifactsExcluded, cytoMax);
}

function computeTopStats(counts, nBins, nTotal, topPct) {
    // (1) Find threshold V*
    nTopTarget = floor(nTotal * topPct / 100.0);
    if (nTopTarget < 1) nTopTarget = 1;
    acc = 0;
    threshold_value = 0;
    for (i = nBins - 1; i >= 0; i--) {
        if (counts[i] > 0) {
            acc += counts[i];
            if (acc >= nTopTarget) { threshold_value = i; break; }
        }
    }

    // (2) Mean and std over the top pool
    nTop = 0; sumV = 0; sumV2 = 0;
    for (i = threshold_value; i < nBins; i++) {
        if (counts[i] > 0) {
            nTop  += counts[i];
            sumV  += i * counts[i];
            sumV2 += i * i * counts[i];
        }
    }
    mean_top = sumV / nTop;
    var_top  = (sumV2 / nTop) - mean_top * mean_top;
    std_top  = sqrt(maxOf(var_top, 0));

    // (3) Median of the top pool
    half = nTop / 2.0;
    a = 0; median_top = threshold_value;
    for (i = threshold_value; i < nBins; i++) {
        if (counts[i] > 0) {
            a += counts[i];
            if (a >= half) { median_top = i; break; }
        }
    }

    // (4) Whole-cytosol percentiles for sanity / outlier detection
    p95     = pctIndex(counts, nBins, nTotal, 95);
    p99     = pctIndex(counts, nBins, nTotal, 99);
    p99_25  = pctIndex(counts, nBins, nTotal, 99.25);
    p99_5   = pctIndex(counts, nBins, nTotal, 99.5);
    p99_9   = pctIndex(counts, nBins, nTotal, 99.9);
    p99_95  = pctIndex(counts, nBins, nTotal, 99.95);
    p99_99  = pctIndex(counts, nBins, nTotal, 99.99);
    p99_995 = pctIndex(counts, nBins, nTotal, 99.995);
    p99_999 = pctIndex(counts, nBins, nTotal, 99.999);
    p99_9995 = pctIndex(counts, nBins, nTotal, 99.9995);
    p99_9999 = pctIndex(counts, nBins, nTotal, 99.9999);
    pmax = pctIndex(counts, nBins, nTotal, 100);

    return newArray(threshold_value, nTop, mean_top, median_top, std_top,
                    p95, p99, p99_25, p99_5, p99_9, p99_95, p99_99, p99_995, p99_999, p99_9995, p99_9999, pmax);
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

function getRawStatisticsCount() {
    getRawStatistics(n);
    return n;
}


// ============== 7. CSV / OUTPUT ============================

function appendCsvRow(csvPath, imgName, cellLine, tp, comboKey, channel,
                      stats, nCyto, nArtifactsExcluded, cytoMax) {
    line = imgName + "," + cellLine + "," + tp + "," + comboKey + "," + channel
        + "," + STAT_METHOD + "," + TOP_PCT
        + "," + MODE + "," + THR_MODE
        + "," + stats[0] + "," + stats[1] + "," + nCyto + "," + nArtifactsExcluded
        + "," + stats[2] + "," + stats[3] + "," + stats[4]
        + "," + stats[5] + "," + stats[6] + "," + stats[7] + "," + stats[8] + "," + stats[9]
        + "," + stats[10] + "," + stats[11] + "," + stats[12] + "," + stats[13] + "," + stats[14] + "," + stats[15] + "," + cytoMax
        + "," + CELL_THR_METHOD + "," + NUC_THR_METHOD + "," + CELL_THR_FACTOR
        + "," + BLUR_SIGMA_CELL + "," + BLUR_SIGMA_NUC
        + "," + ARTIFACT_UPPER_BOUND + "," + MIN_PARTICLE_SIZE
        + "," + MACRO_VERSION + "," + RUN_ID;
    File.append(line, csvPath);
}

function saveMasksTif(imgName) {
    masksDir = MEASURE_DIR_RUN_ID + "masks" + File.separator;
    filename = RUN_ID + "_" + imgName;
    if (isOpen("Cell_Mask"))     { selectWindow("Cell_Mask");     saveAs("Tiff", masksDir + filename + "_cell.tif"); }
    if (isOpen("Nucleus_Mask"))  { selectWindow("Nucleus_Mask");  saveAs("Tiff", masksDir + filename + "_nuc.tif"); }
    if (isOpen("Cytosol_Mask"))  { selectWindow("Cytosol_Mask");  saveAs("Tiff", masksDir + filename + "_cyto.tif"); }
    if (isOpen("Artifact_Mask")) { selectWindow("Artifact_Mask"); saveAs("Tiff", masksDir + filename + "_artifact.tif"); }
}

function saveQcOverlay(imgName, m1) {
    qcDir = MEASURE_DIR_RUN_ID + "qc" + File.separator;
    filename = RUN_ID + "_" + imgName + "_qc.png";
    src = m1 + "_channel";
    if (!isOpen(src)) return;
    selectWindow(src);
    run("Duplicate...", "title=qc_tmp");
    run("Enhance Contrast", "saturated=0.35");
    run("8-bit");
    selectWindow("Cytosol_Mask");
    run("Create Selection");
    selectWindow("qc_tmp");
    run("Restore Selection");
    run("Add Selection...");
    run("Flatten");
    saveAs("PNG", qcDir + filename);
    close();
    if (isOpen("qc_tmp")) { selectWindow("qc_tmp"); close(); }
}

// ------------------------------
// Save an 8-bit RGB composite JPG with all 3 channels merged and
// the cytosol ROI drawn as outline. Lets you visually verify:
//   - Artefacts in EITHER channel are excluded from the ROI
//     (red dot inside ROI = uncaught m1 artefact; green dot inside
//      ROI = uncaught m2 artefact)
//   - Nuclei (blue) are cleanly outside the ROI
//   - The cytosol mask follows the actual cell shape
//
// Conventional channel-to-color mapping for our data:
//   c1=m1_channel  -> RED   (568 marker)
//   c2=m2_channel  -> GREEN (488 marker)
//   c3=DAPI        -> BLUE
//
// Key implementation details:
//  - Merge Channels uses `keep` so the source windows aren't
//    consumed; saveMasksTif and downstream cleanup still see them.
//  - Stack.setDisplayMode("composite") -> all channels overlay
//    in their LUT colors (what we want for visual QC).
//  - Per-channel Enhance Contrast (saturated=0.35) is cosmetic
//    — it adjusts each channel's display LUT so faint signal is
//    visible. Does not change pixel data.
//  - Stack to RGB flattens the composite into one 8-bit RGB image
//    which is what JPG can store.
//  - JPG quality 90 = good balance figure-quality vs file size.
// ------------------------------
function saveCompositeJpg(imgName, m1, m2) {
    qcDir = MEASURE_DIR_RUN_ID + "qc" + File.separator;
    jpgPath = qcDir + RUN_ID + "_" + imgName + "_composite.jpg";

    // Defensive: we need all 3 source channels for the merge.
    if (!isOpen(m1 + "_channel") || !isOpen(m2 + "_channel") || !isOpen("DAPI_channel")) {
        print("  WARN: composite jpg skipped (missing channel window)");
        return;
    }

    // CHANGED V0.5.1: instead of `Merge Channels ... keep` (which leaves
    // the source channels alive but can subtly alter their LUTs / display
    // metadata via Stack.setChannel + Enhance Contrast below), we
    // DUPLICATE each source channel first and merge the duplicates
    // WITHOUT keep. The originals are then guaranteed untouched —
    // no shared state, no LUT leak, no risk of breaking the masks
    // for the next iteration.
    selectWindow(m1 + "_channel"); run("Select None"); run("Duplicate...", "title=c_m1_tmp");
    selectWindow(m2 + "_channel"); run("Select None"); run("Duplicate...", "title=c_m2_tmp");
    selectWindow("DAPI_channel");  run("Select None"); run("Duplicate...", "title=c_dapi_tmp");

    // Merge the DUPLICATES — they get consumed (closed automatically),
    // which is exactly what we want.
    run("Merge Channels...", "c1=c_m1_tmp c2=c_m2_tmp c3=c_dapi_tmp create");
    // Active window is now "Composite"

    // Composite display mode: all channels visible at once in their LUTs.
    Stack.setDisplayMode("composite");

    // Per-channel auto-contrast so dim signal is visible in the JPG.
    Stack.setChannel(1); run("Enhance Contrast", "saturated=0.35");
    Stack.setChannel(2); run("Enhance Contrast", "saturated=0.35");
    Stack.setChannel(3); run("Enhance Contrast", "saturated=0.35");

    // Flatten the composite into a single 8-bit RGB image.
    run("Stack to RGB");
    // New active window: "Composite (RGB)"

    // Overlay cytosol ROI. Use ROI Manager (index 0) — already added
    // in processOneImage and not yet reset.
    if (roiManager("count") > 0) {
        roiManager("Select", 0);
        run("Add Selection...");
    }

    // Burn overlay into pixel data so the JPG actually contains it.
    run("Flatten");

    // Save with high JPEG quality (figure-quality, ~3-5× smaller than PNG).
    run("Options...", "jpeg=90");
    saveAs("Jpeg", jpgPath);
    print("  saved composite jpg -> " + jpgPath);

    // Cleanup. The flattened+saved JPG window is active now → close.
    // Composite (RGB) and Composite (the duplicate-based one) may
    // still be open → close them. The original m1/m2/DAPI channels
    // are completely untouched throughout this function.
    close();
    if (isOpen("Composite (RGB)")) { selectWindow("Composite (RGB)"); close(); }
    if (isOpen("Composite"))       { selectWindow("Composite");       close(); }

    // CRITICAL: re-pin Black Background. Stack to RGB / Merge Channels
    // can flip this option as a side effect on some Fiji builds. If
    // we don't re-pin here, the NEXT image's Convert to Mask might
    // invert (foreground=0, background=255) → broken nucleus mask.
    setOption("BlackBackground", true);
}


// ============== 8. CLEANUP & UTILS =========================

// Save the IJ Log window to a text file alongside the CSVs.
// Call this LAST in MAIN so the file captures every print() up to here,
// including per-image threshold values, artifact counts, particle filter
// counts, and any warnings. Essential for reproducibility and debugging.
function saveLogToFile() {
    if (!isOpen("Log")) return;
    logPath = MEASURE_DIR_RUN_ID + "macro_log_" + RUN_ID + ".txt";
    selectWindow("Log");
    saveAs("Text", logPath);
    print("Saved log: " + logPath);
}

function cleanupBetweenImages() {
    while (nImages > 0) { selectImage(nImages); close(); }
    if (isOpen("ROI Manager")) roiManager("reset");
    if (isOpen("Results"))     run("Clear Results");
    if (selectionType() != -1) run("Select None");
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

function inArray(val, arr) {
    for (i = 0; i < arr.length; i++) if (arr[i] == val) return true;
    return false;
}
