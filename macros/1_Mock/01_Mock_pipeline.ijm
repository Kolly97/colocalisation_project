// ==========================================================
// Mock Top-X% Pipeline  --  V0.3
// Author 	: Kolja Hildenbrand
// Date   	: 2026-05-25
// Status	: OLD
//
// Changes vs V0.2:
// ==========================================================
// Major changes vs v0.2:
//  - NEW:   Added threshold mode selector: automatic vs manual.
//
//  - NEW:   Manual threshold workflow for cell and nucleus masks.
//           The user can inspect and adjust thresholds before
//           mask conversion.
//
//  - NEW:   Added makeCellMaskManual() and
//           makeNucleusMaskManual().
//
//  - NEW:   Added threshold_mode to CSV output.
//
//  - CHANGE: processOneImage() now switches between automatic
//            and manual mask generation based on THR_MODE.
//
//  - CHANGE: CSV files are now written into run-specific output
//            folders and include the selected cell line in the name.
//
//  - REMOVE: No major functionality removed.
// ==========================================================


// ============== 1. CONFIG (globals via `var`) ==============
// In dialog mode the startup dialog OVERWRITES the list-type
// defaults below; in filename mode the defaults are used as-is.
// -----------------------------------------------------------

// Set during main
var INPUT_DIR;
var MEASURE_DIR;
var MEASURE_DIR_RUN_ID;
var RUN_ID;
var MODE;        // "filename" or "dialog"
var THR_MODE; 	// "manual"or "automatic"
var CELL_LINE = newArray("Huh7", "VeroE6");// dialog mode: set at startup; filename mode: per image (overwritten)

// Lists used by CSV setup and per-image dialog
var MARKERS       = newArray("HA568", "HA488", "dsRNA488", "NS4B568");
var ANALYSE_COMBI = newArray("HA568_dsRNA488", "NS4B568_dsRNA488", "NS4B568_HA488");
var TIMEPOINTS    = newArray("12h", "24h");

// Thresholding
var CELL_THR_METHOD = "Li";    // try: Li, Otsu, Triangle, Yen
var NUC_THR_METHOD  = "Otsu";  // try: Otsu, Triangle
var BLUR_SIGMA_CELL = 1;       // in MICROMETERS (Gaussian "scaled")
var BLUR_SIGMA_NUC  = 1;       // in MICROMETERS (Gaussian "scaled")
var NUC_CLOSE_ITER  = 2;       // morphological closing passes for nucleus

// Top-percentile measurement
var TOP_PCT     = 1.0;
var STAT_METHOD = "median_top_hist";

// Output flags
var SAVE_QC    = true;
var SAVE_MASKS = true;

// Reproducibility
var MACRO_VERSION = "0.3.1";

// CSV header (must match appendCsvRow order)
var CSV_HEADER = "image,cell_line,timepoint,combo,channel,"
               + "stat_method,top_pct,"
               + "macro_mode,threshold_mode,"
               + "threshold_value,n_top_pixels,n_cyto_pixels,"
               + "mean_top,median_top,std_top,"
               + "p95,p99,p99_25,p99_5,p99_9,p99_95,p99_99,p99_995,p99_999,"
               + "cell_thr_method,nuc_thr_method,"
               + "blur_sigma_cell,blur_sigma_nuc,"
               + "macro_version,run_id\n";


// ============== MAIN =======================================
chooseInputDir();
askModeAndConfig();   // may overwrite MODE/CELL_LINE/MARKERS/ANALYSE_COMBI/TIMEPOINTS
RUN_ID = makeRunId();
buildOutputDir();
mockFiles = listMockFiles();
initCsvFiles();

print("=== Mock pipeline V" + MACRO_VERSION + " | run_id=" + RUN_ID + " | mode=" + MODE + " ===");
print("Input  : " + INPUT_DIR);
print("Output : " + MEASURE_DIR);
print("Cell line (startup): " + CELL_LINE);
print("Markers   : " + arrToStr(MARKERS));
print("Combos    : " + arrToStr(ANALYSE_COMBI));
print("Timepoints: " + arrToStr(TIMEPOINTS));
print("Script Mode: " + MODE);
print("Threshold Mode: " + THR_MODE);
print("Found " + mockFiles.length + " mock .tif files.");

// No batch mode in v0.3 — dialogs need a visible GUI and you
// want to see the masks anyway while testing.
for (f = 0; f < mockFiles.length; f++) {
    print("[" + (f+1) + "/" + mockFiles.length + "] " + mockFiles[f]);
    processOneImage(mockFiles[f]);
    cleanupBetweenImages();
}
print("=== DONE ===");


// ============== 2. SETUP FUNCTIONS =========================

function chooseInputDir() {
    INPUT_DIR = getDirectory("Choose folder with Mock + MOI .tif images");
    if (INPUT_DIR == "") exit("No folder selected.");
}

// Two-step startup dialog.
//   Step 1: pick mode (filename / dialog) and Threshold mode (manual / automatic)
//   Step 2 (dialog mode only): edit CELL_LINE/MARKERS/COMBOS/TIMEPOINTS.
function askModeAndConfig() {
    // ---- Step 1: mode ----
    // addRadioButtonGroup(label, items[], rows, columns, defaultItem)
    // rows × columns defines the visual layout of the buttons.
    Dialog.create("Pipeline mode");
    Dialog.addMessage("How should the macro know which marker is on which channel\n" + "and how much influence do you want to have?");
    Dialog.addRadioButtonGroup("Mode:",
        newArray("filename", "dialog"), 2, 1, "filename");
    Dialog.addRadioButtonGroup("Threshold Moudus:",
        newArray("manual", "automatic"), 2, 1, "automatic");
    Dialog.addMessage("What cell line are you analysing?");
    Dialog.addRadioButtonGroup("Cell Line:",
        newArray("Huh7", "VeroE6"), 2, 1, "Huh7");
    Dialog.show();
    MODE = Dialog.getRadioButton();
    THR_MODE = Dialog.getRadioButton();
    CELL_LINE = Dialog.getRadioButton();

    // ---- Step 2: config (only in dialog mode) ----
    // In filename mode the var defaults at the top of the file
    // are used. CELL_LINE is taken from each filename's token[1].
    if (MODE == "dialog") {
        Dialog.create("Pipeline setup (dialog mode)");
        Dialog.addMessage("These values define CSV files and the per-image dialog options.\n"
                        + "Lists are comma-separated. Whitespace is trimmed.");
        Dialog.addString("Markers:",      arrToStr(MARKERS),       40);
        Dialog.addString("Combos:",       arrToStr(ANALYSE_COMBI), 40);
        Dialog.addString("Timepoints:",   arrToStr(TIMEPOINTS),    20);
        Dialog.show();

        // get* calls MUST be in the same order as add* calls.
        MARKERS       = parseCsvString(Dialog.getString());
        ANALYSE_COMBI = parseCsvString(Dialog.getString());
        TIMEPOINTS    = parseCsvString(Dialog.getString());
    } 
}

// ISO-ish timestamp without separators.
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

// Create one CSV per (timepoint, combo, marker), header on first
// creation only. Number of CSVs = len(TIMEPOINTS) × len(COMBOS) × 2.
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
//
// Per image we need 4 pieces of metadata:
//   tp        : timepoint (e.g. "12h")
//   cellLine  : cell line (e.g. "Huh7")
//   m1        : marker on channel 1
//   m2        : marker on channel 2
//
// Source of these depends on MODE:
//   filename : parse from tokens; on failure fall back to dialog.
//   dialog   : ask interactively per image (with filename-parsed
//              defaults if parsing succeeded).
// ===========================================================

function processOneImage(fname) {
    open(INPUT_DIR + fname);
    title = getTitle();
    imgName = substring(title, 0, lastIndexOf(title, "."));

    // tryParseFilename returns:
    //   newArray(tp, cellLine, m1, m2)  on success
    //   empty array                     on failure
    parsed = tryParseFilename(imgName);

    // ---- decide metadata source -----------------------------
    tp = ""; cellLine = ""; m1 = ""; m2 = "";

    if (MODE == "filename") {
        if (parsed.length == 0) {
            print("  Parse failed for: " + imgName);
            action = askParseFailureAction(imgName);
            if (action == "skip") { print("  -> skipped"); return; }
            // else: continue into dialog
            meta = askImageMetadata("(parse failed)", "", "", "", "");
            tp = meta[0]; cellLine = meta[1]; m1 = meta[2]; m2 = meta[3];
        } else {
            tp = parsed[0]; cellLine = parsed[1]; m1 = parsed[2]; m2 = parsed[3];
        }
    } else {
        // dialog mode: pre-fill from filename if parse worked, else blanks
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
    }
    
    if (THR_MODE == "manual") {
	    makeCellMaskManual(cellSrc);
	    makeNucleusMaskManual("DAPI_channel");
    }
    
    makeCytosolMask();

    // Guard: empty cytosol = nothing to measure
    selectWindow("Cytosol_Mask");
    getStatistics(area, meanCyto);
    if (meanCyto == 0) {
        print("  SKIP: cytosol mask empty after subtraction");
        return;
    }

    // Make cytosol foreground a re-applicable ROI (no pixel multiplication).
    roiManager("reset");
    selectWindow("Cytosol_Mask");
    run("Create Selection");
    roiManager("Add");
    cytoRoiId = 0;
    nCyto = getRawStatisticsCount();

    // ---- measure each marker ------------------------------
    measureAndWrite(m1, m1 + "_channel", cytoRoiId, nCyto, imgName, cellLine, tp, comboKey);
    measureAndWrite(m2, m2 + "_channel", cytoRoiId, nCyto, imgName, cellLine, tp, comboKey);

    // ---- save artefacts -----------------------------------
    if (SAVE_QC)    saveQcOverlay(imgName, m1);
    if (SAVE_MASKS) saveMasksTif(imgName);
}

// Try to parse filename. Returns:
//   newArray(tp, cellLine, m1, m2) on success
//   newArray()  (length 0)         on failure
//
// Failure conditions:
//   - fewer than 7 tokens
//   - markers not in the MARKERS list (sanity)
//   - timepoint not in TIMEPOINTS    (sanity)
function tryParseFilename(imgName) {
    tokens = split(imgName, "_");
    if (tokens.length < 7) return newArray();
    tp_       = tokens[0];
    cellLine_ = tokens[1];
    // tokens[2] is "Mock" (ignored, we already filtered for it)
    m1_       = tokens[3];
    m2_       = tokens[4];
    if (!inArray(m1_, MARKERS) || !inArray(m2_, MARKERS)) return newArray();
    if (!inArray(tp_, TIMEPOINTS))                       return newArray();
    return newArray(tp_, cellLine_, m1_, m2_);
}

// Per-image dialog. Defaults are pre-filled if known.
// Returns array [tp, cellLine, m1, m2].
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

// Dialog shown when filename parsing fails in filename mode.
// Returns "skip" or "dialog".
function askParseFailureAction(imgName) {
    Dialog.create("Filename parse failed");
    Dialog.addMessage("Cannot parse: " + imgName + "\n\nWhat do you want to do?");
    Dialog.addRadioButtonGroup("Action:",
        newArray("skip", "dialog"), 2, 1, "skip");
    Dialog.show();
    return Dialog.getRadioButton();
}


// ============== 4. CHANNEL HANDLING ========================

function isValidCombo(comboKey) {
    return inArray(comboKey, ANALYSE_COMBI);
}

// C1 = marker1, C2 = marker2, C3 = DAPI (fixed by acquisition).
function splitAndRenameChannels(title, m1, m2) {
    selectWindow(title);
    run("Split Channels");
    selectWindow("C1-" + title); rename(m1 + "_channel");
    selectWindow("C2-" + title); rename(m2 + "_channel");
    selectWindow("C3-" + title); rename("DAPI_channel");
}

// Priority: HA568 > HA488 > NS4B568 > dsRNA488. In Mock HA gives
// the best cytoplasmic autofluorescence coverage.
function pickCellMaskSource(m1, m2) {
    priority = newArray("HA568", "HA488", "NS4B568", "dsRNA488");
    for (i = 0; i < priority.length; i++) {
        cand = priority[i] + "_channel";
        if (isOpen(cand)) return cand;
    }
    return "";
}


// ============== 5. MASKS ===================================
// All masks are 8-bit binary 0/255 (ImageJ native).
// ===========================================================

function makeCellMask(srcTitle) {
    selectWindow(srcTitle);
    run("Duplicate...", "title=Cell_Mask");
    // Enhance Contrast only changes the DISPLAY LUT, NOT pixel
    // values. setAutoThreshold reads raw data, so this is
    // cosmetic — but it does help when you're watching the run.
    run("Enhance Contrast", "saturated=0.35");
    // "scaled" = sigma is in microns (physical units from image
    // metadata). On a 2850×2850 image at 100.57 µm, 1 µm ≈ 28 px,
    // i.e. heavy smoothing — exactly what we need to bridge the
    // gaps between dim cytoplasm and bright spots.
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_CELL + " scaled");
    setAutoThreshold(CELL_THR_METHOD + " dark");
    run("Convert to Mask");
    // Internal holes (dim spots inside cells) get filled.
    run("Fill Holes");
}

function makeCellMaskManual(srcTitle) {
    selectWindow(srcTitle);
    run("Duplicate...", "title=Cell_Mask");
    // Enhance Contrast only changes the DISPLAY LUT, NOT pixel
    // values. setAutoThreshold reads raw data, so this is
    // cosmetic — but it does help when you're watching the run.
    run("Enhance Contrast", "saturated=0.35");
    // "scaled" = sigma is in microns (physical units from image
    // metadata). On a 2850×2850 image at 100.57 µm, 1 µm ≈ 28 px,
    // i.e. heavy smoothing — exactly what we need to bridge the
    // gaps between dim cytoplasm and bright spots.
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_CELL + " scaled");
    run("Threshold...");
    setAutoThreshold(CELL_THR_METHOD + " dark");
    waitForUser("Check and adjust the threshold if necessary, then click OK.");
    run("Convert to Mask");
    // Internal holes (dim spots inside cells) get filled.
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

    // Morphological closing (Dilate × N, then Erode × N) before
    // Fill Holes: closes thin gaps in the nuclear boundary so
    // Fill Holes can actually enclose the holes. Without this,
    // small "bites" stay open because they reach the image edge
    // through a 1-pixel gap.
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

    // Morphological closing (Dilate × N, then Erode × N) before
    // Fill Holes: closes thin gaps in the nuclear boundary so
    // Fill Holes can actually enclose the holes. Without this,
    // small "bites" stay open because they reach the image edge
    // through a 1-pixel gap.
    if (NUC_CLOSE_ITER > 0) {
        run("Options...", "iterations=" + NUC_CLOSE_ITER + " count=1 black do=Nothing");
        run("Dilate");
        run("Fill Holes");
        run("Erode");
    } else {
        run("Fill Holes");
    }
}

// Cytosol = Cell − Nucleus. In 8-bit subtraction is clamped at 0,
// so this gives a clean binary: cell-not-nucleus = 255, else 0.
function makeCytosolMask() {
    imageCalculator("Subtract create", "Cell_Mask", "Nucleus_Mask");
    rename("Cytosol_Mask");
    setThreshold(1, 255);
    run("Convert to Mask");
}


// ============== 6. MEASUREMENT (TOP-X% via histogram) ======
// IMPORTANT — IJM quirk:
//   getHistogram(values, counts, 65536, 0, 65535) on a 16-bit
//   image does NOT populate `values[]` (it stays as scalar 0).
//   This is a memory optimisation since bin index == pixel
//   value for this setup. We therefore use `i` directly as
//   the value in all bin loops, and ignore `values`.
// ===========================================================

function measureAndWrite(marker, chTitle, cytoRoiId, nCyto,
                         imgName, cellLine, tp, comboKey) {
    if (!isOpen(chTitle)) {
        print("  WARN: channel window missing: " + chTitle);
        return;
    }
    selectWindow(chTitle);
    roiManager("Select", cytoRoiId);

    nBins = 65536;
    getHistogram(values, counts, nBins, 0, 65535);  // `values` unreliable!

    nTotal = 0;
    for (i = 0; i < nBins; i++) nTotal += counts[i];
    if (nTotal == 0) {
        print("  WARN: 0 pixels in cytosol selection for " + chTitle);
        return;
    }

    stats = computeTopStats(counts, nBins, nTotal, TOP_PCT);
    // stats = [thr, nTop, meanTop, medianTop, stdTop, p95, p99, p99_25, p99_5, p99_9, p99_95, p99_99, p99_995, p99_999]

    csvPath = MEASURE_DIR_RUN_ID + CELL_LINE + "_mock_" + tp + "_" + marker + "_in_" + comboKey + ".csv";
    appendCsvRow(csvPath, imgName, cellLine, tp, comboKey, marker, stats, nCyto);
}

// Pure number-crunching. No selection, no UI. Uses bin index = value.
function computeTopStats(counts, nBins, nTotal, topPct) {
    // (1) Find threshold V*: lowest pixel value still in the top X%.
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

    // (2) Mean and std over the top pool (value >= V*).
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

    // (3) Median of the top pool: cumulative count to nTop/2.
    half = nTop / 2.0;
    a = 0; median_top = threshold_value;
    for (i = threshold_value; i < nBins; i++) {
        if (counts[i] > 0) {
            a += counts[i];
            if (a >= half) { median_top = i; break; }
        }
    }

    // (4) Whole-cytosol percentiles (sanity / outlier detection).
    p95    = pctIndex(counts, nBins, nTotal, 95);
    p99    = pctIndex(counts, nBins, nTotal, 99);
    p99_25 = pctIndex(counts, nBins, nTotal, 99.25);
    p99_5  = pctIndex(counts, nBins, nTotal, 99.5);
    p99_9  = pctIndex(counts, nBins, nTotal, 99.9);
    p99_95  = pctIndex(counts, nBins, nTotal, 99.95);
    p99_99  = pctIndex(counts, nBins, nTotal, 99.99);
    p99_995  = pctIndex(counts, nBins, nTotal, 99.995);
    p99_999  = pctIndex(counts, nBins, nTotal, 99.999);

    return newArray(threshold_value, nTop, mean_top, median_top, std_top,
                    p95, p99, p99_25, p99_5, p99_9, p99_95, p99_99, p99_995, p99_999);
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
                      stats, nCyto) {
    line = imgName + "," + cellLine + "," + tp + "," + comboKey + "," + channel
        + "," + STAT_METHOD + "," + TOP_PCT
        + "," + MODE + "," + THR_MODE
        + "," + stats[0] + "," + stats[1] + "," + nCyto
        + "," + stats[2] + "," + stats[3] + "," + stats[4]
        + "," + stats[5] + "," + stats[6] + "," + stats[7] + "," + stats[8] + "," + stats[9] + "," + stats[10] + "," + stats[11] + "," + stats[12] + "," + stats[13]
        + "," + CELL_THR_METHOD + "," + NUC_THR_METHOD
        + "," + BLUR_SIGMA_CELL + "," + BLUR_SIGMA_NUC
        + "," + MACRO_VERSION + "," + RUN_ID;
    File.append(line, csvPath);
}

function saveMasksTif(imgName) {
    masksDir = MEASURE_DIR_RUN_ID + "masks" + File.separator;
    filename = RUN_ID + "_" + imgName;
    selectWindow("Cell_Mask");    saveAs("Tiff", masksDir + filename + "_cell.tif");
    selectWindow("Nucleus_Mask"); saveAs("Tiff", masksDir + filename + "_nuc.tif");
    selectWindow("Cytosol_Mask"); saveAs("Tiff", masksDir + filename + "_cyto.tif");
}

function saveQcOverlay(imgName, m1) {
    qcDir = MEASURE_DIR_RUN_ID  + "qc" + File.separator;
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
    close();   // close the flattened (now active after saveAs)
    if (isOpen("qc_tmp")) { selectWindow("qc_tmp"); close(); }
}


// ============== 8. CLEANUP & UTILS =========================

function cleanupBetweenImages() {
    while (nImages > 0) { selectImage(nImages); close(); }
    if (isOpen("ROI Manager")) roiManager("reset");
    if (isOpen("Results"))     run("Clear Results");
    if (selectionType() != -1) run("Select None");
}

// "a, b ,c" -> ["a","b","c"]  (trimmed)
function parseCsvString(s) {
    raw = split(s, ",");
    out = newArray(raw.length);
    for (i = 0; i < raw.length; i++) out[i] = trim(raw[i]);
    return out;
}

// ["a","b","c"] -> "a,b,c"  (for logging / dialog defaults)
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
