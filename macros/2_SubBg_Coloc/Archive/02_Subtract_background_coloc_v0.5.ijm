// ==========================================================
// Background Subtraction + Coloc Prep Pipeline  --  V0.5
// Author : Kolja Hildenbrand
// Date   : 2026-05-28
// Status : OLD (superseded by v0.14.3)
//
// Cell mask = V0.4-style: a single priority-picked channel auto-thresholded
// with Li on the RAW autofluorescence (the v0.3 bg-value-based + combined
// C1+C2 cell-mask experiments were reverted as harder to calibrate).
//
// Changes vs V0.4:
//  - NEW:     Particle-size filter on cytosol mask. Disconnected
//             regions smaller than MIN_PARTICLE_SIZE px are dropped
//             — same logic as Mock pipeline.
//  - NEW:     Cell-line-specific defaults via applyCellLineDefaults
//             (Vero settings differ from Huh7).
//  - NEW:     Dialog-mode startup config (like Mock V5) — markers,
//             combos, timepoints, cell line all overridable.
//  - NEW:     Auto coloc parses the plugin's Log output and writes
//             REAL values to the CSV (not just placeholders).
//  - FIX:     setOption("BlackBackground", true) pinned at MAIN
//             and per iteration — defends against the Convert to
//             Mask inversion bug happened in Mock V5.
//  - FIX:     parseFloat() at bg-value store + retrieve — guards
//             against the "Array.concat on empty array stores numbers
//             as strings" IJM quirk.
//
// Filename schema (binding):
//   timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif
//   C1 = marker1, C2 = marker2, C3 = DAPI (always)
//
// Output (inside INPUT_DIR):
//   <RUN_ID>_bgsub/
//     bgsub_<original>.tif                   multi-ch 16-bit
//     bgsub_<original>.jpg                   8-bit RGB + scale bar (opt.)
//     background_values_used.md              bg values entered
//     coloc_results_<RUN_ID>.csv             coloc rows (real values if auto)
//     macro_log_<RUN_ID>.txt                 full IJ Log
// ==========================================================


// ============== 1. CONFIG ==================================

// Runtime state
var INPUT_DIR;
var OUTPUT_DIR;
var RUN_ID;
var MODE;        // "automatic" or "manual"
var CELL_LINE = "Huh7";

// Domain lists — overridable in dialog mode
var MARKERS       = newArray("HA568", "HA488", "dsRNA488", "NS4B568");
var ANALYSE_COMBI = newArray("HA568_dsRNA488", "NS4B568_dsRNA488", "NS4B568_HA488");
var TIMEPOINTS    = newArray("12h", "24h");

// Background lookup table (filled by askBackgroundValues)
var BG_KEYS   = newArray();
var BG_VALUES = newArray();

// --- Cell mask building (V0.5.1: reverted to V0.4-style auto-threshold) ---
// We build the cell mask from ONE channel (priority HA568 > HA488 > NS4B568 >
// dsRNA488) via setAutoThreshold on the RAW (un-subtracted) signal — the
// autofluorescence carries enough shape information across the whole cell.
// The bg values are NOT used here (they're only for the actual subtraction).
var CELL_THR_METHOD = "Li";       // Li | Otsu | Triangle | Yen
var CELL_THR_FACTOR = 0.5;        // multiplier on auto threshold (<1 = permissive)

// Nucleus mask
var NUC_THR_METHOD  = "Otsu";     // Otsu | Triangle
var NUC_CLOSE_ITER  = 2;

// Smoothing (microns, applied via "scaled" Gaussian)
var BLUR_SIGMA_CELL = 1;
var BLUR_SIGMA_NUC  = 1;

// Particle-size filter on cytosol — drops small disconnected fragments
// outside real cells that survived thresholding by accident.
var MIN_PARTICLE_SIZE = 200;

// Generous dilation of the cell mask AFTER OR + Fill Holes.
// Pads outward so we keep membrane-proximal signal in the ROI.
var CELL_DILATE_ITER = 3;

// Coloc step. Three modes (same as V0.4):
//   "none"   = subtract bg, merge, save — no coloc CSV
//   "manual" = build cytosol ROI, pause per image for user
//   "auto"   = run Colocalisation Threshold automatically per image
var COLOC_MODE = "manual";

// JPG export (optional, same as V0.4)
var SAVE_JPG        = false;
var JPG_SCALEBAR_UM = 20;

// Output naming
var OUT_PREFIX = "bgsub_";

// Reproducibility
var MACRO_VERSION = "0.5.0";


// ============== MAIN =======================================

// CRITICAL: pin Black Background ON. Convert to Mask depends on this;
// without it the mask polarity can flip silently between iterations
// (we hit this exact bug in Mock V0.5).
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

// ------------------------------
// Two-step startup dialog:
//   1) Mode (filename/manual), Pipeline scope, Cell line, JPG, Dialog-config
//   2) (optional) override MARKERS / ANALYSE_COMBI / TIMEPOINTS
// ------------------------------
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
    Dialog.addMessage("Which cell line are you analysing?");
    Dialog.addRadioButtonGroup("Cell line:",
        newArray("Huh7", "VeroE6"), 2, 1, "Huh7");

    Dialog.addMessage("---");
    Dialog.addCheckbox("Also save 8-bit RGB JPG with " + JPG_SCALEBAR_UM + " um scale bar", SAVE_JPG);
    Dialog.addCheckbox("Use dialog config (override markers / combos / timepoints)", false);

    Dialog.show();

    // get* in same order as add*
    MODE          = Dialog.getRadioButton();
    pipelineScope = Dialog.getRadioButton();
    CELL_LINE     = Dialog.getRadioButton();
    SAVE_JPG      = Dialog.getCheckbox();
    dialogConfig  = Dialog.getCheckbox();

    if (pipelineScope == "subtract only")              COLOC_MODE = "none";
    else if (pipelineScope == "subtract + auto coloc") COLOC_MODE = "auto";
    else                                               COLOC_MODE = "manual";

    applyCellLineDefaults(CELL_LINE);

    if (dialogConfig) {
        Dialog.create("Pipeline setup (dialog config)");
        Dialog.addMessage("Edit the lists used by the macro. Comma-separated, whitespace trimmed.");
        Dialog.addString("Markers:",    arrToStr(MARKERS),       40);
        Dialog.addString("Combos:",     arrToStr(ANALYSE_COMBI), 40);
        Dialog.addString("Timepoints:", arrToStr(TIMEPOINTS),    20);
        Dialog.show();
        MARKERS       = parseCsvString(Dialog.getString());
        ANALYSE_COMBI = parseCsvString(Dialog.getString());
        TIMEPOINTS    = parseCsvString(Dialog.getString());
    }
}

// Cell-line-specific mask parameters. Called after CELL_LINE is known.
// Does NOT touch CELL_THR_METHOD or NUC_THR_METHOD — those are global
// strategy choices set in CONFIG. Per-line tuning is limited to the
// "strictness knobs" (blur, factor, particle size).
function applyCellLineDefaults(cellLine) {
    if (cellLine == "VeroE6") {
        // Vero ~4× smaller than Huh7, dimmer periphery
        BLUR_SIGMA_CELL   = 1.5;
        MIN_PARTICLE_SIZE = 100;
        CELL_THR_FACTOR   = 0.3;     // permissive — Li threshold scaled down
    } else {  // Huh7 default
        BLUR_SIGMA_CELL   = 1;
        MIN_PARTICLE_SIZE = 200;
        CELL_THR_FACTOR   = 0.5;
    }
    print("Cell-line tuning for " + cellLine + ":");
    print("  BLUR_SIGMA_CELL   = " + BLUR_SIGMA_CELL);
    print("  MIN_PARTICLE_SIZE = " + MIN_PARTICLE_SIZE);
    print("  CELL_THR_FACTOR   = " + CELL_THR_FACTOR);
}

// ------------------------------
// Single dialog asking for all bg values. Built dynamically from
// TIMEPOINTS x ANALYSE_COMBI; READ loop mirrors ADD loop for safe ordering.
// ------------------------------
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
    header = "image,cell_line,timepoint,combo,channel_1,channel_2,"
           + "Rtotal,m,b,Ch1_thresh,Ch2_thresh,"
           + "Rcoloc,R_below_thresh,M1,M2,tM1,tM2,"
           + "Ncoloc,perc_volume,perc_ch1_vol,perc_ch2_vol,"
           + "macro_version,run_id,status\n";
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
    if (COLOC_MODE == "manual") colocLabel = "MANUAL (pause per image)";
    if (COLOC_MODE == "auto")   colocLabel = "AUTO (plugin runs per image, log parsed)";
    print("coloc step        : " + colocLabel);
    print("cell thr method   : " + CELL_THR_METHOD + " (factor " + CELL_THR_FACTOR + ")");
    jpgLabel = "no";
    if (SAVE_JPG) jpgLabel = "yes (" + JPG_SCALEBAR_UM + " um scale bar)";
    print("save jpg          : " + jpgLabel);
    print("input dir         : " + INPUT_DIR);
    print("output dir        : " + OUTPUT_DIR);
    print("MOI files         : " + nFiles);
    print("---");
}


// ============== 3. PER-IMAGE PIPELINE ======================

function processOneImage(fname) {
    setOption("BlackBackground", true);   // paranoid re-pin per iteration

    open(INPUT_DIR + fname);
    title = getTitle();
    imgName = substring(title, 0, lastIndexOf(title, "."));

    // ---- resolve metadata --------------------------------
    cellLine = "";
    meta = resolveMetadata(imgName);
    if (meta.length == 0) { print("  -> skipped"); return; }
    tp = meta[0]; m1 = meta[1]; m2 = meta[2];
    tokens = split(imgName, "_");
    if (tokens.length >= 2) cellLine = tokens[1]; else cellLine = CELL_LINE;

    // ---- validate ----------------------------------------
    comboKey = m1 + "_" + m2;
    if (!isValidCombo(comboKey)) { print("  SKIP: combo not in ANALYSE_COMBI -> " + comboKey); return; }
    if (!inArray(tp, TIMEPOINTS)) { print("  SKIP: timepoint not in TIMEPOINTS -> " + tp); return; }
    print("  combo=" + comboKey + "  tp=" + tp + "  cellLine=" + cellLine);

    // ---- lookup bg values --------------------------------
    bg1 = getBgValue(m1, comboKey, tp);
    bg2 = getBgValue(m2, comboKey, tp);
    print("  bg(" + m1 + ") = " + bg1 + ", bg(" + m2 + ") = " + bg2);

    // ---- split channels ----------------------------------
    splitAndRenameChannels(title, m1, m2);

    // ---- build masks BEFORE subtraction (only if needed) -
    haveCytosolRoi = false;
    if (COLOC_MODE != "none") {
        // V0.5.1: reverted to V0.4-style single-channel auto-threshold.
        // Single channel via pickCellMaskSource (priority HA568 > HA488 > NS4B568 > dsRNA488),
        // setAutoThreshold(Li) on the RAW signal before bg subtraction.
        cellSrc = pickCellMaskSource(m1, m2);
        if (cellSrc == "") {
            print("  WARN: no usable cell-mask source channel -> no ROI");
        } else {
            print("  Cell-Mask source: " + cellSrc);
            makeCellMask(cellSrc);
            makeNucleusMask("DAPI_channel");
            makeCytosolMask();
            cleanCytosolByParticleSize();
            haveCytosolRoi = extractCytosolRoi();
            closeMaskWindows();
        }
    }

    // ---- NOW subtract background -------------------------
    subtractFromChannel(m1 + "_channel", bg1);
    subtractFromChannel(m2 + "_channel", bg2);

    // ---- colocalisation step (channels still split!) -----
    if (COLOC_MODE != "none") {
        doColocalisationStep(imgName, cellLine, tp, comboKey, m1, m2, haveCytosolRoi);
    }

    // ---- merge back as scrollable Color hyperstack -------
    mergeChannelsBack(m1, m2);
    Stack.setDisplayMode("color");
    Stack.setChannel(1);

    // ---- save TIF ----------------------------------------
    outPath = OUTPUT_DIR + OUT_PREFIX + imgName + ".tif";
    saveAs("Tiff", outPath);
    print("  saved tif -> " + outPath);

    // ---- optional JPG ------------------------------------
    if (SAVE_JPG) {
        jpgPath = OUTPUT_DIR + OUT_PREFIX + imgName + ".jpg";
        saveJpgWithScaleBar(jpgPath);
    }

    // Re-pin black background after merge/save side effects
    setOption("BlackBackground", true);
}

// resolveMetadata, tryParseFilename, askImageMetadata, askParseFailureAction
// are identical to V0.4 — using IJM-safe intermediate variable pattern.
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


// ============== 5. MASKS ===================================
// V0.5.1: reverted to V0.4-style single-channel auto-threshold.
// Cell mask from ONE channel (priority-picked) via setAutoThreshold on
// the RAW (un-subtracted) signal. Autofluorescence carries enough shape
// information; bg-value thresholding turned out to be overkill.
// ===========================================================

// ------------------------------
// Cell-mask source priority: HA568 > HA488 > NS4B568 > dsRNA488.
// HA gives best cytoplasmic coverage; dsRNA only as last resort.
// Returns the window title or "" if no candidate is open.
// ------------------------------
function pickCellMaskSource(m1, m2) {
    priority = newArray("HA568", "HA488", "NS4B568", "dsRNA488");
    for (i = 0; i < priority.length; i++) {
        cand = priority[i] + "_channel";
        if (isOpen(cand)) return cand;
    }
    return "";
}

// ------------------------------
// Build the cell mask from `srcTitle`. Uses setAutoThreshold(CELL_THR_METHOD)
// scaled by CELL_THR_FACTOR (lower = more permissive).
// Generous dilation at the end widens the cytoplasm so the downstream
// subtraction / coloc ROI doesn't crop real signal.
// ------------------------------
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

    if (CELL_DILATE_ITER > 0) {
        run("Options...", "iterations=" + CELL_DILATE_ITER + " count=1 black do=Nothing");
        run("Dilate");
    }
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

function makeCytosolMask() {
    imageCalculator("Subtract create", "Cell_Mask", "Nucleus_Mask");
    rename("Cytosol_Mask");
    setThreshold(1, 255);
    run("Convert to Mask");
}

// ------------------------------
// Drop disconnected regions < MIN_PARTICLE_SIZE px from the cytosol.
// Same pattern as Mock V0.4+: Analyze Particles with show=Masks creates
// a new "Mask of Cytosol_Mask" containing only the large enough regions.
// ------------------------------
function cleanCytosolByParticleSize() {
    nBefore = countMaskForeground("Cytosol_Mask");
    selectWindow("Cytosol_Mask");
    run("Analyze Particles...", "size=" + MIN_PARTICLE_SIZE + "-Infinity show=Masks");
    closeIfOpen("Cytosol_Mask");
    selectWindow("Mask of Cytosol_Mask"); rename("Cytosol_Mask");
    nAfter = countMaskForeground("Cytosol_Mask");
    print("  particle filter (>= " + MIN_PARTICLE_SIZE + " px): " + (nBefore - nAfter) + " px removed");
}

function extractCytosolRoi() {
    selectWindow("Cytosol_Mask");
    getStatistics(area, meanCyto);
    if (meanCyto == 0) { print("  WARN: cytosol mask empty"); return false; }
    run("Create Selection");
    roiManager("reset");
    roiManager("Add");
    roiManager("Select", 0);
    roiManager("Rename", "Cytosol");
    return true;
}

function closeMaskWindows() {
    closeIfOpen("Cell_Mask");
    closeIfOpen("Nucleus_Mask");
    closeIfOpen("Cytosol_Mask");
}

// 8-bit binary mask foreground count via histogram[255].
// IJM quirk: 5-arg getHistogram fails on 8-bit; use the 3-arg form.
function countMaskForeground(maskTitle) {
    selectWindow(maskTitle);
    getHistogram(values, counts, 256);
    return counts[255];
}


// ============== 6. COLOCALISATION ==========================

function doColocalisationStep(imgName, cellLine, tp, comboKey, m1, m2, haveRoi) {
    if (COLOC_MODE == "auto") {
        doColocAuto(imgName, cellLine, tp, comboKey, m1, m2, haveRoi);
    } else {
        doColocManual(imgName, cellLine, tp, comboKey, m1, m2, haveRoi);
    }
}

// ------------------------------
// Manual: 3-step user flow (preview, decision, plugin pause).
// Writes a PLACEHOLDER row (no Log parsing — user may use a different
// plugin and we can't predict its output format).
// ------------------------------
function doColocManual(imgName, cellLine, tp, comboKey, m1, m2, haveRoi) {
    activateChannelsAndRoi(haveRoi, m1);
    roiLineShort  = "(no ROI — full image)";
    useRoiSetting = "<None>";
    if (haveRoi) {
        roiLineShort  = "Cytosol ROI is already applied to " + m1 + "_channel";
        useRoiSetting = "Channel 1";
    }

    waitForUser("Preview  [" + imgName + "]",
        "Three split channels are open:\n"
      + "  C1 : " + m1 + "_channel   (bg subtracted)\n"
      + "  C2 : " + m2 + "_channel   (bg subtracted)\n"
      + "  C3 : DAPI_channel\n\n"
      + "Cytosol: " + roiLineShort + "\n\n"
      + "Click OK when ready to decide.\n(Cancel aborts the WHOLE macro run.)");

    Dialog.create("Decision  [" + imgName + "]");
    Dialog.addMessage("What now?");
    Dialog.addRadioButtonGroup("Action:",
        newArray("run Colocalisation Threshold", "skip coloc (image still saved)"),
        2, 1, "run Colocalisation Threshold");
    Dialog.show();
    choice = Dialog.getRadioButton();

    if (choice != "run Colocalisation Threshold") {
        print("  Coloc skipped by user choice");
        appendColocPlaceholderRow(imgName, cellLine, tp, comboKey, m1, m2, "skipped_by_user");
        return;
    }

    waitForUser("Coloc step  [" + imgName + "]",
        "Run Colocalisation Threshold now:\n"
      + "  Analyze > Colocalisation > Colocalisation Threshold\n\n"
      + "Settings:\n"
      + "  Channel 1            : " + m1 + "_channel\n"
      + "  Channel 2            : " + m2 + "_channel\n"
      + "  Use ROI              : " + useRoiSetting + "\n"
      + "  Channel Combination  : Red : Green\n"
      + "  Include zero-zero pixels in threshold calc : [x]\n\n"
      + "Click OK when done.\n(Cancel aborts the WHOLE macro run.)");

    appendColocPlaceholderRow(imgName, cellLine, tp, comboKey, m1, m2, "ok");
}

// ------------------------------
// Auto: run Colocalisation Threshold via macro, then PARSE the Log
// output to extract numeric values and write them to the CSV.
//
// Approach: snapshot log length BEFORE the plugin runs, take substring
// AFTER → that's the chunk written by this one plugin invocation.
// Then extractLogValue() parses "Label: value" pairs from it.
// ------------------------------
function doColocAuto(imgName, cellLine, tp, comboKey, m1, m2, haveRoi) {
    activateChannelsAndRoi(haveRoi, m1);

    args = "channel_1=" + m1 + "_channel "
         + "channel_2=" + m2 + "_channel "
         + "channel=[Red : Green] "
         + "include";
    if (haveRoi) args = args + " use=[Channel 1]";

    print("  auto coloc: run(\"Colocalization Threshold\", \"" + args + "\")");
    logLenBefore = lengthOf(getInfo("log"));
    run("Colocalization Threshold", args);
    newLog = substring(getInfo("log"), logLenBefore);

    appendColocRowFromLog(imgName, cellLine, tp, comboKey, m1, m2, newLog, "auto_ok");
}

// ------------------------------
// Bring m1 channel to front and apply cytosol ROI as active selection.
// The Colocalisation Threshold "Use ROI: Channel 1" setting reads
// whatever selection is active on the chosen image.
// ------------------------------
function activateChannelsAndRoi(haveRoi, m1) {
    if (!isOpen(m1 + "_channel")) return;
    selectWindow(m1 + "_channel");
    if (haveRoi && roiManager("count") > 0) {
        roiManager("Select", 0);
    }
}

// ------------------------------
// Parse a "Label: value" or "Label = value" pair out of free-form
// log text. Returns "" if the label is not found.
//
// The Colocalisation Threshold plugin uses labels like:
//   "Rtotal", "Slope", "Intercept", "Ch1 thresh", "Ch2 thresh",
//   "Rcoloc", "R<threshold", "M1", "M2", "tM1", "tM2", "Ncoloc",
//   "%Volume", "%Ch1 Vol", "%Ch2 Vol"
// Label names depend on plugin version — verify with Plugins>Macros>Record.
// ------------------------------
function extractLogValue(text, label) {
    idx = indexOf(text, label);
    if (idx < 0) return "";
    after = substring(text, idx + lengthOf(label));
    // Skip the colon / equals / spaces between label and value.
    j = 0;
    while (j < lengthOf(after)) {
        c = substring(after, j, j+1);
        if (c == ":" || c == "=" || c == " " || c == "\t") { j++; continue; }
        break;
    }
    after = substring(after, j);
    // Value ends at first whitespace / newline / opening paren.
    valEnd = lengthOf(after);
    for (i = 0; i < lengthOf(after); i++) {
        c = substring(after, i, i+1);
        if (c == " " || c == "\n" || c == "\t" || c == "(" || c == ",") {
            valEnd = i; break;
        }
    }
    return substring(after, 0, valEnd);
}

function appendColocRowFromLog(imgName, cellLine, tp, comboKey, m1, m2, logChunk, status) {
    rtotal    = extractLogValue(logChunk, "Rtotal");
    slope     = extractLogValue(logChunk, "Slope");
    if (slope == "") slope = extractLogValue(logChunk, "m");
    intercept = extractLogValue(logChunk, "Intercept");
    if (intercept == "") intercept = extractLogValue(logChunk, "b");
    ch1thresh = extractLogValue(logChunk, "Ch1 thresh");
    ch2thresh = extractLogValue(logChunk, "Ch2 thresh");
    rcoloc    = extractLogValue(logChunk, "Rcoloc");
    rbelow    = extractLogValue(logChunk, "R<threshold");
    m1coef    = extractLogValue(logChunk, "M1");
    m2coef    = extractLogValue(logChunk, "M2");
    tm1       = extractLogValue(logChunk, "tM1");
    tm2       = extractLogValue(logChunk, "tM2");
    ncoloc    = extractLogValue(logChunk, "Ncoloc");
    pctvol    = extractLogValue(logChunk, "%Volume");
    pctch1    = extractLogValue(logChunk, "%Ch1 Vol");
    pctch2    = extractLogValue(logChunk, "%Ch2 Vol");

    csvPath = OUTPUT_DIR + "coloc_results_" + RUN_ID + ".csv";
    line = imgName + "," + cellLine + "," + tp + "," + comboKey + ","
         + m1 + "," + m2 + ","
         + rtotal + "," + slope + "," + intercept + ","
         + ch1thresh + "," + ch2thresh + ","
         + rcoloc + "," + rbelow + ","
         + m1coef + "," + m2coef + "," + tm1 + "," + tm2 + ","
         + ncoloc + "," + pctvol + "," + pctch1 + "," + pctch2 + ","
         + MACRO_VERSION + "," + RUN_ID + "," + status;
    File.append(line, csvPath);
    print("  coloc CSV row appended (" + status + ", Rtotal=" + rtotal + ", M1=" + m1coef + ", M2=" + m2coef + ")");
}

function appendColocPlaceholderRow(imgName, cellLine, tp, comboKey, m1, m2, status) {
    csvPath = OUTPUT_DIR + "coloc_results_" + RUN_ID + ".csv";
    line = imgName + "," + cellLine + "," + tp + "," + comboKey + ","
         + m1 + "," + m2 + ","
         + ",,,,,"            // Rtotal, m, b, Ch1_thresh, Ch2_thresh
         + ",,,,,,"           // Rcoloc, R<thr, M1, M2, tM1, tM2
         + ",,,,"             // Ncoloc, %Vol, %Ch1Vol, %Ch2Vol
         + MACRO_VERSION + "," + RUN_ID + "," + status;
    File.append(line, csvPath);
    print("  coloc CSV row appended (" + status + ")");
}


// ============== 7. BACKGROUND LOOKUP TABLE =================

// NB: parseFloat at BOTH boundaries (store + retrieve).
// IJM quirk: Array.concat onto an initially-empty newArray() can
// store numeric values as STRINGS (the array has no type hint from
// its empty state). When you later do `value * factor`, the parser
// fails with "; expected" — it can't multiply a string. parseFloat
// at the storage site canonicalises everything to numeric; the
// parseFloat at the retrieval site is paranoid defense.
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


// ============== 8. OUTPUT / DOCUMENTATION ==================

function writeBgMarkdown(mdPath) {
    if (File.exists(mdPath)) File.delete(mdPath);
    File.append("# Background values used\n\n", mdPath);
    File.append("Generated by `Subtract_background_coloc_v" + MACRO_VERSION + ".ijm`\n\n", mdPath);
    File.append("- Run ID        : " + RUN_ID + "\n", mdPath);
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

// ------------------------------
// 8-bit RGB JPG with scale bar burned in (figure-quality output).
// Composite display mode → all channels overlay → Stack to RGB → flatten.
// ------------------------------
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


// ============== 9. CLEANUP & UTILS =========================

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
