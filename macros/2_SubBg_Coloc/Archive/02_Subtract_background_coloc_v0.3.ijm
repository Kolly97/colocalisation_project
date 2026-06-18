// ==========================================================
// Background Subtraction + Coloc Prep Pipeline  --  V0.3
// Author : Kolja Hildenbrand
// Date   : 2026-05-25
// Status : OLD (superseded by v0.14.3)
//
// Changes vs V0.2:
//  - Version-bump checkpoint only. The code is identical to V0.2 (the Purpose
//    below is unchanged); MACRO_VERSION was bumped 0.2.0 -> 0.3.0 to snapshot
//    this state between the V0.2 prototype and the V0.4 ROI/coloc-bridge work.
//  - (This file was briefly misnamed "..._v0.4.ijm"; renamed to v0.3 to match
//    its header/MACRO_VERSION. A real V0.4 source file is not retained.)
//
// Purpose: For each MOI image,
//   (1) build a generous cytosol ROI from the cell-mask source
//       channel BEFORE subtracting any background (so the
//       threshold sees raw signal),
//   (2) subtract per-channel background values from C1 and C2,
//   (3) re-merge the channels into a scrollable multi-channel
//       hyperstack (Color display, like the original CZI),
//   (4) save the result as 16-bit TIF,
//   (5) pause for the user to run Coloc 2 (or another coloc
//       plugin) on the two non-DAPI channels within the cytosol
//       ROI,
//   (6) append a placeholder row to a per-run colocalisation CSV
//       (the actual Coloc 2 -> CSV bridge will be wired up in V0.4
//       once we know the exact column set we want to keep).
//
// Filename schema (binding):
//   timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif
//   C1 = marker1, C2 = marker2, C3 = DAPI (always)
//
// Output (inside INPUT_DIR):
//   <RUN_ID>_bgsub/
//     bgsub_<original_filename>.tif       multi-ch stack, 16-bit
//     background_values_used.md           bg values entered
//     coloc_results_<RUN_ID>.csv          coloc rows (placeholders for now)
//     macro_log_<RUN_ID>.txt              full IJ Log
// ==========================================================


// ============== 1. CONFIG (globals via `var`) ==============

// Runtime state
var INPUT_DIR;
var OUTPUT_DIR;
var RUN_ID;
var MODE;        // "automatic" or "manual"

// Domain lists
var MARKERS       = newArray("HA568", "HA488", "dsRNA488", "NS4B568");
var ANALYSE_COMBI = newArray("HA568_dsRNA488", "NS4B568_dsRNA488", "NS4B568_HA488");
var TIMEPOINTS    = newArray("12h", "24h");

// Background lookup table (filled by askBackgroundValues)
var BG_KEYS   = newArray();
var BG_VALUES = newArray();

// Mask / ROI parameters
var CELL_THR_METHOD  = "Triangle"; // Li | Otsu | Triangle | Yen (Triangle = permissive)
var NUC_THR_METHOD   = "Otsu";     // Otsu | Triangle
var BLUR_SIGMA_CELL  = 1;          // micrometers (Gaussian "scaled")
var BLUR_SIGMA_NUC   = 1;
var NUC_CLOSE_ITER   = 2;          // morphological closing for nucleus holes
var CELL_DILATE_ITER = 3;          // extra dilation -> "generous" cytoplasm

// Cell-mask threshold factor: the auto-threshold's lower bound is
// MULTIPLIED by this before binarisation. 1.0 = use as-is.
// < 1.0 = MORE permissive (lower threshold = more pixels classified
// as cell). 0.5 means "use half the auto threshold". Useful when
// the auto method (Li/Otsu/Triangle) is too conservative and cuts
// off dim cell periphery.
var CELL_THR_FACTOR  = 0.5;

// Coloc step. Three modes:
//   "none"   = subtract bg, merge, save — no coloc, no masks
//   "manual" = build cytosol ROI, pause per image for the user
//              to run Colocalisation Threshold by hand
//   "auto"   = build cytosol ROI, then RUN Colocalisation Threshold
//              automatically (no pause). The only manual step in
//              the whole run is entering bg values at startup.
var COLOC_MODE = "manual";

// Optional JPG export: in addition to the 16-bit TIF, save an
// 8-bit composite RGB with a scale bar (for slides/figures).
var SAVE_JPG        = false;
var JPG_SCALEBAR_UM = 20;     // scale bar width in micrometers

// Output naming
var OUT_PREFIX = "bgsub_";

// Reproducibility
var MACRO_VERSION = "0.3.0";


// ============== MAIN =======================================
chooseInputDir();
RUN_ID = makeRunId();
askModeAndConfig();
askBackgroundValues();
buildOutputDir();
// Coloc CSV only exists when we actually do coloc — keeps the
// output folder clean in "subtract only" runs.
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

// ------------------------------
// Ask the user for the input folder. Cancel = hard exit.
// ------------------------------
function chooseInputDir() {
    INPUT_DIR = getDirectory("Choose folder with MOI .tif images");
    if (INPUT_DIR == "") exit("No folder selected.");
}

// ------------------------------
// Timestamp tag YYYYMMDD_HHMM. Filesystem-safe.
// ------------------------------
function makeRunId() {
    getDateAndTime(y, m, dw, d, h, mn, s, ms);
    return "" + y + IJ.pad(m+1, 2) + IJ.pad(d, 2)
        + "_" + IJ.pad(h, 2) + IJ.pad(mn, 2);
}

// ------------------------------
// Startup mode dialog — asks THREE independent questions:
//
//   1) Marker source:
//        automatic = parse marker/timepoint from filename
//                    (with dialog fallback on parse failure)
//        manual    = ask per image via dialog
//
//   2) Pipeline scope (controls COLOC_MODE):
//        subtract only       = no masks, no pause, no coloc CSV
//        subtract + manual   = build cytosol ROI, pause per image
//                              for user-driven Colocalisation Threshold
//        subtract + auto     = build cytosol ROI, run Colocalisation
//                              Threshold AUTOMATICALLY (no pause).
//                              Only manual step in the run = bg values.
//
//   3) Save JPG: extra 8-bit RGB JPG (with scale bar) alongside TIF.
//
// All choices update globals used downstream.
// ------------------------------
function askModeAndConfig() {
    Dialog.create("Pipeline mode");

    Dialog.addMessage("How should the macro know which marker is on which channel?");
    Dialog.addRadioButtonGroup("Marker source:",
        newArray("automatic", "manual"), 2, 1, "automatic");
    Dialog.addMessage("  automatic = parse from filename (dialog fallback on parse failure)");
    Dialog.addMessage("  manual    = ask per image via dialog");

    Dialog.addMessage("---");
    Dialog.addMessage("What should the macro do after subtracting the background?");
    Dialog.addRadioButtonGroup("Pipeline:",
        newArray("subtract only",
                 "subtract + manual coloc",
                 "subtract + auto coloc"),
        3, 1, "subtract + manual coloc");
    Dialog.addMessage("  subtract only           = no masks, no pause, just save 16-bit TIFs");
    Dialog.addMessage("  subtract + manual coloc = preview & decision per image, then user runs the plugin");
    Dialog.addMessage("  subtract + auto coloc   = runs Colocalisation Threshold automatically per image");

    Dialog.addMessage("---");
    Dialog.addCheckbox("Also save 8-bit RGB JPG with " + JPG_SCALEBAR_UM + " um scale bar (Lower Right)", SAVE_JPG);

    Dialog.show();

    // get* order must match add* order.
    MODE          = Dialog.getRadioButton();
    pipelineScope = Dialog.getRadioButton();
    SAVE_JPG      = Dialog.getCheckbox();

    // Map the pipeline-scope label to COLOC_MODE
    if (pipelineScope == "subtract only")           COLOC_MODE = "none";
    else if (pipelineScope == "subtract + auto coloc") COLOC_MODE = "auto";
    else                                            COLOC_MODE = "manual";
}

// ------------------------------
// Single dialog for ALL background values. Built dynamically from
// TIMEPOINTS x ANALYSE_COMBI so a new tp/combo in CONFIG is
// reflected automatically (no manual dialog editing). The READ
// loop is identical to the ADD loop — bullet-proof ordering.
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

// ------------------------------
// Output dir = INPUT_DIR/<RUN_ID>_bgsub/. Also writes the
// background_values_used.md record.
// ------------------------------
function buildOutputDir() {
    OUTPUT_DIR = INPUT_DIR + RUN_ID + "_bgsub" + File.separator;
    ensureDir(OUTPUT_DIR);
    writeBgMarkdown(OUTPUT_DIR + "background_values_used.md");
}

// ------------------------------
// Pre-create the coloc CSV with header (idempotent). Header
// mirrors Coloc 2's Costes output so once we wire up the actual
// Results-table extraction in V0.4, columns line up.
// ------------------------------
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

// ------------------------------
// Sorted list of MOI .tif files in INPUT_DIR (case-insensitive).
// ------------------------------
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

// ------------------------------
// Pretty Log header — also gets persisted to the saved log file.
// ------------------------------
function logHeader(nFiles) {
    print("=== Bg-sub + Coloc-prep V" + MACRO_VERSION + " ===");
    print("run_id     : " + RUN_ID);
    print("mode       : " + MODE);
    // Coloc label per mode (IJM has no ternary).
    colocLabel = "disabled (no masks, no coloc CSV)";
    if (COLOC_MODE == "manual") colocLabel = "MANUAL (pause per image)";
    if (COLOC_MODE == "auto")   colocLabel = "AUTO (Colocalisation Threshold runs per image)";
    print("coloc step : " + colocLabel);
    jpgLabel = "no";
    if (SAVE_JPG) jpgLabel = "yes (8-bit RGB JPG + " + JPG_SCALEBAR_UM + " um scale bar)";
    print("save jpg   : " + jpgLabel);
    print("input dir  : " + INPUT_DIR);
    print("output dir : " + OUTPUT_DIR);
    print("MOI files  : " + nFiles);
    print("---");
}


// ============== 3. PER-IMAGE PIPELINE ======================
// Order is deliberate:
//   1) open + resolve metadata + validate
//   2) lookup bg1, bg2  (but do NOT subtract yet)
//   3) split channels + rename windows by marker name
//   4) build masks on the RAW (un-subtracted) channels —
//      thresholds are more reliable here
//   5) extract cytosol -> Selection -> ROI Manager (name "Cytosol")
//   6) close mask windows (ROI persists in manager)
//   7) NOW subtract bg from C1 and C2 (DAPI untouched)
//   8) merge back into a Color-display hyperstack
//   9) save as <OUT_PREFIX><imgName>.tif
//  10) (optional) pause for user-driven coloc, then write a
//      placeholder CSV row
// ===========================================================

function processOneImage(fname) {
    open(INPUT_DIR + fname);
    title = getTitle();
    imgName = substring(title, 0, lastIndexOf(title, "."));

    // ---- resolve metadata --------------------------------
    cellLine = "";
    meta = resolveMetadata(imgName);
    if (meta.length == 0) { print("  -> skipped"); return; }
    tp = meta[0]; m1 = meta[1]; m2 = meta[2];
    // cellLine: from filename when parsable, else "unknown"
    tokens = split(imgName, "_");
    if (tokens.length >= 2) cellLine = tokens[1]; else cellLine = "unknown";

    // ---- validate ----------------------------------------
    comboKey = m1 + "_" + m2;
    if (!isValidCombo(comboKey)) {
        print("  SKIP: combo not in ANALYSE_COMBI -> " + comboKey);
        return;
    }
    if (!inArray(tp, TIMEPOINTS)) {
        print("  SKIP: timepoint not in TIMEPOINTS -> " + tp);
        return;
    }
    print("  combo=" + comboKey + "  tp=" + tp + "  cellLine=" + cellLine);

    // ---- lookup background values ------------------------
    bg1 = getBgValue(m1, comboKey, tp);
    bg2 = getBgValue(m2, comboKey, tp);
    print("  bg(" + m1 + ") = " + bg1);
    print("  bg(" + m2 + ") = " + bg2);

    // ---- split channels ----------------------------------
    splitAndRenameChannels(title, m1, m2);

    // ---- build masks BEFORE subtraction (only if needed) -
    // In COLOC_MODE = "none" we skip mask building entirely:
    // the cytosol ROI would never be used.
    haveCytosolRoi = false;
    if (COLOC_MODE != "none") {
        cellSrc = pickCellMaskSource(m1, m2);
        if (cellSrc == "") {
            print("  WARN: no usable cell-mask source channel -> no ROI");
        } else {
            print("  cell-mask source: " + cellSrc);
            makeCellMask(cellSrc);
            makeNucleusMask("DAPI_channel");
            makeCytosolMask();
            haveCytosolRoi = extractCytosolRoi();   // adds "Cytosol" to ROI Manager
            closeMaskWindows();
        }
    }

    // ---- NOW subtract background -------------------------
    // Channels stay SPLIT (m1_channel, m2_channel, DAPI_channel)
    // because the coloc plugin needs >= 2 separate image windows.
    // We merge + save AFTER the (optional) coloc step.
    subtractFromChannel(m1 + "_channel", bg1);
    subtractFromChannel(m2 + "_channel", bg2);

    // ---- colocalisation step (channels still split!) -----
    // doColocalisationStep dispatches to manual or auto based on
    // COLOC_MODE. When COLOC_MODE == "none" we don't call it at all.
    if (COLOC_MODE != "none") {
        status = doColocalisationStep(imgName, m1, m2, haveCytosolRoi);
        appendColocPlaceholderRow(imgName, cellLine, tp, comboKey, m1, m2, status);
    }

    // ---- merge back as scrollable Color hyperstack -------
    mergeChannelsBack(m1, m2);
    Stack.setDisplayMode("color");
    Stack.setChannel(1);   // start on C1

    // ---- save TIF ----------------------------------------
    outPath = OUTPUT_DIR + OUT_PREFIX + imgName + ".tif";
    saveAs("Tiff", outPath);
    print("  saved tif -> " + outPath);

    // ---- optional JPG with scale bar ---------------------
    if (SAVE_JPG) {
        jpgPath = OUTPUT_DIR + OUT_PREFIX + imgName + ".jpg";
        saveJpgWithScaleBar(jpgPath);
    }
}

// ------------------------------
// Wraps the filename-vs-dialog logic. Returns [tp, m1, m2] or
// empty array (= skip this image).
//
// IMPORTANT — same IJM quirk as the coloc dispatcher:
//   `return askImageMetadata(...)` (direct return of a user-function
//   call that returns an array) confuses IJM's type inferencer
//   and the array gets read back as `0` in the caller. Always go
//   through an intermediate variable.
// ------------------------------
function resolveMetadata(imgName) {
    parsed = tryParseFilename(imgName);

    // We collect the answer into `result` and return ONCE at the end.
    // This is the IJM-safe pattern for any function that returns
    // arrays produced by other user functions.
    result = newArray();   // default: skip

    if (MODE == "automatic") {
        if (parsed.length > 0) {
            result = parsed;
        } else {
            print("  Parse failed for: " + imgName);
            action = askParseFailureAction(imgName);
            if (action == "skip") {
                result = newArray();
            } else {
                result = askImageMetadata("(parse failed)", "", "", "");
            }
        }
    } else {
        // manual mode
        defTp = ""; defM1 = ""; defM2 = "";
        if (parsed.length > 0) { defTp = parsed[0]; defM1 = parsed[1]; defM2 = parsed[2]; }
        result = askImageMetadata(imgName, defTp, defM1, defM2);
    }

    return result;
}


// ============== 4. CHANNEL HANDLING ========================

// ------------------------------
// Split into per-channel windows and rename them by marker name.
// C1 -> <m1>_channel, C2 -> <m2>_channel, C3 -> DAPI_channel.
// Downstream code references markers by name, not by C-number.
// ------------------------------
function splitAndRenameChannels(title, m1, m2) {
    selectWindow(title);
    run("Split Channels");
    selectWindow("C1-" + title); rename(m1 + "_channel");
    selectWindow("C2-" + title); rename(m2 + "_channel");
    selectWindow("C3-" + title); rename("DAPI_channel");
}

// ------------------------------
// Subtract a scalar from a 16-bit channel. Subtract clamps at 0
// (no negatives), keeps 16-bit dtype.
// ------------------------------
function subtractFromChannel(chTitle, bgValue) {
    if (!isOpen(chTitle)) {
        print("  WARN: channel window missing: " + chTitle);
        return;
    }
    selectWindow(chTitle);
    // Defensive: a stray active selection would scope the subtract
    // to that selection only. Clear it so we always operate on the
    // full image.
    run("Select None");
    run("Subtract...", "value=" + bgValue);
}

// ------------------------------
// Merge the three renamed channels into a Composite hyperstack.
// `create` PRESERVES 16-bit data — without it you get 8-bit RGB
// and lose your dynamic range silently.
// Source windows are consumed (closed) by Merge Channels.
// Caller should then setDisplayMode("color") to get the
// scrollable per-channel view.
// ------------------------------
function mergeChannelsBack(m1, m2) {
    run("Merge Channels...",
        "c1=" + m1 + "_channel "
      + "c2=" + m2 + "_channel "
      + "c3=DAPI_channel create");
}


// ============== 5. MASKS ===================================
// Same family of functions as the Mock measurement pipeline V0.2 —
// keep them in sync if you tune one, tune the other.
// All masks 8-bit binary 0/255 (ImageJ standard).
// ===========================================================

// ------------------------------
// Cell-mask source priority: HA568 > HA488 > NS4B568 > dsRNA488.
// HA gives best cytoplasmic coverage even in MOI; dsRNA only as
// last resort (very sparse signal makes a bad cell outline).
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
// Build the cell mask from `srcTitle`. Generous on purpose:
// extra dilation steps at the end widen the cytoplasm so the
// downstream subtraction / coloc ROI doesn't crop real signal.
// ------------------------------
function makeCellMask(srcTitle) {
    selectWindow(srcTitle);
    run("Duplicate...", "title=Cell_Mask");
    run("Enhance Contrast", "saturated=0.35");
    // "scaled" => sigma in microns. On 2850 px / 100 um images,
    // sigma=1 um ~= 28 px = strong smoothing => connected cytoplasm.
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_CELL + " scaled");

    // Auto threshold ONCE to get the algorithm's pick, then scale
    // the lower bound by CELL_THR_FACTOR so we can be more permissive
    // (smaller factor = more cell area kept) without changing the
    // method. This is the cleanest way to widen the mask: the method
    // still adapts per image, we just shift the cutoff.
    setAutoThreshold(CELL_THR_METHOD + " dark");
    getThreshold(lo, hi);
    loNew = lo * CELL_THR_FACTOR;
    if (loNew < 0) loNew = 0;
    setThreshold(loNew, hi);
    print("  cell thr (" + CELL_THR_METHOD + "): auto=" + lo
        + ", scaled by " + CELL_THR_FACTOR + " -> " + loNew);

    run("Convert to Mask");
    run("Fill Holes");
    // Generous dilation: pad the mask outward so we don't lose
    // membrane-proximal signal in the ROI.
    if (CELL_DILATE_ITER > 0) {
        run("Options...", "iterations=" + CELL_DILATE_ITER + " count=1 black do=Nothing");
        run("Dilate");
    }
}

// ------------------------------
// Nucleus mask from DAPI. Morphological closing (Dilate->Fill->
// Erode) before Fill Holes handles boundary gaps that would
// otherwise leave little "bites" in the nucleus shape.
// ------------------------------
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

// ------------------------------
// Cytosol = Cell - Nucleus. 8-bit subtraction clamps at 0, so
// 255-255=0, 255-0=255, 0-anything=0 gives a clean binary.
// Defensive re-binarise ensures 0/255 even if upstream left a
// stray intermediate value.
// ------------------------------
function makeCytosolMask() {
    imageCalculator("Subtract create", "Cell_Mask", "Nucleus_Mask");
    rename("Cytosol_Mask");
    setThreshold(1, 255);
    run("Convert to Mask");
}

// ------------------------------
// Turn the Cytosol_Mask's foreground into an ROI and add to the
// ROI Manager as "Cytosol". Returns true on success, false if
// the cytosol mask is empty (no ROI added).
// ------------------------------
function extractCytosolRoi() {
    selectWindow("Cytosol_Mask");
    getStatistics(area, meanCyto);
    if (meanCyto == 0) {
        print("  WARN: cytosol mask empty - no ROI created");
        return false;
    }
    run("Create Selection");
    roiManager("reset");
    roiManager("Add");
    roiManager("Select", 0);
    roiManager("Rename", "Cytosol");
    return true;
}

// ------------------------------
// Close mask windows (Cell, Nucleus, Cytosol) — we already
// captured the ROI; raw masks are no longer needed and would
// clutter the workspace for the coloc step.
// ------------------------------
function closeMaskWindows() {
    closeIfOpen("Cell_Mask");
    closeIfOpen("Nucleus_Mask");
    closeIfOpen("Cytosol_Mask");
}


// ============== 6. METADATA HELPERS ========================

// ------------------------------
// Strict positional parse. Returns [tp, m1, m2] or [].
// Validates each token against MARKERS / TIMEPOINTS so typos in
// filenames are caught upfront instead of poisoning later steps.
// ------------------------------
function tryParseFilename(imgName) {
    tokens = split(imgName, "_");
    if (tokens.length < 7) return newArray();
    tp_ = tokens[0];
    // tokens[1] = cellLine, tokens[2] = condition (MOI1 / MOI5)
    m1_ = tokens[3];
    m2_ = tokens[4];
    if (!inArray(m1_, MARKERS) || !inArray(m2_, MARKERS)) return newArray();
    if (!inArray(tp_, TIMEPOINTS))                        return newArray();
    return newArray(tp_, m1_, m2_);
}

// ------------------------------
// Per-image dialog with safe defaults. Returns [tp, m1, m2].
// ------------------------------
function askImageMetadata(imgLabel, defTp, defM1, defM2) {
    if (defTp == "" || !inArray(defTp, TIMEPOINTS)) defTp = TIMEPOINTS[0];
    if (defM1 == "" || !inArray(defM1, MARKERS))    defM1 = MARKERS[0];
    if (defM2 == "" || !inArray(defM2, MARKERS))    defM2 = MARKERS[1 % MARKERS.length];

    Dialog.create("Image metadata");
    Dialog.addMessage("Image: " + imgLabel + "\nC3 is always DAPI.");
    Dialog.addRadioButtonGroup("Timepoint:",    TIMEPOINTS, 1, TIMEPOINTS.length, defTp);
    Dialog.addRadioButtonGroup("Channel 1 (C1):", MARKERS,  1, MARKERS.length,    defM1);
    Dialog.addRadioButtonGroup("Channel 2 (C2):", MARKERS,  1, MARKERS.length,    defM2);
    Dialog.show();

    tp_ = Dialog.getRadioButton();
    m1_ = Dialog.getRadioButton();
    m2_ = Dialog.getRadioButton();

    if (m1_ == m2_) {
        showMessage("Error", "C1 and C2 must be different markers.");
        exit("Aborted: identical markers for C1 and C2.");
    }
    return newArray(tp_, m1_, m2_);
}

// ------------------------------
// Shown in automatic mode when filename parsing fails.
// Returns "skip" or "manual".
// ------------------------------
function askParseFailureAction(imgName) {
    Dialog.create("Filename parse failed");
    Dialog.addMessage("Cannot parse: " + imgName + "\n\nWhat do you want to do?");
    Dialog.addRadioButtonGroup("Action:",
        newArray("skip", "manual"), 2, 1, "skip");
    Dialog.show();
    return Dialog.getRadioButton();
}

// ------------------------------
// Whitelist check for combo.
// NB: intermediate variable on purpose — see the IJM quirk note
// in doColocalisationStep / resolveMetadata. Never `return userFunc(...)`.
// ------------------------------
function isValidCombo(comboKey) {
    result = inArray(comboKey, ANALYSE_COMBI);
    return result;
}


// ============== 7. BACKGROUND LOOKUP TABLE =================
// Parallel-arrays "associative table" — see V0.2 commentary.
// ===========================================================

function setBgValue(marker, combo, tp, value) {
    key = marker + "_in_" + combo + "_" + tp;
    BG_KEYS   = Array.concat(BG_KEYS,   key);
    BG_VALUES = Array.concat(BG_VALUES, value);
}

function getBgValue(marker, combo, tp) {
    key = marker + "_in_" + combo + "_" + tp;
    for (i = 0; i < BG_KEYS.length; i++)
        if (BG_KEYS[i] == key) return BG_VALUES[i];
    exit("Background value missing for key: " + key
       + "\nDid you forget to enter it in the startup dialog?");
}


// ============== 8. COLOCALISATION (semi-automatic) =========
// Coloc 2 CAN be scripted (`run("Coloc 2", "...")`) but Costes'
// iterative threshold is slow and non-deterministic, so V0.3
// keeps the user in the loop. The macro:
//   - selects the merged image and the Cytosol ROI
//   - shows a dialog with the exact recommended Coloc 2 settings
//   - waits for the user to finish
//   - writes a CSV row (PLACEHOLDER until we wire up the actual
//     Results table extraction in V0.4)
// ===========================================================

// ------------------------------
// Drive the coloc pause for one image. Activates the merged
// hyperstack and the cytosol ROI so the user can immediately
// run Coloc 2 with the right inputs.
// ------------------------------
// ------------------------------
// Dispatcher: pick manual vs auto coloc based on COLOC_MODE.
// Returns the status string to put in the CSV row.
//
// IMPORTANT — IJM quirk:
//   `return userFunction(...)` confuses the IJM type inferencer
//   and produces a bogus "Numeric return value expected" error.
//   Always go through an intermediate variable when the value
//   returned IS another user-function call. The single-return
//   pattern below is the IJM-safe form.
// ------------------------------
function doColocalisationStep(imgName, m1, m2, haveRoi) {
    if (COLOC_MODE == "auto") {
        result = doColocAuto(imgName, m1, m2, haveRoi);
    } else {
        result = doColocManual(imgName, m1, m2, haveRoi);
    }
    return result;
}

// ------------------------------
// MANUAL coloc: three-step user flow per image:
//   1) PREVIEW (waitForUser, non-modal) — user inspects channels
//   2) DECISION (Dialog, modal) — run plugin or skip?
//   3) COLOC PAUSE (waitForUser, non-modal) — user runs plugin
// Returns "ok" or "skipped_by_user".
// ------------------------------
function doColocManual(imgName, m1, m2, haveRoi) {
    activateChannelsAndRoi(haveRoi, m1);

    roiLineShort  = "(no ROI — full image)";
    useRoiSetting = "<None>";
    if (haveRoi) {
        roiLineShort  = "Cytosol ROI is already applied to " + m1 + "_channel";
        useRoiSetting = "Channel 1";
    }

    // ----- 1. PREVIEW -----
    waitForUser("Preview  [" + imgName + "]",
        "Three split channels are open:\n"
      + "  C1 : " + m1 + "_channel   (bg subtracted)\n"
      + "  C2 : " + m2 + "_channel   (bg subtracted)\n"
      + "  C3 : DAPI_channel\n\n"
      + "Cytosol: " + roiLineShort + "\n\n"
      + "Click OK when you're ready to decide what to do.\n"
      + "(Cancel aborts the WHOLE macro run.)");

    // ----- 2. DECISION -----
    Dialog.create("Decision  [" + imgName + "]");
    Dialog.addMessage("What now?");
    Dialog.addRadioButtonGroup("Action:",
        newArray("run Colocalisation Threshold", "skip coloc (image will still be saved)"),
        2, 1, "run Colocalisation Threshold");
    Dialog.show();
    choice = Dialog.getRadioButton();

    if (choice != "run Colocalisation Threshold") {
        print("  Coloc skipped by user choice");
        return "skipped_by_user";
    }

    // ----- 3. COLOC PAUSE -----
    waitForUser("Coloc step  [" + imgName + "]",
        "Run Colocalisation Threshold now:\n"
      + "  Analyze > Colocalisation > Colocalisation Threshold\n\n"
      + "Settings:\n"
      + "  Channel 1            : " + m1 + "_channel\n"
      + "  Channel 2            : " + m2 + "_channel\n"
      + "  Use ROI              : " + useRoiSetting + "\n"
      + "  Channel Combination  : Red : Green\n"
      + "  Show Scatter plot                          : [x]\n"
      + "  Include zero-zero pixels in threshold calc : [x]\n"
      + "  (other checkboxes: leave unchecked)\n\n"
      + "Click OK when the analysis is complete.\n"
      + "(Cancel aborts the WHOLE macro run.)");

    return "ok";
}

// ------------------------------
// AUTO coloc: run Colocalisation Threshold programmatically.
// No pause — only the bg-value dialog at startup is manual.
//
// Args format derived DIRECTLY from Plugins > Macros > Record...:
//   run("Colocalization Threshold",
//       "channel_1=<title1> channel_2=<title2> use=[Channel 1] channel=[Red : Green] include");
//
// Critical Fiji inconsistencies to remember:
//   - Plugin name in run() is "Colocali**Z**ation" (American),
//     even though the Fiji menu shows "Colocali**S**ation"
//     (British). Always match the recorder, not the menu.
//   - The `use` value is "[Channel 1]" — square brackets are
//     required because of the space inside, and the literal text
//     is "Channel 1" with a space (NOT "Channel_1").
//   - Channel titles like "HA568_channel" need no brackets
//     because they contain no spaces.
//   - `include` = include zero-zero pixels in threshold calc
//     (matches supervisor's workflow).
//   - Omit `use=...` entirely when we have no ROI -> full image.
//
// Returns "auto_ok" — the actual numeric values land in the Log
// window. CSV parsing of those values is a V0.5+ feature.
// ------------------------------
function doColocAuto(imgName, m1, m2, haveRoi) {
    activateChannelsAndRoi(haveRoi, m1);

    args = "channel_1=" + m1 + "_channel "
         + "channel_2=" + m2 + "_channel "
         + "channel=[Red : Green] "
         + "include";
    if (haveRoi) args = args + " use=[Channel 1]";

    print("  auto coloc: run(\"Colocalization Threshold\", \"" + args + "\")");
    run("Colocalization Threshold", args);
    return "auto_ok";
}

// ------------------------------
// Prepare the workspace for the Colocalisation Threshold plugin:
//   - bring the first marker channel (m1) to the front
//   - apply the cytosol ROI to it as the ACTIVE SELECTION
// Why: Colocalisation Threshold's "Use ROI" dropdown only offers
// "Channel 1" or "Channel 2" (NOT ROI Manager). It then uses
// whatever selection is currently active on that channel. By
// pre-loading our cytosol ROI on m1_channel, the user just picks
// "Use ROI: Channel 1" and we're done — no manual ROI drawing.
// ------------------------------
function activateChannelsAndRoi(haveRoi, m1) {
    if (!isOpen(m1 + "_channel")) return;
    selectWindow(m1 + "_channel");
    if (haveRoi && roiManager("count") > 0) {
        roiManager("Select", 0);   // applies ROI 0 to the active window
    }
}

// ------------------------------
// Append a placeholder CSV row for this image. Real Coloc 2
// values will be filled in V0.4 by reading from the Results
// table (getResult("Rtotal", n-1) etc.) — for now, empty fields
// + a `status` column (ok | skipped | error) so we can later
// distinguish missing data from zero values.
// ------------------------------
function appendColocPlaceholderRow(imgName, cellLine, tp, comboKey, m1, m2, status) {
    csvPath = OUTPUT_DIR + "coloc_results_" + RUN_ID + ".csv";
    line = imgName + "," + cellLine + "," + tp + "," + comboKey + ","
         + m1 + "," + m2 + ","
         + ",,,,,"           // Rtotal, m, b, Ch1_thresh, Ch2_thresh
         + ",,,,,,"          // Rcoloc, R<thr, M1, M2, tM1, tM2
         + ",,,,"            // Ncoloc, %Vol, %Ch1Vol, %Ch2Vol
         + MACRO_VERSION + "," + RUN_ID + "," + status;
    File.append(line, csvPath);
    print("  coloc CSV row appended (" + status + ")");
}


// ------------------------------
// Save an 8-bit RGB JPG of the currently active composite, with a
// scale bar burned in. The composite must be the active window.
//
// Steps:
//   1) Switch display to "composite" (all channels overlayed in
//      their LUT colors). Without this, "Stack to RGB" would only
//      flatten the currently displayed single channel.
//   2) Stack to RGB → creates a new 8-bit RGB window where every
//      channel's contribution is flattened in.
//   3) Add a Scale Bar OVERLAY (non-destructive) sized in microns.
//      The image must have its pixel-to-micron calibration set
//      (CZIs from Bio-Formats always do).
//   4) Flatten the overlay → it becomes part of the pixel data so
//      it gets saved in the JPG.
//   5) Set JPEG quality and save.
//
// JPEG quality is set to 90 (a good compromise for figure-quality
// output without huge files). It's a global Fiji setting; that's
// fine since this macro is the only thing changing it.
// ------------------------------
function saveJpgWithScaleBar(jpgPath) {
    // Composite display so all 3 channels are visible.
    Stack.setDisplayMode("composite");

    // Create an 8-bit RGB flatten of the composite.
    run("Stack to RGB");
    // Active window is now "<title> (RGB)".

    // Add scale bar as overlay in the bottom-right corner.
    // width=N is in calibrated units (microns since the image is calibrated).
    run("Scale Bar...",
        "width=" + JPG_SCALEBAR_UM + " "
      + "height=10 thickness=20 font=100 "
      + "color=White background=None "
      + "location=[Lower Right] overlay");

    // Burn the overlay into the pixel data.
    run("Flatten");

    // Save. Set quality first so the JPG isn't ugly.
    run("Options...", "jpeg=90");
    saveAs("Jpeg", jpgPath);
    print("  saved jpg -> " + jpgPath);
}


// ============== 9. OUTPUT / DOCUMENTATION ==================

// ------------------------------
// Human-readable markdown record of all bg values entered.
// Generated dynamically from CONFIG -> guaranteed in sync with
// what was actually subtracted.
// ------------------------------
function writeBgMarkdown(mdPath) {
    if (File.exists(mdPath)) File.delete(mdPath);

    File.append("# Background values used\n\n", mdPath);
    File.append("Generated by `01_Subtract_background_v" + MACRO_VERSION + ".ijm`\n\n", mdPath);
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
// Persist the IJ Log to a text file. Call LAST so it captures
// everything (including coloc step messages).
// ------------------------------
function saveLogToFile() {
    if (!isOpen("Log")) return;
    logPath = OUTPUT_DIR + "macro_log_" + RUN_ID + ".txt";
    selectWindow("Log");
    saveAs("Text", logPath);
    print("Saved log: " + logPath);
}


// ============== 10. CLEANUP & UTILS ========================

// ------------------------------
// Close all open images, reset ROI/Results, drop selection.
// Do NOT close the Log — saved at end of run.
// ------------------------------
function cleanupBetweenImages() {
    while (nImages > 0) { selectImage(nImages); close(); }
    if (isOpen("ROI Manager")) roiManager("reset");
    if (isOpen("Results"))     run("Clear Results");
}

// ------------------------------
// File.makeDirectory is NOT recursive — parents must exist first.
// ------------------------------
function ensureDir(path) {
    if (!File.exists(path)) File.makeDirectory(path);
}

// ------------------------------
// Close a window by title if it's open. No-op otherwise.
// ------------------------------
function closeIfOpen(title) {
    if (isOpen(title)) { selectWindow(title); close(); }
}

// ------------------------------
// True if `val` is in `arr` (string equality).
// ------------------------------
function inArray(val, arr) {
    for (i = 0; i < arr.length; i++) if (arr[i] == val) return true;
    return false;
}
