// ==========================================================
// Background Subtraction + Coloc Prep Pipeline  --  V0.6
// Author : Kolja Hildenbrand
// Date   : 2026-05-28
// Status : OLD (superseded by v0.14.3)
//
// Changes vs V0.5:
//  - CHANGE:  Subtract background FIRST, THEN build the masks. On the
//             bg-subtracted image the histogram is a clean 0-peak + signal
//             tail, which thresholds more predictably than raw autofluorescence.
//             (This reverses V0.4's ROI-FIRST order; V0.7 reverts to ROI-FIRST
//             again for an unbiased boundary.)
//  - CHANGE:  Cell mask = OR of BOTH channels' signal masks
//             (makeChannelSignalMask per channel -> imageCalculator OR ->
//             fill holes), replacing V0.5's single priority channel. Captures
//             signal present in EITHER marker. Tuned for bg-subtracted data:
//             CELL_THR_METHOD = Triangle, CELL_THR_FACTOR = 1.0,
//             BLUR_SIGMA_CELL = 2, MIN_PARTICLE_SIZE = 20 (keep small puncta).
//  - CHANGE:  Coloc Log parser made ROBUST: tryLabels() tries several
//             alternative label spellings per value (the "Colocalization
//             Threshold" plugin's wording varies); when nothing parses it
//             dumps the raw Log chunk so the labels can be read and added.
//
// KNOWN BIAS: the V0.6 ROI is the full-image signal-OR cytosol — it measures
// colocalisation only where there IS signal, which INFLATES Pearson's R
// (methodologically circular). This is why V0.7 switches to a single,
// morphology-defined representative-cell ROI.
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

// --- Cell mask building (V0.5.2: bg-subtracted + combined C1∪C2) ---
// Strategy:
//   1. Background is subtracted FIRST (bg-subtracted image has most pixels
//      at 0 and signal as bright dots above 0 — a clean bimodal histogram).
//   2. Each marker channel gets its own threshold (Triangle on bimodal
//      → finds the inflection point, ideal for "0 vs signal" data).
//   3. The two per-channel masks are OR'd → cell mask captures areas with
//      signal in EITHER channel (the actual "real signal" regions).
// Triangle is the right method here because after bg subtraction the
// histogram is heavily skewed (huge 0-peak + tail). Li/Otsu can pick
// non-sensical values on such histograms.
var CELL_THR_METHOD = "Triangle"; // Triangle | Li | Otsu | Yen
var CELL_THR_FACTOR = 1.0;        // multiplier on auto threshold (<1 = permissive)

// Nucleus mask
var NUC_THR_METHOD  = "Otsu";     // Otsu | Triangle
var NUC_CLOSE_ITER  = 2;

// Smoothing (microns, applied via "scaled" Gaussian)
var BLUR_SIGMA_CELL = 2;
var BLUR_SIGMA_NUC  = 1;

// Particle-size filter on cytosol — drops only the smallest disconnected
// fragments. Set LOW (20 px) so we don't lose real punctate signal — virus
// replication sites can be very small after bg subtraction.
var MIN_PARTICLE_SIZE = 20;

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
var MACRO_VERSION = "0.6.0";


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
        // Vero ~4× smaller than Huh7, dimmer periphery → smaller smoothing.
        BLUR_SIGMA_CELL   = 0.5;
        MIN_PARTICLE_SIZE = 20;
        CELL_THR_FACTOR   = 1.0;     // Triangle on bg-subtracted = sensible
    } else {  // Huh7 default
        BLUR_SIGMA_CELL   = 0.5;
        MIN_PARTICLE_SIZE = 20;
        CELL_THR_FACTOR   = 1.0;
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

    // ---- SUBTRACT FIRST (V0.5.2 order change) ------------
    // The cell mask is built AFTER subtraction so the threshold sees a
    // clean "0 vs signal" bimodal histogram instead of bg-dominated noise.
    // This makes the auto-threshold (Triangle) pick a sensible cutoff
    // even when the signal is small punctate dots.
    subtractFromChannel(m1 + "_channel", bg1);
    subtractFromChannel(m2 + "_channel", bg2);

    // ---- build cell mask from BOTH bg-subtracted channels ----
    haveCytosolRoi = false;
    if (COLOC_MODE != "none") {
        makeCellMaskFromChannels(m1, m2);   // OR of m1_mask + m2_mask
        makeNucleusMask("DAPI_channel");    // DAPI is NOT bg-subtracted
        makeCytosolMask();
        cleanCytosolByParticleSize();
        haveCytosolRoi = extractCytosolRoi();
        closeMaskWindows();
    }

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
// V0.5.2: combined cell mask from BOTH bg-subtracted channels.
// Each marker channel is thresholded separately (Triangle on its
// bimodal "0 vs signal" histogram), then the per-channel masks are
// OR'd → cell mask = where EITHER marker has real signal.
// ===========================================================

// ------------------------------
// Build a binary mask of "above-threshold" pixels in one channel.
// Used twice (m1, m2), then OR'd by makeCellMaskFromChannels.
// Operates on the bg-SUBTRACTED channel — threshold is meaningful
// only because we know background is at 0 after subtraction.
// ------------------------------
function makeChannelSignalMask(srcTitle, maskName) {
    selectWindow(srcTitle);
    run("Duplicate...", "title=" + maskName);
    // Mild blur (0.5 µm default) connects punctate signal into small
    // blobs without smearing dots into bg. NOT scaled too aggressively.
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_CELL + " scaled");

    setAutoThreshold(CELL_THR_METHOD + " dark");
    getThreshold(lo, hi);
    loNew = lo * CELL_THR_FACTOR;
    if (loNew < 1) loNew = 1;   // never threshold at 0 (whole image)
    setThreshold(loNew, hi);
    print("  " + maskName + " thr (" + CELL_THR_METHOD + " ×" + CELL_THR_FACTOR + "): auto=" + lo + " -> " + loNew);

    run("Convert to Mask");
    run("Fill Holes");
}

// ------------------------------
// Combined cell mask = M1_signal_mask OR M2_signal_mask.
// Pixel is "cell" if EITHER channel shows real signal there.
// On bg-subtracted data this captures all "real biology" locations.
// ------------------------------
function makeCellMaskFromChannels(m1, m2) {
    makeChannelSignalMask(m1 + "_channel", "Cell_M1");
    makeChannelSignalMask(m2 + "_channel", "Cell_M2");
    imageCalculator("OR create", "Cell_M1", "Cell_M2");
    rename("Cell_Mask");
    closeIfOpen("Cell_M1");
    closeIfOpen("Cell_M2");

    // Fill holes on the combined mask (OR can leave new internal holes).
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

// Try several alternative labels for one value (Coloc Threshold plugin
// uses different strings across versions). Returns the first non-empty.
function tryLabels(text, labels) {
    for (i = 0; i < labels.length; i++) {
        v = extractLogValue(text, labels[i]);
        if (v != "") return v;
    }
    return "";
}

function appendColocRowFromLog(imgName, cellLine, tp, comboKey, m1, m2, logChunk, status) {
    // Each row tries multiple label variants — match the user's CSV columns:
    // Rtotal, m, b, Ch1_thresh, Ch2_thresh, Rcoloc, R_below_thresh,
    // M1, M2, tM1, tM2, Ncoloc, perc_volume, perc_ch1_vol, perc_ch2_vol
    rtotal    = tryLabels(logChunk, newArray("Rtotal", "Pearson's R total", "R total"));
    slope     = tryLabels(logChunk, newArray("Slope", "m"));
    intercept = tryLabels(logChunk, newArray("Intercept", "b"));
    ch1thresh = tryLabels(logChunk, newArray("Ch1 thresh", "Threshold image 1", "Image1 Min Threshold", "Image 1 Min Threshold"));
    ch2thresh = tryLabels(logChunk, newArray("Ch2 thresh", "Threshold image 2", "Image2 Min Threshold", "Image 2 Min Threshold"));
    rcoloc    = tryLabels(logChunk, newArray("Rcoloc", "Pearson's R coloc", "Pearson R coloc", "R coloc"));
    rbelow    = tryLabels(logChunk, newArray("R<threshold", "R<Threshold", "R below threshold"));
    m1coef    = tryLabels(logChunk, newArray("M1:", "M1 =", "Manders' M1", "Image1 Manders"));
    m2coef    = tryLabels(logChunk, newArray("M2:", "M2 =", "Manders' M2", "Image2 Manders"));
    tm1       = tryLabels(logChunk, newArray("tM1", "thresholded M1", "Image1 thresholded Manders"));
    tm2       = tryLabels(logChunk, newArray("tM2", "thresholded M2", "Image2 thresholded Manders"));
    ncoloc    = tryLabels(logChunk, newArray("Ncoloc", "Number of coloc pixels", "Number of pixels above thresholds", "N coloc"));
    pctvol    = tryLabels(logChunk, newArray("%Volume", "% volume above thresholds", "% Volume above thresholds"));
    pctch1    = tryLabels(logChunk, newArray("%Ch1 Vol", "% above threshold in Ch1", "% Ch1 Vol above thresholds"));
    pctch2    = tryLabels(logChunk, newArray("%Ch2 Vol", "% above threshold in Ch2", "% Ch2 Vol above thresholds"));

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
    print("  coloc CSV row appended (" + status + "): Rtotal=" + rtotal
        + " Rcoloc=" + rcoloc + " M1=" + m1coef + " M2=" + m2coef);

    // If most values came back empty, the labels probably don't match this
    // plugin version → print the raw log chunk so the user can see the
    // actual label strings and adjust tryLabels() accordingly.
    if (rtotal == "" && m1coef == "") {
        print("  WARN: no coloc labels parsed from log — raw chunk below:");
        print("  ---vvv---");
        print(logChunk);
        print("  ---^^^---");
    }
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
