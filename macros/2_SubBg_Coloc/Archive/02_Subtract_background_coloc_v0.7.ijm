// ==========================================================
// Background Subtraction + Coloc Prep Pipeline  --  V0.7
// Author : Kolja Hildenbrand
// Date   : 2026-06-01
// Status : OLD (superseded by v0.14.3)
//
// Changes vs V0.6:
//  - NEW:     Central-cell ROI for the colocalisation step. Instead of
//             the V0.6 signal-OR full-image cytosol mask (which is BIASED
//             — it only measures coloc where there is signal, inflating
//             Pearson R), the ROI is now built around ONE representative
//             cell in the middle of the image.
//             Two strategies, chosen at startup (ROI_MODE):
//               * auto_central : pick the nucleus closest to image centre
//                                whose dilated territory does NOT touch
//                                the image border (= fully-imaged cell),
//                                build cytosol = territory − nucleus.
//					Problem: Sometimes the middle cell is not infected.
//
//               * manual_draw  : pause per image, user draws a freehand
//                                ROI around the central cell (à la Thomas).
//  - CHANGE:  ROI is built on the RAW (un-subtracted) channels / DAPI,
//             BEFORE background subtraction — the V0.4 ordering. This
//             keeps the ROI signal-INDEPENDENT (unbiased) and, in manual
//             mode, leaves the whole cell visible via autofluorescence
//             (a bg-subtracted image is mostly zeros + dots, hard to trace).
//             Per-image order: split → build ROI → subtract bg → coloc
//             → merge → save.
//  - NEW:     Per-image coloc QC JPG (saveColocQc): composite of all three
//             channels with the ACTUAL ROI burnt in as an outline, so the
//             chosen / drawn cell is visually verifiable.
//  - CHANGE:  Domain-list dialog (markers / combos / timepoints) is now
//             ALWAYS shown at startup, pre-filled with the workflow
//             standards (Feature A, like Mock V0.6). The "Use dialog
//             config" checkbox is gone.
//  - NEW:     CSV provenance columns roi_mode, nuc_dilate_um.
//  - DROP:    V0.6 signal-OR cell mask (makeChannelSignalMask /
//             makeCellMaskFromChannels) and its particle filter — replaced
//             by the central-cell logic.
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
//     coloc_results_<RUN_ID>.csv             coloc rows (real values if auto)
//     macro_log_<RUN_ID>.txt                 full IJ Log
// ==========================================================


// ============== 1. CONFIG ==================================

// Runtime state
var INPUT_DIR;
var OUTPUT_DIR;
var RUN_ID;
var MODE;        // "automatic" or "manual" (per-image metadata source)
var CELL_LINE = "Huh7";

// Domain lists — always editable at startup (pre-filled with these)
var MARKERS       = newArray("HA568", "HA488", "dsRNA488", "NS4B568");
var ANALYSE_COMBI = newArray("HA568_dsRNA488", "NS4B568_dsRNA488", "NS4B568_HA488");
var TIMEPOINTS    = newArray("12h", "24h");

// Background lookup table (filled by askBackgroundValues)
var BG_KEYS   = newArray();
var BG_VALUES = newArray();

// --- ROI strategy for the colocalisation step ---
// "auto_central" : nucleus nearest image centre with a fully-imaged
//                  (non-border-touching) territory → cytosol = territory − nucleus.
// "manual_draw"  : user draws a freehand ROI per image (Thomas workflow).
// Chosen at startup; default below.
var ROI_MODE         = "auto_central";

// DAPI-dilation distance: how far the nucleus is grown to approximate the
// cell territory. In CALIBRATED micrometers (images are calibrated), so
// run("Enlarge...") handles the px conversion — no manual maths, no
// hundreds of Dilate iterations.
var NUC_DILATE_UM    = 8;

// A territory whose bounding box comes within BORDER_MARGIN_PX of any image
// edge is treated as a partially-imaged cell and skipped in auto_central.
var BORDER_MARGIN_PX = 0;

// Minimum nucleus area (pixels) in Analyze Particles — drops DAPI debris /
// tiny speckles so they are not mistaken for a nucleus.
var NUC_MIN_SIZE     = 1000;

// Nucleus mask
var NUC_THR_METHOD  = "Otsu";     // Otsu | Triangle
var NUC_CLOSE_ITER  = 2;
var BLUR_SIGMA_NUC  = 1;          // microns (scaled Gaussian)

// Coloc step. Three modes:
//   "none"   = subtract bg, merge, save — no coloc CSV
//   "manual" = build ROI, pause per image for the user to run the plugin
//   "auto"   = run Colocalisation Threshold automatically per image
var COLOC_MODE = "manual";

// QC / figure exports
var SAVE_COLOC_QC   = true;       // per-image composite + ROI outline (JPG)
var SAVE_JPG        = false;      // bg-subtracted composite + scale bar (JPG)
var JPG_SCALEBAR_UM = 20;

// Output naming
var OUT_PREFIX = "bgsub_";

// Reproducibility
var MACRO_VERSION = "0.7.0";


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
//   1) Marker source, Pipeline scope, ROI strategy, Cell line, JPG
//   2) ALWAYS: edit MARKERS / ANALYSE_COMBI / TIMEPOINTS (Feature A)
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
    Dialog.addMessage("How should the colocalisation ROI be built?");
    Dialog.addRadioButtonGroup("ROI strategy:",
        newArray("automatic central cell", "manual draw"), 2, 1, "automatic central cell");

    Dialog.addMessage("---");
    Dialog.addMessage("Which cell line are you analysing?");
    Dialog.addRadioButtonGroup("Cell line:",
        newArray("Huh7", "VeroE6"), 2, 1, "Huh7");

    Dialog.addMessage("---");
    Dialog.addCheckbox("Also save 8-bit RGB JPG with " + JPG_SCALEBAR_UM + " um scale bar", SAVE_JPG);

    Dialog.show();

    // get* in same order as add*
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

    // ALWAYS show the domain-list dialog (Feature A): pre-filled with the
    // workflow standards so a plain OK keeps the defaults.
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

// Cell-line-specific tuning (strictness knobs only — NOT the strategy).
function applyCellLineDefaults(cellLine) {
    if (cellLine == "VeroE6") {
        // Vero ~4× smaller than Huh7 → smaller nuclei, smaller territory.
        BLUR_SIGMA_NUC = 1;
        NUC_MIN_SIZE   = 500;
        NUC_DILATE_UM  = 5;
    } else {  // Huh7 default
        BLUR_SIGMA_NUC = 1;
        NUC_MIN_SIZE   = 1000;
        NUC_DILATE_UM  = 8;
    }
    print("Cell-line tuning for " + cellLine + ":");
    print("  BLUR_SIGMA_NUC = " + BLUR_SIGMA_NUC);
    print("  NUC_MIN_SIZE   = " + NUC_MIN_SIZE);
    print("  NUC_DILATE_UM  = " + NUC_DILATE_UM);
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
           + "roi_mode,nuc_dilate_um,"
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
    print("ROI strategy      : " + ROI_MODE + " (dilate " + NUC_DILATE_UM + " um)");
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
// Order (V0.4-style, ROI-first):
//   1) open + resolve metadata + validate
//   2) lookup bg1, bg2  (do NOT subtract yet)
//   3) split + rename channels
//   4) build the coloc ROI on the RAW channels / DAPI (unbiased, visible)
//   5) save coloc QC (composite + ROI) from the raw channels
//   6) NOW subtract bg from C1 and C2 (DAPI untouched)
//   7) colocalisation step (channels still split, ROI in manager)
//   8) merge back → save 16-bit TIF
//   9) optional bg-subtracted JPG
// ===========================================================

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

    // ---- lookup bg values (NOT subtracted yet) -----------
    bg1 = getBgValue(m1, comboKey, tp);
    bg2 = getBgValue(m2, comboKey, tp);
    print("  bg(" + m1 + ") = " + bg1 + ", bg(" + m2 + ") = " + bg2);

    // ---- split channels ----------------------------------
    splitAndRenameChannels(title, m1, m2);

    // ---- build ROI on RAW channels (BEFORE subtraction) --
    // Signal-independent ROI (unbiased coloc) and, in manual mode, the
    // whole cell is visible via autofluorescence.
    haveCytosolRoi = false;
    if (COLOC_MODE != "none") {
        haveCytosolRoi = buildColocRoi(m1);
        if (SAVE_COLOC_QC) saveColocQc(imgName, m1, m2);
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

// resolveMetadata / tryParseFilename / askImageMetadata / askParseFailureAction
// use the IJM-safe intermediate-variable pattern.
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


// ============== 5. ROI CONSTRUCTION ========================
// The coloc ROI ends up at ROI Manager index 0 in every path, so the
// downstream coloc step (activateChannelsAndRoi → use=[Channel 1]) is
// strategy-agnostic. Returns true if a usable ROI was created.
// ===========================================================

function buildColocRoi(m1) {
    if (ROI_MODE == "manual_draw") {
        r = buildManualRoi(m1);
        return r;
    }
    r = buildCentralCellRoi();
    return r;
}

// ------------------------------
// Nucleus mask from DAPI (raw). Morphological closing before Fill Holes
// handles boundary gaps. Same recipe as Mock / V0.4.
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
// AUTOMATIC central-cell ROI.
//   1) nucleus mask → all nuclei as ROIs (Analyze Particles "add").
//   2) for each nucleus: grow by NUC_DILATE_UM (calibrated µm) to a
//      territory; SKIP it if the territory's bbox touches the image
//      border (partially-imaged cell — user requirement).
//   3) among the fully-imaged candidates, choose the one whose nucleus
//      centre is closest to the image centre.
//   4) cytosol = territory XOR nucleus  (= territory ∖ nucleus, since
//      nucleus ⊂ territory), stored as the sole ROI at index 0.
// Returns false (→ full-image ROI fallback) if no nucleus is found.
// ------------------------------
function buildCentralCellRoi() {
    makeNucleusMask("DAPI_channel");

    roiManager("reset");
    selectWindow("Nucleus_Mask");
    run("Analyze Particles...", "size=" + NUC_MIN_SIZE + "-Infinity pixel add");
    nNuc = roiManager("count");
    if (nNuc == 0) {
        print("  WARN: no nuclei found — coloc ROI will be the full image");
        closeIfOpen("Nucleus_Mask");
        return false;
    }

    imgW = getWidth();
    imgH = getHeight();
    cx = imgW / 2.0;
    cy = imgH / 2.0;

    // Dilation distance in PIXELS (unambiguous): NUC_DILATE_UM is calibrated
    // µm; convert via the pixel width. If the image is uncalibrated, pw=1 and
    // the value is used as pixels directly. The "pixel" keyword on Enlarge
    // then forces pixel interpretation regardless of calibration.
    getPixelSize(unit, pw, ph);
    dilatePx = NUC_DILATE_UM;
    if (pw > 0) dilatePx = round(NUC_DILATE_UM / pw);
    if (dilatePx < 1) dilatePx = 1;
    print("  nucleus dilation: " + NUC_DILATE_UM + " " + unit + " -> " + dilatePx + " px");

    bestComplete = -1; bestCompleteDist = 1e18;   // fully-imaged candidates
    bestAny      = -1; bestAnyDist      = 1e18;    // fallback (ignore border)

    for (i = 0; i < nNuc; i = i + 1) {
        roiManager("Select", i);
        getSelectionBounds(nx, ny, nw, nh);        // nucleus bbox (px)
        ncx = nx + nw / 2.0;
        ncy = ny + nh / 2.0;
        dx = ncx - cx; dy = ncy - cy;
        dist = sqrt(dx*dx + dy*dy);

        // Grow to territory and test border touch (bbox approximation).
        run("Enlarge...", "enlarge=" + dilatePx + " pixel");
        getSelectionBounds(tx, ty, tw, th);        // territory bbox (px)
        touches = (tx <= BORDER_MARGIN_PX)
               || (ty <= BORDER_MARGIN_PX)
               || (tx + tw >= imgW - BORDER_MARGIN_PX)
               || (ty + th >= imgH - BORDER_MARGIN_PX);

        if (dist < bestAnyDist) { bestAnyDist = dist; bestAny = i; }
        if (!touches && dist < bestCompleteDist) { bestCompleteDist = dist; bestComplete = i; }
    }

    chosen = bestComplete; chosenDist = bestCompleteDist;
    if (chosen < 0) {
        print("  WARN: no fully-imaged cell — falling back to nearest-centre cell (touches border)");
        chosen = bestAny; chosenDist = bestAnyDist;
    }
    if (chosen < 0) { closeIfOpen("Nucleus_Mask"); return false; }
    print("  central cell = nucleus #" + chosen + " (dist to centre " + d2s(chosenDist, 1) + " px)");

    // territory ROI for the chosen nucleus
    roiManager("Select", chosen);
    run("Enlarge...", "enlarge=" + dilatePx + " pixel");
    roiManager("Add");
    terrIdx = roiManager("count") - 1;

    // cytosol = territory XOR nucleus → active image selection
    roiManager("Select", newArray(terrIdx, chosen));
    roiManager("XOR");

    // Keep ONLY the cytosol: reset the manager (image selection survives)
    // and re-add it so it lands at index 0 for the downstream coloc step.
    roiManager("reset");
    roiManager("Add");
    roiManager("Select", 0);
    roiManager("Rename", "Cytosol_central");

    closeIfOpen("Nucleus_Mask");
    return true;
}

// ------------------------------
// MANUAL ROI: user draws a freehand ROI around the central cell on the
// RAW m1 channel (autofluorescence makes the cell body visible), exactly
// like the supervisor's example. Returns false if nothing was drawn.
// ------------------------------
function buildManualRoi(m1) {
    ch = m1 + "_channel";
    if (!isOpen(ch)) { print("  WARN: " + ch + " missing for manual ROI"); return false; }
    selectWindow(ch);
    run("Enhance Contrast", "saturated=0.35");   // display only, not pixel data
    run("Select None");
    roiManager("reset");
    setTool("freehand");
    waitForUser("Draw ROI  [" + ch + "]",
        "Draw a FREEHAND ROI around the central cell on\n"
      + "  " + ch + "  (raw autofluorescence — whole cell visible)\n\n"
      + "Then click OK.\n(Cancel aborts the WHOLE macro run.)");
    if (selectionType() == -1) {
        print("  WARN: no ROI drawn — proceeding without ROI (full image)");
        return false;
    }
    roiManager("Add");
    roiManager("Select", 0);
    roiManager("Rename", "Cytosol_manual");
    return true;
}

// ------------------------------
// Per-image coloc QC: composite of all three RAW channels (m1 red,
// m2 green, DAPI blue) with the ACTUAL ROI burnt in as an outline.
// Built from DUPLICATES so the original channel windows stay intact for
// the subtraction / coloc / merge that follow (Mock pitfall: Merge
// Channels consumes its sources).
// ------------------------------
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
    if (roiManager("count") > 0) {
        roiManager("Select", 0);
        run("Add Selection...");
    }
    run("Flatten");

    run("Options...", "jpeg=90");
    qcPath = OUTPUT_DIR + "qc_" + imgName + "_roi.jpg";
    saveAs("Jpeg", qcPath);
    print("  saved coloc QC -> " + qcPath);

    close();
    closeIfOpen("Composite (RGB)");
    closeIfOpen("Composite");
    setOption("BlackBackground", true);   // Stack to RGB can flip the option
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
// Manual: preview → decision → plugin pause. Writes a PLACEHOLDER row.
// ------------------------------
function doColocManual(imgName, cellLine, tp, comboKey, m1, m2, haveRoi) {
    activateChannelsAndRoi(haveRoi, m1);
    roiLineShort  = "(no ROI — full image)";
    useRoiSetting = "<None>";
    if (haveRoi) {
        roiLineShort  = "Cell ROI is already applied to " + m1 + "_channel";
        useRoiSetting = "Channel 1";
    }

    waitForUser("Preview  [" + imgName + "]",
        "Three split channels are open:\n"
      + "  C1 : " + m1 + "_channel   (bg subtracted)\n"
      + "  C2 : " + m2 + "_channel   (bg subtracted)\n"
      + "  C3 : DAPI_channel\n\n"
      + "ROI: " + roiLineShort + "\n\n"
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
// output into the CSV. Snapshot log length BEFORE → substring AFTER is
// this invocation's chunk.
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
// Apply the cytosol ROI to m1_channel as the active selection. The plugin's
// "Use ROI: Channel 1" reads whatever selection is active on that image.
// ------------------------------
function activateChannelsAndRoi(haveRoi, m1) {
    if (!isOpen(m1 + "_channel")) return;
    selectWindow(m1 + "_channel");
    if (haveRoi && roiManager("count") > 0) {
        roiManager("Select", 0);
    }
}

// ------------------------------
// Parse a "Label: value" / "Label = value" pair from free-form log text.
// Returns "" if the label is not found.
// ------------------------------
function extractLogValue(text, label) {
    idx = indexOf(text, label);
    if (idx < 0) return "";
    after = substring(text, idx + lengthOf(label));
    j = 0;
    while (j < lengthOf(after)) {
        c = substring(after, j, j+1);
        if (c == ":" || c == "=" || c == " " || c == "\t") { j++; continue; }
        break;
    }
    after = substring(after, j);
    valEnd = lengthOf(after);
    for (i = 0; i < lengthOf(after); i++) {
        c = substring(after, i, i+1);
        if (c == " " || c == "\n" || c == "\t" || c == "(" || c == ",") {
            valEnd = i; break;
        }
    }
    return substring(after, 0, valEnd);
}

function tryLabels(text, labels) {
    for (i = 0; i < labels.length; i++) {
        v = extractLogValue(text, labels[i]);
        if (v != "") return v;
    }
    return "";
}

function appendColocRowFromLog(imgName, cellLine, tp, comboKey, m1, m2, logChunk, status) {
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
         + ROI_MODE + "," + NUC_DILATE_UM + ","
         + MACRO_VERSION + "," + RUN_ID + "," + status;
    File.append(line, csvPath);
    print("  coloc CSV row appended (" + status + "): Rtotal=" + rtotal
        + " Rcoloc=" + rcoloc + " M1=" + m1coef + " M2=" + m2coef);

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
         + ROI_MODE + "," + NUC_DILATE_UM + ","
         + MACRO_VERSION + "," + RUN_ID + "," + status;
    File.append(line, csvPath);
    print("  coloc CSV row appended (" + status + ")");
}


// ============== 7. BACKGROUND LOOKUP TABLE =================
// parseFloat at BOTH store + retrieve (IJM "Array.concat on empty array
// stores numbers as strings" quirk — HANDOFF §7).

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
    File.append("- ROI strategy  : " + ROI_MODE + " (dilate " + NUC_DILATE_UM + " um)\n", mdPath);
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
// 8-bit RGB JPG of the bg-subtracted composite with a scale bar burned in.
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
