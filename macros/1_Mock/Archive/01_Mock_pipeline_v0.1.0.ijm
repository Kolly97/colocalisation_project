// ==========================================================
// Mock Top-X% Pipeline  --  V0.1 (Production Mode)
// Author 	: Kolja Hildenbrand
// Date   	: 2026-05-07
// Status	: old
// Purpose: For each Mock image in a folder, compute the
//          Top-X% pixel intensity statistics inside the
//          cytosol of each marker channel. Write per-image
//          rows into per-(timepoint, combo, marker) CSVs.
// Filename schema (binding):
//     timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif
//     e.g.  12h_Huh7_Mock_HA568_dsRNA488_CS1_1.tif
//     C1 = marker1, C2 = marker2, C3 = DAPI (always)
// ==========================================================


// ---------- 1. CONFIG (var = global) -----------------------
// Why `var`? In IJM, variables defined inside a function are
// LOCAL to that function. To share state across functions
// (paths, parameters), declare them with `var` at top level.
// We deliberately keep CONFIG separate from logic — change
// thresholds/percentages here, never inside functions.
// -----------------------------------------------------------
var INPUT_DIR;
var MEASURE_DIR;
var RUN_ID;

// what we analyse
var ANALYSE_COMBI = newArray("HA568_dsRNA488", "NS4B568_dsRNA488", "NS4B568_HA488");
var TIMEPOINTS    = newArray("12h", "24h");

// thresholding
var CELL_THR_METHOD = "Li";    // tried later: Li, Otsu, Triangle, Yen
var NUC_THR_METHOD  = "Otsu";  // tried later: Otsu, Triangle
var BLUR_SIGMA_CELL = 1;
var BLUR_SIGMA_NUC  = 1;

// top-percentile measurement
var TOP_PCT     = 1.0;                 // top 1% of cytosol pixels
var STAT_METHOD = "median_top_hist";   // bookkeeping label only

// outputs
var SAVE_QC    = true;
var SAVE_MASKS = true;

// reproducibility
var MACRO_VERSION = "0.1.0";

// CSV header. Order must match appendCsvRow() exactly.
var CSV_HEADER = "image,cell_line,timepoint,combo,channel,"
               + "stat_method,top_pct,"
               + "threshold_value,n_top_pixels,n_cyto_pixels,"
               + "mean_top,median_top,std_top,"
               + "p95,p99,p99_25,p99_5,p99_9,"
               + "cell_thr_method,nuc_thr_method,"
               + "blur_sigma_cell,blur_sigma_nuc,"
               + "macro_version,run_id\n";


// ============== MAIN =======================================
// Every real step is a named function. Easy to comment out, easy to test.
// ===========================================================
chooseInputDir();
RUN_ID = makeRunId();
buildOutputDir();
initCsvFiles();

mockFiles = listMockFiles();
print("=== Mock pipeline V" + MACRO_VERSION + " | run_id=" + RUN_ID + " ===");
print("Input  : " + INPUT_DIR);
print("Output : " + MEASURE_DIR);
print("Found " + mockFiles.length + " mock .tif files.");

// setBatchMode("hide") = images don't draw, but everything
// else (selectWindow, ROI, getHistogram) still works. ~5-10x
// speedup on big batches. Toggle to false for debugging.
// setBatchMode("hide");
for (f = 0; f < mockFiles.length; f++) {
    print("[" + (f+1) + "/" + mockFiles.length + "] " + mockFiles[f]);
    processOneImage(mockFiles[f]);
    cleanupBetweenImages();
}
setBatchMode("exit and display");

print("=== DONE ===");


// ============== 2. SETUP FUNCTIONS =========================

function chooseInputDir() {
    INPUT_DIR = getDirectory("Choose folder with Mock + MOI .tif images");
    if (INPUT_DIR == "") exit("No folder selected.");
}

// ISO-ish timestamp without separators — safe for filesystems.
// IJ.pad(n,2) zero-pads to width 2.
function makeRunId() {
    getDateAndTime(y, m, dw, d, h, mn, s, ms);
    return "" + y + IJ.pad(m+1, 2) + IJ.pad(d, 2)
        + "_" + IJ.pad(h, 2) + IJ.pad(mn, 2);
}

function buildOutputDir() {
    MEASURE_DIR = INPUT_DIR + "measure_mock" + File.separator;
    if (!File.exists(MEASURE_DIR)) File.makeDirectory(MEASURE_DIR);
    if (SAVE_MASKS && !File.exists(MEASURE_DIR + "masks"))
        File.makeDirectory(MEASURE_DIR + "masks");
    if (SAVE_QC && !File.exists(MEASURE_DIR + "qc"))
        File.makeDirectory(MEASURE_DIR + "qc");
}

// Create one CSV per (timepoint, combo, marker). 12 files total
// (2 tp × 3 combos × 2 markers). Header only on first creation —
// otherwise we'd append a header line to existing data.
function initCsvFiles() {
    for (i = 0; i < ANALYSE_COMBI.length; i++) {
        combo = ANALYSE_COMBI[i];
        parts = split(combo, "_");
        m1 = parts[0]; m2 = parts[1];
        for (t = 0; t < TIMEPOINTS.length; t++) {
            tp = TIMEPOINTS[t];
            csv1 = MEASURE_DIR + "Mock_" + tp + "_" + m1 + "_in_" + combo + ".csv";
            csv2 = MEASURE_DIR + "Mock_" + tp + "_" + m2 + "_in_" + combo + ".csv";
            if (!File.exists(csv1)) File.append(CSV_HEADER, csv1);
            if (!File.exists(csv2)) File.append(CSV_HEADER, csv2);
        }
    }
}

// Filter: only `.tif` files containing "mock" (case-insensitive).
// Sort to make order reproducible across OS file-listing quirks.
function listMockFiles() {
    files = getFileList(INPUT_DIR);
    out = newArray();
    for (i = 0; i < files.length; i++) {
        n = files[i];
        ln = toLowerCase(n);
        if (endsWith(ln, ".tif") && indexOf(ln, "mock") >= 0)
            out = Array.concat(out, n);
    }
    out = Array.sort(out);
    return out;
}


// ============== 3. PER-IMAGE PIPELINE ======================

function processOneImage(fname) {
    open(INPUT_DIR + fname);
    title = getTitle();
    imgName = substring(title, 0, lastIndexOf(title, "."));

    // ---- parse filename (positional, schema is binding) ----
    tokens = split(imgName, "_");
    if (tokens.length < 7) {
        print("  SKIP: bad filename schema (need 7 tokens): " + imgName);
        return;
    }
    tp       = tokens[0];
    cellLine = tokens[1];
    cond     = tokens[2];
    m1       = tokens[3];
    m2       = tokens[4];
    comboKey = m1 + "_" + m2;

    if (!isValidCombo(comboKey)) {
        print("  SKIP: combo not in ANALYSE_COMBI -> " + comboKey);
        return;
    }
    print("  combo=" + comboKey + " tp=" + tp + " cellLine=" + cellLine);

    // ---- split + rename channels by marker name -----------
    splitAndRenameChannels(title, m1, m2);

    // ---- masks --------------------------------------------
    cellSrc = pickCellMaskSource(m1, m2);
    if (cellSrc == "") {
        print("  SKIP: no marker channel suitable as Cell-Mask source");
        return;
    }
    print("  Cell-Mask source: " + cellSrc);

    makeCellMask(cellSrc);
    makeNucleusMask("DAPI_channel");
    makeCytosolMask();

    // Guard: empty cytosol = nothing to measure
    selectWindow("Cytosol_Mask");
    getStatistics(area, meanCyto);
    if (meanCyto == 0) {
        print("  SKIP: cytosol mask empty after subtraction");
        return;
    }

    // Convert cytosol foreground to a real ROI we can re-apply
    // to ANY channel — cleaner than multiplying images.
    roiManager("reset");
    selectWindow("Cytosol_Mask");
    run("Create Selection");
    roiManager("Add");
    cytoRoiId = 0;
    nCyto = getRawStatisticsCount();   // count of foreground pixels

    // ---- measure each marker channel ----------------------
    measureAndWrite(m1, "" + m1 + "_channel", cytoRoiId, nCyto,
                    imgName, cellLine, tp, comboKey);
    measureAndWrite(m2, "" + m2 + "_channel", cytoRoiId, nCyto,
                    imgName, cellLine, tp, comboKey);

    // ---- save artefacts -----------------------------------
    if (SAVE_QC)    saveQcOverlay(imgName, m1);
    if (SAVE_MASKS) saveMasksTif(imgName);
}


// ============== 4. CHANNEL HANDLING ========================

function isValidCombo(comboKey) {
    for (i = 0; i < ANALYSE_COMBI.length; i++)
        if (ANALYSE_COMBI[i] == comboKey) return true;
    return false;
}

// Channel-1 = marker1, Channel-2 = marker2, Channel-3 = DAPI.
// We rename windows by marker name, so downstream code stays
// readable: "HA568_channel" not "C1-foo.tif".
function splitAndRenameChannels(title, m1, m2) {
    selectWindow(title);
    run("Split Channels");
    selectWindow("C1-" + title); rename(m1 + "_channel");
    selectWindow("C2-" + title); rename(m2 + "_channel");
    selectWindow("C3-" + title); rename("DAPI_channel");
}

// Cell-Mask source priority (your spec): HA568 > HA488 > NS4B568 > dsRNA488.
// Reason: HA gives the best autofluorescence/background coverage of the
// whole cell shape in Mock, where viral markers are absent.
function pickCellMaskSource(m1, m2) {
    priority = newArray("HA568", "HA488", "NS4B568", "dsRNA488");
    for (i = 0; i < priority.length; i++) {
        cand = priority[i] + "_channel";
        if (isOpen(cand)) return cand;
    }
    return "";
}


// ============== 5. MASKS ===================================
// All masks are 8-bit binary 0/255 — the ImageJ standard.
// Reason vs 0/1: native binary tools (Fill Holes, Create
// Selection, Watershed, ...) all expect 0/255. Using 0/1
// breaks them silently. We don't need 0/1 for measurement
// because we use ROIs, not multiplication.
// ===========================================================

function makeCellMask(srcTitle) {
    selectWindow(srcTitle);
    run("Duplicate...", "title=Cell_Mask");
    // Gaussian smoothing before global threshold reduces salt-and-pepper
    // artefacts; sigma=1 px is a conservative default for IF data.
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_CELL);
    setAutoThreshold(CELL_THR_METHOD + " dark");
    run("Convert to Mask");
    // Fill Holes: cell interior often dimmer than rim → small holes
    // get closed. Crucial for clean cytosol later.
    run("Fill Holes");
}

function makeNucleusMask(dapiTitle) {
    selectWindow(dapiTitle);
    run("Duplicate...", "title=Nucleus_Mask");
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA_NUC);
    setAutoThreshold(NUC_THR_METHOD + " dark");
    run("Convert to Mask");
    run("Fill Holes");
}

// Cytosol = Cell AND NOT Nucleus.
// Done as `Cell - Nucleus` because in 8-bit binary:
//   255 - 255 = 0  (nucleus area dropped)
//   255 - 0   = 255 (cell area kept)
//     0 - *   = 0  (outside cell stays out)
// 8-bit subtraction clamps at 0, so this is safe & simple.
function makeCytosolMask() {
    imageCalculator("Subtract create", "Cell_Mask", "Nucleus_Mask");
    rename("Cytosol_Mask");
    // Defensive re-binarise: guarantees clean 0/255 even if any
    // upstream step left intermediate values somewhere.
    setThreshold(1, 255);
    run("Convert to Mask");
}


// ============== 6. MEASUREMENT (TOP-X% via histogram) ======
// Why histogram: with a 16-bit image and a possibly large
// cytosol selection, sorting all pixel values in IJM is slow
// and memory-heavy. The histogram has fixed size 65536 and
// gives us EXACT per-integer counts — every percentile,
// mean, std of the top pool falls out by walking the bins.
// Active ROI selection limits the histogram to cytosol only.
// ===========================================================

function measureAndWrite(marker, chTitle, cytoRoiId, nCyto,
                         imgName, cellLine, tp, comboKey) {
    if (!isOpen(chTitle)) {
        print("  WARN: channel window missing: " + chTitle);
        return;
    }
    selectWindow(chTitle);
    roiManager("Select", cytoRoiId);

    // Full integer histogram for 16-bit: 65536 bins covering 0..65535.
    // values[i] = i, counts[i] = #pixels in selection with value i.
    nBins = 65536;
    getHistogram(values, counts, nBins, 0, 65535);

    nTotal = 0;
    for (i = 0; i < nBins; i++) nTotal += counts[i];
    if (nTotal == 0) {
        print("  WARN: 0 pixels in cytosol selection for " + chTitle);
        return;
    }

    stats = computeTopStats(values, counts, nBins, nTotal, TOP_PCT);
    // stats = [thr, nTop, meanTop, medianTop, stdTop, p95, p99, p99_25, p99_5, p99_9]

    csvPath = MEASURE_DIR + "Mock_" + tp + "_" + marker + "_in_" + comboKey + ".csv";
    appendCsvRow(csvPath, imgName, cellLine, tp, comboKey, marker,
                 stats, nCyto);
}

// Pure number-crunching, no UI / no globals beyond config.
// Easier to test and to lift into the tuning dispatcher later.
function computeTopStats(values, counts, nBins, nTotal, topPct) {
    // 1) Find the threshold value V* such that
    //    count(pixels with value >= V*) >= topPct% of nTotal.
    nTopTarget = floor(nTotal * topPct / 100.0);
    if (nTopTarget < 1) nTopTarget = 1;
    acc = 0;
    threshold_value = 0;
    for (i = nBins - 1; i >= 0; i--) {
        if (counts[i] > 0) {
            acc += counts[i];
            if (acc >= nTopTarget) { threshold_value = values[i]; break; }
        }
    }

    // 2) Stats over the top pool (value >= V*). We use weighted
    //    sums over histogram bins — exact for integer-valued data.
    nTop = 0; sumV = 0; sumV2 = 0;
    for (i = 0; i < nBins; i++) {
        if (values[i] >= threshold_value) {
            nTop  += counts[i];
            sumV  += values[i] * counts[i];
            sumV2 += values[i] * values[i] * counts[i];
        }
    }
    mean_top = sumV / nTop;
    var_top  = (sumV2 / nTop) - mean_top * mean_top;
    std_top  = sqrt(maxOf(var_top, 0));

    // 3) Median of the top pool: walk bins from V* upward,
    //    pick the bin where cumulative count crosses nTop/2.
    half = nTop / 2.0;
    a = 0; median_top = threshold_value;
    for (i = 0; i < nBins; i++) {
        if (values[i] >= threshold_value) {
            a += counts[i];
            if (a >= half) { median_top = values[i]; break; }
        }
    }

    // 4) Whole-cytosol percentiles (NOT just top pool) — useful
    //    later for outlier-image detection / sanity checks.
    p95    = pct(values, counts, nBins, nTotal, 95);
    p99    = pct(values, counts, nBins, nTotal, 99);
    p99_25 = pct(values, counts, nBins, nTotal, 99.25);
    p99_5  = pct(values, counts, nBins, nTotal, 99.5);
    p99_9  = pct(values, counts, nBins, nTotal, 99.9);

    return newArray(threshold_value, nTop, mean_top, median_top, std_top,
                    p95, p99, p99_25, p99_5, p99_9);
}

function pct(values, counts, nBins, nTotal, p) {
    target = nTotal * p / 100.0;
    a = 0;
    for (i = 0; i < nBins; i++) {
        a += counts[i];
        if (a >= target) return values[i];
    }
    return values[nBins - 1];
}

// Helper: count of pixels in the active selection.
function getRawStatisticsCount() {
    getRawStatistics(n);
    return n;
}


// ============== 7. CSV / OUTPUT ============================
// CSV append: must match CSV_HEADER order exactly. If you add
// a column, change BOTH this function and CSV_HEADER.
function appendCsvRow(csvPath, imgName, cellLine, tp, comboKey, channel,
                      stats, nCyto) {
    line = imgName + "," + cellLine + "," + tp + "," + comboKey + "," + channel
        + "," + STAT_METHOD + "," + TOP_PCT
        + "," + stats[0] + "," + stats[1] + "," + nCyto
        + "," + stats[2] + "," + stats[3] + "," + stats[4]
        + "," + stats[5] + "," + stats[6] + "," + stats[7] + "," + stats[8] + "," + stats[9]
        + "," + CELL_THR_METHOD + "," + NUC_THR_METHOD
        + "," + BLUR_SIGMA_CELL + "," + BLUR_SIGMA_NUC
        + "," + MACRO_VERSION + "," + RUN_ID;
    File.append(line, csvPath);
}

// Save the three masks as separate binary TIFs. Note: saveAs
// renames the active window to the saved filename → we close
// at end of iteration anyway, so renaming side effects don't
// matter. We do it AFTER measurements for that reason.
function saveMasksTif(imgName) {
    masksDir = MEASURE_DIR + "masks" + File.separator;
    selectWindow("Cell_Mask");    saveAs("Tiff", masksDir + imgName + "_cell.tif");
    selectWindow("Nucleus_Mask"); saveAs("Tiff", masksDir + imgName + "_nuc.tif");
    selectWindow("Cytosol_Mask"); saveAs("Tiff", masksDir + imgName + "_cyto.tif");
}

// QC PNG: marker1 channel auto-contrasted to 8-bit, with the
// cytosol selection drawn as outline. Quick visual sanity
// check — if the outline doesn't follow the cells, your
// thresholding parameters are wrong.
function saveQcOverlay(imgName, m1) {
    qcDir = MEASURE_DIR + "qc" + File.separator;
    src = m1 + "_channel";
    if (!isOpen(src)) return;

    selectWindow(src);
    run("Duplicate...", "title=qc_tmp");
    run("Enhance Contrast", "saturated=0.35");  // display only
    run("8-bit");

    selectWindow("Cytosol_Mask");
    run("Create Selection");
    selectWindow("qc_tmp");
    run("Restore Selection");
    run("Add Selection...");
    run("Flatten");
    saveAs("PNG", qcDir + imgName + "_qc.png");

    // close the flattened image (now active after saveAs)
    close();
    if (isOpen("qc_tmp")) { selectWindow("qc_tmp"); close(); }
}


// ============== 8. CLEANUP =================================

function cleanupBetweenImages() {
    // Close every open image. Don't touch Log / Results window.
    while (nImages > 0) { selectImage(nImages); close(); }
    if (isOpen("ROI Manager")) roiManager("reset");
    if (isOpen("Results"))     run("Clear Results");
    if (selectionType() != -1) run("Select None");
}