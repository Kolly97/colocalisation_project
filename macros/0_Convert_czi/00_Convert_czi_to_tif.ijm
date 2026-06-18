// ==========================================================
// CZI -> TIF Batch Converter  (Stage 0 / pre-processing)  --  V0.1.0
// Author : Kolja Hildenbrand
// Date   : 2026-06-11
// Status : Current
//
// Purpose: Convert a folder of confocal .czi files to multi-channel .tif so the
//          Mock (Stage 1) and SubBg+Coloc (Stage 2) macros can read them.
//          Handles multi-SERIES .czi (one tif per real series), drops Zeiss
//          label/overview thumbnails, optionally collapses a z-STACK, and can
//          rename to the binding filename schema.
//
//          WHY a separate tool: the analysis macros do ONE job (measure). They
//          assume a clean single-plane, 3-channel tif (C1=marker1, C2=marker2,
//          C3=DAPI). Microscope .czi can carry series / z-stacks / label images
//          / odd names. Keeping the format-wrangling HERE keeps the analysis
//          macros simple and robust (one job per tool).
//
// Requires: Bio-Formats (ships with Fiji).
//
// Filename schema (Stage 1/2 parse this):
//   timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif
//   C3 = DAPI ALWAYS. The converter PRESERVES the czi channel order — so verify
//   in your acquisition that channels are marker1, marker2, DAPI in that order.
//
// Output: <input>/tif_converted/
//           <name>.tif                 multi-channel 16-bit (one per real series)
//           name_map.csv               original .czi -> output .tif (provenance)
//           convert_log_<RUN_ID>.txt   full IJ Log
// ==========================================================


// ============== 1. CONFIG =================================

var INPUT_DIR;
var OUTPUT_DIR;
var RUN_ID;

// z handling for z-STACK czi: "as_is" | "single_mid" | "max"
//  - as_is      : keep the plane(s) as they are. Use ONLY if data is already 2D.
//  - single_mid : keep the middle z-plane. RECOMMENDED for colocalisation — a max
//                 projection merges signal from different depths and can INFLATE
//                 the apparent coloc value.
//  - max        : maximum-intensity projection over z.
var Z_MODE = "as_is";

// naming: "keep" the original base name | "template" build the schema name.
var NAME_MODE = "keep";
// Template fields (used only if NAME_MODE=="template"). Assumes the whole folder
// is ONE batch = one (timepoint x cellLine x condition x marker combo).
var T_TP   = "24h";
var T_CELL = "Huh7";
var T_COND = "MOI1";
var T_M1   = "HA568";
var T_M2   = "dsRNA488";

// Series smaller than this (px, width OR height) are treated as Zeiss
// label/overview thumbnails and skipped (not real acquisition data).
var MIN_DIM = 256;

var MACRO_VERSION = "0.1.0";


// ============== MAIN ======================================

chooseInputDir();
RUN_ID = makeRunId();
askConfig();
buildOutputDir();
files = listCziFiles();

print("=== CZI->TIF converter V" + MACRO_VERSION + " | run " + RUN_ID + " ===");
print("input  : " + INPUT_DIR);
print("output : " + OUTPUT_DIR);
print("z mode : " + Z_MODE + "    naming: " + NAME_MODE);
print(".czi files: " + files.length);
print("---");
setBatchMode(true);

nOut = 0;
for (f = 0; f < files.length; f = f + 1) {
    print("[" + (f+1) + "/" + files.length + "] " + files[f]);
    nOut = nOut + convertOne(files[f]);
    cleanup();
}
print("=== DONE: wrote " + nOut + " tif(s) ===");
saveLog();


// ============== 2. SETUP ==================================

function chooseInputDir() {
    INPUT_DIR = getDirectory("Choose folder with .czi images");
    if (INPUT_DIR == "") exit("No folder selected.");
}

function makeRunId() {
    getDateAndTime(y, mo, dw, d, h, mn, s, ms);
    return "" + y + IJ.pad(mo+1, 2) + IJ.pad(d, 2) + "_" + IJ.pad(h, 2) + IJ.pad(mn, 2);
}

function askConfig() {
    Dialog.create("CZI -> TIF converter");
    Dialog.addMessage("Z handling (only matters for z-stacks):");
    Dialog.addRadioButtonGroup("z:",
        newArray("as_is (already 2D)", "single middle plane", "max projection"),
        3, 1, "as_is (already 2D)");
    Dialog.addMessage("---");
    Dialog.addMessage("File naming:");
    Dialog.addRadioButtonGroup("name:",
        newArray("keep original name", "build schema name (template)"),
        2, 1, "keep original name");
    Dialog.show();

    zc = Dialog.getRadioButton();
    if (zc == "max projection")            Z_MODE = "max";
    else if (zc == "single middle plane")  Z_MODE = "single_mid";
    else                                   Z_MODE = "as_is";

    nc = Dialog.getRadioButton();
    if (nc == "build schema name (template)") NAME_MODE = "template";
    else                                      NAME_MODE = "keep";

    // Second dialog only when a template name is requested.
    if (NAME_MODE == "template") {
        Dialog.create("Schema template (whole folder = one batch)");
        Dialog.addMessage("Output name = tp_cellLine_condition_marker1_marker2_<orig>.tif\n"
                        + "condition MUST contain 'mock' or 'moi' so Stage 1/2 pick the right files.");
        Dialog.addString("timepoint:",     T_TP);
        Dialog.addString("cell line:",      T_CELL);
        Dialog.addString("condition:",      T_COND);
        Dialog.addString("marker1 (C1):",   T_M1);
        Dialog.addString("marker2 (C2):",   T_M2);
        Dialog.show();
        T_TP   = trim(Dialog.getString());
        T_CELL = trim(Dialog.getString());
        T_COND = trim(Dialog.getString());
        T_M1   = trim(Dialog.getString());
        T_M2   = trim(Dialog.getString());
    }
}

function buildOutputDir() {
    OUTPUT_DIR = INPUT_DIR + "tif_converted" + File.separator;
    if (!File.exists(OUTPUT_DIR)) File.makeDirectory(OUTPUT_DIR);
    // Provenance map: which .czi (and which series) became which .tif.
    mapPath = OUTPUT_DIR + "name_map.csv";
    if (!File.exists(mapPath))
        File.append("original_czi,series_index,output_tif,z_mode,macro_version,run_id", mapPath);
}

function listCziFiles() {
    all = getFileList(INPUT_DIR);
    out = newArray();
    for (i = 0; i < all.length; i = i + 1)
        if (endsWith(toLowerCase(all[i]), ".czi")) out = Array.concat(out, all[i]);
    out = Array.sort(out);
    return out;
}


// ============== 3. CONVERT ONE FILE =======================
// Opens every series of one .czi via Bio-Formats, drops label/overview
// thumbnails, applies the z choice, saves each real series as a tif.
// Returns how many tif files were written.

function convertOne(fname) {
    base = fname;
    if (lastIndexOf(base, ".") > 0) base = substring(base, 0, lastIndexOf(base, "."));

    before = getList("image.titles");
    // Bio-Formats: open ALL series as separate images (multi-position czi -> many;
    // single-position -> one + maybe label/overview thumbnails).
    // NOTE: confirm this option string via Plugins>Macros>Record on your build.
    run("Bio-Formats Importer",
        "open=[" + INPUT_DIR + fname + "] color_mode=Default view=Hyperstack "
      + "stack_order=XYCZT open_all_series");
    opened = newImagesSince(before);
    if (opened.length == 0) { print("  WARN: nothing opened for " + fname); return 0; }

    written = 0;
    realIdx = 0;                                    // counts only the kept (real) series
    for (s = 0; s < opened.length; s = s + 1) {
        selectWindow(opened[s]);
        // Drop Zeiss label / overview thumbnails (much smaller than acquisition).
        if (getWidth() < MIN_DIM || getHeight() < MIN_DIM) {
            print("  skip series '" + opened[s] + "' (" + getWidth() + "x" + getHeight()
                + " < " + MIN_DIM + " px) — looks like a label/overview");
            close();
            continue;
        }
        realIdx = realIdx + 1;
        applyZ();                                   // collapse z if needed (active image)
        outName = buildOutName(base, realIdx);
        saveAs("Tiff", OUTPUT_DIR + outName + ".tif");
        print("  -> " + outName + ".tif");
        File.append(fname + "," + realIdx + "," + outName + ".tif," + Z_MODE
                  + "," + MACRO_VERSION + "," + RUN_ID, OUTPUT_DIR + "name_map.csv");
        written = written + 1;
    }
    return written;
}

// Collapse the z dimension of the ACTIVE image per Z_MODE (no-op if single plane).
function applyZ() {
    Stack.getDimensions(w, h, c, z, t);
    if (z <= 1) return;                             // already 2D — nothing to do
    if (Z_MODE == "max") {
        run("Z Project...", "projection=[Max Intensity]");   // -> MAX_<title>, channels kept
    } else if (Z_MODE == "single_mid") {
        mid = round(z / 2); if (mid < 1) mid = 1;
        run("Duplicate...", "title=zpick duplicate slices=" + mid);   // all channels at z=mid
    } else {
        print("  WARN: z-stack with " + z + " planes kept AS-IS — Stage 1/2 expect ONE plane");
    }
}

// Build the output base name (no extension). naming = keep | template.
function buildOutName(origBase, realIdx) {
    clean = sanitize(origBase);
    if (NAME_MODE == "template")
        name = T_TP + "_" + T_CELL + "_" + T_COND + "_" + T_M1 + "_" + T_M2 + "_" + clean;
    else
        name = clean;
    if (realIdx > 1) name = name + "_p" + realIdx;  // disambiguate multi-series czi
    return name;
}

// Replace characters that would break a filename or the underscore token scheme.
function sanitize(s) {
    s = replace(s, " ", "-");
    s = replace(s, "/", "-");
    return s;
}


// ============== 4. UTILS ==================================

// Titles of images opened since the `before` snapshot (Bio-Formats may open
// several series at once; the active window alone is not reliable).
function newImagesSince(before) {
    cur = getList("image.titles");
    out = newArray();
    for (i = 0; i < cur.length; i = i + 1)
        if (!inArray(cur[i], before)) out = Array.concat(out, cur[i]);
    return out;
}

function inArray(v, arr) {
    for (i = 0; i < arr.length; i = i + 1) if (arr[i] == v) return true;
    return false;
}

function cleanup() {
    while (nImages > 0) { selectImage(nImages); close(); }
}

function saveLog() {
    if (!isOpen("Log")) return;
    selectWindow("Log");
    saveAs("Text", OUTPUT_DIR + "convert_log_" + RUN_ID + ".txt");
}
