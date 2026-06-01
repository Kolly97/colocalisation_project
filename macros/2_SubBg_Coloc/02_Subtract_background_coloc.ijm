// ==========================================================
// Background Subtraction (+ Coloc prep)  --  V0.0.1
// Author : Kolja Hildenbrand
// Status : FIRST CONCEPTUAL DRAFT ("vibe coding") — NOT functional.
//          Contains syntax errors and incomplete logic. Kept as the
//          historical starting point of the pipeline. Do NOT run this;
//          use the latest version instead.
//
// Purpose:
//   Sketch how infected (MOI) images should be background-corrected
//   before colocalisation analysis. Background values are measured by
//   the Mock pipeline (separate script) and entered here by hand.
//
// Expected filename scheme (binding):
//   timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif
//   (C1 = marker1, C2 = marker2, C3 = DAPI)
//
// Intended order of operations:
//   1. chooseInputDir()        - pick folder with MOI .tif images
//   2. makeRunId()             - timestamp YYYYMMDD_HHMM for the run
//   3. askModeAndConfig()      - automatic vs. manual metadata mode
//   4. get_background_values() - one dialog for all bg values
//                                (per marker x combo x timepoint)
//   5. buildOutputDir()        - <RUN_ID>_substracted_images/ +
//                                background_values_used.md (provenance)
//   6. listMOIFiles()          - all .tif files containing "moi"
//   7. loop over images -> processOneImage():
//        open -> parse filename (tryParseFilename / askImageMetadata)
//        -> validate combo (isValidCombo)
//        -> split + rename channels (splitAndRenameChannels)
//        -> subtract matching bg value from each marker channel
//   8. cleanupBetweenImages()  - close windows between images
//
// What it can do (as concepts in this draft):
//   - folder selection + reproducible run ID
//   - manual entry of all background values in one dialog
//   - markdown documentation of the entered values
//   - filename parsing with a manual-entry fallback
//   - marker-combination validation
//   - channel splitting + marker-based renaming
//   - batch loop over all MOI images
//
// Not yet implemented (added in later versions):
//   - the actual background subtraction body
//   - cytosol-ROI generation, colocalisation step, image saving
//   - scalable bg-value lookup table (uses 12 named globals here)
// ==========================================================

// Set during Main

var HA568_in_HA568_dsRNA488_12h;
var dsRNA488_in_HA568_dsRNA488_12h;

var NS4B568_in_NS4B568_dsRNA488_12h;
var dsRNA488_in_NS4B568_dsRNA488_12h;

var NS4B568_in_NS4B568_HA488_12h;
var HA488_in_NS4B568_HA488_12h;

var HA568_in_HA568_dsRNA488_24h;
var dsRNA488_in_HA568_dsRNA488_24h;

var NS4B568_in_NS4B568_dsRNA488_24h;
var dsRNA488_in_NS4B568_dsRNA488_24h;

var NS4B568_in_NS4B568_HA488_24h;
var HA488_in_NS4B568_HA488_24h;

var M1_bg
var M2_bg

var MARKERS       = newArray("HA568", "HA488", "dsRNA488", "NS4B568");
var ANALYSE_COMBI = newArray("HA568_dsRNA488", "NS4B568_dsRNA488", "NS4B568_HA488");
var TIMEPOINT = newArray("12 h", "24 h");

var INPUT_DIR;
var OUTPUT_DIR;
var MODE
var RUN_ID
var MACRO_VERSION = "0.0.1";


// ============== MAIN =======================================
chooseInputDir(); 
RUN_ID = makeRunId();
askModeAndConfig ()
get_background_values();
buildOutputDir();
ImageFiles = listMOIFiles();

print("=== Mock pipeline V" + MACRO_VERSION + " | run_id=" + RUN_ID + " ===");
print("Input  : " + INPUT_DIR);
print("Output : " + OUTPUT_DIR);
print("Found " + ImageFiles.length + " MOI .tif files.");
print("=== RUN BOY RUN ===");

for (f = 0; f < ImageFiles.length; f++) {
    print("[" + (f+1) + "/" + ImageFiles.length + "] " + ImageFiles[f]);
    processOneImage(ImageFiles[f]);
    cleanupBetweenImages();
}
print("=== DONE ===");


// ============== FUNCTIONS =======================================

// ------------------------------
// Function for mode selection
// ------------------------------

function askModeAndConfig() {
    Dialog.create("Pipeline mode");
    Dialog.addMessage("How should the macro run?");
    Dialog.addRadioButtonGroup("Mode:",
        newArray("automatic", "manual"), 2, 1, "automatic");
    Dialog.show();
    MODE = Dialog.getRadioButton();
}

// ------------------------------
// Function for Image Handling
// ------------------------------

function processOneImage(fname) {
    open(INPUT_DIR + fname);
    title = getTitle();
    imgName = substring(title, 0, lastIndexOf(title, "."));

    // tryParseFilename returns:
    //   newArray(tp, cellLine, m1, m2)  on success
    //   empty array                     on failure
    parsed = tryParseFilename(imgName);

    // ---- decide metadata source -----------------------------
    tp = ""; m1 = ""; m2 = "";

    if (MODE == "manual") {
        if (parsed.length == 0) {
            print("  Parse failed for: " + imgName);
            action = askParseFailureAction(imgName);
            if (action == "skip") { print("  -> skipped"); return; }
            // else: continue into manual
            meta = askImageMetadata("(parse failed)", "", "", "", "");
            tp = meta[0], m1 = meta[0]; m2 = meta[1];
        } else {
        	tp = parsed[0], m1 = parsed[1]; m2 = parsed[2];
        }
        }
    } else {
        // manual mode: pre-fill from filename if parse worked, else blanks
        defM1 = ""; defM2 = "";
        if (parsed.length > 0) { defM1 = parsed[0]; defM2 = parsed[1]; }
        meta = askImageMetadata(imgName, defTP, defM1, defM2);
        m1 = meta[0]; m2 = meta[1];
    }

    // ---- validate -----------------------------------------
    comboKey = m1 + "_" + m2;
    if (!isValidCombo(comboKey)) {
        print("  SKIP: combo not in ANALYSE_COMBI -> " + comboKey);
        return;
    }
	/// Hier will ich printen, welche wert von welchem channel abgezogen wird
    print("  combo=" + comboKey + " tp=" + tp;

    // ---- channels -----------------------------------------
    splitAndRenameChannels(title, m1, m2);
    
    substactBackground(m1, m2);
    
// ------------------------------
// Substract background
// ------------------------------
if ANALYSE_COMBI == comboKey
	

// ============== 4. CHANNEL HANDLING ========================

// ------------------------------
// Function to split and rename channels
// ------------------------------

function splitAndRenameChannels(title, m1, m2) {
    selectWindow(title);
    run("Split Channels");
    selectWindow("C1-" + title); rename(m1 + "_channel");
    selectWindow("C2-" + title); rename(m2 + "_channel");
    selectWindow("C3-" + title); rename("DAPI_channel");
}

// ------------------------------
// Check if combi is valid
// ------------------------------

function isValidCombo(comboKey) {
    return inArray(comboKey, ANALYSE_COMBI);
}

// ------------------------------
// Function for getting file name
// ------------------------------

function tryParseFilename(imgName) {
    tokens = split(imgName, "_");
    if (tokens.length < 7) return newArray();
    tp_ 	= tokens[0];
    m1_     = tokens[3];
    m2_     = tokens[4];
    if (!inArray(m1_, MARKERS) || !inArray(m2_, MARKERS)) return newArray();
    if (!inArray(tp_, TIMEPOINTS))							return newArray();
    return newArray(tp_, m1_, m2_);
}


// ------------------------------
// Function if parsing failes
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
// Function to parse meta data
// ------------------------------


function askImageMetadata(imgLabel, defTP, defM1, defM2) {
    if (defTP == "" || !inArray(defTP, TIMEPOINT))    defTP = TIMEPOINT[0];
    if (defM1 == "" || !inArray(defM1, MARKERS))    defM1 = MARKERS[0];
    if (defM2 == "" || !inArray(defM2, MARKERS))    defM2 = MARKERS[1 % MARKERS.length];

    Dialog.create("Image metadata");
    Dialog.addMessage("Image: " + imgLabel + "\nC3 is always DAPI.");
    Dialog.addRadioButtonGroup("Timepoint:", TIMEPOINT, 1, TIMEPOINT.length, defTP);
    Dialog.addRadioButtonGroup("Channel 1 (C1):", MARKERS, 1, MARKERS.length, defM1);
    Dialog.addRadioButtonGroup("Channel 2 (C2):", MARKERS, 1, MARKERS.length, defM2);
    Dialog.show();
	
	tp_		= Dialog.getRadioButton();
    m1_     = Dialog.getRadioButton();
    m2_     = Dialog.getRadioButton();

    if (m1_ == m2_) {
        showMessage("Error", "C1 and C2 must be different markers.");
        exit("Aborted: identical markers for C1 and C2.");
    }
    return newArray(tp_, m1_, m2_);
}




// ------------------------------
// Function for Background input dialog
// ------------------------------

function get_background_values() {
	Dialog.create("Background Values");
	Dialog.addMessage("Enter measured background values", 14);
	Dialog.addMessage("--------------------");
	Dialog.addMessage("12 h samples", 13);
	Dialog.addMessage("--------------------");
	
	Dialog.addMessage("HA568 + dsRNA488");
	Dialog.addNumber("HA568:", 0);
	Dialog.addNumber("dsRNA488:", 0);
	
	Dialog.addMessage("NS4B568 + dsRNA488");
	Dialog.addNumber("NS4B568:", 0);
	Dialog.addNumber("dsRNA488:", 0);
	
	Dialog.addMessage("NS4B568 + HA488");
	Dialog.addNumber("NS4B568:", 0);
	Dialog.addNumber("HA488:", 0);
	
	Dialog.addMessage("--------------------");
	Dialog.addMessage("24 h samples", 13);
	Dialog.addMessage("--------------------");
	
	Dialog.addMessage("HA568 + dsRNA488");
	Dialog.addNumber("HA568:", 0);
	Dialog.addNumber("dsRNA488:", 0);
	
	Dialog.addMessage("NS4B568 + dsRNA488");
	Dialog.addNumber("NS4B568:", 0);
	Dialog.addNumber("dsRNA488:", 0);
	
	Dialog.addMessage("NS4B568 + HA488");
	Dialog.addNumber("NS4B568:", 0);
	Dialog.addNumber("HA488:", 0);
	
	Dialog.show();
	
	
	// IMPORTANT: getNumber() exactly in same order as addNumber()
	
	HA568_in_HA568_dsRNA488_12h      = Dialog.getNumber();
	dsRNA488_in_HA568_dsRNA488_12h   = Dialog.getNumber();
	
	NS4B568_in_NS4B568_dsRNA488_12h  = Dialog.getNumber();
	dsRNA488_in_NS4B568_dsRNA488_12h = Dialog.getNumber();
	
	NS4B568_in_NS4B568_HA488_12h     = Dialog.getNumber();
	HA488_in_NS4B568_HA488_12h       = Dialog.getNumber();
	
	HA568_in_HA568_dsRNA488_24h      = Dialog.getNumber();
	dsRNA488_in_HA568_dsRNA488_24h   = Dialog.getNumber();
	
	NS4B568_in_NS4B568_dsRNA488_24h  = Dialog.getNumber();
	dsRNA488_in_NS4B568_dsRNA488_24h = Dialog.getNumber();
	
	NS4B568_in_NS4B568_HA488_24h     = Dialog.getNumber();
	HA488_in_NS4B568_HA488_24h       = Dialog.getNumber();

}

// ------------------------------
// Function to get Input_dir
// ------------------------------

function chooseInputDir() {
    INPUT_DIR = getDirectory("Choose folder with Mock + MOI .tif images");
    if (INPUT_DIR == "") exit("No folder selected.");
}

// ------------------------------
// Function to create Run_ID
// ------------------------------
function makeRunId() {
    getDateAndTime(y, m, dw, d, h, mn, s, ms);
    return "" + y + IJ.pad(m+1, 2) + IJ.pad(d, 2)
        + "_" + IJ.pad(h, 2) + IJ.pad(mn, 2);
        
// ------------------------------
// Function to create output folders
// ------------------------------

function buildOutputDir() {
    OUTPUT_DIR = INPUT_DIR + RUN_ID + "_substracted_images" + File.separator;
    if (!File.exists(OUTPUT_DIR)) File.makeDirectory(OUTPUT_DIR);
	mdPath = OUTPUT_DIR + File.separator + "background_values_used.md";
	createMDFile(mdPath);
}

// ------------------------------
// Function to create md file for values
// ------------------------------

function createMDFile (mdPath) {
	// alte Datei löschen, falls sie existiert
	if (File.exists(mdPath)) {
	    File.delete(mdPath);
	}
	
	// Markdown schreiben
	File.append("# Background values used\n\n", mdPath);
	File.append("Generated by Fiji macro.\n\n", mdPath);
	File.append("Macro versio used : " + MACRO_VERSION + "\n\n", mdPath);
	File.append("Run ID : " + RUN_ID + "\n\n", mdPath);
	
	File.append("## 12 h samples\n\n", mdPath);
	
	File.append("### HA568 + dsRNA488\n", mdPath);
	File.append("- HA568: " + HA568_in_HA568_dsRNA488_12h + "\n", mdPath);
	File.append("- dsRNA488: " + dsRNA488_in_HA568_dsRNA488_12h + "\n\n", mdPath);
	
	File.append("### NS4B568 + dsRNA488\n", mdPath);
	File.append("- NS4B568: " + NS4B568_in_NS4B568_dsRNA488_12h + "\n", mdPath);
	File.append("- dsRNA488: " + dsRNA488_in_NS4B568_dsRNA488_12h + "\n\n", mdPath);
	
	File.append("### NS4B568 + HA488\n", mdPath);
	File.append("- NS4B568: " + NS4B568_in_NS4B568_HA488_12h + "\n", mdPath);
	File.append("- HA488: " + HA488_in_NS4B568_HA488_12h + "\n\n", mdPath);
	
	
	File.append("## 24 h samples\n\n", mdPath);
	
	File.append("### HA568 + dsRNA488\n", mdPath);
	File.append("- HA568: " + HA568_in_HA568_dsRNA488_24h + "\n", mdPath);
	File.append("- dsRNA488: " + dsRNA488_in_HA568_dsRNA488_24h + "\n\n", mdPath);
	
	File.append("### NS4B568 + dsRNA488\n", mdPath);
	File.append("- NS4B568: " + NS4B568_in_NS4B568_dsRNA488_24h + "\n", mdPath);
	File.append("- dsRNA488: " + dsRNA488_in_NS4B568_dsRNA488_24h + "\n\n", mdPath);
	
	File.append("### NS4B568 + HA488\n", mdPath);
	File.append("- NS4B568: " + NS4B568_in_NS4B568_HA488_24h + "\n", mdPath);
	File.append("- HA488: " + HA488_in_NS4B568_HA488_24h + "\n\n", mdPath);
	
	print("Saved background documentation to:");
	print(mdPath);
}

// ------------------------------
// Function to create list of all images
// ------------------------------

function listMoOIFiles() {
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

// ============== 8. CLEANUP & UTILS =========================

function cleanupBetweenImages() {
    while (nImages > 0) { selectImage(nImages); close(); }
    if (isOpen("ROI Manager")) roiManager("reset");
    if (isOpen("Results"))     run("Clear Results");
    if (selectionType() != -1) run("Select None");
}

function inArray(val, arr) {
    for (i = 0; i < arr.length; i++) if (arr[i] == val) return true;
    return false;
}
