# Project Handoff — Colocalisation Pipeline (Fiji/IJM)

> Purpose of this file: hand off the full state of this project to a new chat
> session so work can continue without re-deriving context. Read this top to
> bottom before touching code.

---

## 0. Who I am / how I want to be taught

I'm Kolja, MSc student doing an internship in the **Bartenschlager lab**
(infectious diseases / virology, Heidelberg). I'm building Fiji macros for an
automated **colocalisation imaging pipeline**. The macros will later be used by
my **supervisor for a publication**, so this is research software:
reproducibility, documentation, and a clean handover matter.

**I communicate in German.** Respond in English.

**My learning goal:** I don't just want working code — I want to *become more
of an expert* in automated image analysis. When you change or write code, act as
an **image-analysis + Fiji + virology expert acting as my tutor**. Use evidence-
based teaching methods:
- Explain not just *what* but *why* and *what the mental model / thought process
  of an image analyst is*.
- Make reasoning explicit (decompose problems, name trade-offs, anticipate
  failure modes).
- Connect each decision to the biology and to the downstream analysis.
- Be didactic but token-efficient — high signal, scientific-standard quality.
- After non-trivial changes, give a short "lesson" section that generalises the
  insight so I can reuse it.

---

## 1. The biology & experiment

- Cells: **Huh7** (large) and **VeroE6** (~4× smaller, dimmer, less cytoplasm).
- Conditions: **Mock** (uninfected control), **MOI1**, **MOI5** (infected).
- Timepoints: **12h**, **24h**.
- Markers (immunofluorescence, confocal, 16-bit TIFs, single z-plane):
  `HA568`, `HA488`, `dsRNA488`, `NS4B568`. **DAPI is always channel 3.**
- Marker combos analysed (order matters: marker1=C1, marker2=C2):
  `HA568_dsRNA488`, `NS4B568_dsRNA488`, `NS4B568_HA488`.
- Signal intensities (raw, before bg subtraction):
  - Huh7: max values ~4500–20000 depending on marker.
    dsRNA strongest; NS4B & HA488 fairly strong; **HA568 weakest**.
  - Vero: similar, sometimes higher.
- Signal morphology: punctate (virus replication sites), perinuclear, small dots.
- ~30 Mock images per (cell line × timepoint), ~10 per marker combo. ~60+ MOI.

**Filename schema (binding, positional, underscore-separated):**
```
<timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<coverslip>_<imgIndex>.tif
e.g.  12h_Huh7_Mock_HA568_dsRNA488_CS1_1.tif
```
Token positions: 0=tp, 1=cellLine, 2=condition, 3=marker1(→C1), 4=marker2(→C2),
5=coverslip(CS#), 6=imgIndex.

---

## 2. The 3-script pipeline (big picture)

```
.czi → Bio-Formats → .tif (Mock + MOI mixed in folders)
   │
   ├── SCRIPT 1 (Mock): measure Top-X% pixel stats in cytosol of Mock images
   │        → CSVs of per-image background-estimate values
   │        → aggregate (median over images) → ONE bg value per
   │          (marker × combo × timepoint × cell line)
   │
   ├── SCRIPT 2 (SubBg+Coloc): for each MOI image, subtract those bg values,
   │        build a cytosol ROI, run Colocalisation Threshold, write coloc CSV
   │
   └── (downstream) aggregate coloc CSVs → stats / figures (not built yet)
```

The **manual link between Script 1 and 2**: after Script 1, I look at the
output CSVs, take the **median of the `p99_9995` column** (the recommended
percentile for the background value — high enough to sit above autofluorescence
noise, below real artefacts) per marker/combo/timepoint, and **type those values
into Script 2's startup dialog**.

---

## 3. Current file layout

```
macros/
├── HANDOFF.md                  ← this file
├── 1_Mock/
│   ├── 01_Mock_pipeline.ijm    ← CURRENT Mock script, MACRO_VERSION "0.6.0"
│   ├── CHANGELOG.md            ← version history (also on GitHub)
│   ├── ToDo.md
│   └── README.md               ← overview readme (overview only; USAGE comes later)
└── 2_SubBg_Coloc/
    ├── 01_Subtract_background_v0.1.ijm   (old)
    ├── 01_Subtract_background_v0.2.ijm   (old)
    ├── 02_Subtract_background_coloc_v0.4.ijm  (old, but the V0.4 ROI-FIRST logic was "good")
    ├── 05_Subtract_background_coloc_v0.5.ijm  (old, header V0.6 — subtract-FIRST + signal-OR ROI)
    └── 06_Subtract_background_coloc_v0.7.ijm  ← CURRENT (MACRO_VERSION "0.7.0")
```

GitHub repo: `Kollybook/colocalisation_project` (default branch `main`).
Versions are tracked by separate filenames + CHANGELOG.md. Data/output are NOT
committed (only code + docs).

---

## 4. SCRIPT 1 — Mock pipeline (`1_Mock/01_Mock_pipeline.ijm`, v0.6.0)

> v0.6.0 change: the domain-list dialog (markers / combos / timepoints) is now
> ALWAYS shown at startup, pre-filled with the workflow standards (plain OK =
> defaults). `MODE` (filename/dialog) now only governs the per-image metadata
> source, not list editability. Everything below still applies.


**Goal:** per Mock image, measure the Top-X% brightest pixel statistics inside
the cytosol of each marker channel → CSV rows for downstream bg estimation.

**Per-image pipeline:**
1. Parse filename (or dialog mode) → tp, cellLine, m1, m2.
2. Split channels, rename by marker (`HA568_channel`, etc.), DAPI=`DAPI_channel`.
3. Build masks (on RAW signal):
   - Cell mask: single priority channel (HA568>HA488>NS4B568>dsRNA488),
     `setAutoThreshold(Li dark)` scaled by `CELL_THR_FACTOR`, Fill Holes.
   - Nucleus mask: DAPI, Otsu, morphological closing (Dilate→Fill Holes→Erode).
   - Cytosol = Cell − Nucleus (8-bit subtraction clamps at 0).
4. Artefact cleanup (Mock-specific!):
   - Upper-bound threshold per channel (pixel > `ARTIFACT_UPPER_BOUND` = artefact),
     OR'd across both channels, subtracted from cytosol. Catches bright dust.
   - Particle filter: drop cytosol regions < `MIN_PARTICLE_SIZE` px (coverslip dirt).
5. Cytosol → ROI (Create Selection → ROI Manager).
6. Per marker: `getHistogram` inside ROI → `computeTopStats` → CSV row.
7. Save QC PNG (marker1 + cytosol outline), composite JPG (all channels + ROI),
   binary mask TIFs.

**Top-X% measurement (the statistical core, `computeTopStats`):**
- Uses the 65536-bin histogram (bin index == pixel value for 16-bit).
- Pass 1: walk from brightest down, accumulate until ≥ topPct% of pixels → that
  bin = `threshold_value` (where the top pool starts).
- Pass 2: mean + std of the top pool via weighted bin sums (exact for ints).
- Pass 3: median of top pool (cumulative to half).
- Pass 4: whole-cytosol percentiles p95 … p99_9999 (outlier detection).

**Key CONFIG vars (top of file):**
`MARKERS`, `ANALYSE_COMBI`, `TIMEPOINTS`, `CELL_THR_METHOD="Li"`,
`NUC_THR_METHOD="Otsu"`, `BLUR_SIGMA_CELL`, `CELL_THR_FACTOR=0.5`,
`ARTIFACT_UPPER_BOUND` (Huh7=2000, Vero=1000 via `applyCellLineDefaults`),
`MIN_PARTICLE_SIZE` (Huh7=200, Vero=100), `TOP_PCT=1.0`,
`SAVE_QC/SAVE_MASKS/SAVE_COMPOSITE_JPG`.

**Cell-line defaults** are set in `applyCellLineDefaults(cellLine)` AFTER the
startup dialog (so they react to the user's cell-line choice). Does NOT override
threshold *method*, only the "strictness knobs".

**Modes (startup dialog):** Marker source = filename | dialog;
Threshold mode = automatic | manual (manual pauses with Threshold slider).

**Output:** `<inputDir>/measure_mock/<RUN_ID>/{csv files, masks/, qc/, macro_log.txt}`.
12 CSVs = 2 tp × 3 combos × 2 markers. CSV columns include provenance
(`macro_version`, `run_id`, threshold params).

**Status:** working. ROI for Mock is fine because we only care about the bright
top pool (which is inside any reasonable cytosol mask). Vero needed tuning
(smaller particle size, Li not Triangle — Triangle gave near-zero thresholds on
dim Vero data).

---

## 5. SCRIPT 2 — SubBg + Coloc (`2_SubBg_Coloc/06_..._v0.7.ijm`, V0.7.0)

> **V0.7 is the CURRENT script.** It replaces the V0.6 full-image signal-OR ROI
> with a **single-cell ROI** and reverts to the **ROI-FIRST ordering** of V0.4
> (build ROI on raw channels BEFORE bg subtraction → ROI is signal-independent
> = unbiased; in manual mode the cell stays visible via autofluorescence).
> Per-image order: `split → build ROI → subtract bg → coloc → merge → save`.
>
> **ROI_MODE (startup radio):**
> - `auto_central` (default): nucleus mask (DAPI) → all nuclei as ROIs → for each,
>   grow by `NUC_DILATE_UM` (→ pixels via pixel size, `Enlarge ... pixel`) to a
>   territory; **skip territories touching the image border** (partially-imaged
>   cells); among the rest pick the nucleus nearest the image centre; cytosol =
>   territory **XOR** nucleus → ROI Manager index 0 ("Cytosol_central").
> - `manual_draw`: pause per image, user draws a freehand ROI on the raw m1
>   channel (Thomas' workflow) → ROI Manager index 0 ("Cytosol_manual").
> - Both end at ROI index 0, so the coloc step (`use=[Channel 1]`) is unchanged.
>
> **New CONFIG:** `ROI_MODE`, `NUC_DILATE_UM` (Huh7 8 / Vero 5 µm), `BORDER_MARGIN_PX`,
> `NUC_MIN_SIZE` (Huh7 1000 / Vero 500 px), `SAVE_COLOC_QC`.
> **New output:** `qc_<img>_roi.jpg` (composite of all channels + ROI outline).
> **New CSV columns:** `roi_mode`, `nuc_dilate_um`.
> **Feature A:** domain-list dialog always shown at startup (pre-filled standards).
>
> The text below describes the older V0.6 (file 05) for reference.

**Goal:** for each MOI image: subtract per-channel bg, build a cytosol ROI, run
Colocalisation Threshold, parse its output into a CSV.

**Background values:** entered ONCE at startup via a dynamic dialog
(TIMEPOINTS × ANALYSE_COMBI × 2 markers). Stored in a parallel-array lookup
table (`BG_KEYS`/`BG_VALUES`, key = `<marker>_in_<combo>_<tp>`). Also written to
`background_values_used.md` for provenance.

**Per-image pipeline (CURRENT, V0.6 — being actively tuned):**
1. Split + rename channels.
2. **Subtract bg FIRST** (m1, m2; DAPI untouched). ← changed in V0.6 from "mask
   first" because thresholds were weird; bg-subtracted histogram is clean
   bimodal (0-peak + signal tail).
3. Build cell mask from BOTH bg-subtracted channels:
   `makeChannelSignalMask` (Triangle, blur 0.5–2 µm) on each → OR'd →
   Fill Holes → Dilate. (`makeCellMaskFromChannels`)
4. Nucleus mask from DAPI; Cytosol = Cell − Nucleus.
5. Particle filter `MIN_PARTICLE_SIZE=20` (LOW, to keep small punctate signal).
6. Cytosol → ROI Manager (named "Cytosol").
7. Coloc step (modes: none | manual | auto):
   - manual: preview → decision → pause for user to run plugin.
   - auto: `run("Colocalization Threshold", ...)` then parse Log.
8. Merge channels back (Color hyperstack), save 16-bit TIF, optional JPG.

**Coloc plugin = "Colocalization Threshold"** (Tony Collins / WCIF), NOT Coloc 2.
- Macro call (from recorder): `run("Colocalization Threshold",
  "channel_1=<t1> channel_2=<t2> use=[Channel 1] channel=[Red : Green] include")`.
- NOTE the spelling: menu shows "Colocali**s**ation" (British) but `run()` needs
  "Colocali**z**ation" (American).
- `use=[Channel 1]` reads the ACTIVE selection on channel-1 window → we pre-apply
  the cytosol ROI to m1_channel so the user just picks "Use ROI: Channel 1".
- Output goes to the **Log window** (not Results table). We snapshot log length
  before/after and parse the diff (`extractLogValue` + `tryLabels` with multiple
  fallback label strings).
- 15 values to parse into CSV: `Rtotal, m, b, Ch1_thresh, Ch2_thresh, Rcoloc,
  R_below_thresh, M1, M2, tM1, tM2, Ncoloc, perc_volume, perc_ch1_vol,
  perc_ch2_vol`.

**Output:** `<inputDir>/<RUN_ID>_bgsub/{bgsub_*.tif, *.jpg, background_values_used.md,
coloc_results_<RUN_ID>.csv, macro_log.txt}`.

---

## 6. THE ROI QUESTION — RESOLVED in V0.7

> **DECISION (implemented in V0.7):** single representative cell, ROI built
> signal-independently. `auto_central` = nucleus-nearest-centre with a
> non-border-touching DAPI-dilated territory (cytosol = territory − nucleus);
> `manual_draw` = freehand around the central cell (matches the supervisor's
> example image). ROI is built BEFORE bg subtraction (unbiased + visible).
> Per-cell granularity = exactly one cell per image (no full-image bias).
> **Still TODO:** confirm the definition with Thomas; verify auto-coloc log labels.

The original discussion (three approaches) is kept below for context.

**How should the ROI for Colocalisation Threshold be built?** Coloc values
currently look "weird". Discussed three approaches:

1. **Signal-based** (current V0.6): threshold both channels, OR → ROI = where
   signal is. **Problem: biased** — measuring coloc only where there's signal
   inflates Pearson R; methodically circular for publication values.
2. **DAPI-dilation** (RECOMMENDED): Nucleus mask → dilate ~8–10 µm → cell
   territory; cytosol = territory − nucleus. **Unbiased** (ROI independent of the
   measured signal), robust to weak HA568, consistent across combos. This is the
   standard, defensible approach for publication; what the supervisor likely
   expects. Reviewer-proof answer: "nuclear stain dilated ~10 µm, nucleus
   subtracted".
3. **Per-cell Watershed/Voronoi**: highest resolution (per-cell stats), good for
   dense layers (the images ARE dense), but complex and 10× more CSV rows.

**My recommendation was Option 2 (DAPI-dilation) as default, with signal_or as a
fallback for comparison runs, and optional Voronoi for dense cultures.**
Proposed config: `CELL_MASK_STRATEGY = "dapi_dilation" | "voronoi" | "signal_or"`,
`NUC_DILATE_UM = 8`. **NOT yet implemented.** Also: I want to discuss the ROI
question with Thomas (supervisor / "mit thomas besprechen wie ROI aussehen muss").

Next steps for Script 2:
- [ ] Implement DAPI-dilation ROI strategy (the proposed V0.7).
- [ ] Verify the auto-coloc Log parsing labels against the actual plugin output
      (run once, check the `WARN: no coloc labels parsed` raw dump if empty).
- [ ] Decide per-cell vs per-image granularity.

---

## 7. IJM PITFALLS LEARNED (hard-won — respect all of these)

| Pitfall | Symptom | Fix |
|---|---|---|
| Ternary `a ? b : c` | `')' expected`, `<?>` in Debug | use `if/else` |
| `i++`, `x += 1` | syntax error | `i = i + 1` |
| `return userFunc(...)` (returning another user fn's result directly) | "Array/Numeric expected", value reads as `0` | intermediate var: `r = f(); return r;` |
| `replace(s, regex)` | nothing replaced (no regex!) | `substring(s,0,lastIndexOf(s,"."))` |
| `var` inside a function | shadows the global | declare `var` only at top |
| `File.makeDirectory` not recursive | "file not found" on saveAs | create each level top-down |
| `getHistogram(v,c,65536,0,65535)` on 16-bit | `values[]` stays scalar `0` | use bin index `i` as the value |
| `getHistogram(v,c,n,min,max)` on 8-bit | "16 or 32-bit required" | use 3-arg form `getHistogram(v,c,256)` |
| `run("Subtract...")` | clamps at 0 (no negatives) | intended for bg-sub; be aware |
| **Black Background** global flips | mask inverts (fg=0) after iteration 1 | `setOption("BlackBackground", true)` at MAIN + per iteration + after Merge/RGB |
| Plugin name UI vs macro | "command not found" | use Plugins>Macros>Record, copy verbatim |
| `run()` arg with space | "missing parameter" | bracket it: `use=[Channel 1]` |
| `Dialog.create` is MODAL | Fiji GUI fully blocked | use `waitForUser` for non-modal pauses |
| `Array.concat` onto empty `newArray()` | numbers stored as STRINGS → `"1570"*1` fails "; expected" | `parseFloat()` at store AND retrieve |
| `Merge Channels` without `create` | 8-bit RGB, dynamic range lost | always add `create` for quantitative data |
| `Merge Channels` consumes sources | windows gone for later steps | add `keep`, OR duplicate sources first |
| `Dialog.getX` order | values silently swapped | get* order MUST mirror add* order (use identical loops) |

**Debugging method that works:** when a macro errors, the Debug window shows the
exact token (`<?>`, `<*>`, `<]>`) where the parser choked, plus all variable
values. Always read the variable types (string `"5"` vs number `5`) — many bugs
are type-inference issues.

---

## 8. CODE STYLE / ARCHITECTURE CONVENTIONS (keep consistent)

- **CONFIG block at top** (all `var` globals). Never tune inside functions.
- **MAIN is short**: pure orchestration, one line per phase.
- **One function = one job**, each with a `// ----` docstring (what / inputs /
  returns / assumptions).
- **Domain lists** (`MARKERS`, `ANALYSE_COMBI`, `TIMEPOINTS`) are the single
  source of truth — CSV creation, dialogs, validation all derive from them.
- **Provenance in every CSV row**: `macro_version`, `run_id`, threshold params.
- **`run_id`** = `YYYYMMDD_HHMM`, also names the output subfolder → parallel runs
  never overwrite; pandas can filter by `run_id`.
- **Defensive validation**: skip+print (not exit) on bad combo/timepoint so one
  bad image doesn't kill the batch.
- **Cell-line-specific params** via `applyCellLineDefaults()` (strictness knobs
  only, not strategy).
- **Filter composition**: artefacts decomposed by discriminator
  (intensity → upper-bound threshold; geometry → particle filter). Order matters:
  subtract artefacts BEFORE particle filter (subtraction can fragment regions).
- **Mock vs MOI difference**: Mock uses an upper-bound artefact cut (only
  autofluorescence, bright = dirt). MOI must NOT (real specific signal is bright).

---

## 9. KNOWN TUNING KNOBS & current values

Script 1 (Mock): `CELL_THR_METHOD=Li`, `CELL_THR_FACTOR=0.5`, `TOP_PCT=1.0`,
`ARTIFACT_UPPER_BOUND` 2000(Huh7)/1000(Vero), `MIN_PARTICLE_SIZE` 200/100.

Script 2 (MOI, V0.7): `ROI_MODE=auto_central` (or `manual_draw`),
`NUC_DILATE_UM` Huh7 8 / Vero 5, `BORDER_MARGIN_PX=0`, `NUC_MIN_SIZE` Huh7 1000 /
Vero 500, `NUC_THR_METHOD=Otsu`, `NUC_CLOSE_ITER=2`, `COLOC_MODE=manual`,
`SAVE_COLOC_QC=true`. (V0.6's `CELL_THR_*` / signal-OR knobs are gone.)

Recommended bg percentile to take from Script 1 output: **p99_9995**, median over
images per marker/combo/tp.

---

## 10. Immediate ToDo (from ToDo.md + our discussion)
- [x] Single-cell ROI in Script 2 (V0.7): `auto_central` (DAPI-dilation, border-excluded,
      nearest-centre) + `manual_draw`. ROI built BEFORE subtraction. DONE.
- [x] Startup parameter entry always shown in both scripts (Feature A). DONE
      (Mock v0.6.0, Coloc v0.7.0).
- [ ] **Test V0.7 in Fiji** (auto_central + manual_draw + border-exclusion + QC) — not yet run.
- [ ] Discuss the single-cell ROI definition with Thomas (supervisor).
- [ ] Verify/fix coloc Log-parsing labels against real plugin output.
- [ ] Annotate both scripts thoroughly.
- [ ] Write an overview README for both scripts (Mock README exists; SubBg needs one).
- [ ] Later: detailed USAGE.md (pitfalls, step-by-step) — separate from overview README.
```
