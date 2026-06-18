# CHANGELOG for Subtract_background_coloc.ijm

Background-subtraction + colocalisation-preparation pipeline for FIJI/ImageJ.
Takes the per-channel background values measured by the Mock pipeline, subtracts
them from infected (MOI) images, builds a cytosol ROI, and prepares / runs the
colocalisation analysis.

---

| Version | Short description                                            |
| ------- | ------------------------------------------------------------ |
| v0.0.1. | First conceptual draft (vibe-coding). Sketches the intended workflow: choose folder, enter background values per marker/combo/timepoint, list MOI images, split channels, subtract background. NOT functional — contains syntax errors and incomplete logic; serves as the initial blueprint. |
| v0.1.0  | First functional version. Background subtraction now actually works end-to-end: scalable background lookup table, dynamic value-entry dialog, per-channel subtraction, channel re-merge into a 16-bit composite, saving of bg-subtracted TIFs, plus markdown + full IJ-log provenance per run. |
| v0.2.0  | First **coloc-prep** version. Builds a generous cytosol ROI from the cell-mask channel on the RAW signal **before** subtracting background, subtracts bg from C1/C2, re-merges into a scrollable multi-channel hyperstack, saves a 16-bit TIF, then **pauses for the user to run Coloc 2** on the two non-DAPI channels inside the ROI; appends a placeholder row to a per-run coloc CSV (the real CSV bridge planned for v0.4). |
| v0.3.0  | **Version-bump checkpoint** — code identical to v0.2 with `MACRO_VERSION` → `0.3.0`; the cytosol ROI + manual Coloc-2 hand-off + placeholder CSV are unchanged. *(Was briefly on disk as `..._v0.4.ijm`; renamed to `..._v0.3.ijm` to match its header/version.)* |
| v0.4.0  | *(no source file retained — reconstructed from later headers.)* Established the **"ROI-FIRST" ordering** that all later versions keep (build the cytosol ROI on the raw channels **before** background subtraction, so the ROI is signal-independent = unbiased), a **single priority-picked channel** cell mask auto-thresholded with `Li` on the raw autofluorescence, and the **first working Coloc → CSV bridge**. Per v0.5's header, the **bg-value-based** and **combined C1+C2** cell-mask experiments were tried and reverted around here. |
| v0.5.0  | Added the **particle-size filter** on the cytosol mask (drop fragments < `MIN_PARTICLE_SIZE`), **cell-line-specific defaults** (`applyCellLineDefaults`, Vero ≠ Huh7), **dialog-mode startup config** (markers/combos/timepoints/cell line overridable), and **auto-coloc that parses the plugin's Log** into real CSV values. Fixes: `BlackBackground` pinned at MAIN + per iteration (mask-inversion bug); `parseFloat()` at bg store + retrieve (the "Array.concat stores numbers as strings" IJM quirk). |
| v0.6.0  | **Subtract background FIRST, then build the mask** (the bg-subtracted histogram is a clean 0-peak + signal tail). Cell mask is now the **OR of both channels' signal masks** (`makeChannelSignalMask` per channel, `Triangle`, blur 2 µm, `MIN_PARTICLE_SIZE` 20) instead of a single autofluorescence channel — captures signal in either marker. Coloc Log parser made **robust** (`tryLabels` tries many alternative label spellings per value, and dumps the raw Log chunk when nothing parses). *Note: this signal-OR ROI is **biased** — it measures coloc only where there is signal — which motivated the v0.7 single-cell ROI.* |
| v0.7.0  | First **single representative-cell ROI** (replaces the biased v0.6 signal-OR full-image ROI). Two strategies chosen at startup (`ROI_MODE`): **`auto_central`** (the nucleus nearest the image centre whose dilated, non-border-touching territory − nucleus = cytosol) or **`manual_draw`** (freehand around the central cell). ROI is built on the **raw** channels **before** subtraction (back to the v0.4 ROI-FIRST order → unbiased + cell visible while drawing). Adds a **per-image coloc QC JPG** (all channels + the actual ROI burnt in). *Known problem: the central cell is sometimes uninfected (fixed by the v0.8 infection gate).* |
| v0.8.0  | Auto cytosol redesigned (infection gate + autofluorescence cell mask ∩ nucleus-Voronoi split − nucleus); manual ROI overhaul (composite drawing, multiple ROIs/image, auto nucleus exclusion); coloc values now read from the plugin's Results table instead of the IJ Log. |
| v0.9.0  | Automatic ROI is now MULTI-CELL via MorphoLibJ marker-controlled watershed (nuclei seed one labelled region per cell — fixes fragmentation/random-area). Keeps every cell ≥ CELL_MIN_SIZE px AND infected (p99 both channels), measures them all; optional nucleus subtraction; MorphoLibJ existence+version gate. Flexible filename tokens (askTokenMapping) in both scripts. |
| v0.10.0 | Fixes: manual mode shows already-drawn ROIs; token dialog uses radio buttons (NaN bug fixed); tif/jpg/qc in their own sub-folders; coloc results no longer parsed (write provenance row only). |
| v0.11.0 | Dropped MorphoLibJ: the marker-controlled watershed got a BINARY marker and merged touching cells into one ROI. Replaced by a **pure-IJM nucleus-Voronoi split** of a solid cell mask (cuts only between cells; one ROI per cell, no fragmentation). New keep rule: area > CELL_MIN_SIZE px AND raw-channel p99 ≥ SIGNAL_BG_MULT (=2) × entered bg, both channels. |
| v0.12.0 | Auto ROI rewritten as **"grow each nucleus inside the cell mask"**: Voronoi cut (between-cell boundaries) + **flood-fill (doWand) from each nucleus centroid**. v0.11 ran Analyze Particles on the cut pieces → many tiny non-cell fragments on real data. v0.12 anchors every ROI to exactly one nucleus (no fragments), isolated nuclei grab their whole blob, touching cells stop at the ridge. Keep rule unchanged. |
| v0.13.0 | **Branched from v0.10** (MorphoLibJ marker-controlled watershed, the preferred segmentation) — v0.11/v0.12 Voronoi kept only as history. Single behavioural change: the **infection gate is now relative and nucleus-centred**. Each nucleus is grown ~`GATE_RING_UM` (4 µm) into the cytoplasm, clipped to its own watershed territory; the cell is kept if the **raw p99.9 in that ring ≥ `SIGNAL_BG_MULT` (2) × the entered background in BOTH channels**. Replaces the absolute `INFECT_MIN_P99 = 100` test on bg-subtracted data (which dropped the weak HA568 channel → empty QC). ROI boundary unchanged (still pure morphology). |
| v0.14.0 | **Issue B fixed — touching cells now separate.** The watershed was fed the **binary** nucleus mask as its marker, so MorphoLibJ merged neighbouring cells into one region. The marker is now a **label image** (Connected Components Labeling → one integer per nucleus), so each nucleus seeds its own basin. **QC ROIs are numbered** (number == `roi_index` in the CSV) for both auto + manual, so a coloc row can be matched to its cell. **Manual UX:** the "exclude nucleus?" question is asked **once at startup**, not per image. |
| v0.14.1 | Field tuning (Fiji): `CELL_MIN_SIZE` 200 → **80000 px²** (real Huh7/Vero cells are large; drops watershed crumbs), gate percentile `GATE_PCTL` 99.9 → **99.99**. |
| v0.14.2 | Field tuning (Fiji): `ROI_MODE` renamed `auto_central` → **`auto_watershed`** (dialog + internals), `CELL_MIN_SIZE` 80000 → **120000 px²**, **token-mapping dialog read-order fixed** (timepoint/marker1/marker2 were read in the wrong order), and dropped cells now **log their area vs the minimum**. |
| v1.0.0 (current) | Field tuning (Fiji): **merged `AUTO_EXCLUDE_NUC` + `MANUAL_EXCLUDE_NUC` → one `EXCLUDE_NUC`** (one checkbox governs both modes), `NUC_MIN_AREA` 100 → **5000 px** (only large, real nuclei seed a cell — drops debris and tiny mitotic fragments), `SIGNAL_BG_MULT` 2 → **1.5** (keep weak-but-real HA568 cells), `GATE_RING_UM` 4 → **8 µm** (ring scales with the large cells). *(The shipped file initially carried the v0.14.2 header/version string by copy-paste; corrected to 0.14.3.)* |

---

### v1.0.0 (current)
- **Merged `AUTO_EXCLUDE_NUC` and `MANUAL_EXCLUDE_NUC` into a single `EXCLUDE_NUC`.** Having two
  separate "exclude nucleus?" flags (one for auto, one set in the once-only manual dialog) was
  redundant and could disagree; now one global + one checkbox label "Exclude nucleus from each
  cell ROI" governs both ROI modes.
- **`NUC_MIN_AREA` 100 → 5000 px.** Only nuclei ≥ 5000 px seed a watershed cell. At 100 px every
  speck of DAPI debris and every fragment of a mitotic/pyknotic nucleus seeded a spurious basin
  (→ extra fragmented ROIs). 5000 px keeps only whole interphase nuclei. *Re-scale with pixel
  size if magnification changes.*
- **`SIGNAL_BG_MULT` 2 → 1.5**, **`GATE_RING_UM` 4 → 8 µm** (the ring is grown further into the
  now-large cells so the gate samples cytoplasmic puncta, not just the perinuclear rim).
- **Provenance fix:** the file was created by copy from v0.14.2 and still reported
  `MACRO_VERSION = "0.14.2"` (so a run's CSV `macro_version` column would have been wrong);
  corrected the header title and `MACRO_VERSION` to **0.14.3**.

---

### Current auto-ROI defaults (v1.0.0)
| Knob | Value | Meaning |
|---|---|---|
| `ROI_MODE` | `auto_watershed` | nucleus-seeded MorphoLibJ watershed |
| `NUC_MIN_AREA` | 5000 px | min nucleus area to seed a cell |
| `CELL_MIN_SIZE` | 120000 px² | min kept-ROI area |
| `SIGNAL_BG_MULT` | 1.5 | gate: ring p99.99 ≥ 1.5 × bg, both channels |
| `GATE_RING_UM` | 8 µm | nucleus growth for the gate ring |
| `GATE_PCTL` | 99.99 | percentile measured in the ring |
| `EXCLUDE_NUC` | true | subtract the nucleus from each ROI |

### To verify / watch (in Fiji)
- Pixel-based knobs (`CELL_MIN_SIZE`, `NUC_MIN_AREA`) are tuned for the **current magnification**
  — recompute them if the pixel size changes (they are raw px / px², not µm).
- If too few cells survive: lower `CELL_MIN_SIZE` or `SIGNAL_BG_MULT`. If watershed crumbs or
  uninfected cells slip through: raise them. Read the per-cell `kept/drop` Log lines to decide
  which filter is firing.

---

## v0.14.1 – v0.14.3  (field tuning, edited live in Fiji)

> [!NOTE]
> These are **parameter / naming refinements** made while running v0.14 on real Huh7 and
> VeroE6 data — no algorithm change. The segmentation (labelled-marker watershed), the
> relative nucleus-ring infection gate, the numbered QC and the once-only manual question
> are all unchanged from v0.14.0. Current file: `02_Subtract_background_coloc_v0.14.3.ijm`.

### Why these knobs moved (the mental model)
Two parameters drive almost all of the auto-ROI behaviour, and v0.14.0's defaults were
inherited from earlier, smaller-ROI logic:

- **`CELL_MIN_SIZE`** (the area below which a cell ROI is discarded). v0.14.0 shipped `200 px²`,
  a leftover from when ROIs were small fragments. A real confluent Huh7/VeroE6 cytosol is
  **tens to hundreds of thousands of px²**, so `200` kept everything — including watershed
  crumbs and partial cells. Raised in two steps (80000 → **120000 px²**) to keep only
  whole, well-segmented cells. *This is the single most important auto-ROI knob — re-check it
  for any new magnification or cell line (it is in raw pixels, not µm²).*
- **`SIGNAL_BG_MULT`** (the infection gate: keep a cell iff its nucleus-ring `GATE_PCTL`
  ≥ this × the entered background, both channels). Lowered `2 → 1.5` because the weak **HA568**
  channel produced real-but-dim infected cells that a 2× cut dropped. 1.5 = "signal ≥ 50 %
  above the Mock background".

---

### v0.14.2
- **Renamed `ROI_MODE` value `auto_central` → `auto_watershed`** everywhere (config default,
  the "ROI strategy" dialog label `automatic watershed ROI`, all internal `if (ROI_MODE …)`
  checks, the `writeBgMarkdown` provenance). Purely a clarity rename — the mode IS the
  watershed now, not the old "central cell".
- `CELL_MIN_SIZE` 80000 → **120000 px²**.
- **Fixed the token-mapping dialog read order.** `askTokenMapping()` added the radio groups in
  the order timepoint → marker1 → marker2 but read them back in a different order, so the three
  `Dialog.getRadioButton()` results could be misassigned. The get-order now mirrors the
  add-order (HANDOFF §7: "Dialog.getX order must mirror addX order"). Defaults are still
  marker1 = first token, marker2 = last token, timepoint = first token.
- **Dropped cells now log their area**: `drop cell #i: area=… (>=CELL_MIN_SIZE) …` so you can
  see from the Log whether the area filter or the infection gate rejected a cell.

---

### v0.14.1
- `CELL_MIN_SIZE` 200 → **80000 px²**.
- `GATE_PCTL` 99.9 → **99.99** (sample the very brightest puncta in the ring; with the larger
  ring a slightly higher percentile still lands on real signal, not noise).

---

## v0.14.0

> [!NOTE]
> Current version (`02_Subtract_background_coloc_v0.14.ijm`). Same MorphoLibJ
> dependency as v0.13. Code written, not yet validated in Fiji.

### Issue B — touching cells now SEPARATE (the main fix)
v0.10–v0.13 fed `watershedCellsFromNuclei()` the **binary** `Nuc_Mask` as the watershed
marker. With a binary marker MorphoLibJ does not run a true marker-controlled watershed; it
collapses to **connected components of the cell mask**, so two cells that touch (one blob in
the autofluorescence mask) become **one region** — exactly the "connected cells not
separated" failure. Fix:
- `Connected Components Labeling` (connectivity 8, 16-bit) turns `Nuc_Mask` into a **label
  image** `Nuc_Labels` (one distinct integer per nucleus).
- The watershed is now seeded by `marker=Nuc_Labels`, so **each nucleus floods its own basin**
  over the distance-to-nucleus surface and neighbours are split at the watershed line.
- Log prints `nucleus seeds (labelled): N`.

### QC ROIs are numbered
`saveColocQc()` now draws each ROI's index on the flattened QC JPG (yellow, after the outline
overlay). The number equals `roi_index` in `coloc_results_*.csv` (manager position + 1), so a
colocalisation row can be matched back to its cell — useful for spotting / excluding an
uninfected cell after the fact. Applies to both automatic and manual ROIs.

### Manual mode UX
The "Exclude nucleus from drawn ROIs?" question is asked **once at startup**
(`askManualRoiOptions()` → global `MANUAL_EXCLUDE_NUC`), called from MAIN before the image
loop, instead of popping up for every image.

### Added
- `askManualRoiOptions()`, global `MANUAL_EXCLUDE_NUC`.
- Nucleus labelling in `watershedCellsFromNuclei()` (`Nuc_Labels`, closed in
  `closeAutoMaskWindows`).
- ROI-number drawing block in `saveColocQc()`.

### To verify (in Fiji)
- `Connected Components Labeling` and `Marker-controlled Watershed` option strings match your
  MorphoLibJ build (confirm with `Plugins ▸ Macros ▸ Record`).
- Touching infected cells now produce **separate** numbered ROIs (not one merged outline).
- The QC numbers line up with the `roi_index` column in the CSV.
- The manual nucleus question appears once, before the first image.

---

## v0.13.0

> [!NOTE]
> Current version‚ (`02_Subtract_background_coloc_v0.13.ijm`). **Branched from v0.10**,
> not v0.12: it uses the **MorphoLibJ marker-controlled watershed** segmentation
> (preferred over the v0.11/v0.12 pure-IJM Voronoi). Requires the **IJPB-plugins**
> update site. Only the infection gate changed; code written, not yet validated in Fiji.

### Why branch from v0.10?
On real data the watershed (v0.10) gives cleaner per-cell territories than the Voronoi
cut, *when it fires*. Its two field failures — "no ROI shown" and "touching cells not
separated" — were traced to (a) an over-strict absolute infection gate dropping every
cell, and (b) the watershed receiving a **binary** marker. This version fixes (a); (b)
(labelled markers so touching cells split) is the **next** planned change (issue B), not
included here.

### Major change — relative, nucleus-centred infection gate
v0.10 kept a cell iff the **p99 of the background-subtracted channel** inside the whole
cell ROI was ≥ an **absolute** `INFECT_MIN_P99` (= 100), in both channels. Problems:
an absolute cut does not transfer between markers / cell lines (the weak **HA568** was
wrongly dropped → empty QC, "no ROI"), and measuring inside the whole watershed cell
couples the gate to the segmentation quality.

Now the gate is **relative** and **nucleus-centred**:
- For each nucleus, grow it by `GATE_RING_UM` (≈ 4 µm) into the cytoplasm and **clip the
  disk to that cell's own watershed territory** (a neighbour's bright puncta cannot leak in).
- In that per-cell region, measure the **raw** (pre-subtraction) `GATE_PCTL` (= 99.9)
  percentile of **both** channels.
- **Keep** the cell iff both ≥ `SIGNAL_BG_MULT` (= 2) × the background value entered for
  that channel — "real signal sits ≥ 2× above the Mock background you measured".
- Log prints, per nucleus, `kept`/`drop` with the measured p99.9 vs each threshold.

### Added
- CONFIG `SIGNAL_BG_MULT`, `GATE_RING_UM`, `GATE_PCTL`.
- `ringUmToPx()` (µm → dilation px via image calibration), `buildNucleusRing()`
  (per-cell nucleus ring, clipped to territory → `Gate_Mask`), `activeSelectionPercentile()`
  (16-bit percentile on the active selection of a channel).

### Removed
- `INFECT_MIN_P99`; `makeBgSubDuplicate()` and the `m1_sub` / `m2_sub` duplicates (the gate
  now reads the raw channels directly); `roiPercentile()` (superseded by
  `activeSelectionPercentile()`).

### To verify (in Fiji)
- MorphoLibJ installed (the `ensureMorphoLibJ()` gate enforces it).
- Images that previously showed **no ROI** now keep their infected cells; the Log shows
  per-nucleus `kept`/`drop` p99.9 values — tune `SIGNAL_BG_MULT` (weak HA568) and
  `GATE_RING_UM` from those numbers.
- **Known open issue (next change):** touching cells may still share one ROI until the
  watershed marker is made *labelled* instead of binary.

---

## v0.12.0

> [!NOTE]
> Current version (`02_Subtract_background_coloc_v0.12.ijm`). Pure ImageJ macro, no plugin.

### Major change — "grow each nucleus inside the cell mask"
v0.11 cut the solid cell mask by the nucleus-Voronoi ridges and then ran `Analyze Particles` on
the pieces. On real Vero data (patchy autofluorescence, confluent sheets) this produced **many
tiny fragments that are not cells**, and most images kept 0 (see the test log: "0 of 15 cells"
etc.). The split itself was fine; the per-piece labelling was the problem.

v0.12 keeps the Voronoi cut (the between-cell boundaries) but changes how territories are read out:
- **Cell mask is OR'd with the nuclei** so every nucleus is a valid seed.
- For **each nucleus**, the macro **flood-fills its cell territory** with `doWand(centroid, …,
  "8-connected")` on `Cells_Split`. So:
  - every ROI is anchored to exactly **one nucleus** → no nucleus-less fragments;
  - an **isolated** nucleus grabs its **whole** mask blob (no spurious cut);
  - **touching** cells stop at the Voronoi ridge between them — the hand-drawn-style separation.
- Then nucleus optionally subtracted; keep iff area ≥ `CELL_MIN_SIZE` AND raw p99 ≥
  `SIGNAL_BG_MULT × bg` in both channels. Log: `kept nucleus #i: area=… p99(m1)=… p99(m2)=…`.

### Changed
- `buildAutoCytosolRoi` rewritten (flood-fill per nucleus instead of Analyze-Particles per piece).
- `closeAutoMaskWindows` also closes `cell_tmp`.

### To verify (in Fiji)
- No more tiny non-cell ROIs; each kept ROI is one nucleus's cytosol.
- Isolated infected cells = whole-cell ROI; touching infected cells = separated at the ridge.
- `auto ROIs kept: K of N nuclei` in the log; tune `SIGNAL_BG_MULT` for the weak HA568 channel.

---

## v0.11.0

> [!NOTE]
> Current version (`02_Subtract_background_coloc_v0.11.ijm`). **No plugin dependency** — the auto
> mode is pure ImageJ macro again. Code written, not yet validated in Fiji.

### Major change — automatic per-cell ROI without MorphoLibJ
v0.9/v0.10 used MorphoLibJ marker-controlled watershed. In the field it passed a **binary** marker,
so the plugin collapsed to connected-components and **merged touching cells into one ROI** (QC
showed one outline over several cells). Replaced by a **pure-IJM nucleus-seeded split**:
- `makeCellMask` builds a **solid** cell mask (blur → threshold → morphological close → fill holes).
- `splitCellsByNucleusVoronoi()` takes the **Voronoi of the nuclei** — its ridges are the
  boundaries between neighbouring cells — and subtracts (dilated) ridges from the cell mask,
  **cutting only between cells**. Each cell stays one connected blob → **one ROI per cell, no
  within-cell fragmentation** (the v0.8 failure mode is avoided by the solid mask + cut-between).
- `Analyze Particles` → one ROI per cell; per cell the nucleus is optionally subtracted and the
  cytosol captured as a **single** selection (a ring stays one ROI).

### New keep rule (replaces the bg-subtracted p99 ≥ 100 gate)
A cell is measured iff **area > `CELL_MIN_SIZE` px** AND its **raw-channel 99th percentile ≥
`SIGNAL_BG_MULT` (default 2) × the entered background**, in BOTH channels. Because channels are
still raw when the ROI is built, this is literally "signal ≥ 2× the bg you typed". Log prints
`kept cell #i: p99(m1)=… (>=…) p99(m2)=… (>=…)`.

### Removed
- `watershedCellsFromNuclei`, `ensureMorphoLibJ`, `MORPHOLIBJ_TESTED_VERSION`, `newestImageNotIn`,
  `makeBgSubDuplicate`, the `m1_sub`/`m2_sub` duplicates, and `INFECT_MIN_P99` (→ `SIGNAL_BG_MULT`).

### To verify (in Fiji)
- Touching cells now split into separate ROIs; each infected cell → its own ROI + CSV row.
- Images that previously showed no ROI now show their infected cells; tune `SIGNAL_BG_MULT` if the
  weak HA568 channel wrongly drops cells.

---

## v0.9.0

> [!NOTE]
> Current version. Adds a dependency: **MorphoLibJ (IJPB-plugins update site)**, required only for
> automatic ROI mode. Code written, NOT yet validated in Fiji.

### Major changes
- **Automatic ROI is now MULTI-CELL and watershed-based.** The v0.8 approach (autofluorescence cell
  mask cut by nucleus-Voronoi, then *connected components*, keep one central cell) fragmented cells,
  picked tiny/random areas, and only ever measured one cell. v0.9 replaces it with a **nucleus-seeded
  MorphoLibJ marker-controlled watershed**, which yields **exactly one labelled region per nucleus**
  (no fragmentation), then keeps and measures **every** qualifying cell.
- **Flexible filename tokens** in both scripts (see Mock CHANGELOG v0.8.1 for the Mock side).

### Why (the v0.8 failures, from `_issues/`)
ROIs not matching cell morphology, "small area selected as cell", "middle cell not selected /
random area", and "13 fragmented ROIs for 2 infected cells". Root cause: a patchy cell mask +
*connected-components* cannot guarantee one-cell-one-region. A **labelled watershed** can.

### Added
- `ensureMorphoLibJ()` — runs only when automatic ROI is chosen: a hard existence check
  (`List.setCommands` / `Marker-controlled Watershed`) that exits with install instructions if
  missing, plus a version-confirmation dialog (`MORPHOLIBJ_TESTED_VERSION`).
- `makeSeedNuclei()` — nucleus mask filtered to nuclei ≥ `NUC_MIN_AREA` (100 px) → `Nuc_Mask` seeds.
- `watershedCellsFromNuclei()` — distance-to-nucleus EDM as input, `Nuc_Mask` markers, `Cell_Mask`
  mask → `Cell_Labels` (one integer label per cell). *(Exact MorphoLibJ option string to confirm
  via the macro recorder on first run.)*
- `buildAutoCytosolRoi` rewritten: one ROI per label → optional nucleus subtraction
  (`AUTO_EXCLUDE_NUC`) → keep if area ≥ `CELL_MIN_SIZE` **and** infected (p99 of both bg-subtracted
  channels ≥ `INFECT_MIN_P99`). Returns N kept ROIs; the existing per-ROI coloc loop measures each.
- `askTokenMapping()` (+ `tokenIndexOf`, `labelAtOr`) — startup dialog to map which filename token
  holds marker1 / marker2 / timepoint; `tryParseFilename` uses `TOK_M1/M2/TP`. Cell line is the
  dialog choice, not parsed.
- `makeCellMask` now morphologically **closes** + fills holes → solid per-cell blobs (less leakage).

### Changed
- `CELL_MIN_SIZE` is now the **per-kept-cell** minimum area (default 200 px), not a cell-mask
  particle floor. New config: `NUC_MIN_AREA`, `AUTO_EXCLUDE_NUC`, `MORPHOLIBJ_TESTED_VERSION`,
  `TOK_M1/M2/TP`.

### Removed
- `splitTouchingCellsVoronoi()` (incl. its stray `waitForUser("can I close everything?")` debug
  pause) and the single-cell centrality selection.

### To verify on first run (in Fiji)
- MorphoLibJ installed (the gate enforces this); recorded `Marker-controlled Watershed` options
  match the macro call; tune the `input` surface (distance-to-nucleus vs intensity) if cells
  leak/over-split.
- The "13 fragments / 2 infected" image now yields ~2 clean per-cell ROIs; uninfected/tiny cells
  dropped; one CSV row per kept cell.

---

## v0.8.0

> [!NOTE]
> Status: current version. The single-cell colocalisation workflow is considered
> feature-complete; further changes are expected to be fixes, not major rewrites.
> **Code written, not yet validated in Fiji** (see "To verify on first run" below).

### Recap — what v0.7.0 did

v0.7.0 was the first **single-cell** version. For each MOI image it built ONE cytosol
ROI and ran the colocalisation step on it. The ROI was built **before** background
subtraction (so it is independent of the measured signal). Two ROI modes:

- `auto_central`: take the DAPI nucleus closest to the image centre, **grow it by a
  fixed distance** (`NUC_DILATE_UM`, via `Enlarge`) to approximate the cell territory,
  then subtract the nucleus → cytosol.
- `manual_draw`: pause per image and let the user draw a freehand ROI on channel 1.

Automatic mode (`auto`) tried to capture the colocalisation numbers by **parsing the
IJ Log** with `extractLogValue()` / `tryLabels()`.

### Problems in v0.7.0 (and how v0.8.0 fixes them)

| # | Problem in v0.7.0 | Fix in v0.8.0 |
|---|---|---|
| **A1** | "Nucleus nearest the centre" often selected an **uninfected** cell (no marker signal) or a **mitotic** cell → the colocalisation value was meaningless. | **Infection gate**: a candidate cell is accepted only if the **99th percentile of BOTH background-subtracted marker channels** inside it is ≥ `INFECT_MIN_P99` (default 100). p99 (not max) is robust to single hot pixels; on bg-subtracted data 100 cleanly separates real puncta from ~0 background. Selection priority: infected + fully-imaged → infected → complete → any. |
| **A2** | The **fixed circular dilation** of the nucleus did not match the real cytoplasm (nuclei are round, cytoplasm is elongated/asymmetric and differs between cell sizes), so the ROI missed most of the true cytosol. | The territory is now derived from the **actual cell shape**: a cell mask from the 568-channel autofluorescence, **split between touching cells by a nucleus-seeded Voronoi tessellation**, minus the nucleus. Pure ImageJ macro — no extra plugin. Adapts to different cell sizes because the mask follows the signal, not a fixed radius. |
| **M1** | In manual mode the nucleus had to be drawn around by hand to exclude it — slow and error-prone. | An up-front checkbox **"Exclude nucleus from drawn ROIs"** automatically subtracts the DAPI nucleus mask from every drawn ROI (image-based subtraction, so it never adds nucleus area outside the ROI). |
| **M2** | The ROI could only be drawn on **channel 1**, where the other markers / nucleus were not visible. | Drawing now happens on a **merged RGB composite** of all three channels (`Draw_RGB`). |
| **M3** | Only **one ROI per image** could be measured, even when several good cells were present. | Manual mode now collects **multiple ROIs per image** first; the colocalisation step then runs on each ROI in turn (one CSV row per `roi_index`). |
| **G1** | Automatic colocalisation wrote **empty CSV rows** and spammed the Log with `WARN: no coloc labels parsed`. Cause: the "Colocalization Threshold" plugin writes to its own **Results table**, not the IJ Log, so the Log parser found nothing. | The Log parser (`extractLogValue` / `tryLabels` / `appendColocRowFromLog`) is **deleted**. Values are now read from the plugin's **Results table** via `getResult()` (for both `auto` and `manual` modes). The first read dumps `Table.headings` to the Log so the exact column names can be confirmed. |

### Added

- **Infection gate** (`roiPercentile()` + p99 check on bg-subtracted duplicates) so only
  infected cells are analysed in `auto_central`.
- **Real-shape cytosol** for `auto_central`: `makeCellMask()` (autofluorescence, ported
  from the Mock pipeline), `splitTouchingCellsVoronoi()` (nucleus-seeded Voronoi cut),
  and per-cell `cell − nucleus` cytosol.
- **Manual multi-ROI workflow** (`buildManualRois()`): draw on a composite (`makeDrawComposite()`),
  add several ROIs, optional automatic nucleus exclusion (`excludeNucleusFromAllRois()`).
- **Per-ROI colocalisation loop**: `processOneImage()` iterates over every ROI; one CSV
  row per ROI with a new **`roi_index`** column.
- **Results-table reader** (`appendColocRowFromTable()` + `safeResult()`), with a one-time
  `Table.headings` dump for header verification.
- QC composite (`saveColocQc()`) now overlays **all** ROIs of the image.
- New CONFIG: `CELL_THR_METHOD`, `CELL_THR_FACTOR`, `BLUR_SIGMA_CELL`, `CELL_MIN_SIZE`,
  `INFECT_MIN_P99`, `BORDER_MARGIN_PX`.

### Changed

- `auto_central` ROI is now signal-aware for **cell selection** but still
  morphology-defined for the **ROI boundary** (so the colocalisation value stays unbiased).
- CSV header gains `roi_index`; `roi_mode` retained; `nuc_dilate_um` removed.
- `applyCellLineDefaults()` now tunes the cell-mask / cell-size knobs
  (Huh7 `CELL_MIN_SIZE` 3000, VeroE6 1000).

### Removed

- The IJ-Log colocalisation parser: `extractLogValue()`, `tryLabels()`,
  `appendColocRowFromLog()`, and the `getInfo("log")` diff in `doColocAuto()`.
- The fixed DAPI-dilation territory logic (`NUC_DILATE_UM`, `Enlarge`-based territory).

### To verify on first run (in Fiji)

- The colocalisation reader assumes the "Colocalization Threshold" plugin writes to the
  **standard ImageJ Results table**. If the CSV columns come back blank with
  `..._noresult` status, the plugin uses a custom results window — switch the reader to
  `Table.get(col, row, "<window title>")` using the title printed alongside the headings.
- `auto_central` chooses an infected cell; the QC JPG shows a real-shape, neighbour-split
  cytosol; manual multi-ROI produces one CSV row per ROI.

---

## v0.7.0

> [!NOTE]
> Status: OLD (`02_..._v0.7.ijm`). First **single representative-cell** ROI; ROI built on the
> raw channels **before** subtraction (back to the v0.4 ROI-FIRST order).

### Changed — from full-image to a single central cell
Replaces v0.6's biased signal-OR full-image ROI with a ROI around **one representative cell**.
Two strategies, chosen at startup via `ROI_MODE`:
- **`auto_central`** — pick the nucleus closest to the image centre whose **dilated territory
  does not touch the image border** (= a fully-imaged cell); cytosol = territory − nucleus.
- **`manual_draw`** — pause per image; the user draws a freehand ROI around the central cell
  (matches the supervisor's manual workflow).

### Changed — ordering and QC
- ROI is built on the **raw** (un-subtracted) channels/DAPI, **before** background subtraction
  → the boundary is **signal-independent (unbiased)** and, in manual mode, the whole cell is
  still visible via autofluorescence. Per-image order: split → build ROI → subtract bg → coloc
  → merge → save.
- Adds a **per-image coloc QC JPG** (`saveColocQc`): a composite of all three channels with the
  **actual ROI burnt in** as an outline, so the chosen/drawn cell is visually verifiable.

> [!WARNING]
> Known problem: "the nucleus nearest the centre" is sometimes an **uninfected** (or mitotic)
> cell, so the coloc value is meaningless. v0.8 fixes this with an **infection gate**.

---

## v0.6.0

> [!NOTE]
> Status: OLD (`02_..._v0.6.ijm`). *(The file's own header was accidentally a copy of v0.5's;
> the real v0.6 changes — recovered from the code diff v0.5 → v0.6 — are below, and the
> header has been corrected.)*

### Changed — order of operations and cell mask
- **Subtract background FIRST, then build the masks.** On the bg-subtracted image the
  histogram is a clean **0-peak + signal tail**, which thresholds more predictably than the
  raw autofluorescence. *(This reverses v0.4's ROI-FIRST order — and is itself reverted again
  in v0.7, which goes back to ROI-FIRST for an unbiased boundary.)*
- **Cell mask = OR of BOTH channels' signal masks** (`makeChannelSignalMask` per channel,
  then `imageCalculator("OR")` + fill holes), replacing the single priority channel. Captures
  signal present in **either** marker. Tuned for bg-subtracted data: `CELL_THR_METHOD =
  Triangle`, `CELL_THR_FACTOR = 1.0`, `BLUR_SIGMA_CELL = 2`, `MIN_PARTICLE_SIZE = 20` (keeps
  small punctate signal).

### Changed — coloc Log parsing made robust
- `tryLabels()` tries **several alternative label spellings** per value (e.g. `Rtotal` /
  `Pearson's R total`; `Ch1 thresh` / `Image1 Min Threshold`; `M1:` / `Manders' M1`), because
  the "Colocalization Threshold" plugin's Log wording varies. When nothing parses, it now
  **dumps the raw Log chunk** so the labels can be read and added.

> [!WARNING]
> The v0.6 ROI is the **full-image signal-OR cytosol** — it measures colocalisation only
> where there is signal, which **inflates Pearson's R** (methodologically circular). This bias
> is the reason v0.7 switched to a single, morphology-defined representative-cell ROI.

---

## v0.5.0

> [!NOTE]
> Status: OLD (`02_..._v0.5.ijm`). Cell mask = V0.4-style single priority channel,
> `setAutoThreshold(Li)` on the raw autofluorescence.

### Added
- **Particle-size filter** on the cytosol mask — disconnected regions smaller than
  `MIN_PARTICLE_SIZE` px are dropped (same logic as the Mock pipeline; removes coverslip
  dirt that thresholded in by accident).
- **Cell-line-specific defaults** via `applyCellLineDefaults()` (VeroE6 settings differ from
  Huh7 — smaller, dimmer cells need a smaller particle floor and a more permissive factor).
- **Dialog-mode startup config** (like Mock V5): markers, combos, timepoints and cell line
  are all overridable at startup.
- **Auto-coloc that parses the plugin's Log output** and writes **real** values to the CSV
  (no longer just placeholders).

### Fixed
- `setOption("BlackBackground", true)` pinned at MAIN **and** per iteration — defends against
  the `Convert to Mask` inversion bug seen in Mock V5.
- `parseFloat()` at background-value **store and retrieve** — guards against the IJM quirk
  where `Array.concat` onto an empty array stores numbers as **strings**.

---

## v0.4.0

> [!NOTE]
> Status: OLD — **no source file is retained** for v0.4.0 (the file that looked like v0.4 turned
> out to be v0.3, above). Reconstructed from later headers (`v0.5`/`v0.7` refer back to "the
> V0.4 ordering / V0.4-style cell mask").

The version that established two conventions every later version keeps:

- **ROI-FIRST ordering** — build the cytosol ROI on the **raw** channels **before**
  subtracting the background. This makes the ROI **independent of the measured signal**
  (unbiased) and, in manual workflows, keeps the whole cell visible (a bg-subtracted image
  is mostly zeros + dots and is hard to trace).
- **Single priority-picked cell-mask channel** auto-thresholded with **`Li`** on the raw
  autofluorescence.

It also wired up the first working **Coloc → CSV bridge** (real values, the column set that
v0.5 then parsed from the plugin Log). Per v0.5's header, the **bg-value-based** and
**combined C1 + C2** cell-mask experiments were tried and **reverted** around here, in favour
of the simpler single-channel `Li` auto-threshold.

---

## v0.3.0

> [!NOTE]
> Status: OLD. Source file = `02_..._v0.3.ijm` (header "V0.3", `MACRO_VERSION "0.3.0"`). A diff
> against v0.2 shows the **only** code change is the version constant. *(The file was briefly
> on disk misnamed `..._v0.4.ijm`; renamed to v0.3 to match its header.)*

A **version-bump checkpoint**: the code is identical to v0.2 (cytosol ROI built on the raw
channels before subtraction → subtract bg → re-merge → manual Coloc-2 hand-off → placeholder
CSV row), with `MACRO_VERSION` bumped `0.2.0 → 0.3.0`. No functional change — a saved snapshot
between the v0.2 prototype and the v0.4 work.

> [!WARNING]
> **A real v0.4 source file is not retained** (the file that looked like v0.4 was this v0.3).
> The v0.4.0 entry below is reconstructed from references in later headers.

---

## v0.2.0

> [!NOTE]
> Status: OLD. First **colocalisation-prep** version of Script 2 — it adds the cytosol
> ROI + the manual Coloc-2 hand-off on top of v0.1.0's pure background subtraction.
> *(An earlier revision of this file mistakenly held the Mock pipeline's v0.2 notes; the
> entry below describes the actual `02_..._v0.2.ijm` background/coloc-prep macro.)*

### Description

For each MOI image this version:
1. builds a **generous cytosol ROI** from the cell-mask source channel on the **raw**
   signal — **before** any background subtraction, so the threshold sees real signal;
2. subtracts the per-channel background values from C1 and C2 (DAPI untouched);
3. re-merges the channels into a scrollable multi-channel **Color hyperstack** (like the
   original CZI) and saves it as a 16-bit TIF;
4. **pauses** for the user to run **Coloc 2** (or another coloc plugin) on the two
   non-DAPI channels inside the cytosol ROI;
5. appends a **placeholder row** to a per-run colocalisation CSV.

### Added

- Cytosol-ROI construction from the cell-mask channel (cell mask − nucleus mask), built on
  the raw channels before subtraction.
- Channel re-merge into a Color hyperstack + 16-bit TIF save.
- Per-run colocalisation CSV with a placeholder row per image (manual Coloc-2 step).

### Not yet implemented

- An automated Coloc → CSV bridge (real values instead of placeholders) — planned for v0.4
  once the exact column set was decided.
- Particle-size filtering, cell-line-specific defaults, QC JPG export, single-cell ROI.

---

## v0.1.0

**Major changes**

- First **functional** version — the whole background-subtraction path runs
  end-to-end without errors.
- **Scalable background lookup table** replacing the 12 named globals.
- **Background subtraction + channel re-merge + saving** actually implemented.
- Coloc Threshold is not possible yet

### Description

For each MOI `.tif` in the chosen folder, this version subtracts the per-channel
background values (measured by the Mock pipeline) from C1 and C2, leaves DAPI
(C3) untouched, re-merges the channels into a 16-bit multi-channel TIF, and saves
it into a run-specific output folder together with a markdown record of the
values used and the full IJ Log of the run.

### Added

- **Background lookup table** (`BG_KEYS` / `BG_VALUES`) with `setBgValue()` /
  `getBgValue()`. Keyed by `<marker>_in_<combo>_<tp>`. Scales to new combos /
  timepoints with zero extra code.
- **`getBgValue()` hard-exits if a key is missing** — fail loudly instead of
  silently subtracting 0 after a typo.
- **Dynamic background-value dialog** (`askBackgroundValues()`): the entry fields
  are generated by looping over `TIMEPOINTS × ANALYSE_COMBI`. The read-back loop
  mirrors the add loop exactly, so the value order can never drift.
- **Actual background subtraction** (`subtractFromChannel()`): `run("Subtract...")`
  on each marker channel; clamps at 0, keeps 16-bit.
- **Channel re-merge** (`mergeChannelsBack()`) using `Merge Channels ... create`
  to preserve the 16-bit data (without `create` you would get 8-bit RGB and lose
  dynamic range).
- **Saving** of each bg-subtracted image as `bgsub_<original>.tif`.
- **Timepoint validation** in addition to the combo validation.
- **`logHeader()`** for a clean run header in the Log.
- **`saveLogToFile()`** to persist the full IJ Log as a text file per run.
- **`ensureDir()`** helper (non-recursive `File.makeDirectory` guard).
- **`writeBgMarkdown()`** generates the provenance markdown dynamically from the
  same domain lists that drive everything else (cannot drift out of sync).
- Per-image Log lines documenting which bg value was subtracted from which
  channel (audit trail).

### Changed

- Reworked `processOneImage()` into a clean, complete sequence:
  open → resolve metadata → validate → look up bg → split → subtract → merge → save.
- Metadata resolution now properly supports automatic (filename) and manual
  (dialog) modes, with a parse-failure fallback dialog in automatic mode.
- `tryParseFilename()` returns `[tp, m1, m2]` and validates tokens against
  `MARKERS` / `TIMEPOINTS` up front.
- Output folder renamed to `<RUN_ID>_bgsub/` (was `<RUN_ID>_substracted_images/`).

### Fixed

- Fixed the broken `if ANALYSE_COMBI == comboKey` statement and missing
  semicolons from the draft.
- Fixed the function-name typo (`listMoOIFiles` → `listMoiFiles`, called
  correctly).
- Fixed `makeRunId()` missing its closing brace (which had nested
  `buildOutputDir()` inside it).
- Completed `processOneImage()` (the draft referenced an undefined
  `substactBackground()` and had no merge/save).

### Removed

- Removed the 12 separate named background globals (replaced by the lookup table).
- Removed the hard-coded `get_background_values()` / `createMDFile()` dialogs
  (replaced by the dynamic, loop-generated versions).

### Not yet implemented (comes in later versions)

- Cytosol-ROI generation, the colocalisation step, particle-size filtering,
  cell-line-specific defaults, QC JPG export.

---

## v0.0.1

> [!NOTE]
> **Status: first conceptual draft ("vibe coding"), NOT a functional release.**
> This version was written to sketch out the *intended* structure and workflow
> of the background-subtraction pipeline. It contains syntax errors, incomplete
> branches, and a broken `processOneImage()` / background-subtraction step. It is
> kept in the repository as the historical starting point — to show where the
> pipeline began and how it evolved. Do not run this version; use the latest one.

### Idea / intended workflow

The goal of this first draft was to outline how infected (MOI) images should be
background-corrected before colocalisation. The planned order of operations was:

> [!Important]
>
> Expected filename scheme:
> timepoint_cellLine_condition_marker1_marker2_CS#_imgIdx.tif

1. **Choose the input folder** containing the MOI `.tif` images (`chooseInputDir()`).
2. **Create a run ID** (timestamp `YYYYMMDD_HHMM`) so every run is uniquely
   labelled and reproducible (`makeRunId()`).
3. **Ask for the run mode** (automatic vs. manual metadata handling)
   (`askModeAndConfig()`).
4. **Collect all background values** in one dialog — one value per
   (marker × combo × timepoint) — entered by hand from the Mock-pipeline results
   (`get_background_values()`). In this draft each value was stored in its own
   named global variable (e.g. `HA568_in_HA568_dsRNA488_12h`).
5. **Build the output folder** `<RUN_ID>_substracted_images/` and write a
   `background_values_used.md` file documenting exactly which value was entered
   for which channel/combo/timepoint (`buildOutputDir()`, `createMDFile()`).
6. **List all MOI images** in the folder (filename must contain `moi`)
   (`listMOIFiles()`).
7. **Per image:** open it, parse the filename into timepoint / marker1 / marker2
   (`tryParseFilename()`, with `askImageMetadata()` / `askParseFailureAction()`
   as manual fallbacks), validate the marker combination against the allowed list
   (`isValidCombo()`), split the channels and rename them by marker
   (`splitAndRenameChannels()`), and **subtract the matching background value**
   from each marker channel.
8. **Clean up** between images (`cleanupBetweenImages()`).

### What this draft already contained (as concepts)

- Folder selection + run-ID generation.
- A single large dialog for entering all background values by hand.
- Markdown documentation of the entered background values (provenance).
- Filename parsing with a manual-entry fallback dialog.
- Marker-combination validation against `ANALYSE_COMBI`.
- Channel splitting + marker-based channel renaming.
- A batch loop over all MOI images.
- Utility helpers (`inArray()`, `cleanupBetweenImages()`).

### Known limitations of this draft (fixed in later versions)

- Missing semicolons and a broken `if ANALYSE_COMBI == comboKey` statement.
- `processOneImage()` is incomplete — the actual subtraction (`substactBackground()`)
  is referenced but never defined; the merge/save steps are absent.
- Function name typo: `listMoOIFiles()` defined but `listMOIFiles()` called.
- `makeRunId()` is missing its closing brace, so `buildOutputDir()` is nested
  inside it.
- Background values stored as 12 separate named globals instead of a scalable
  lookup table.
- No cytosol-ROI generation, no colocalisation step, no image saving yet.