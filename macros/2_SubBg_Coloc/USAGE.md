# USAGE — Stage 2: Background subtraction + Colocalisation

Step-by-step guide for `02_Subtract_background_coloc_v0.14.3.ijm` (always run the **highest
version number** in this folder). Detailed companion to [`README.md`](README.md). It covers
every startup dialog, both ROI modes, all three coloc modes, and every tunable setting.

> **What this macro does, per MOI image:** subtract the per-channel Mock background → build
> one **cytosol ROI per infected cell** (automatically, or by hand) → run the **Colocalisation
> Threshold** analysis per cell → write one CSV row per ROI + a numbered QC image.
>
> **Order matters (and is deliberate):** the ROI is built on the **raw** channels *before*
> background subtraction, so the ROI boundary is **independent of the signal you measure**
> (unbiased), and the whole cell is visible while drawing.

---

## 0. Before you run — preparation

1. **Finish Stage 1 first.** You need one background value per *(cell line × timepoint ×
   combo × marker)* — the **median of the `p99_9995`** column from the Mock CSVs (see the
   [Mock USAGE](../1_Mock/USAGE.md) §6). Have that little table in front of you.
2. **Install plugins:**
   - **Colocalization Threshold** — ships with Fiji (`Analyze ▸ Colocalisation ▸
     Colocalisation Threshold`). Needed for the coloc step.
   - **MorphoLibJ** (IJPB-plugins update site) — needed **only for automatic ROI mode**.
     `Help ▸ Update… ▸ Manage update sites ▸ tick IJPB-plugins ▸ Apply ▸ restart`. The macro
     hard-checks for it and aborts with instructions if missing. (Manual ROI mode needs no
     extra plugin.)
3. **Calibrated 16-bit `.tif`**, C1=marker1, C2=marker2, **C3=DAPI**. Calibration matters: the
   infection-gate ring is specified in **µm** and converted to pixels via the image calibration.
4. Only files whose name contains **`moi`** (case-insensitive) are processed.
5. **Test on a few images first** (one per combo, ideally a dense field with touching cells)
   and check the numbered QC before running the whole batch.

---

## 1. The run, end to end

```
choose MOI folder
  → Dialog 1  "Pipeline mode"           (§2.1)  marker source · pipeline · ROI strategy · cell line · 2 checkboxes
  → [auto ROI] MorphoLibJ check         (§2.2)
  → [manual ROI] "Manual ROI options"   (§2.3)  once-only: exclude nucleus?
  → Dialog "Pipeline setup"             (§2.4)  markers / combos / timepoints lists
  → Dialog "Background values"          (§2.5)  type the Stage-1 medians
  → [auto marker source] token mapping  (§2.6)
  → loop over every MOI image:
        split → build ROI(s) on RAW channels → save numbered QC
              → subtract bg → coloc per ROI → merge → save TIF
```

When the Log prints `=== DONE ===` it has finished.

---

## 2. Startup dialogs (all the modes)

### 2.1 Dialog 1 — "Pipeline mode"

| Field | Options | Meaning |
|---|---|---|
| **Marker source** | `automatic` / `manual` | Where per-image metadata (timepoint + the two markers) comes from. `automatic` parses the filename (token mapping §2.6); `manual` asks per image. Cell line is always the dialog choice below. |
| **Pipeline** | `subtract only` · `subtract + manual coloc` · `subtract + auto coloc` | How far the macro goes. **`subtract only`** = just write bg-subtracted TIFs (no ROI, no coloc). **`+ manual coloc`** = build ROIs, then pause so you run the plugin by hand (§4). **`+ auto coloc`** = build ROIs and run the plugin automatically (§4). |
| **ROI strategy** | `automatic watershed ROI` / `manual draw` | How cytosol ROIs are made: the nucleus-seeded watershed (§3.1) or hand-drawn (§3.2). Ignored if Pipeline = `subtract only`. |
| **Cell line** | `Huh7` / `VeroE6` | Recorded in every CSV row + the bg-provenance file. |
| ☑ **Exclude nucleus from each cell ROI** | on/off (`EXCLUDE_NUC`) | If on, the DAPI nucleus is subtracted from every ROI (auto **and** manual) → a cytosol ring. If off, the whole cell (incl. nucleus) is the ROI. |
| ☑ **Also save 8-bit RGB JPG with scale bar** | on/off (`SAVE_JPG`) | Optional publication-style composite of the bg-subtracted image with a 20 µm scale bar. |

### 2.2 MorphoLibJ check (only if ROI strategy = automatic)

Two quick dialogs: a hard existence check on the `Marker-controlled Watershed` command (aborts
with install instructions if missing), then a version-confirmation prompt (the macro was
written against MorphoLibJ `1.6.x`). Click **Yes** to continue.

### 2.3 "Manual ROI options" (only if ROI strategy = manual draw)

Asked **once for the whole run** (not per image): a single checkbox **"Exclude nucleus from
drawn ROIs"**. If on, the DAPI nucleus is auto-subtracted from each ROI you draw, so you don't
have to trace around it.

### 2.4 "Pipeline setup" — the domain lists

Markers / Combos / Timepoints, comma-separated, pre-filled with the workflow standards —
**click OK to keep**. They drive validation and the background dialog; an image whose combo or
timepoint isn't listed is skipped (Log note, not an error).

### 2.5 "Background values" — the Stage-1 link

One numeric field per **(timepoint × combo × marker)**, grouped by timepoint. Type the
**median `p99_9995`** value you computed from the Mock CSVs into each field. These are
subtracted from the matching channel, and they also define the infection-gate threshold
(`SIGNAL_BG_MULT × bg`, §3.1). They are written to `background_values_used.md` for provenance.

### 2.6 "Filename token mapping" (only if Marker source = automatic)

Radio groups to pick **which token holds timepoint / marker1 / marker2** (e.g. `0: 24h`,
`3: NS4B568`, `4: dsRNA488`), shown once using your first MOI file. **Cell line is not
parsed** — it is the §2.1 choice.

---

## 3. The two ROI strategies

The ROI(s) end up in the ROI Manager (index 0…n-1); the coloc step then runs on each. Every
ROI is **outlined and numbered** in the QC JPG, and the **number equals the `roi_index`
column** in the CSV, so you can match any coloc row back to its cell.

### 3.1 `automatic watershed ROI` — one ROI per infected cell

Pipeline (all on the **raw** channels, before subtraction):

1. **Seed nuclei** — DAPI → Otsu mask, keep only nuclei **≥ `NUC_MIN_AREA` (5000 px)** so
   debris and small mitotic fragments don't seed spurious cells.
2. **Cell mask** — autofluorescence of the 568 channel → blur → `Li` × `CELL_THR_FACTOR`
   (0.5) → morphological **close** → fill holes → one solid blob per cell.
3. **Label the nuclei** — Connected Components Labeling → one integer per nucleus (a
   **labelled** marker; this is what makes touching cells separate — a binary marker merges
   them).
4. **Marker-controlled watershed** — each labelled nucleus floods its own basin over the
   distance-to-nucleus surface, bounded by the cell mask → **one region per cell, neighbours
   split at the watershed line**.
5. **Per cell** — one ROI per region, nucleus optionally subtracted (`EXCLUDE_NUC`).
6. **Keep filter (which cells are measured):** keep a cell iff
   - **area ≥ `CELL_MIN_SIZE` (120000 px²)** — drops crumbs / partial cells, **and**
   - **infected** — grow its nucleus by **`GATE_RING_UM` (8 µm)** into the cytoplasm, clip
     that disk to the cell's own territory, and require the **raw `GATE_PCTL` (99.99th)
     percentile ≥ `SIGNAL_BG_MULT` (1.5) × the entered background** in **both** channels.

   The Log prints, per cell, `kept`/`drop` with the area and the measured p99.99 vs each
   threshold, then `auto ROIs kept: K of N cells`. **Read these lines to tune** — they tell
   you whether the area filter or the infection gate rejected a cell.

> **Why this is unbiased:** the signal only chooses *which* cells are infected enough to
> measure (an uninfected cell carries no biology); the ROI *boundary* is pure morphology
> (watershed cell mask − nucleus), and the gate is sampled in a nucleus ring, not in the data
> you later correlate. Three independent criteria → a defensible coloc value.

### 3.2 `manual draw` — you draw the ROIs

For each image the macro shows a **merged RGB composite** of all three channels and pauses:

- Draw a **freehand** ROI around a cell → OK. Already-drawn ROIs stay outlined and numbered on
  the image, so you can see what's done.
- Choose **"draw another"** to add more cells, or **"done"**. (Cancel aborts the whole run.)
- If you ticked "exclude nucleus" (§2.3), the DAPI nucleus is subtracted from each ROI
  automatically.

Coloc then runs on each drawn ROI in turn (one CSV row per ROI). This matches the supervisor's
hand-picked-cell workflow.

---

## 4. The coloc step (set by "Pipeline")

- **`subtract only`** — no ROI, no coloc; the macro just writes bg-subtracted TIFs.
- **`subtract + manual coloc`** — the macro applies ROI *i* to `marker1_channel` and pauses
  with a `waitForUser` box listing the exact settings. Run it by hand:
  `Analyze ▸ Colocalisation ▸ Colocalisation Threshold`, then
  - **Channel 1** = `<marker1>_channel`, **Channel 2** = `<marker2>_channel`
  - **Use ROI** = `Channel 1` (reads the active selection the macro applied)
  - **Channel** = `Red : Green`, **Include zero-zero pixels** = ☑
  - Read the values off the plugin window, click OK.
- **`subtract + auto coloc`** — the macro calls the plugin itself for each ROI.

> **Note (spelling):** the menu reads "Colocali**s**ation" (British) but the macro's `run()`
> call uses "Colocali**z**ation" (American) — that's intentional, not a typo.

### Getting the numbers out

The coloc **values are NOT auto-parsed** — the plugin's output could not be read reliably, so
the macro writes one **provenance row per ROI** with the metadata, `roi_index`, `roi_mode`,
`macro_version`, `run_id`, `status`, and **empty value columns**. You read/export Rtotal /
Manders M1,M2 / thresholds from the plugin's own results window, and use the **QC ROI number =
`roi_index`** to know which row is which cell (and to drop any uninfected ones afterwards).

CSV columns: `image, cell_line, timepoint, combo, channel_1, channel_2, roi_index, Rtotal, m,
b, Ch1_thresh, Ch2_thresh, Rcoloc, R_below_thresh, M1, M2, tM1, tM2, Ncoloc, perc_volume,
perc_ch1_vol, perc_ch2_vol, roi_mode, macro_version, run_id, status`.

---

## 5. The settings (CONFIG block at the top of the `.ijm`)

Most are set by the dialogs; edit the file to change defaults. Current v0.14.3 values:

| Variable | Default | Controls |
|---|---|---|
| `ROI_MODE` | `auto_watershed` | also set by §2.1 (`manual_draw` otherwise) |
| `COLOC_MODE` | `manual` | also set by §2.1 (`none` / `auto`) |
| **`NUC_MIN_AREA`** | `5000` px | min nucleus area to seed a cell. **Raw px — rescale with magnification.** |
| **`CELL_MIN_SIZE`** | `120000` px² | min kept-ROI area. **The single most important auto knob; raw px².** |
| **`SIGNAL_BG_MULT`** | `1.5` | infection gate: ring p99.99 must be ≥ this × entered bg, both channels |
| **`GATE_RING_UM`** | `8` µm | how far the nucleus is grown to sample cytoplasm for the gate |
| **`GATE_PCTL`** | `99.99` | percentile measured in the ring (robust to hot pixels) |
| `EXCLUDE_NUC` | `true` | subtract nucleus from each ROI (auto + manual); also the §2.1/§2.3 checkbox |
| `CELL_THR_METHOD` / `CELL_THR_FACTOR` | `Li` / `0.5` | autofluorescence cell-mask threshold / permissiveness (<1 = more cytoplasm) |
| `NUC_THR_METHOD` / `NUC_CLOSE_ITER` | `Otsu` / `2` | nucleus-mask threshold / closing iterations |
| `BLUR_SIGMA_CELL` / `BLUR_SIGMA_NUC` | `1` / `1` µm | Gaussian pre-blur for the masks |
| `BORDER_MARGIN_PX` | `0` | de-prioritise cells whose bbox is within this many px of an edge |
| `SAVE_COLOC_QC` / `SAVE_JPG` | `true` / `false` | numbered QC composite / scale-bar JPG |
| `JPG_SCALEBAR_UM` | `20` | scale-bar length for the optional JPG |
| `MORPHOLIBJ_TESTED_VERSION` | `1.6.x` | version named in the §2.2 confirmation |

> **Pixel vs micron knobs:** `NUC_MIN_AREA` and `CELL_MIN_SIZE` are in **raw pixels** and were
> tuned for the current objective. If pixel size changes, recompute them (area scales with the
> square of the pixel size). `GATE_RING_UM` is in microns and auto-converts, so it is
> magnification-independent.

---

## 6. Output

`<input>/<RUN_ID>_bgsub/` (`RUN_ID = YYYYMMDD_HHMM`):

```
<RUN_ID>_bgsub/
├── tif/  bgsub_<img>.tif            multi-channel 16-bit, background-subtracted
├── qc/   qc_<img>_roi.jpg          all channels + ROI outlines, NUMBERED (number == roi_index)
├── jpg/  bgsub_<img>.jpg           (optional) 8-bit RGB + scale bar
├── coloc_results_<RUN_ID>.csv      one provenance row per ROI (values read from the plugin)
├── background_values_used.md       the bg values entered, + gate description
└── macro_log_<RUN_ID>.txt          full IJ Log (per-cell kept/drop lines, thresholds)
```

---

## 7. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| **`MorphoLibJ … not found`** abort | Install IJPB-plugins (§0.2) or choose `manual draw`. |
| Watershed/CCL command **errors** | Option strings differ between MorphoLibJ builds. `Plugins ▸ Macros ▸ Record`, run the command once by hand, paste the exact string into the macro (`watershedCellsFromNuclei`). |
| **No ROI kept** (`auto ROIs kept: 0 of N`) | Read the `drop` Log lines. Area too big a bar → lower `CELL_MIN_SIZE`. Gate too strict (weak HA568) → lower `SIGNAL_BG_MULT`. Nuclei not seeding → lower `NUC_MIN_AREA`. |
| **Touching cells share one ROI** | Should be fixed (labelled markers). If it recurs, the watershed got a binary marker — verify `Connected Components Labeling` ran (Log: `nucleus seeds (labelled): N`) and the option strings. |
| Too many junk ROIs / crumbs | Raise `CELL_MIN_SIZE`; raise `NUC_MIN_AREA`; raise `CELL_THR_FACTOR` toward 1 (tighter cell mask). |
| Coloc CSV value columns are empty | **Expected** — read the numbers from the plugin window; match by QC number. |
| QC shows no outline | No ROI passed the keep filter (auto) or none drawn (manual) — see the "No ROI kept" row. |

For the design rationale, the watershed binary-vs-label fix, the relative-gate reasoning, and
the IJM pitfalls, see [`../HANDOFF.md`](../HANDOFF.md) and [`CHANGELOG.md`](CHANGELOG.md).
