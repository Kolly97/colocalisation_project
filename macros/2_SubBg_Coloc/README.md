# Subtract Background + Colocalisation Prep (Script 2)

Fiji/ImageJ macro for the **second stage** of the colocalisation pipeline. For every
infected (**MOI**) image it subtracts the per-channel background, builds **per-cell cytosol
ROIs** (one per infected cell, or hand-drawn), and runs (or prepares) the colocalisation
analysis, writing one CSV row per ROI.

> **Current version:** `02_Subtract_background_coloc_v0.14.3.ijm` (`MACRO_VERSION 0.14.3`).
> Older `02_..._v0.x.ijm` files are kept for history; always run the highest version.
> See [`CHANGELOG.md`](CHANGELOG.md) for the version history.
> **`automatic watershed ROI` mode requires the MorphoLibJ plugin** (IJPB-plugins update
> site) for the nucleus-seeded marker-controlled watershed; `manual draw` mode does not.

---

## Where this fits in the pipeline

```
.czi → Bio-Formats → .tif  (Mock + MOI images)
        │
        │  SCRIPT 1  (../1_Mock/01_Mock_pipeline.ijm)
        │    measures Top-1% cytosol intensity in MOCK images
        │    → per-image CSVs of background-estimate values
        │    → you take the MEDIAN of the `p99_9995` column per
        │      (marker × combo × timepoint × cell line)
        ▼
   ┌─────────────────────────────────────────────────────────┐
   │  SCRIPT 2  (this macro)                                  │
   │   • subtract those background values from each MOI image │
   │   • build one cytosol ROI per infected cell (or by hand) │
   │   • run Colocalisation Threshold per ROI → coloc CSV     │
   └─────────────────────────────────────────────────────────┘
        ▼
   (downstream) aggregate coloc CSVs → stats / figures
```

The **link between the two scripts is manual**: after running Script 1, read its output
CSVs, take the median `p99_9995` value per marker/combo/timepoint, and type those numbers
into this macro's background-values dialog at startup.

---

## Requirements

- **Fiji** (ImageJ ≥ 1.53).
- The **"Colocalization Threshold"** plugin
  (`Analyze ▸ Colocalisation ▸ Colocalisation Threshold`; Tony Collins / WCIF). Ships with
  standard Fiji.
- **MorphoLibJ** (IJPB-plugins update site) — only for `automatic watershed ROI` ROI mode
  (marker-controlled watershed). Install via `Help ▸ Update… ▸ Manage update sites ▸ tick
  IJPB-plugins ▸ Apply ▸ restart`. The macro hard-checks for it and aborts with these
  instructions if missing. Not needed for `manual draw`.
- Input images must be **calibrated** confocal `.tif` (Bio-Formats export keeps calibration).
  Calibration also sets the nucleus-ring width (`GATE_RING_UM` µm → px).

---

## Input expectations

**Filename schema (binding, positional, underscore-separated):**

```
<timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<coverslip>_<imgIndex>.tif
e.g.  24h_Huh7_MOI1_HA568_dsRNA488_CS1_2.tif
```

- Token 0 = timepoint, 1 = cell line, 2 = condition, 3 = marker1 (**C1**),
  4 = marker2 (**C2**), 5 = coverslip, 6 = image index.
- **Channel order is fixed:** C1 = marker1, C2 = marker2, **C3 = DAPI (always)**.
- Only files whose name contains `moi` (case-insensitive) are processed.
- If a filename cannot be parsed, the macro falls back to a per-image metadata dialog
  (or you can choose manual mode at startup).

---

## How to run

1. In Fiji: `File ▸ Open…` the `.ijm`, or drag it onto the toolbar, then **Run** (or
   `Macros ▸ Run`).
2. **Choose the input folder** with the MOI `.tif` images.
3. **Pipeline-mode dialog:**
   - *Marker source* — `automatic` (parse from filename) or `manual` (ask per image).
   - *Pipeline* — `subtract only` / `subtract + manual coloc` / `subtract + auto coloc`.
   - *ROI strategy* — `automatic watershed ROI` or `manual draw` (see below).
   - *Cell line* — `Huh7` or `VeroE6` (sets the strictness defaults).
   - *Save JPG* — optional 8-bit RGB composite with a scale bar.
4. **Pipeline-setup dialog** (always shown) — the marker / combo / timepoint lists,
   pre-filled with the workflow standards. Just click OK to keep them.
5. **Background-values dialog** — one value per *(marker × combo × timepoint)*. Enter the
   medians from Script 1.
6. The macro then loops over every MOI image.

---

## The two ROI strategies

The colocalisation ROI is always built on the **raw** channels, *before* background
subtraction, so the ROI boundary is independent of the measured colocalisation signal
(unbiased) and the whole cell is visible while drawing.

### `automatic watershed ROI` (multi-cell, pure IJM)

Builds **one cytosol ROI per infected cell** automatically:

1. **Cell mask** from the 568-channel autofluorescence — blur → threshold → morphological
   **close** → fill holes, so each cell is one **solid** blob.
2. **Nucleus mask** from DAPI; keep nuclei ≥ `NUC_MIN_AREA` (**5000 px**) as cell seeds (drops
   debris and small mitotic fragments).
3. **Marker-controlled watershed** (MorphoLibJ): the nuclei are turned into a **label image**
   (one integer per nucleus) and seed one region per cell, bounded by the cell mask, so
   touching cells are split into separate territories.
4. One ROI per cell (area ≥ `CELL_MIN_SIZE`, **120000 px²**), nucleus optionally subtracted
   (`EXCLUDE_NUC`).
5. **Keep rule (infection gate):** for each nucleus, grow it ≈ `GATE_RING_UM` (**8 µm**) into
   the cytoplasm, **clip the disk to that cell's own territory**, and measure the **raw**
   `GATE_PCTL` (**99.99th**) percentile of both channels there. The cell is kept only if **both
   channels** reach **≥ `SIGNAL_BG_MULT` (default 1.5) × the entered background** — i.e. it
   shows real signal at least 50 % above the bg you typed. Uninfected cells are dropped.
6. The colocalisation step then runs on **each kept cell** (one CSV row per cell); the QC JPG
   outlines all of them.

> Choosing *which* cells to measure uses the signal (an uninfected cell carries no biology); the
> ROI *boundary* is still pure morphology, so the colocalisation value stays unbiased. The gate
> is measured in a ring around the nucleus (not the whole cell), so it is robust to how exactly
> the watershed drew the cell edge.

### `manual draw`

You draw the ROI(s) yourself, matching the supervisor's manual workflow:

- Drawing happens on a **merged composite** of all three channels.
- You can add **several ROIs per image** (draw all good cells first; coloc then runs on
  each in turn).
- The **"exclude nucleus from drawn ROIs?"** choice is asked **once at startup** (not per
  image); if on, the DAPI nucleus is subtracted from every drawn ROI automatically.

---

## Colocalisation step

Set by the *Pipeline* choice:

- **subtract only** — no ROI, no coloc; just background-subtracted TIFs.
- **subtract + manual coloc** — the macro applies the ROI and pauses; you run
  `Colocalisation Threshold` by hand with the settings shown in the dialog.
- **subtract + auto coloc** — the macro runs `Colocalization Threshold` automatically per
  ROI.

The coloc **numbers are NOT auto-parsed** (the plugin's output could not be read reliably). The
macro writes one **provenance row per ROI** (metadata + `roi_index` + empty value columns +
status); read/export the Rtotal / Manders / threshold values from the plugin's own results window.

Recommended plugin settings (also shown in the manual dialog):
`Channel 1 = marker1`, `Channel 2 = marker2`, `Use ROI = Channel 1`,
`Channel = Red : Green`, *Include zero-zero pixels in threshold calc* = on.

---

## Outputs

Everything lands in `<inputDir>/<RUN_ID>_bgsub/` (`RUN_ID = YYYYMMDD_HHMM`). Images go in
**sub-folders**; CSV / markdown / log stay at the root:

| File | Contents |
|---|---|
| `tif/bgsub_<original>.tif` | Background-subtracted multi-channel **16-bit** image. |
| `qc/qc_<original>_roi.jpg` | Composite of all channels with the **ROI(s) outlined and numbered** (the number == `roi_index` in the CSV) — for checking cell choice and matching a coloc row to its cell. |
| `jpg/bgsub_<original>.jpg` | *(optional)* 8-bit RGB composite with a scale bar. |
| `coloc_results_<RUN_ID>.csv` | One row per ROI: metadata + `roi_index` + (empty) coloc value columns + provenance. |
| `background_values_used.md` | The background values entered for this run (provenance). |
| `macro_log_<RUN_ID>.txt` | Full IJ Log of the run. |

Every CSV row carries `roi_mode`, `macro_version`, and `run_id`, so runs never collide and
can be filtered downstream (e.g. with pandas).

---

## Key configuration (top of the `.ijm`)

| Variable | Meaning | Default |
|---|---|---|
| `ROI_MODE` | `auto_watershed` or `manual_draw` (also set in the dialog) | `auto_watershed` |
| `SIGNAL_BG_MULT` | Infection gate: raw p99.99 in the nucleus ring must be ≥ this × entered bg, **both** channels | `1.5` |
| `GATE_RING_UM` | How far the nucleus is grown into the cytoplasm to sample signal (µm → px) | `8` µm |
| `GATE_PCTL` | Percentile measured in the ring (high = robust to hot pixels) | `99.99` |
| `NUC_MIN_AREA` | Min nucleus area (px) to seed a cell — drops debris/mitotic fragments | `5000` |
| `CELL_MIN_SIZE` | Minimum kept-ROI footprint (**raw px²** — rescale with magnification) | `120000` |
| `EXCLUDE_NUC` | Subtract the nucleus from each ROI (auto + manual) | `true` |
| `CELL_THR_METHOD` / `CELL_THR_FACTOR` | Autofluorescence cell-mask threshold method / permissiveness | `Li` / `0.5` |
| `BLUR_SIGMA_CELL` / `BLUR_SIGMA_NUC` | Gaussian smoothing (µm) for the masks | `1` / `1` |
| `NUC_THR_METHOD` / `NUC_CLOSE_ITER` | Nucleus-mask threshold / morphological closing | `Otsu` / `2` |
| `BORDER_MARGIN_PX` | Cells whose territory comes within this many px of an edge are de-prioritised | `0` |
| `COLOC_MODE` | `none` / `manual` / `auto` (also set in the dialog) | `manual` |
| `SAVE_COLOC_QC` / `SAVE_JPG` | QC composite / scale-bar JPG export | `true` / `false` |

> See [`USAGE.md`](USAGE.md) for the full step-by-step (every dialog and setting) and
> [`CHANGELOG.md`](CHANGELOG.md) for why the v0.14.x defaults were tuned to these values.

Cell-line-specific values are applied by `applyCellLineDefaults()` after you pick the cell
line. Tune the strictness knobs there, not inside the per-image functions.

---

## Known limitations / things to watch

- **Coloc numbers are not auto-captured:** the macro runs / prepares the plugin and writes a
  provenance row per ROI with **empty** value columns (the plugin's output could not be read
  reliably). Read/export Rtotal / Manders / thresholds from the plugin's own window.
- **Cytosol shape** is approximated from autofluorescence + a nucleus-seeded watershed (there
  is no membrane stain); it is not a true cell boundary.
- **MorphoLibJ option strings** (`Connected Components Labeling`, `Marker-controlled
  Watershed`) can differ between builds; confirm them with `Plugins ▸ Macros ▸ Record` on
  first run if the watershed errors.
- The ROI definition for publication should still be **confirmed with the supervisor**.

---

## See also

- [`CHANGELOG.md`](CHANGELOG.md) — version history and the v0.7 → v0.8 fixes.
- [`../1_Mock/`](../1_Mock/) — Script 1, which produces the background values used here. 
