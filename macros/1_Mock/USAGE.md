# USAGE — Stage 1: Mock background pipeline

Step-by-step guide for `01_Mock_pipeline.ijm` (v0.8.1). This is the detailed companion to
the shorter [`README.md`](README.md). It covers every startup mode, every tunable setting,
and what to check in the output.

> **What this macro is for:** measure, on **Mock (uninfected)** images, how bright the
> cytosol autofluorescence is in each marker channel — expressed as high **percentiles** of
> the cytosol pixel histogram. You then aggregate those into one **background value per
> (cell line × timepoint × combo × marker)** and type it into Stage 2. Mock images have *no
> specific signal*, so whatever brightness sits in their cytosol is the background floor.

---

## 0. Before you run — preparation

1. **Install Fiji** (≥ 1.54). No extra plugins are needed for Stage 1.
2. **Calibrated 16-bit `.tif`** images, single z-plane, channel order **C1 = marker1,
   C2 = marker2, C3 = DAPI**. Bio-Formats `.czi → .tif` preserves the µm/px calibration.
3. **Filenames** follow the schema (underscore-separated):
   ```
   <timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<coverslip>_<imgIndex>.tif
   e.g. 12h_Huh7_Mock_HA568_dsRNA488_CS1_1.tif
   ```
   Only files whose name contains **`mock`** (case-insensitive) are processed; everything
   else in the folder is skipped. The token *positions* are remappable at startup, so other
   layouts work too (§2.3).
4. **Make a small test folder first.** Copy **~2 images per condition** into a separate
   folder, and **deliberately include images that contain bright artefacts/dust** — you need
   them to calibrate `ARTIFACT_UPPER_BOUND` (§4). Run the whole macro on the test folder,
   look at the QC, tune the config, *then* run the full set.

---

## 1. Launch

`Plugins ▸ Macros ▸ Edit…` → open `01_Mock_pipeline.ijm` → **Run** (`Ctrl/Cmd + R`).
A folder chooser appears → pick the folder with your `.tif` images.

The macro then shows **two startup dialogs** (§2), loops over every Mock image, and writes a
timestamped output folder (§5). When the Log prints `=== DONE ===` it is finished.

---

## 2. Startup dialogs (the modes)

### 2.1 Dialog 1 — "Pipeline mode" (3 choices)

| Field | Options | What it does |
|---|---|---|
| **Mode** | `filename` / `dialog` | **Where per-image metadata comes from.** `filename` parses timepoint + the two markers from the file name (fast, for well-named batches). `dialog` asks you per image (for messy/mixed names). |
| **Threshold mode** | `automatic` / `manual` | **How the cell & nucleus masks are thresholded.** `automatic` applies the configured auto-method silently. `manual` pauses on each mask with the Threshold slider open so you can check/adjust, then click OK (§3). |
| **Cell line** | `Huh7` / `VeroE6` / `Other cell line` | Recorded in every CSV row and used as the output filename prefix. Hook for cell-line-specific defaults (`applyCellLineDefaults`, §4). |

> **Tip — start in `automatic` + `filename`.** Switch a single image to `manual` only if the
> QC shows the mask is wrong; switch to `dialog` only if filenames don't parse.

### 2.2 Dialog 2 — "Pipeline setup" (the domain lists)

Three comma-separated lists, pre-filled with the workflow standards — **just click OK to keep
them**:

- **Markers** — `HA568, HA488, dsRNA488, NS4B568`
- **Combos** — `HA568_dsRNA488, NS4B568_dsRNA488, NS4B568_HA488` (order = `marker1_marker2`)
- **Timepoints** — `12h, 24h`

These lists are the single source of truth: they define which CSV files are created, the
options in the per-image dialog, and which images are validated/skipped. Add a new marker or
timepoint here and everything downstream follows. *(An image whose parsed combo or timepoint
is not in these lists is skipped with a Log note — not an error.)*

### 2.3 Dialog 3 — "Filename token mapping" (only in `filename` mode)

Shown once, using your first Mock file as the example. Radio-button groups let you pick **which
underscore token holds timepoint / marker1 / marker2** (e.g. `0: 12h`, `3: HA568`, `4:
dsRNA488`). Use this when your naming layout differs from the standard positions.
**Cell line is the dialog choice from §2.1 — it is never parsed from the name.**

---

## 3. Manual threshold mode (what the pauses are)

If you chose **Threshold mode = manual**, the macro pauses **twice per image**:

1. on the **cell mask** (built from the brightest available marker channel), and
2. on the **nucleus mask** (built from DAPI).

Each time, the **Threshold** window is open with the configured auto-method pre-applied.
Adjust the lower slider if the mask under-/over-segments, then click **OK** to continue. This
is the slow, high-control path — use it on a few representative images to find good settings,
then move those into the config and run `automatic` for the batch.

---

## 4. The settings (CONFIG block at the top of the `.ijm`)

All tunables live in the `// CONFIG` block. Edit them in the file (not mid-run). The ones you
will actually touch:

| Variable | Default | What it controls / when to change |
|---|---|---|
| `TOP_PCT` | `0.01` | The "Top-X %" pool size for `mean_top/median_top/std_top` (here 0.01 % of cytosol pixels). The downstream background uses the **percentile** columns, so this mostly affects the `*_top` diagnostics. |
| `CELL_THR_METHOD` | `Li` | Auto-threshold for the cell mask. **Use `Li` for dim VeroE6** — `Triangle` collapses to near-zero on dim data. |
| `NUC_THR_METHOD` | `Otsu` | Auto-threshold for the DAPI nucleus mask. |
| `CELL_THR_FACTOR` | `0.5` | Multiplies the auto threshold's lower bound. **< 1 = more permissive** (keeps dim cytoplasm — important for Vero periphery). Raise toward 1 if the mask leaks into background. |
| `BLUR_SIGMA_CELL` / `BLUR_SIGMA_NUC` | `1` / `1` | Gaussian pre-blur (µm, calibrated) for the masks. Larger = smoother/blobbier masks. |
| `NUC_CLOSE_ITER` | `2` | Morphological closing of the nucleus mask (Dilate→Fill→Erode). |
| `ARTIFACT_UPPER_BOUND` | `4000` | **Artefact cut:** a cytosol pixel brighter than this (any channel) is removed before measuring. Set **deliberately high** so it only removes extreme dust. **Calibrate it** from the logged `cyto_max_raw` (see below). |
| `ARTIFACT_DILATE_ITER` | `20` | Grows each artefact slightly to catch its bright rim. (Candidate to lower — 20 is aggressive.) |
| `MIN_PARTICLE_SIZE` | `200` | After artefact removal, cytosol fragments smaller than this (px) are dropped (coverslip dirt). Use ~`100` for the smaller VeroE6. |
| `SAVE_QC` / `SAVE_MASKS` / `SAVE_COMPOSITE_JPG` | `true` | Toggle the QC PNG / mask TIFs / composite JPG outputs. |

### Calibrating `ARTIFACT_UPPER_BOUND` (why the test folder needs artefacts)

For every image the Log prints, per channel, the **brightest raw cytosol pixel** before
artefact removal, and the same value is written to the CSV as **`cyto_max_raw`**:

```
max intensity   C1 HA568 : max=18432  (cytosol; ARTIFACT_UPPER_BOUND=4000)
```

Look across your test images: **clean** cells top out at a few thousand; an image **with
dust** spikes to ~15–20 k. Set `ARTIFACT_UPPER_BOUND` *above* the clean-cell maxima but
*below* the dust spikes. The philosophy is intentionally loose: this cut only needs to remove
the worst outliers, because the **particle-size filter** and the **cross-image median** (your
aggregation step) absorb the rest. Do **not** chase per-image perfection.

---

## 5. Output

Everything lands in `<input>/measure_mock/<RUN_ID>/` (`RUN_ID = YYYYMMDD_HHMM`, so reruns
never overwrite):

```
measure_mock/<RUN_ID>/
├── <CellLine>_mock_<tp>_<marker>_in_<combo>.csv   ← one CSV per (timepoint × combo × marker)
├── masks/   <RUN_ID>_<img>_{cell,nuc,cyto,artifact}.tif
├── qc/      <RUN_ID>_<img>_qc.png         (marker1 + cytosol outline)
│            <RUN_ID>_<img>_composite.jpg  (all 3 channels + cytosol outline)
└── macro_log_<RUN_ID>.txt                 (full IJ Log — provenance)
```

**CSV columns** (per image × marker): metadata (`image, cell_line, timepoint, combo,
channel`), the run mode, `threshold_value`, pixel counts (`n_top_pixels, n_cyto_pixels,
n_artifacts_excluded`), the top-pool stats (`mean_top, median_top, std_top`), the
**percentile ladder** `p95 … p99_9999`, `cyto_max_raw`, all the threshold parameters, and
`macro_version, run_id`. The percentile ladder is the important part — it is what you
aggregate next.

### QC checklist (look at a handful before trusting any CSV)

- **Cytosol outline follows the real cell shape** (not leaking into black background, not
  cutting into the cell).
- **Nuclei (blue) are outside the outline** — the cytosol = cell − nucleus.
- **No bright dust inside the outline** (a red/green dot inside the ROI = an artefact that
  slipped the cut → lower `ARTIFACT_UPPER_BOUND` or check `MIN_PARTICLE_SIZE`).

---

## 6. The aggregation step → background values for Stage 2

The CSVs are **per image**. To get the single background value Stage 2 needs:

1. Open all CSVs for one `(cell line × timepoint × combo × marker)` group (Excel / pandas / R).
2. Take the **`p99_9995`** column. *(Why this percentile: it is the upper tail of Mock
   autofluorescence — high enough to sit above the noise bulk, but below the very last few
   outlier pixels that `p99_9999`/`max` would chase. It is a robust "background ceiling".)*
3. Compute the **median of `p99_9995` across all images** in the group. *(Median, not mean,
   so one image with a stray bright spot can't move the value.)*
4. You end up with ~12 numbers (2 timepoints × 3 combos × 2 markers). These are exactly the
   fields in Stage 2's background dialog.

➡ Continue with [the Stage 2 USAGE](../2_SubBg_Coloc/USAGE.md).

---

## 7. Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| Image skipped: `combo not in ANALYSE_COMBI` / `timepoint not in TIMEPOINTS` | The parsed combo/timepoint isn't in your §2.2 lists, or token mapping (§2.3) is wrong. Fix the lists or the mapping. |
| Near-zero thresholds on VeroE6 | Don't use `Triangle`; keep `CELL_THR_METHOD = Li`, and lower `CELL_THR_FACTOR` / `MIN_PARTICLE_SIZE` for the smaller, dimmer cells. |
| Masks look **inverted** (foreground black) | A `BlackBackground` flip. The macro re-pins `setOption("BlackBackground", true)` per image; if you edited the mask code, keep that pin (HANDOFF §7). |
| `cytosol mask empty after cleanup` | Threshold too strict or artefact cut too low for this image — check the QC, raise `CELL_THR_FACTOR` toward 1 or `ARTIFACT_UPPER_BOUND`. |
| Background looks too high | An artefact slipped into the cytosol → it inflates the high percentiles. Lower `ARTIFACT_UPPER_BOUND` using `cyto_max_raw`, or rely on the cross-image median to absorb it. |

For the design rationale and the hard-won IJM pitfalls, see [`../HANDOFF.md`](../HANDOFF.md).
