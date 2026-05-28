# Mock Top-X% Pipeline

Fiji/ImageJ macro that measures the **Top X% brightest pixel statistics inside the
cytosol** of each marker channel in Mock immunofluorescence images. The per-image
values feed into the downstream background-subtraction script (`2_SubBg/`).

> **Script:** `03_Mock_top1_pipeline_v0.3.ijm`
> **Version:** `0.3.1` · **Status:** in active use

---

## What it does

For each `Mock` `.tif` in an input folder, the macro:

1. Parses metadata from the filename (or asks via dialog).
2. Builds a **cytosol mask** = (cell mask − nucleus mask) from one cytoplasmic
   marker and DAPI.
3. Inside that cytosol, takes a **histogram** of each marker channel and computes
   statistics over the top X % brightest pixels (mean, median, std, several
   percentiles, plus the threshold value at which the top-X% pool starts).
4. Writes one row per (image × marker) into a per-(timepoint × combo) CSV.
5. Optionally saves the three binary masks (TIF) and a QC overlay (PNG) for each
   image, so the segmentation can be reviewed visually.

These CSVs are aggregated (median over images) downstream to derive per-channel
background values used in `2_SubBg/`.

---

## Requirements

| Tool | Version | Notes |
|---|---|---|
| Fiji | ≥ 2.16.0 | bundled ImageJ ≥ 1.54p |

No additional plugins, no update sites required.

---

## Input format

### Folder

A folder containing Mock (and optionally MOI) `.tif` images. Only files whose
filename contains `mock` (case-insensitive) are processed; everything else is
skipped automatically.

### Filename schema (positional, strict)

```
<timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<coverslip>_<imgIndex>.tif
```

Example:
```
12h_Huh7_Mock_HA568_dsRNA488_CS1_1.tif
```

Token positions are fixed:

| Position | Token | Allowed values (default config) |
|---|---|---|
| 0 | timepoint | `12h`, `24h` |
| 1 | cell line | `Huh7`, `VeroE6` |
| 2 | condition | `Mock`, `MOI1`, `MOI5` |
| 3 | marker 1 (→ C1) | `HA568`, `HA488`, `dsRNA488`, `NS4B568` |
| 4 | marker 2 (→ C2) | same set as marker 1, must differ |
| 5 | coverslip | starts with `CS` |
| 6 | image index | integer |

**C3 is always DAPI** by acquisition convention.

---

## Step-by-step usage

1. Open Fiji.
2. **Plugins → Macros → Edit…** → open
   `macros/1_Mock/03_Mock_top1_pipeline_v0.3.ijm`.
3. Click **Run** (or `Cmd/Ctrl + R`).
4. **Folder dialog**: pick the folder with your `.tif` images.
5. **Pipeline mode dialog** asks three things:

   | Question | Options | Meaning |
   |---|---|---|
   | Mode | `filename` / `dialog` | parse from filename, or ask per image |
   | Threshold modus | `automatic` / `manual` | use auto-threshold, or pause per image to adjust the slider |
   | Cell line | `Huh7` / `VeroE6` | applied to all images in the run |

6. **(dialog mode only)** A second dialog lets you edit the marker/combo/timepoint
   lists for this specific run. Comma-separated, whitespace is trimmed.

7. The macro loops through every Mock image, prints progress to the **Log**
   window, and writes outputs into a freshly created subfolder of the input
   folder.

8. When done, the Log shows `=== DONE ===`. Open any of the resulting CSV files
   in Excel / pandas / R to inspect.

### Manual threshold mode — what to expect

When `Threshold modus = manual`, the macro pauses twice per image:
- once on the cell mask, with the Threshold dialog open (adjust slider, click
  OK on the "Check and adjust…" dialog),
- once on the nucleus mask (same).

After your OK click, the macro applies the chosen threshold and continues.
Use this mode for a few test images, then switch to `automatic` for the full run.

---

## Output

Inside the chosen input folder:

```
measure_mock/
└── <RUN_ID>/                              # e.g. 20260526_1535
    ├── <CellLine>_mock_<tp>_<marker>_in_<combo>.csv
    │   ...  (12 CSV files = 2 timepoints × 3 combos × 2 markers)
    ├── masks/
    │   ├── <RUN_ID>_<imgname>_cell.tif
    │   ├── <RUN_ID>_<imgname>_nuc.tif
    │   └── <RUN_ID>_<imgname>_cyto.tif
    └── qc/
        └── <RUN_ID>_<imgname>_qc.png      # marker1 + cytosol outline
```

### CSV columns

| Column | Meaning |
|---|---|
| `image`, `cell_line`, `timepoint`, `combo`, `channel` | metadata of the row |
| `stat_method`, `top_pct` | bookkeeping (`median_top_hist`, `1.0`) |
| `macro_mode`, `threshold_mode` | what mode was used for this run |
| `threshold_value` | pixel intensity at which the top-X% pool starts |
| `n_top_pixels`, `n_cyto_pixels` | size of the top pool and the cytosol |
| `mean_top`, `median_top`, `std_top` | statistics of the top pool |
| `p95 … p99_999` | whole-cytosol percentiles (sanity / outlier detection) |
| `cell_thr_method`, `nuc_thr_method`, `blur_sigma_*` | mask parameters |
| `macro_version`, `run_id` | provenance |

The `run_id` column lets you mix outputs of multiple runs in the same CSV and
filter them later in pandas:
```python
df[df.run_id == "20260526_1535"]
```

### QC PNG — what to check

The PNG shows the marker-1 channel auto-contrasted, with a yellow outline of the
cytosol ROI. Verify visually:

1. The outline **follows the cells** (no big chunks missed, no empty background).
2. **Nuclei are excluded** (no closed loops around the dark central holes).
3. **No obvious debris** is inside the outline.

If the outline is consistently wrong → tune the threshold / blur parameters
(see the next section).

---

## Configuration: what to change when samples change

All tunable values live in the `// ============== 1. CONFIG` block at the top of
the script. Change them **there** (never inside functions).

| Change | Variable to edit | Example |
|---|---|---|
| **New cell line** (e.g. A549) | `Cell Line:` radio in `askModeAndConfig` (≈ line 123) | add `"A549"` to the array, adjust the `rows, cols` of the radio group |
| **New timepoint** (e.g. 48h) | `TIMEPOINTS` | `newArray("12h", "24h", "48h")` |
| **New marker** (e.g. NS5A647) | `MARKERS` **and** `ANALYSE_COMBI` | add the marker to MARKERS, add new combos (e.g. `"NS5A647_dsRNA488"`) to ANALYSE_COMBI |
| **Different top-X%** | `TOP_PCT` | `0.5` for top 0.5%, `5.0` for top 5% |
| **Different cell threshold method** | `CELL_THR_METHOD` | `"Triangle"` (more permissive), `"Yen"` (more conservative) |
| **Different nucleus threshold method** | `NUC_THR_METHOD` | `"Triangle"` if Otsu cuts off |
| **More / less smoothing on cell mask** | `BLUR_SIGMA_CELL` | `2` for more smoothing (in µm), `0.5` for less |
| **Stronger / weaker nucleus closing** | `NUC_CLOSE_ITER` | `3` for stubborn boundary gaps, `0` to disable |
| **Don't save masks / QC** | `SAVE_MASKS`, `SAVE_QC` | `false` for fast batch runs |

After any change, bump `MACRO_VERSION` (e.g. `"0.3.2"`) so the provenance column
distinguishes runs made with the new settings.

### Adding a new channel order

If a future dataset has C1 = DAPI instead of C3 = DAPI, edit
`splitAndRenameChannels` (≈ line 372). Only that one function — the rest of the
script references markers by name, not channel number.

### Adding columns to the CSV

Three coordinated changes, all in this file:
1. `CSV_HEADER` (line 62) — add the new column name in the right position.
2. `computeTopStats` — compute and return the new value as part of the stats
   array.
3. `appendCsvRow` — concatenate the new value at the matching position.

Out-of-sync changes silently misalign columns. After editing, run a test image
and open the resulting CSV in Excel to verify the header / data align.

---

## Common issues

| Symptom | Likely cause | Fix |
|---|---|---|
| Macro silently skips every image | Filenames don't contain `mock` (case-insensitive) | rename to follow the schema |
| `SKIP: combo not in ANALYSE_COMBI` | A `marker1_marker2` combo isn't whitelisted | add it to `ANALYSE_COMBI` or rename the file |
| `SKIP: cytosol mask empty` | Cell mask and nucleus mask had identical foreground | use Triangle for cell mask, or check thresholding visually (manual mode) |
| Cell mask too small (cells cut off) | Threshold too strict | switch `CELL_THR_METHOD` to `Triangle`, or lower `BLUR_SIGMA_CELL` |
| Cell mask too large (background included) | Threshold too lenient | switch `CELL_THR_METHOD` to `Yen` or `Otsu` |
| Nuclei have "bites" (open holes) | Boundary gaps prevent Fill Holes | increase `NUC_CLOSE_ITER` to `3` or `4` |
| `Array expected` error in Log | IJM type quirk (rare) — usually `return userFunction(...)` somewhere | report; fix is intermediate variable pattern |

When in doubt, switch to `Threshold modus = manual` for a few images and watch
the masks build live — the Log window shows the cell-mask source, threshold
values used, and skip reasons line by line.

---

## Reproducibility note

Every CSV row records `macro_version` and `run_id`, plus the mask parameters in
effect (`cell_thr_method`, `blur_sigma_cell`, etc.). To reproduce an analysis
six months from now:

1. Note the `run_id` and `macro_version` from the CSV row of interest.
2. Check out the matching version of the script from git.
3. Re-run on the same images — the values should be identical (deterministic
   pipeline).

The `RUN_ID` is also part of the output folder name, so multiple parallel runs
on the same input never overwrite each other.

---

## File locations in this repo

```
macros/
└── 1_Mock/
    ├── 03_Mock_top1_pipeline_v0.3.ijm        # the macro
    └── README.md                              # this file
```

