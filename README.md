# README for 01_Mock_background_pipeline_ijm

Fiji/ImageJ macro that measures the **Top X% brightest pixel statistics inside the
cytosol** of each marker channel in Mock immunofluorescence images. The per‑image
percentile values are aggregated downstream into the background values used by
the second script (`2_SubBg_Coloc/`).

> **Script (current):** `01_Mock_pipeline.ijm`
> **Version:** `0.5.0` · **Status:** active
> **Version history:** see [`CHANGELOG.md`](CHANGELOG.md) for all previous
> versions and the changes between them. Older script files
> (Version 0.1 - 0.4.1.) are kept in the repository on GitHub for traceability.

---

## Where this script sits in the workflow

```
   raw .tif (Mock + MOI mixed)
            │
            │  ── only Mock images are processed ──
            ▼
   ┌──────────────────────────────────┐
   │ THIS SCRIPT  (01_Mock_pipeline)  │  ← top‑X% statistics in cytosol
   └──────────────────────────────────┘
            │
            │  per‑image CSVs (one row per image × marker)
            ▼
   manual aggregation (see "What to do with the output" below)
            │
            │  one background value per (cell line × tp × combo × marker)
            ▼
   ┌──────────────────────────────────┐
   │ 2_SubBg_Coloc/ — subtract bg     │  ← values typed into its startup dialog
   └──────────────────────────────────┘
```

---

## What it does (in one paragraph)

For each `Mock`‑containing `.tif` in the chosen input folder, the macro splits
the channels, builds **cell**, **nucleus** and **cytosol** masks, cleans the
cytosol mask of bright artefacts and small disconnected fragments, and then
measures inside that cytosol the **histogram of each marker channel**. From
the histogram it computes the top‑X% statistics (mean, median, std, threshold
value at which the top‑X% pool starts) plus a series of whole‑cytosol
percentiles (p95 … p99.9999) used for sanity / outlier detection. One row per
(image × marker) is appended to a per‑(timepoint × combo) CSV.

For visual review it also saves the three binary masks (cell / nucleus /
cytosol / artefact) as TIFs, a single‑channel **QC PNG** with the cytosol
outline overlay, and an **8‑bit RGB composite JPG** with all three channels +
cytosol outline for cross‑channel verification.

---

## Requirements

| Tool | Version | Notes |
|---|---|---|
| Fiji | ≥ 2.16.0 | bundled ImageJ ≥ 1.54p |

No additional plugins or update sites required.

---

## Input format

A folder with `.tif` images. Only files whose filename contains `mock`
(case‑insensitive) are processed; everything else is skipped automatically.

Filename schema (positional, strict):
```
<timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<coverslip>_<imgIndex>.tif
```
Example: `12h_Huh7_Mock_HA568_dsRNA488_CS1_1.tif`

| Position | Token | Allowed defaults |
|---|---|---|
| 0 | timepoint | `12h`, `24h` |
| 1 | cell line | `Huh7`, `VeroE6` |
| 2 | condition | `Mock`, `MOI1`, `MOI5` |
| 3 | marker 1 (→ C1) | `HA568`, `HA488`, `dsRNA488`, `NS4B568` |
| 4 | marker 2 (→ C2) | same set, must differ from marker 1 |
| 5 | coverslip | starts with `CS` |
| 6 | image index | integer |

**C3 is always DAPI** by acquisition convention.

---

## Quick start

1. Open Fiji.
2. **Plugins → Macros → Edit…** → open `macros/1_Mock/01_Mock_pipeline.ijm`.
3. **Run** (`Cmd/Ctrl + R`).
4. Folder dialog → pick the folder with your `.tif` images.
5. Mode dialog → choose:
   - **Mode**: `filename` (parse from name) or `dialog` (ask per image)
   - **Threshold mode**: `automatic` or `manual`
   - **Cell line**: `Huh7` or `VeroE6` (applies tuned defaults)
6. Macro loops over all Mock images and writes outputs into
   `<input>/measure_mock/<RUN_ID>/` (one subfolder per run, timestamped).
7. When the Log shows `=== DONE ===`, the CSVs, masks, QC PNGs, composite JPGs
   and a copy of the Log are all in that subfolder.

---

## What to do with the output (manual aggregation step)

The CSVs are **per‑image** values. To use them as a background value in the next
script, you have to aggregate them per `(cell line × timepoint × combo × marker)`
group. Recommended procedure:

1. Open the relevant CSV in Excel / pandas / R. Each row is one Mock image.
2. **Pick the `p99_9995` column** as the background value source.
   - This is the 99.9995% percentile of the cytosol pixel intensities — the
     upper tail of Mock autofluorescence, but robust to the very last few
     outlier pixels (which `p99.9999` or `max` would be sensitive to).
   - Below this percentile is "typical Mock background"; above it is mostly
     artefacts and the very last noise spikes. Subtracting this from infected
     samples gives a conservative, defensible background floor.
3. Compute the **median of `p99_9995` across all images** in the group.
   - Median (not mean) because a single bad image with a stray bright spot
     should not move the threshold.
4. Note that median value for each `(cell line, timepoint, combo, marker)`.

You'll end up with one table — typically 12 entries for our current setup
(2 timepoints × 3 combos × 2 markers per combo).

5. Start the second script (`02_background_subtraction_coloc_pipeline.ijm`). When the **Background values** dialog opens, type these median values into the matching fields.
   That dialog is built dynamically from the same `TIMEPOINTS × ANALYSE_COMBI`
   lists, so the field labels match the CSV groups one‑to‑one.

---

## Output structure (brief)

```
measure_mock/
└── <RUN_ID>/                                  # e.g. 20260527_2243
    ├── <CellLine>_mock_<tp>_<marker>_in_<combo>.csv
    │   ...   (one CSV per timepoint × combo × marker)
    ├── masks/
    │   ├── <RUN_ID>_<imgname>_cell.tif
    │   ├── <RUN_ID>_<imgname>_nuc.tif
    │   ├── <RUN_ID>_<imgname>_cyto.tif
    │   └── <RUN_ID>_<imgname>_artifact.tif
    ├── qc/
    │   ├── <RUN_ID>_<imgname>_qc.png          # single channel + cytosol outline
    │   └── <RUN_ID>_<imgname>_composite.jpg   # 3‑channel RGB + cytosol outline
    └── macro_log_<RUN_ID>.txt
    └── mock_analysis/
    		└── <CellLine>_mock_background_analysis.csv
```

The QC files exist to verify visually that the cytosol mask is correct — that
artefacts are excluded, that nuclei are excluded, that the outline follows the
real cell shape. **Always look at a handful** before trusting a CSV.

---

## Configuration: what to change when samples change

All tunable parameters live in the `// CONFIG` block at the top of the script.

| Change | Variable to edit |
|---|---|
| New cell line (e.g. A549) | add option to `Cell line` radio in `askModeAndConfig`, add an `else if` branch in `applyCellLineDefaults` |
| New timepoint (e.g. 48h) | `TIMEPOINTS` |
| New marker | `MARKERS` **and** `ANALYSE_COMBI` (add at least one combo using it) |
| Different top‑X% | `TOP_PCT` |
| Different threshold method | `CELL_THR_METHOD`, `NUC_THR_METHOD` |
| Stricter / more permissive cell mask | `CELL_THR_FACTOR` (lower = more permissive) |
| Artefact upper cutoff | `ARTIFACT_UPPER_BOUND` (per‑cell‑line in `applyCellLineDefaults`) |
| Min disconnected region size | `MIN_PARTICLE_SIZE` |

After any change, bump `MACRO_VERSION` so the provenance column in the CSV
distinguishes runs made with the new settings, and add the change to
`CHANGELOG.md`.

---

## Note on documentation

This README is a **high‑level overview**. A separate `USAGE.md` is planned with
detailed step‑by‑step instructions. For further questions contact Kolja Hildenbrand

> [!IMPORTANT]
>
> Kolja.Hildenbrand@gmail.com

