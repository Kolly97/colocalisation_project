# Usage Guide

Step‑by‑step instructions for running the `colocalisation_pipeline` macros on your data.
Aimed at a first‑time user (e.g. a colleague or supervisor receiving the macro).

---

## 1. Prerequisites

- **Fiji** ≥ 2.16.0 (bundled ImageJ ≥ 1.54p). Download: https://imagej.net/software/fiji/downloads
- **Disk space**: ~2× your image folder for intermediate masks/QC.
- **No additional plugins** required for `v0.1` (Mock pipeline).

To check your version: in Fiji, **Help → About ImageJ**.

## 2. Recommended local folder layout

Keep your images, the macro repo, and the outputs side‑by‑side, but **only the repo
goes to GitHub**:

```
~/Bartenschlager_project/                    ← parent (NOT a git repo)
├── colocalisation_pipeline/                 ← the cloned git repo
│   ├── macros/01_mock_top_pct.ijm
│   └── ...
├── data/                                    ← raw TIFs (never committed)
│   ├── rep1/Huh7_rep1/converted_tif/*.tif
│   └── rep2/...
└── output/                                  ← results land here
    └── 2026-05-07_Huh7_12h/
        ├── measure_mock/*.csv
        ├── measure_mock/qc/*.png
        └── measure_mock/masks/*.tif
```

Avoid placing this tree inside iCloud Drive — Git's internals can be evicted from
the local cache and corrupt the repo.

## 3. Filename schema (strictly enforced)

Each `.tif` must be named:

```
<timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<CS>_<imgIdx>.tif
```

| Token | Allowed values | Example |
|---|---|---|
| `timepoint` | `12h`, `24h`, ... | `12h` |
| `cellLine` | `Huh7`, `VeroE6`, ... | `Huh7` |
| `condition` | `Mock`, `MOI1`, `MOI5` | `Mock` |
| `marker1` | e.g. `HA568`, `dsRNA488`, `NS4B568`, `HA488` | `HA568` |
| `marker2` | (same domain as marker1) | `dsRNA488` |
| `CS` | coverslip ID, must start with `CS` | `CS1` |
| `imgIdx` | running integer | `1` |

Full example: `24h_Huh7_Mock_HA568_dsRNA488_CS1_2.tif`

**The order of channels in the file matters**: the macro assumes
`C1 = marker1`, `C2 = marker2`, `C3 = DAPI`. Verify in Fiji once
(File → Open → check channel order) before running on a new dataset.

## 4. Running the Mock pipeline (`v0.1`)

### A. From the GUI

1. Launch Fiji.
2. **File → New → Script…** (or `Plugins → Macros → Edit…`).
3. Open `colocalisation_pipeline/macros/01_mock_top_pct.ijm`.
4. Click **Run** (or `Ctrl/Cmd + R`).
5. A folder dialog appears — choose the folder containing your Mock + MOI TIFs.
6. The macro filters for `Mock`‑containing filenames automatically.
7. Watch the **Log** window for progress (`[3/30] 12h_Huh7_Mock_...`).

### B. Headless (advanced, optional)

```bash
/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx \
  --headless \
  --console \
  -macro /path/to/macros/01_mock_top_pct.ijm \
  /path/to/data/folder
```

(Headless requires reading `INPUT_DIR` from a macro argument — currently the script
uses `getDirectory()`. Convert is on the roadmap for `v0.3`.)

## 5. Output

Inside the chosen image folder, a `measure_mock/` subdirectory is created:

```
measure_mock/
├── Mock_12h_HA568_in_dsRNA488_HA568.csv     ← one CSV per (timepoint, combo, marker)
├── Mock_12h_dsRNA488_in_dsRNA488_HA568.csv
├── Mock_24h_HA568_in_dsRNA488_HA568.csv
├── ... (12 CSVs total: 2 timepoints × 3 combos × 2 markers)
├── masks/
│   ├── <imgname>_cell.tif                    ← binary cell mask (8-bit, 0/255)
│   ├── <imgname>_nuc.tif                     ← binary nucleus mask
│   └── <imgname>_cyto.tif                    ← binary cytosol mask
└── qc/
    └── <imgname>_qc.png                      ← marker1 channel + cytosol outline
```

### CSV columns

| Column | Meaning |
|---|---|
| `image` | image filename without extension |
| `cell_line`, `timepoint`, `combo`, `channel` | parsed from filename |
| `stat_method`, `top_pct` | bookkeeping (e.g. `median_top_hist`, `1.0`) |
| `threshold_value` | pixel intensity at which the top‑X% pool starts (P_(100−X)) |
| `n_top_pixels` | number of pixels in the top‑X% pool |
| `n_cyto_pixels` | total pixels in the cytosol mask |
| `mean_top`, `median_top`, `std_top` | statistics of the top pool |
| `p95, p99, p99_25, p99_5, p99_9` | whole‑cytosol percentiles (sanity / outlier detection) |
| `cell_thr_method`, `nuc_thr_method`, `blur_sigma_*` | mask parameters |
| `macro_version`, `run_id` | provenance — never edit by hand |

Open in Excel / pandas / R. Provenance columns let you mix outputs of different
runs in one CSV and filter by `run_id` later.

### QC PNG — what to look for

The PNG shows the marker‑1 channel (auto‑contrasted) with the **cytosol selection
drawn as a yellow outline**. Three things to verify visually:

1. The outline **follows the cells** (no large empty regions covered, no obvious cells missed).
2. **Nuclei are excluded** (outline does not enclose dark central holes).
3. **No big out‑of‑focus / debris** is included.

If the outline is consistently wrong → adjust `BLUR_SIGMA_CELL` or
`CELL_THR_METHOD` in the macro CONFIG block (see §6).

## 6. Tuning parameters

All tunable values live in the `// CONFIG` block at the top of the macro:

```ijm
var CELL_THR_METHOD = "Li";        // try: Li, Otsu, Triangle, Yen
var NUC_THR_METHOD  = "Otsu";      // try: Otsu, Triangle
var BLUR_SIGMA_CELL = 1;           // 0.5–3 typical
var BLUR_SIGMA_NUC  = 1;
var TOP_PCT         = 1.0;         // try 0.1, 0.5, 1.0, 5.0
```

**Heuristic guide:**

| Symptom | Likely cause | First fix |
|---|---|---|
| Cell mask too small (cells cut off) | threshold too strict | try `Triangle` or lower `BLUR_SIGMA_CELL` |
| Cell mask too large (background included) | threshold too lenient | try `Yen` or `Otsu`, raise `BLUR_SIGMA_CELL` |
| Nuclei not removed cleanly | DAPI threshold too lenient | try `Triangle` for nucleus, check DAPI staining |
| `n_top_pixels < 10` warnings | cytosol mask too small | smaller cells (Vero?) → reduce expected min cell area, or check Cell mask |

A systematic tuning mode (grid over methods × `top_pct`, separate output folder per
combination) is planned for `v0.2`.

## 7. Common errors

**`Bad filename schema: ...`**
The file does not split into ≥ 7 underscore‑separated tokens. Rename your files.

**`SKIP: combo not in ANALYSE_COMBI`**
The `marker1_marker2` combo is not in the configured list. Either rename the file
or extend `ANALYSE_COMBI` at the top of the macro.

**`SKIP: cytosol mask empty`**
Cell mask and nucleus mask had identical foreground — check thresholding visually
(set `setBatchMode(false)` to see intermediate windows).

**Macro silently does nothing**
Ensure files actually contain `mock` (case‑insensitive) in their name.

## 8. Reproducibility

Every CSV row carries `macro_version` and `run_id`, so you can:

- Re‑run with new parameters and append to the same CSVs without losing history —
  filter by `run_id` in pandas: `df[df.run_id == "20260507_1422"]`.
- Cite a specific Zenodo DOI of the macro version used in your paper (see README).
- Share CSVs with collaborators who can verify the analysis without re‑running it.

## 9. Getting help

- Check the **Log** window in Fiji — every skip/warn is logged there.
- Open an issue on GitHub: https://github.com/<your-user>/colocalisation_pipeline/issues
- Include: macro version, Fiji version, an example filename, and the Log output.
