# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **Fiji/ImageJ Macro Language (IJM)** pipeline for fluorescence microscopy image analysis, used in virology research (Bartenschlager Lab). It quantifies background signal and colocalisation between fluorescent markers in confocal images of virus-infected cells. There is no traditional build system — the macros run directly inside the Fiji desktop application.

## Running the Macros

Both macros are run interactively inside Fiji:
1. Open Fiji (≥ 2.16.0, bundled ImageJ ≥ 1.54p)
2. `Plugins → Macros → Edit…` → open the target `.ijm` file
3. Press `Ctrl+R` (or click **Run**) — a series of GUI dialogs will guide the user through configuration

**Script 1 — Mock Background Estimation** (`macros/1_Mock/01_Mock_pipeline.ijm`, v0.6.1):
- Processes only files whose name contains `"mock"`
- Dialogs collect: pipeline mode (filename-parsed vs dialog), threshold mode, cell line (Huh7 | VeroE6)
- Then prompts for a folder of `.tif` images
- Outputs per-image CSVs, QC PNGs, composite JPGs, and mask TIFs to `<INPUT_DIR>/measure_mock/<RUN_ID>/`

**Script 2 — Background Subtraction + Colocalisation** (`macros/2_SubBg_Coloc/02_Subtract_background_coloc.ijm`, v0.1.0):
- Processes only files whose name contains `"moi"`
- Dialogs collect: ROI mode (auto_central | manual_draw), background values for each marker/combo/timepoint combination
- Then prompts for a folder of `.tif` images
- Outputs background-subtracted TIFs, QC JPGs, colocalisation CSV, and `background_values_used.md` to `<INPUT_DIR>/<RUN_ID>_bgsub/`

**Optional post-analysis** (Python/Jupyter):
```bash
cd analysis
jupyter notebook mock_bg_comparison.ipynb
```

## Two-Script Pipeline Architecture

The pipeline is sequential — Script 2 depends on outputs from Script 1:

```
Raw .czi confocal images
        ↓  (Bio-Formats → 16-bit 3-channel .tif)
Script 1: Mock Background Estimation
        ↓  Per-image CSVs with intensity percentiles (p95–p99.9999)
Manual step: Aggregate CSVs → derive one background value per marker/combo/timepoint
        ↓  User enters these background values into Script 2 dialogs
Script 2: Background Subtraction + Colocalisation
        ↓  Colocalisation CSVs + background-subtracted images
Jupyter notebook: Aggregate statistics and figures
```

Background values from Script 1 are **manually aggregated** (median of p99.9995 percentile across replicates) before being entered into Script 2 — there is no automated handoff between the two scripts.

## Image & Domain Conventions

**Input format:** 16-bit, 3-channel, single z-plane `.tif` (confocal, exported via Bio-Formats from `.czi`)

**Channel assignment:** C1 = marker 1, C2 = marker 2, C3 = DAPI (nucleus; never background-subtracted)

**Filename schema (positional, underscore-separated):**
```
<timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<coverslip>_<imgIndex>.tif
# e.g.: 12h_Huh7_Mock_HA568_dsRNA488_CS1_1.tif
```
Scripts parse this schema to look up per-image metadata. Files that don't match are skipped with a logged warning (batch-robust; one bad file never aborts the run).

**Supported cell lines:** Huh7 (large cells) and VeroE6 (~4× smaller, dimmer)
**Markers:** HA568, HA488, dsRNA488, NS4B568
**Analysed combos:** HA568_dsRNA488, NS4B568_dsRNA488, NS4B568_HA488
**Timepoints:** 12h, 24h
**Conditions:** Mock (Script 1), MOI1 / MOI5 (Script 2)

Cell-line-specific parameter tuning is done through `applyCellLineDefaults(cellLine)` — this changes threshold strictness and mask sizes, not strategy.

## Code Structure & Conventions

Each macro follows this layout:
1. **CONFIG block** at the top — all tunable parameters as `var` globals (single source of truth)
2. **MAIN** — short orchestrator (one call per phase)
3. **Named functions** — each does one job, preceded by a `// ----` docstring header (what it does, inputs, returns, assumptions)
4. **Domain lists** (MARKERS, ANALYSE_COMBI, TIMEPOINTS) — defined once in CONFIG, reused everywhere

**Reproducibility columns in every CSV row:** `macro_version`, `run_id` (YYYYMMDD_HHMM), threshold parameters used.

**Output directory naming:**
- Script 1: `<INPUT_DIR>/measure_mock/<RUN_ID>/` with subdirectories `masks/` and `qc/`
- Script 2: `<INPUT_DIR>/<RUN_ID>_bgsub/`

## IJM Language Pitfalls

IJM has several non-obvious quirks — violating these causes silent wrong results or crashes:

- **No `i++` or `+=`** — use `i = i + 1`, `x = x + delta`
- **No ternary operator** — use `if / else`
- **Return from user functions via intermediate variable** — `return userFunc(args)` triggers a type-inference bug; assign to a temp var first, then return it
- **`var` only at top level** — declaring `var` inside a function shadows the global silently
- **`File.makeDirectory` is not recursive** — create parent directories explicitly before children
- **16-bit histogram indexing** — bin index `i` equals pixel value directly; do not scale
- **8-bit histogram** — use the 3-argument form: `getHistogram(values, counts, 256)`
- **"Black Background" Fiji option** — can flip between runs; pin it explicitly with `setOption("BlackBackground", true/false)` whenever binary operations are used
- **Plugin name spelling** — use `Plugins > Macros > Record` to get the exact string; colocalisation plugin is spelled `"Colocalization Threshold"` (American English) in macro calls
- **`Dialog.create` is modal** — it blocks the Fiji GUI; use `waitForUser` for non-modal pause points
- **`Array.concat` stores numbers as strings** — call `parseFloat()` when numeric operations follow

## Current Status & Known TODOs

See `macros/HANDOFF.md` for the full context. Key open items at time of writing:

- Script 2 (v0.1.0) is the **first functional version but has not yet been tested in Fiji** with real data
- Coloc plugin Log-parsing labels need verification against actual plugin output
- Single-cell ROI strategy (`auto_central` vs `manual_draw`) is pending discussion with supervisor (Thomas)
- USAGE.md documentation for both scripts is incomplete

The authoritative reference for architecture decisions, design rationale, and known edge cases is **`macros/HANDOFF.md`**.
