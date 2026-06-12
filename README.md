# Colocalisation Pipeline — Fiji/ImageJ macros

Automated image-analysis pipeline for **confocal immunofluorescence colocalisation**
in virus-infected cells, built in the **Bartenschlager lab** (Dept. of Infectious
Diseases, Heidelberg). The macros measure a per-channel autofluorescence background on
**Mock** (uninfected) images, subtract it from **infected (MOI)** images, build per-cell
cytosol ROIs, and run the **Colocalisation Threshold** analysis (Manders / Pearson) per
cell.

This is research software intended for a publication: every run is timestamped and every
output row carries the macro version + run ID, so results are reproducible and auditable.

---

## The biology in one paragraph

Cells (**Huh7**, large; **VeroE6**, ~4× smaller and dimmer) are either **Mock**
(uninfected control), **MOI1** or **MOI5** (infected), fixed at **12 h** or **24 h**, and
stained for two viral/host markers plus DAPI. Markers: `HA568`, `HA488`, `dsRNA488`,
`NS4B568` (punctate replication-site signal; **dsRNA strongest, HA568 weakest**). **DAPI
is always channel 3.** Three marker pairs are analysed (order = `marker1`=C1, `marker2`=C2):
`HA568_dsRNA488`, `NS4B568_dsRNA488`, `NS4B568_HA488`. The question is *how much the two
markers colocalise* in the cytosol of infected cells — which requires first removing the
cell's intrinsic autofluorescence (measured on Mock) and restricting the measurement to a
clean, unbiased cytosol ROI.

---

## The three-stage pipeline

```
 .czi  ──Bio-Formats──▶  .tif   (Mock + MOI images mixed in one folder, single z-plane, 16-bit)
                          │
   ┌──────────────────────┴───────────────────────────────────────────────┐
   │ STAGE 1 — Mock background   (macros/1_Mock/01_Mock_pipeline.ijm)       │
   │   For each MOCK image: cytosol mask → Top-X% pixel statistics per      │
   │   marker channel → per-image CSVs.                                     │
   └──────────────────────┬───────────────────────────────────────────────┘
                          │   MANUAL aggregation (you, in Excel/pandas/R):
                          │   median of the `p99_9995` column per
                          │   (cell line × timepoint × combo × marker)
                          ▼
   ┌──────────────────────────────────────────────────────────────────────┐
   │ STAGE 2 — Subtract bg + Coloc  (macros/2_SubBg_Coloc/02_…_v0.14.3.ijm) │
   │   For each MOI image: subtract those bg values → build per-cell        │
   │   cytosol ROIs (auto watershed OR manual draw) → run Colocalisation    │
   │   Threshold per cell → one CSV row per ROI.                            │
   └──────────────────────┬───────────────────────────────────────────────┘
                          ▼
            STAGE 3 — downstream stats / figures  (analysis/, not part of the macros)
```

**The link between Stage 1 and Stage 2 is manual and deliberate:** you read the Stage 1
CSVs, take the **median of the `p99_9995` percentile** per group, and type those numbers
into Stage 2's background dialog. (Why `p99_9995` and median: see the Mock
[USAGE](macros/1_Mock/USAGE.md) — it is the robust upper tail of autofluorescence,
insensitive to a few stray bright pixels or one bad image.)

---

## Repository layout

```
colocalisation_project/
├── README.md                     ← you are here (project overview)
├── CITATION.cff · LICENSE.txt
├── macros/
│   ├── HANDOFF.md                ← full project context / state (start here when resuming dev)
│   ├── USAGE.md                  ← index → the two per-macro USAGE guides
│   ├── ToDo.md
│   ├── 1_Mock/
│   │   ├── 01_Mock_pipeline.ijm  ← STAGE 1 (current, v0.8.1; edited in place)
│   │   ├── README.md · USAGE.md · CHANGELOG.md
│   └── 2_SubBg_Coloc/
│       ├── 02_Subtract_background_coloc_v0.14.3.ijm  ← STAGE 2 (current — highest N)
│       ├── 02_…_v0.2 … v0.14.2.ijm                   ← history (one file per version)
│       ├── README.md · USAGE.md · CHANGELOG.md
├── analysis/                     ← downstream notebooks (Stage 3)
├── images/ · example_thomas/     ← data + supervisor reference (not the macros)
└── _issues/ · _templates/ · _tutorials/
```

**Versioning conventions** (please keep):
- **Stage 1 (Mock)** is edited **in place** — bump `MACRO_VERSION` in the header and add a
  `CHANGELOG.md` entry.
- **Stage 2 (Coloc)** keeps **one file per version** (`02_…_v0.N.ijm`); always run/edit the
  **highest N**, and update `CHANGELOG.md` (+ `README.md`/`USAGE.md` when behaviour changes).
- Every CSV row carries `macro_version` + `run_id`; outputs go in a `run_id`-named folder, so
  parallel runs never overwrite and pandas can filter by run.

---

## The two macros at a glance

| | Stage 1 — Mock | Stage 2 — SubBg + Coloc |
|---|---|---|
| **File** | `1_Mock/01_Mock_pipeline.ijm` (v0.8.1) | `2_SubBg_Coloc/02_…_v0.14.3.ijm` |
| **Input** | Mock `.tif` (filename contains `mock`) | MOI `.tif` (filename contains `moi`) |
| **Does** | Top-X% cytosol intensity stats per marker | subtract bg, per-cell ROIs, Colocalisation Threshold |
| **Output** | per-image CSVs + masks + QC | bg-subtracted TIFs + numbered QC + per-ROI coloc CSV |
| **Plugins** | none | **MorphoLibJ** (auto ROI) + **Colocalization Threshold** |
| **Guide** | [README](macros/1_Mock/README.md) · [USAGE](macros/1_Mock/USAGE.md) · [CHANGELOG](macros/1_Mock/CHANGELOG.md) | [README](macros/2_SubBg_Coloc/README.md) · [USAGE](macros/2_SubBg_Coloc/USAGE.md) · [CHANGELOG](macros/2_SubBg_Coloc/CHANGELOG.md) |

---

## Requirements

- **Fiji** (ImageJ ≥ 1.54). Download: <https://fiji.sc>.
- **Colocalization Threshold** plugin (Tony Collins / WCIF) — ships with standard Fiji,
  under `Analyze ▸ Colocalisation ▸ Colocalisation Threshold`. *(Needed only for Stage 2's
  coloc step.)*
- **MorphoLibJ** (IJPB-plugins update site) — needed only for Stage 2's **automatic** ROI
  mode (marker-controlled watershed). Install: `Help ▸ Update… ▸ Manage update sites ▸ tick
  IJPB-plugins ▸ Apply ▸ restart`. Stage 2 hard-checks for it and aborts with instructions if
  missing; **manual ROI mode and Stage 1 need no extra plugins.**
- Input images **calibrated** (Bio-Formats `.czi → .tif` keeps the µm/px calibration, which
  the watershed ring and scale bars rely on).

---

## Filename schema (binding, positional, underscore-separated)

```
<timepoint>_<cellLine>_<condition>_<marker1>_<marker2>_<coverslip>_<imgIndex>.tif
e.g.  12h_Huh7_Mock_HA568_dsRNA488_CS1_1.tif      24h_VeroE6_MOI5_NS4B568_HA488_CS2_3.tif
```

| token | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| field | timepoint | cell line | condition | marker1 → **C1** | marker2 → **C2** | coverslip | index |

Both macros let you **remap which token holds marker1 / marker2 / timepoint** at startup, so
other naming layouts work. **Cell line is chosen in a dialog, not parsed.** **C3 is always
DAPI.** Stage 1 only processes files containing `mock`; Stage 2 only files containing `moi`.

---

## Quick start

1. Install Fiji (+ MorphoLibJ if you'll use Stage 2 auto mode).
2. **Stage 1:** open `macros/1_Mock/01_Mock_pipeline.ijm` in Fiji → Run → point it at your
   image folder. Follow [the Mock USAGE](macros/1_Mock/USAGE.md). Aggregate the CSVs
   (median `p99_9995` per group).
3. **Stage 2:** open `macros/2_SubBg_Coloc/02_…_v0.14.3.ijm` → Run → type the Stage 1 medians
   into the background dialog. Follow [the Coloc USAGE](macros/2_SubBg_Coloc/USAGE.md).

Both macros are **batch** (loop over a whole folder) and write a timestamped output folder
plus a copy of the IJ Log for provenance.

---

## Status & contact

The macros are feature-complete and in the **optimisation + documentation** phase; the
automatic Stage-2 ROI parameters are still being tuned on real data (see
`macros/2_SubBg_Coloc/CHANGELOG.md` v0.14.x and `macros/HANDOFF.md`). The coloc **numbers
are not auto-parsed** from the plugin — the macro writes one provenance row per ROI and you
read/export the values from the plugin window (the numbered QC tells you which row is which
cell).

- **Author:** Kolja Hildenbrand — Kolja.Hildenbrand@gmail.com
- **Licence:** see [`LICENSE.txt`](LICENSE.txt) · **Citation:** see [`CITATION.cff`](CITATION.cff)
- **Dev context / handover:** [`macros/HANDOFF.md`](macros/HANDOFF.md)
