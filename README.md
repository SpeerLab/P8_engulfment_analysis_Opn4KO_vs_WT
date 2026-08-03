# Microglial Engulfment Analysis — WT vs *Opn4* KO

Code accompanying the manuscript's analysis of microglial engulfment of presynaptic
input (VGluT2) and anterograde tracer (CTB), imaged by confocal microscopy. This
repository contains the full image-analysis and plotting pipeline, plus a
reproducible notebook that generates the SuperPlots, including the breakdown by
sex provided in the peer-review response.

## Blinding

Analysis was performed by an independent observer with **no knowledge of the
experimental grouping**. Genotypes were coded as **A** and **B** throughout the
analysis, and imaging channels were kept separate and unlabelled during
processing. Decoded assignments:

| Coded label | Actual group |
|---|---|
| A | Wild-type (WT) |
| B | *Opn4* knockout (*Opn4* KO) |

| Channel | Marker | Role |
|---|---|---|
| red | VGluT2 | engulfed presynaptic puncta |
| green | microglia | engulfing volume (engulfment reference) |
| blue | CTB | engulfed anterograde tracer |

Engulfment = fraction of microglial (green) volume that also contains VGluT2 (red) or
CTB (blue) signal.

## Contents

```
.
├── README.md                        This file
├── Opn4_Engulfment_Analysis.ipynb   Notebook: SuperPlots from CSV (WT vs KO + sex)
├── engulfment_plots.py              Standalone Python module behind the notebook
├── requirements.txt                 Python dependencies
├── data/
│   └── engulfment_results.csv       Per-imaging-field batch output (analysis result)
├── figures/                         Generated figures (PNG + vector PDF + EPS)
└── scripts/
    ├── SimplePunctaAnalysis_v5.m    Image analysis: tune + batch (MATLAB)
    └── SuperPlot_FromCSV_v3.m       SuperPlots, WT vs KO (MATLAB)
```

## The pipeline at a glance

```
raw confocal stacks (red / green / blue TIFF series per imaging field)
        │
        │  SimplePunctaAnalysis_v5.m   (MODE = 'tune'  -> set parameters)
        │  SimplePunctaAnalysis_v5.m   (MODE = 'batch' -> process everything)
        v
engulfment_results.csv   +   QC images, binary masks, MIPs, montages
        │
        ├──>  SuperPlot_FromCSV_v3.m           MATLAB SuperPlots (WT vs KO)
        │
        └──>  Opn4_Engulfment_Analysis.ipynb   Python SuperPlots (WT vs KO, plus
                                               a descriptive breakdown by sex)
```

---

## 1. Image analysis - `scripts/SimplePunctaAnalysis_v5.m`

MATLAB tool with two modes.

**Tune mode** (`MODE = 'tune'`): loads one representative imaging field and opens
an interactive figure with sliders for the per-channel threshold percentiles and
object-size limits. Adjust until the segmentation overlays look right, preview the
engulfment with *Calculate 3D Stats*, then *Save & Close* to write parameters to
`tuned_parameters.mat`.

**Batch mode** (`MODE = 'batch'`): recursively finds every imaging field under the
data root, applies the saved parameters, and writes results plus quality-control
outputs.

### Processing steps per field
1. **Load** three channels - `red/` (VGluT2), `green/` (microglia), `blue/` (CTB) -
   each a numbered TIFF z-series.
2. **Normalize**: per-image linear contrast stretch between the
   `stretch_low`/`stretch_high` percentiles of non-zero voxels -> 0-255.
3. **Blur**: 3-D Gaussian, `blur_sigma` pixels.
4. **Threshold**: adaptive per-channel percentile (`*_percentile`). Fixed and Otsu
   methods are also available.
5. **Size-filter**: 26-connected components in 3-D, kept within
   `*_min_size_um3`/`*_max_size_um3`.
6. **Quantify**: engulfment = fraction of microglial (green) volume also positive in
   VGluT2 (red) or CTB (blue), plus object counts, engulfed volumes, mean object
   volume, and object density.

### Key parameters (edit at the top of the script)
| Parameter | Meaning |
|---|---|
| `stretch_low`, `stretch_high` | normalization percentiles |
| `blur_sigma` | Gaussian blur (pixels) |
| `green_percentile`, `green_min_size_um3` | microglia threshold + minimum object volume |
| `red_percentile`, `red_min_size_um3`, `red_max_size_um3` | VGluT2 threshold + size gate |
| `ctb_percentile`, `ctb_min_size_um3`, `ctb_max_size_um3` | CTB threshold + size gate |

Voxel dimensions are read from CZI metadata via Bio-Formats when available,
otherwise from TIFF resolution tags. Set `BFMATLAB_DIR` to your Bio-Formats path.

### Expected input folder structure
```
<root>/
├── Sample_A/                     (coded group A = WT)
│   ├── Sample_1A/                (animal / biological replicate)
│   │   ├── <imaging_field_1>/{red,green,blue}/*.tif
│   │   └── <imaging_field_2>/...
│   └── Sample_2A/ ...
└── Sample_B/                     (coded group B = Opn4 KO)
    └── Sample_1B/ ...
```
The first folder level under `Sample_A`/`Sample_B` (e.g. `Sample_1A`) is the
**animal** (biological replicate); folders below it are imaging fields.

### Outputs (batch)
- `engulfment_results.csv` - one row per imaging field (the file in `data/`).
- `replicate_averages.csv` - per-animal means.
- `summary_report.txt`, `final_statistics.txt`, `methods_section.txt` - linear
  mixed model (`Response ~ Group + (1|Animal)`) and a drop-in methods paragraph.
- `QC/` overlays and histograms; per-field masks and MIPs; per-sample montages.

### `engulfment_results.csv` column dictionary

One row per imaging field. Channel columns are named `color_marker_*` so both the
acquisition channel and the biological marker are explicit: **red = VGluT2**,
**green = microglia** (the engulfing reference), **blue = CTB**.

| Column | Meaning |
|---|---|
| `voxel_vol_um3` | voxel volume (\u00b5m\u00b3) |
| `image_volume_um3` | total imaged volume of the field (\u00b5m\u00b3) |
| `green_microglia_threshold_used` | intensity threshold applied to the microglia channel |
| `red_VGluT2_threshold_used` | intensity threshold applied to the VGluT2 channel |
| `blue_CTB_threshold_used` | intensity threshold applied to the CTB channel |
| `green_microglia_volume_um3` | segmented microglial volume (\u00b5m\u00b3) |
| `green_microglia_voxels` | segmented microglial volume (voxels) |
| `red_VGluT2_total_objects` | count of all VGluT2 objects in the field |
| `red_VGluT2_total_volume_um3` | total VGluT2 volume (\u00b5m\u00b3) |
| `red_VGluT2_engulfed_objects` | VGluT2 objects overlapping microglia |
| `red_VGluT2_engulfed_volume_um3` | VGluT2 volume inside microglia (\u00b5m\u00b3) |
| `red_VGluT2_engulfment_pct` | engulfed VGluT2 as % of microglial volume |
| `red_VGluT2_mean_object_volume_um3` | mean volume per engulfed VGluT2 object |
| `red_VGluT2_density_per_1000um3` | VGluT2 object density (objects / 1000 \u00b5m\u00b3) |
| `blue_CTB_total_volume_um3` | total CTB volume (\u00b5m\u00b3) |
| `blue_CTB_engulfed_objects` | CTB objects overlapping microglia |
| `blue_CTB_engulfed_volume_um3` | CTB volume inside microglia (\u00b5m\u00b3) |
| `blue_CTB_engulfment_pct` | engulfed CTB as % of microglial volume |
| `blue_CTB_mean_object_volume_um3` | mean volume per engulfed CTB object |
| `dataset` | imaging-field folder name |
| `group` | blinded group code (A = WT, B = *Opn4* KO) |
| `replicate` | animal ID (biological replicate) |
| `path` | source folder path |

---

## 2. SuperPlots - `scripts/SuperPlot_FromCSV_v3.m`

Reads `engulfment_results.csv` and renders the SuperPlots comparing WT (A) vs
*Opn4* KO (B), with individual imaging fields, per-animal means, a box plot of
animal means, and mean +/- SD. Vector EPS/PDF plus PNG.

---

## 3. Notebook - `Opn4_Engulfment_Analysis.ipynb`

A self-contained Python notebook that reads the same CSV and regenerates clean
SuperPlots. It first shows the primary **WT vs *Opn4* KO** comparison for every
metric, then, in a dedicated section, a **breakdown by sex** within each genotype.

The sex breakdown is included because limited cohort availability meant several
female mice were used in the final analysis, and the notebook makes that
composition explicit. Those plots are **descriptive only** - the female subgroups
are small (WT: 2 females; *Opn4* KO: 1 female), so no formal test of a sex
difference is performed. Animal and field counts are printed on every column.

### Running it
```bash
pip install -r requirements.txt
jupyter notebook Opn4_Engulfment_Analysis.ipynb
```
Then *Run All*. Genotype and sex assignments live in editable dictionaries
(`GROUP_LABELS`, `SEX_MAP`) at the top of `engulfment_plots.py`; the loader asserts
that every animal in the CSV is mapped, so nothing drops silently. An interactive
cell lets you switch metric/layout via dropdowns, and a final cell exports every
figure to `figures/` as 300-dpi PNG plus vector PDF and EPS.

> **Vector format note:** EPS is provided for consistency with the MATLAB
> workflow, but EPS cannot store transparency, so the semi-transparent field
> points and faint field-to-mean connecting lines are rendered opaque in the EPS
> files. For figure assembly (e.g. Illustrator), prefer the **PDF**, which
> preserves transparency; use EPS only if your workflow specifically requires it.

`engulfment_plots.py` contains the same logic as an importable module if you'd
rather script it than use the notebook.

---

## Requirements

**MATLAB** R2025a with Image Processing and Statistics & Machine Learning
Toolboxes; Bio-Formats for MATLAB (reading voxel sizes from CZI). Optional:
Computer Vision Toolbox (montage text labels).

**Python** 3.9+ - see `requirements.txt` (`pandas`, `numpy`, `matplotlib`,
`nbformat`, `jupyter`, optional `ipywidgets`).

## Reproducibility

Given `data/engulfment_results.csv`, the notebook regenerates all figures
deterministically (the beeswarm layout is computed, not randomized). From raw
images, run the MATLAB tune -> batch stages to regenerate the CSV first.

## Citation / contact

If you use this code, please cite the accompanying manuscript. Questions:
Colenso M. Speer (cspeer@umd.edu).
