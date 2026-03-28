
# IC Trajectory

## Overview

This repository contains the R scripts used to generate the main-text
figures. All input data are summary-level results provided in the
Supplementary Tables (Excel workbook) accompanying the manuscript.

## Requirements

R >= 4.2 with the following packages:

``` r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "forcats", "scales",
  "patchwork", "ggvenn", "ggbreak", "readr",
  "survival", "survminer", "corrplot"
))
```

## Usage

1.  Place the supplementary data CSVs in `data/` (see individual scripts
    for expected filenames, or export the relevant sheets from the
    Supplementary Tables workbook).
2.  Source the helpers first, then run each figure script:

``` r
source("00_helpers.R")
source("fig1_trajectories.R")
source("fig2_forest_outcomes.R")
source("fig3_molecular.R")
source("fig4_mediation_cindex.R")
```

Figures are saved as SVG + PNG in `output/`.

## File Manifest

| Script | Figure |
|------------------------------------|-------------------------------------|
| `00_helpers.R` | Shared constants, palettes, themes, utilities |
| `fig1_trajectories.R` | **Figure 1** — Trajectories, KM survival, domain correlation |
| `fig2_forest_outcomes.R` | **Figure 2** — Forest plots |
| `fig3_molecular.R` | **Figure 3** — Venn diagrams, feature counts, GSEA enrichment |
| `fig4_mediation_cindex.R` | **Figure 4** — Mediation + C-index heatmaps |

## Data

All input files are summary-level results exported from the
Supplementary Tables. Each script header lists the expected CSV
filenames and columns.
