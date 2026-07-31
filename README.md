
# Intrinsic Capacity Figure Code

Author: Ting Zhai

## Overview

This repository contains R scripts for the main and supporting figures,
including the external-validation analyses.

Input data are supplied separately and are not stored in this repository.

## Requirements

R >= 4.2 with the following packages:

``` r
install.packages(c(
  "aplot", "corrplot", "dplyr", "ggbreak", "ggplot2", "ggplotify",
  "ggvenn", "lme4", "mgcv", "patchwork", "readr", "readxl", "scales",
  "stringr", "survival", "survminer", "svglite", "tibble", "tidyr"
))
```

## Inputs

The supplementary workbook and source data are not distributed in this
repository. The scripts run when users supply the workbook and CSV files
listed below. Runtime checks stop with an explicit message when a required
file or column is absent.

Set the supplementary workbook path:

``` sh
export UKB_IC_SUPP_XLSX="/path/to/Supp_Tables_R1.xlsx"
```

Set the directories containing the original figure inputs and the added
validation inputs:

``` sh
export UKB_IC_DATA_DIR="/path/to/original_figure_inputs"
export UKB_IC_DERIVED_DIR="/path/to/validation_figure_inputs"
```

If these variables are not set, both locations default to `data/`.

### Supplementary workbook tables

| Script | Sheet and required columns |
|---|---|
| `fig2.R` | ST3: `outcome_class`, `endpoint`, `exposure`, `model`, `stratum`, `hr`, `hr_lo`, `hr_hi`, `p_fdr`; ST6: `Model`, `Source`, `Lifestyle factor`, `IC domain`, `beta`, `se`, `p_fdr` |
| `fig3_validation.R` | ST8: `Outcome`, `Model 2 HR (95% CI)`; ST9: `Outcome`, `White`, `Black`, `Hispanic`, `Chinese`; ST10: `Outcome`, `IC`, `Frailty index`, `Fried` |
| `fig4a_omic_breadth_depth.R` | ST12: `IC Domain`, `platform`, `N features (reduced)`, `delta_r2` |
| `fig5_mediation_cindex.R` | ST14: `effect_code`, `estimate`, `lower`, `upper`, `pval`, `outcome`, `outcome_label`, `outcome_group`, `outcome_order`, `score_group`; ST15: `model`, `c_index_cv`, `n`, `events`, `outcome`, `outcome_label`, `outcome_group`, `outcome_order`, `score_group`, `platform` |
| `fig6_omics_validation.R` | ST17: `Domain`, `UKB-discovered BH-FDR<0.05 (n)`, `Replication rate (%)`, `Sign-concordance (%)` |
| `supp_fig3_validation.R` | ST18: `Outcome`, `IC level HR (95% CI)`, `Prior-slope HR (95% CI)`; ST10: columns beginning `DNAm clock` and `IC adds over clock` |

The validation scripts use the fixed table ranges specified in their source
files. Estimate-and-confidence-interval fields are read in the displayed
`estimate (lower–upper)` form.

### CSV files

| Script | File | Required columns or structure |
|---|---|---|
| `fig1.R` | `ic_trajectories.csv` | `participant_id`, `domain`, `score`, `age`, `sex`, `visit` |
|  | `ic_survival.csv` | `ic_score`, `time_years`, `event` |
|  | `ic_domain_corr.csv` | Square, finite numeric correlation matrix with column headers |
| `fig3_validation.R` | `validation_longitudinal.csv` | `participant_id`, `exam`, `age`, `sex`, `domain`, `score` |
|  | `validation_participant_outcomes.csv` | `participant_id`, `baseline_ic`, `followup_years`, `death`, `baseline_condition_count`, `incident_any` |
| `fig4_molecular.R` | `omic_features.csv` | `platform`, `domain`, `feature` |
|  | `gsea_results.csv` | `pathway`, `domain`, `database`, `NES`, `p_adjust` |
| `fig5_mediation_cindex.R` | `nmr_mediation_results.csv` | `effect_code`, `estimate`, `lower`, `upper`, `pval`, `outcome`, `outcome_label`, `outcome_group`, `outcome_order`, `score_group` |
| `fig6_omics_validation.R` | `validation_olink_concordance.csv` | `domain`, `ukb_effect`, `validation_effect` |
|  | `validation_metabolite_concordance.csv` | `domain`, `ukb_effect`, `validation_effect` |
| `supp_fig3_validation.R` | `validation_domain_correlations.csv` | `domain_x`, `domain_y`, `correlation` |
|  | `validation_benchmark_associations.csv` | `domain`, `benchmark`, `family`, `estimate`, `significant` |
|  | `validation_attrition.csv` | `exam`, `age_stratum`, `retained_percent` |

Numeric estimates, standard errors, confidence limits, correlations, and
percentages must be stored as numeric values. Binary event fields use `0` and
`1`; `significant` is logical. Sex is labelled `Female` or `Male`. Scripts
accept `Overall IC` or `IC` for the composite and `Cognition` for the cognitive
domain where normalization is needed.

The four participant-level trajectory and outcome files are intended to be
prepared by authorized users in their approved analysis setting. The
repository contains no participant records or simulated replacements.

## Run

Run any script from the repository root; it will load `00_helpers.R`
automatically.

``` sh
Rscript fig1.R
Rscript fig2.R
Rscript fig3_validation.R
Rscript fig4a_omic_breadth_depth.R
Rscript fig4_molecular.R
Rscript fig5_mediation_cindex.R
Rscript fig6_omics_validation.R
Rscript supp_fig3_validation.R
```

Outputs are written to `output/` in the formats listed in each script header.
Input files and rendered figures are not included in the repository.

## File manifest

| Script | Figure |
|---|---|
| `00_helpers.R` | Shared paths, palettes, themes, parsers, and input checks |
| `fig1.R` | Figure 1b, 1c, and 1e code; Supplementary Figure 1f |
| `fig2.R` | Figure 2a–c |
| `fig3_validation.R` | Figure 3 |
| `fig4a_omic_breadth_depth.R` | Figure 4a |
| `fig4_molecular.R` | Selected Figure 4 molecular panels |
| `fig5_mediation_cindex.R` | Figure 5 |
| `fig6_omics_validation.R` | Figure 6 |
| `supp_fig3_validation.R` | Supplementary Figure 3 |

Figure 1a and the final multipanel figures are assembled separately.
