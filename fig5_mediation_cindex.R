# ============================================================================
# fig5_mediation_cindex.R — Figure 5
#
# Inputs:
#   Supplementary Tables workbook, sheets ST14 and ST15
#   <derived_dir>/nmr_mediation_results.csv
#     effect_code, estimate, lower, upper, pval, outcome, outcome_label,
#     outcome_group, outcome_order, score_group
#
# Output: output/fig5_mediation_cindex.svg/.png
# ============================================================================

if (!exists("read_supp_table")) source("00_helpers.R")

nmr_mediation_path <- file.path(derived_dir, "nmr_mediation_results.csv")
require_input(nmr_mediation_path, "NMR mediation input")

mediation_columns <- c(
  "effect_code", "estimate", "lower", "upper", "pval", "outcome",
  "outcome_label", "outcome_group", "outcome_order", "score_group"
)
cindex_columns <- c(
  "model", "c_index_cv", "n", "events", "outcome", "outcome_label",
  "outcome_group", "outcome_order", "score_group", "platform"
)

olink_mediation <- read_supp_table("ST14", header_row = 2)
nmr_mediation <- read_csv(nmr_mediation_path, show_col_types = FALSE)
cindex_results <- read_supp_table("ST15", header_row = 2)

require_columns(olink_mediation, mediation_columns, "ST14")
require_columns(
  nmr_mediation, mediation_columns, "nmr_mediation_results.csv"
)
require_columns(cindex_results, cindex_columns, "ST15")

disease_exclude <- c(
  "neurodegenerative_other", "cvd_athero_other", "eye_retinal_other",
  "muscle_soft_tissue", "eye_retinal_vascular", "psych_cmd",
  "diabetes_any", "resp_chronic_lower"
)

system_order <- c(
  "Neurological", "Psychiatric", "Musculoskeletal", "Sensory",
  "Respiratory", "Cardiovascular", "Metabolic", "Renal/Hepatic"
)

prepare_mediation <- function(med_res) {
  pm <- med_res %>%
    filter(
      effect_code == "pm",
      score_group == "overall",
      !outcome %in% disease_exclude,
      !is.na(outcome)
    ) %>%
    transmute(
      outcome, outcome_label, outcome_group,
      outcome_order = as.numeric(outcome_order),
      pm = as.numeric(estimate),
      pm_pval = as.numeric(pval)
    )

  te <- med_res %>%
    filter(
      effect_code == "Rte",
      score_group == "overall",
      !outcome %in% disease_exclude,
      !is.na(outcome)
    ) %>%
    transmute(
      outcome,
      te_lower = as.numeric(lower),
      te_upper = as.numeric(upper)
    )

  if (anyDuplicated(pm$outcome) || anyDuplicated(te$outcome)) {
    stop("Mediation input contains duplicate overall-score outcomes.", call. = FALSE)
  }

  label_levels <- pm %>%
    arrange(outcome_order) %>%
    pull(outcome_label)

  out <- pm %>%
    left_join(te, by = "outcome") %>%
    mutate(
      system = factor(outcome_group, levels = system_order),
      label = factor(outcome_label, levels = rev(label_levels)),
      pm_sig = pm_pval < 0.05,
      te_near_null = te_lower <= 1 & te_upper >= 1
    ) %>%
    filter(!is.na(system), is.finite(pm))

  if (n_distinct(out$outcome) != 25) {
    stop(
      "Expected 25 curated outcomes in the mediation input; found ",
      n_distinct(out$outcome), ".",
      call. = FALSE
    )
  }

  out
}

prepare_cindex <- function(ci_res, platform_name) {
  model_levels <- c(
    "Age + sex" = "age_sex",
    "Phenotype score" = "phenotype",
    "Surrogate score" = "surrogate",
    "Combined" = "combined"
  )

  d <- ci_res %>%
    filter(
      platform == platform_name,
      score_group == "overall",
      !outcome %in% disease_exclude,
      !is.na(outcome)
    ) %>%
    mutate(
      model_key = recode(model, !!!model_levels, .default = NA_character_),
      c_index_cv = as.numeric(c_index_cv),
      outcome_order = as.numeric(outcome_order)
    ) %>%
    filter(!is.na(model_key), is.finite(c_index_cv))

  key_counts <- d %>% count(outcome, model_key)
  if (any(key_counts$n != 1)) {
    stop(
      platform_name,
      " C-index input contains duplicate outcome/model rows.",
      call. = FALSE
    )
  }

  expected_keys <- tidyr::expand_grid(
    outcome = unique(d$outcome),
    model_key = unname(model_levels)
  )
  missing_keys <- anti_join(
    expected_keys,
    distinct(d, outcome, model_key),
    by = c("outcome", "model_key")
  )
  if (nrow(missing_keys)) {
    stop(
      platform_name,
      " C-index input is missing one or more required models.",
      call. = FALSE
    )
  }

  outcome_meta <- d %>%
    arrange(outcome_order) %>%
    distinct(
      outcome, outcome_label, outcome_group, outcome_order, n, events
    )

  wide <- d %>%
    select(outcome, model_key, c_index_cv) %>%
    pivot_wider(
      names_from = model_key,
      values_from = c_index_cv,
      names_prefix = "cv_"
    )

  require_columns(
    wide,
    c(
      "cv_age_sex", "cv_phenotype", "cv_surrogate", "cv_combined"
    ),
    paste(platform_name, "C-index input")
  )

  label_levels <- outcome_meta %>%
    arrange(outcome_order) %>%
    pull(outcome_label)

  outcome_meta %>%
    left_join(wide, by = "outcome") %>%
    mutate(
      delta_phenotype = cv_phenotype - cv_age_sex,
      delta_surrogate = cv_surrogate - cv_age_sex,
      delta_combined = cv_combined - cv_age_sex
    ) %>%
    select(
      outcome, outcome_label, outcome_group, outcome_order, n, events,
      starts_with("delta_")
    ) %>%
    pivot_longer(
      cols = starts_with("delta_"),
      names_to = "comparison",
      values_to = "delta_cv"
    ) %>%
    mutate(
      comparison = recode(
        comparison,
        "delta_phenotype" = "Phenotype",
        "delta_surrogate" = "Surrogate",
        "delta_combined" = "Combined"
      ),
      comparison = factor(
        comparison,
        levels = c("Phenotype", "Surrogate", "Combined")
      ),
      system = factor(outcome_group, levels = system_order),
      label = factor(outcome_label, levels = rev(label_levels))
    ) %>%
    filter(!is.na(system), is.finite(delta_cv))
}

plot_mediation <- function(dat) {
  ggplot(dat, aes(x = "Mediation", y = label, fill = pm)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_point(
      data = filter(dat, pm_sig),
      shape = 21, size = 2.5, stroke = 0.55,
      fill = NA, colour = "black"
    ) +
    geom_text(
      data = filter(dat, te_near_null),
      aes(label = "\u00d7"),
      size = 4.2, colour = "black", alpha = 0.8
    ) +
    facet_grid(system ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(
      low = "#F2D479", high = "#C44A3A",
      limits = c(0, 1), oob = squish,
      name = "Proportion\nmediated"
    ) +
    labs(x = NULL, y = NULL) +
    theme_heatmap() +
    theme(
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key.height = unit(0.6, "cm"),
      legend.position = "bottom"
    )
}

plot_cindex <- function(dat) {
  ggplot(dat, aes(x = comparison, y = label, fill = delta_cv)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(system ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(
      low = "white", high = "#C44A3A",
      limits = c(0, 0.18), oob = squish,
      name = expression(Delta * " CV C-index")
    ) +
    labs(x = NULL, y = NULL) +
    theme_heatmap() +
    theme(
      axis.text.y = element_blank(),
      legend.key.height = unit(0.6, "cm"),
      legend.position = "bottom"
    )
}

make_platform_panel <- function(
  med_res, ci_res, platform_name, title, panel_tag
) {
  med_dat <- prepare_mediation(med_res)
  ci_dat <- prepare_cindex(ci_res, platform_name)

  title_grob <- grid::grobTree(
    grid::textGrob(
      panel_tag, x = grid::unit(0, "npc"), just = "left",
      gp = grid::gpar(fontsize = 15, fontface = "bold")
    ),
    grid::textGrob(
      title, x = grid::unit(0.5, "npc"),
      gp = grid::gpar(fontsize = 13, fontface = "bold")
    )
  )

  patchwork::wrap_plots(
    patchwork::wrap_elements(full = title_grob),
    plot_mediation(med_dat) + plot_cindex(ci_dat) +
      plot_layout(widths = c(1, 2.5)),
    ncol = 1,
    heights = c(0.06, 1)
  )
}

panel_olink <- make_platform_panel(
  olink_mediation, cindex_results, "Olink",
  "Protein (Olink) surrogate performance", "a"
)
panel_nmr <- make_platform_panel(
  nmr_mediation, cindex_results, "NMR",
  "Metabolite (NMR) surrogate performance", "b"
)

fig5 <- (panel_olink | panel_nmr) +
  plot_layout(widths = c(1, 1), guides = "collect")

save_fig(fig5, "fig5_mediation_cindex", width = 15, height = 9.5)
