# ============================================================================
# fig4_mediation_cindex.R — Figure 4
#
# Input:  data/mediation_results.csv, data/omic_cindex.csv,
#         data/outcome_meta.csv
# Output: output/fig4_mediation_cindex.svg/.png
#
# Prerequisites: source("00_helpers.R")
# ============================================================================

# ── Load outcome metadata ───────────────────────────────────────────────────
outcome_meta <- read_csv(file.path(data_dir, "outcome_meta.csv"),
                         show_col_types = FALSE) %>%
  mutate(
    system = factor(system, levels = unique(system)),
    label  = factor(label, levels = rev(label))
  )

disease_exclude <- read_csv(file.path(data_dir, "disease_exclude.csv"),
                            show_col_types = FALSE)$outcome

# ── Load results ─────────────────────────────────────────────────────────────
med_results    <- read_csv(file.path(data_dir, "mediation_results.csv"),
                           show_col_types = FALSE)
cindex_results <- read_csv(file.path(data_dir, "omic_cindex.csv"),
                           show_col_types = FALSE)

# ── Prepare mediation ────────────────────────────────────────────────────────
prepare_mediation <- function(med_res) {
  pm <- med_res %>%
    filter(effect_code == "pm", score_group == "overall",
           !outcome %in% disease_exclude) %>%
    select(outcome, estimate, lower, upper, pval) %>%
    rename(pm = estimate, pm_lower = lower, pm_upper = upper, pm_pval = pval)

  te <- med_res %>%
    filter(effect_code == "Rte", score_group == "overall",
           !outcome %in% disease_exclude) %>%
    transmute(outcome, te = estimate, te_lower = lower, te_upper = upper)

  pm %>%
    left_join(te, by = "outcome") %>%
    inner_join(outcome_meta, by = "outcome") %>%
    mutate(te_near_null = te_lower <= 1 & te_upper >= 1,
           pm_sig = pm_pval < 0.05)
}

# ── Panel A: Mediation heatmap ───────────────────────────────────────────────
med_dat <- prepare_mediation(med_results)

fig4a <- ggplot(med_dat, aes(x = "Overall IC", y = label, fill = pm)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_point(data = subset(med_dat, pm_sig),
             shape = 21, size = 2.5, stroke = 0.5, fill = NA, colour = "black") +
  geom_text(data = subset(med_dat, te_near_null),
            aes(label = "\u00d7"), size = 4.5, colour = "black", alpha = 0.8) +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#F2D479", high = "#C44A3A",
                      limits = c(0, 1), oob = squish,
                      name = "Proportion\nmediated") +
  labs(x = NULL, y = NULL) +
  theme_heatmap() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        legend.key.height = unit(0.8, "cm"))

# ── Prepare C-index ──────────────────────────────────────────────────────────
prepare_cindex <- function(ci_res) {
  d <- ci_res %>%
    filter(score_group == "overall", !outcome %in% disease_exclude) %>%
    select(outcome, model, c_index_cv, n, events) %>%
    pivot_wider(names_from = model, values_from = c_index_cv, names_prefix = "cv_")

  if ("cv_phenotype_score" %in% names(d))
    d$delta_phenotype <- d$cv_phenotype_score - d$cv_age_sex
  if ("cv_surrogate_score" %in% names(d))
    d$delta_surrogate <- d$cv_surrogate_score - d$cv_age_sex
  if ("cv_combined" %in% names(d))
    d$delta_combined  <- d$cv_combined - d$cv_age_sex

  d %>%
    select(outcome, n, events, starts_with("delta_")) %>%
    pivot_longer(cols = starts_with("delta_"), names_to = "comparison",
                 values_to = "delta_cv") %>%
    mutate(comparison = recode(comparison,
                               delta_phenotype = "Phenotype",
                               delta_surrogate = "Surrogate",
                               delta_combined  = "Combined"),
           comparison = factor(comparison,
                               levels = c("Phenotype", "Surrogate", "Combined"))) %>%
    inner_join(outcome_meta, by = "outcome")
}

# ── Panel B: C-index heatmap ────────────────────────────────────────────────
ci_dat <- prepare_cindex(cindex_results)

fig4b <- ggplot(ci_dat, aes(x = comparison, y = label, fill = delta_cv)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#DCEAF4", mid = "white", high = "#C44A3A",
                       midpoint = 0,
                       name = expression(Delta * " CV C-index")) +
  labs(x = NULL, y = NULL) +
  theme_heatmap() +
  theme(legend.key.height = unit(0.8, "cm"))

# ── Combine ──────────────────────────────────────────────────────────────────
fig4 <- fig4a + fig4b + plot_layout(widths = c(1, 2.2))

save_fig(fig4, "fig4_mediation_cindex", width = 12, height = 10)
