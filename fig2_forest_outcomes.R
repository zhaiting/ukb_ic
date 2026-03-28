# ============================================================================
# fig2_forest_outcomes.R — Figure 2
#
# Input:  data/res_disease.csv, data/res_mortality.csv,
#         data/disease_meta.csv, data/mort_meta.csv
# Output: output/fig2_forest_outcomes.svg/.png
#
# Prerequisites: source("00_helpers.R")
# ============================================================================

# ── Load ─────────────────────────────────────────────────────────────────────
dis_overall  <- read_csv(file.path(data_dir, "res_disease.csv"),
                         show_col_types = FALSE)
mort_overall <- read_csv(file.path(data_dir, "res_mortality.csv"),
                         show_col_types = FALSE)

disease_meta <- load_disease_meta()
mort_meta    <- load_mort_meta()

# ── Tidy ─────────────────────────────────────────────────────────────────────
dis_df <- tidy_dis(
  dis_overall %>% filter(analysis_set == "overall"), "Overall", disease_meta
)

outcome_lev <- disease_meta %>% arrange(system, label) %>% pull(label)
system_order <- levels(disease_meta$system)

dis_plot <- dis_df %>%
  mutate(
    score_label   = factor(score_label, levels = rev(score_order)),
    outcome_label = factor(label, levels = rev(outcome_lev)),
    system        = factor(system, levels = system_order),
    sig           = p.val < 0.05,
    ring_alpha    = as.numeric(as.character(score_label) == primary_domain)
  )

# ── Helper: two-panel forest ────────────────────────────────────────────────
make_forest_pair <- function(df, panel_title) {

  p_left <- ggplot(
    df %>% filter(as.character(score_label) == "IC composite"),
    aes(x = hr, y = outcome_label)
  ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_point(aes(alpha = sig), size = 2.2,
               color = score_colors[["IC composite"]]) +
    scale_x_continuous(
      trans  = compose_trans(log10_trans(), reverse_trans()),
      labels = label_number(accuracy = 0.01),
      breaks = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.18), guide = "none") +
    facet_grid(system ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "IC composite HR", y = NULL, title = panel_title) +
    theme_ic() + theme(legend.position = "none")

  domain_levels <- rev(score_order[score_order != "IC composite"])
  df_dom <- df %>%
    filter(as.character(score_label) != "IC composite") %>%
    mutate(score_label = factor(as.character(score_label), levels = domain_levels))

  p_right <- ggplot(df_dom, aes(x = hr, y = outcome_label, color = score_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_point(aes(alpha = sig), size = 2,
               position = position_dodgev(height = 0.05)) +
    geom_point(aes(alpha = I(ring_alpha)),
               shape = 21, size = 2.5, stroke = 0.85,
               color = "black", fill = NA,
               position = position_dodgev(height = 0.05)) +
    scale_x_log10(labels = label_number(accuracy = 0.01),
                  breaks = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0)) +
    scale_color_manual(
      values = score_colors[names(score_colors) != "IC composite"],
      name = "Domain score",
      guide = guide_legend(reverse = TRUE, nrow = 1)
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.18), guide = "none") +
    facet_grid(system ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Hazard ratio (log scale)", y = NULL) +
    theme_ic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y = element_blank(), legend.position = "bottom")

  p_left + p_right + plot_layout(widths = c(1, 4))
}

# ── Panels ───────────────────────────────────────────────────────────────────
fig2a <- make_forest_pair(dis_plot, "A")

mort_df <- tidy_mort(mort_overall, "Overall", mort_meta) %>%
  mutate(
    score_label   = factor(score_label, levels = rev(score_order)),
    outcome_label = factor(label, levels = rev(unique(label))),
    system        = domain_group,
    sig           = p.val < 0.05,
    ring_alpha    = as.numeric(as.character(score_label) == primary_domain)
  )

fig2b <- make_forest_pair(mort_df, "B")

fig2 <- fig2a / fig2b

save_fig(fig2, "fig2_forest_outcomes", width = 8, height = 18)
