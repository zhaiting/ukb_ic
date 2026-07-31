# ============================================================================
# fig6_omics_validation.R — Cross-cohort omics validation
#
# Inputs:
#   Supplementary Tables workbook, sheet ST17
#   <derived_dir>/validation_olink_concordance.csv
#   <derived_dir>/validation_metabolite_concordance.csv
#     Both files: domain, ukb_effect, validation_effect
#
# Output: output/fig6_omics_validation.svg/.png
# ============================================================================

if (!exists("read_supp_range")) source("00_helpers.R")

olink_path <- file.path(derived_dir, "validation_olink_concordance.csv")
metabolite_path <- file.path(derived_dir, "validation_metabolite_concordance.csv")
require_input(olink_path, "Validation Olink concordance input")
require_input(metabolite_path, "Validation metabolite concordance input")

olink <- read_csv(olink_path, show_col_types = FALSE)
metabolite <- read_csv(metabolite_path, show_col_types = FALSE)

required_concordance <- c("domain", "ukb_effect", "validation_effect")
require_columns(olink, required_concordance, "validation_olink_concordance.csv")
require_columns(
  metabolite, required_concordance, "validation_metabolite_concordance.csv"
)

validation_domain_order <- c(
  "IC composite", "Vitality", "Locomotion",
  "Cognitive", "Psychological", "Sensory"
)

clean_concordance <- function(df, standardize = FALSE) {
  out <- df %>%
    transmute(
      domain = recode(
        tolower(domain),
        "ic" = "IC composite",
        "overall ic" = "IC composite",
        "vitality" = "Vitality",
        "locomotion" = "Locomotion",
        "cognition" = "Cognitive",
        "cognitive" = "Cognitive",
        "psychology" = "Psychological",
        "psychological" = "Psychological",
        "sensory" = "Sensory"
      ),
      ukb_effect = as.numeric(ukb_effect),
      validation_effect = as.numeric(validation_effect)
    ) %>%
    filter(
      domain %in% validation_domain_order,
      is.finite(ukb_effect), is.finite(validation_effect)
    )

  if (standardize) {
    out <- out %>%
      group_by(domain) %>%
      mutate(
        ukb_effect = ukb_effect / sd(ukb_effect),
        validation_effect = validation_effect / sd(validation_effect)
      ) %>%
      ungroup()
  }

  out %>%
    mutate(domain = factor(domain, levels = validation_domain_order))
}

make_concordance_plot <- function(df, point_size, point_alpha, axis_label) {
  rho_df <- df %>%
    group_by(domain) %>%
    summarise(
      rho = cor(ukb_effect, validation_effect, method = "spearman"),
      .groups = "drop"
    )

  ggplot(df, aes(ukb_effect, validation_effect, color = domain)) +
    geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey88", linewidth = 0.3) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "grey60", linewidth = 0.35
    ) +
    geom_point(alpha = point_alpha, size = point_size, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.75) +
    geom_text(
      data = rho_df,
      aes(x = -Inf, y = Inf, label = sprintf("\u03c1 = %.2f", rho)),
      inherit.aes = FALSE,
      hjust = -0.12, vjust = 1.3, size = 3.4
    ) +
    facet_wrap(~ domain, scales = "free", nrow = 1) +
    scale_color_manual(
      values = score_colors[validation_domain_order], guide = "none"
    ) +
    labs(
      x = paste(axis_label, "effect estimate"),
      y = "Validation-cohort effect estimate"
    ) +
    theme_ic() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "plain")
    )
}

olink_plot_df <- clean_concordance(olink, standardize = FALSE)
metabolite_plot_df <- clean_concordance(metabolite, standardize = TRUE)

p_olink <- make_concordance_plot(
  olink_plot_df, point_size = 0.45, point_alpha = 0.28,
  axis_label = "UK Biobank"
)

# ── Aggregate discovery-to-validation summary from ST17 ──────────────────────
st17 <- read_supp_range("ST17", "A2:F8")
require_columns(
  st17,
  c("Domain", "UKB-discovered BH-FDR<0.05 (n)",
    "Replication rate (%)", "Sign-concordance (%)"),
  "ST17"
)

replication_df <- st17 %>%
  transmute(
    domain = recode(
      Domain,
      "IC composite" = "IC composite",
      "Cognition" = "Cognitive"
    ),
    discovered = readr::parse_number(`UKB-discovered BH-FDR<0.05 (n)`),
    replication_rate = as.numeric(`Replication rate (%)`),
    sign_concordance = as.numeric(`Sign-concordance (%)`)
  ) %>%
  mutate(domain = factor(domain, levels = validation_domain_order))

p_replication <- ggplot(
  replication_df,
  aes(domain, replication_rate, fill = domain)
) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.1f%%", replication_rate)),
    vjust = -0.35, size = 3.2
  ) +
  scale_fill_manual(values = score_colors[validation_domain_order]) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 60, 20),
    labels = label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = NULL, y = "Replication rate") +
  theme_ic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_metabolite <- make_concordance_plot(
  metabolite_plot_df, point_size = 1.7, point_alpha = 0.72,
  axis_label = "UK Biobank standardized"
) +
  labs(y = "Validation-cohort standardized effect estimate")

fig6_validation <- p_olink / p_replication / p_metabolite +
  plot_layout(heights = c(1.1, 0.8, 1.1))

save_fig(fig6_validation, "fig6_omics_validation", width = 15, height = 10.8)
