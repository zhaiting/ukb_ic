# ============================================================================
# fig4a_omic_breadth_depth.R — Figure 4a
#
# Reduced-model feature counts (top) and incremental model performance
# (bottom).
#
# Input:  Supplementary Tables workbook, sheet ST12
#           required columns: IC Domain, platform, N features (reduced), delta_r2
# Output: output/fig4a_omic_breadth_depth.svg/.png
# ============================================================================

if (!exists("read_supp_table")) source("00_helpers.R")

library(ggbreak)

st12 <- read_supp_table("ST12", header_row = 2)
require_columns(
  st12,
  c("IC Domain", "delta_r2", "N features (reduced)", "platform"),
  "ST12"
)

# Platform display names, in legend order.
platform_labels <- c(
  "Biochemistry"     = "Biochem",
  "CBC"              = "CBC",
  "NMR metabolomics" = "NMR",
  "Olink proteomics" = "Olink"
)
platform_order <- unname(platform_labels)

platform_colors <- c(
  Biochem = "#C9A227",
  CBC     = "#6FA86B",
  NMR     = "#2E9E8F",
  Olink   = "#8E6FB0"
)

# Figure 4a reports the five IC domains; the composite is not shown.
domain_labels <- c(
  "Cognition"     = "Cognitive",
  "Locomotion"    = "Loco",
  "Psychological" = "Psych",
  "Sensory"       = "Sensory",
  "Vitality"      = "Vitality"
)
domain_order <- unname(domain_labels)

omic_df <- st12 %>%
  filter(
    !is.na(platform),
    `IC Domain` %in% names(domain_labels),
    platform %in% names(platform_labels)
  ) %>%
  transmute(
    domain = factor(unname(domain_labels[`IC Domain`]), levels = domain_order),
    platform = factor(unname(platform_labels[platform]), levels = platform_order),
    feature_count = as.numeric(`N features (reduced)`),
    delta_r2_pct = as.numeric(delta_r2) * 100
  )

if (nrow(omic_df) != length(domain_order) * length(platform_order)) {
  stop(
    "Expected ", length(domain_order) * length(platform_order),
    " domain x platform cells in ST12, found ", nrow(omic_df), ".",
    call. = FALSE
  )
}

dodge <- position_dodge(width = 0.8)

p_count <- ggplot(omic_df, aes(domain, feature_count, fill = platform)) +
  geom_col(position = dodge, width = 0.72) +
  geom_text(
    aes(
      y = if_else(feature_count >= 700, 725, feature_count),
      label = feature_count,
      vjust = if_else(feature_count >= 700, 0.5, -0.35),
      color = if_else(feature_count >= 700, "white", "black"),
      group = platform
    ),
    position = dodge, size = 2.7, show.legend = FALSE
  ) +
  scale_color_identity() +
  scale_fill_manual(values = platform_colors, name = "OmicType") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_y_break(c(220, 700), scales = 0.35, space = 0.15) +
  labs(
    title = "Number of Retained Omic Features per IC Domain",
    x = NULL, y = "Features in reduced model"
  ) +
  theme_ic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
  )

p_delta <- ggplot(omic_df, aes(domain, delta_r2_pct, fill = platform)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_col(position = dodge, width = 0.72) +
  geom_text(
    aes(label = sprintf("%.1f", delta_r2_pct)),
    position = dodge, vjust = -0.35, size = 2.6
  ) +
  scale_fill_manual(values = platform_colors, name = "OmicType") +
  scale_y_continuous(
    breaks = seq(0, 20, 5),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(x = "IC Domain", y = expression("Incremental " * Delta * R^2 * " (%)")) +
  theme_ic() +
  theme(legend.position = "bottom")

fig4a <- aplot::plot_list(p_count, p_delta, ncol = 1, heights = c(1.15, 1))

fig4a_w <- 8.4
fig4a_h <- 8.2
svglite::svglite(file.path(fig_dir, "fig4a_omic_breadth_depth.svg"),
                 width = fig4a_w, height = fig4a_h, bg = "white")
print(fig4a)
dev.off()
png(file.path(fig_dir, "fig4a_omic_breadth_depth.png"),
    width = fig4a_w, height = fig4a_h, units = "in", res = 300, bg = "white")
print(fig4a)
dev.off()
message("Saved: fig4a_omic_breadth_depth (.svg + .png)")
