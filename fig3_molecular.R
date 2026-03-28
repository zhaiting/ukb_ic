# ============================================================================
# fig3_molecular.R — Figure 3
#
# Input:  data/omic_features.csv, data/gsea_results.csv
# Output: output/fig3_molecular.svg/.png
#
# Prerequisites: source("00_helpers.R"); install.packages("ggvenn")
#                install.packages("ggbreak")
# ============================================================================

library(ggvenn)
library(ggbreak)

omic_feat <- read_csv(file.path(data_dir, "omic_features.csv"),
                      show_col_types = FALSE)
gsea_res  <- read_csv(file.path(data_dir, "gsea_results.csv"),
                      show_col_types = FALSE)

# ── Panel theme (matches original plot.format) ─────────────────────────────
plot.format <- theme(
  plot.background  = element_blank(),
  panel.background = element_blank(),
  panel.border     = element_rect(color = "black", linewidth = 0.5, fill = NA),
  panel.grid       = element_blank(),
  strip.background = element_rect(color = NA, fill = NA, linewidth = 0.5),
  strip.text       = element_text(color = "black", size = 10),
  axis.line        = element_blank(),
  axis.ticks       = element_line(color = "black", linewidth = 0.5),
  axis.ticks.length = unit(-0.1, "cm"),
  axis.text        = element_text(color = "black", size = 9),
  axis.title       = element_text(color = "black", size = 10),
  plot.title       = element_text(color = "black", size = 10, face = "bold",
                                  hjust = 0.5),
  legend.background = element_blank(),
  legend.key       = element_blank(),
  legend.text      = element_text(color = "black", size = 9),
  legend.title     = element_text(color = "black", size = 10, face = "bold"),
  legend.position  = "right",
  panel.spacing    = unit(0.1, "lines")
)

# ── Venn fill colours (positional, matching domain order) ──────────────────
venn_domain_order  <- c("Vitality", "Psychological", "Cognitive",
                        "Locomotion", "Sensory")
venn_fill_colors   <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# ── Panel A: Venn diagrams per platform ──────────────────────────────────────
platforms <- unique(omic_feat$platform)

make_venn <- function(df, platform_name) {
  df_plat <- df %>% filter(platform == platform_name)
  coef_list <- split(df_plat$feature, df_plat$domain)
  coef_list <- coef_list[intersect(venn_domain_order, names(coef_list))]
  if (length(coef_list) < 2) return(NULL)
  idx <- match(names(coef_list), venn_domain_order)
  ggvenn(coef_list, fill_color = venn_fill_colors[idx]) +
    labs(title = platform_name) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

venn_plots <- Filter(Negate(is.null), lapply(platforms, function(p) make_venn(omic_feat, p)))
fig3a <- wrap_plots(venn_plots, ncol = 2)

# ── Panel B: Feature counts ─────────────────────────────────────────────────
domain_counts <- omic_feat %>%
  group_by(domain, platform) %>%
  summarise(count = n_distinct(feature), .groups = "drop")

fig3b <- ggplot(domain_counts, aes(x = domain, y = count, fill = platform)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9),
            vjust = -0.3, size = 4) +
  labs(x = NULL, y = "Count of Significant Features") +
  plot.format + scale_y_break(c(220, 700)) +
  theme(legend.position = "bottom")

# ── Panel C: GSEA cross-domain ───────────────────────────────────────────────
q_cut <- 0.25
domain_gsea_order <- c("Vitality", "Psychological", "Cognitive",
                       "Locomotion", "Sensory")

gsea_plot_df <- gsea_res %>%
  mutate(
    pathway_clean = clean_pathway_label(pathway),
    significant   = p_adjust < q_cut,
    neg_log_q     = pmin(-log10(p_adjust), 3),
    domain        = factor(domain, levels = domain_gsea_order)
  )

pathway_order <- gsea_plot_df %>%
  group_by(pathway_clean) %>%
  summarise(n_sig = sum(significant), n_dom = n(), .groups = "drop") %>%
  arrange(desc(n_sig), desc(n_dom)) %>%
  pull(pathway_clean)

gsea_plot_df <- gsea_plot_df %>%
  mutate(pathway_clean = factor(pathway_clean, levels = rev(pathway_order)))

fig3c <- ggplot(gsea_plot_df,
                aes(x = domain, y = pathway_clean,
                    color = NES, size = neg_log_q, shape = significant)) +
  geom_point(alpha = 0.88) +
  facet_wrap(~ database, scales = "free_y") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "q < 0.25", "FALSE" = "nominal"),
                     name   = "Significance") +
  scale_color_gradient2(low = "steelblue3", mid = "grey88", high = "firebrick",
                        midpoint = 0, name = "NES",
                        limits = c(-2, 2), oob = squish) +
  scale_size_continuous(name = expression(-log[10](q)), range = c(1.5, 7)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 8),
        strip.text  = element_text(face = "bold"),
        legend.position = "right")

# ── Combine ──────────────────────────────────────────────────────────────────
fig3 <- (fig3a | fig3b) / fig3c + plot_layout(heights = c(1, 1.5))

save_fig(fig3, "fig3_molecular", width = 14, height = 16)
