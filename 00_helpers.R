# ============================================================================
# 00_helpers.R
# Shared constants, colour palettes, labels, and utility functions.
#
# Source this file before running any figure script:
#   source("00_helpers.R")
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(scales)
library(patchwork)
library(readr)

# ── Paths ────────────────────────────────────────────────────────────────────
data_dir <- file.path("data")
fig_dir  <- file.path("output")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── IC score labels & colours ────────────────────────────────────────────────
score_labels <- c(
  ic_score_i0            = "IC composite",
  vitality_score_i0      = "Vitality",
  psychological_score_i0 = "Psychological",
  cog_comp_i0            = "Cognitive",
  locomotion_score_i0    = "Locomotion",
  sensory_score_i0       = "Sensory"
)

score_colors <- c(
  "IC composite"  = "#333333",
  "Vitality"      = "#E41A1C",
  "Psychological" = "#377EB8",
  "Cognitive"     = "#4DAF4A",
  "Locomotion"    = "#984EA3",
  "Sensory"       = "#C4A000"
)

score_order <- c(
  "IC composite", "Vitality", "Psychological",
  "Cognitive", "Locomotion", "Sensory"
)

# ── Domain colour palette ────────────────────────────────────────────────────
domain_colors <- c(
  Vitality      = "#E41A1C",
  Psychological = "#377EB8",
  Cognitive     = "#4DAF4A",
  Locomotion    = "#984EA3",
  Sensory       = "#FF7F00"
)

# ── Sex colours ──────────────────────────────────────────────────────────────
sex_colors <- c(Female = "#C0392B", Male = "#2E86C1")

# ── Outcome metadata (loaded from data/) ─────────────────────────────────────
# These are read at runtime from the supplementary CSVs.
# Columns expected: outcome, label, system, primary_domain
load_disease_meta <- function(path = file.path(data_dir, "disease_meta.csv")) {
  if (file.exists(path)) {
    read_csv(path, show_col_types = FALSE) %>%
      mutate(system = factor(system, levels = unique(system)))
  } else {
    message("disease_meta.csv not found — some scripts may need it in data/")
    NULL
  }
}

load_mort_meta <- function(path = file.path(data_dir, "mort_meta.csv")) {
  if (file.exists(path)) {
    read_csv(path, show_col_types = FALSE)
  } else {
    message("mort_meta.csv not found — some scripts may need it in data/")
    NULL
  }
}

# ── Shared base theme ────────────────────────────────────────────────────────
theme_ic <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
      strip.background   = element_rect(fill = "grey92", color = NA),
      panel.grid.major.y = element_blank(),
      plot.title         = element_text(face = "bold")
    )
}

# ── Heatmap theme ────────────────────────────────────────────────────────────
theme_heatmap <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text.y     = element_text(face = "bold", size = 9, angle = 0, hjust = 0),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y      = element_text(size = 9, colour = "grey20"),
      axis.ticks       = element_blank(),
      legend.position  = "bottom",
      plot.title       = element_text(size = 12, face = "bold"),
      plot.subtitle    = element_text(size = 9, colour = "grey40")
    )
}

# ── Save helper ──────────────────────────────────────────────────────────────
save_fig <- function(p, name, width = 8, height = 10) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p,
         width = width, height = height)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300)
  message("Saved: ", name, " (.svg + .png)")
}

# ── Tidy helpers ─────────────────────────────────────────────────────────────

tidy_mort <- function(df, sex_set, mort_meta) {
  domain_map <- setNames(mort_meta$primary_domain, mort_meta$label)
  group_map  <- setNames(mort_meta$group, mort_meta$label)
  group_order <- unique(mort_meta$group)

  df %>%
    mutate(p.val = .data[["p"]]) %>%
    mutate(
      score_label    = recode(exposure, !!!score_labels),
      domain_group   = recode(label, !!!group_map),
      domain_group   = factor(domain_group, levels = group_order),
      primary_domain = recode(label, !!!domain_map),
      sex_set        = sex_set
    ) %>%
    filter(!is.na(score_label))
}

tidy_dis <- function(df, sex_set, disease_meta) {
  df %>%
    mutate(
      lo      = .data[["conf.low"]],
      hi      = .data[["conf.high"]],
      p.val   = .data[["p.value"]],
      score_label = recode(score, !!!score_labels),
      sex_set = sex_set
    ) %>%
    left_join(disease_meta, by = "outcome") %>%
    filter(!is.na(score_label), !is.na(label))
}

# ── GSEA helper ──────────────────────────────────────────────────────────────
clean_pathway_label <- function(x) {
  gsub("_", " ", gsub("^HALLMARK_|^REACTOME_", "", x))
}

message("Helpers loaded. Place data CSVs in '", data_dir, "/' and run figure scripts.")
