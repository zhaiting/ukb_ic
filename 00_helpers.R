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
library(scales)
library(patchwork)
library(readr)

# ── Paths ────────────────────────────────────────────────────────────────────
data_dir <- Sys.getenv("UKB_IC_DATA_DIR", unset = "data")
fig_dir  <- file.path("output")
supp_xlsx <- Sys.getenv(
  "UKB_IC_SUPP_XLSX",
  unset = file.path(data_dir, "Supp_Tables_R1.xlsx")
)
derived_dir <- Sys.getenv("UKB_IC_DERIVED_DIR", unset = data_dir)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── IC score colours ─────────────────────────────────────────────────────────
score_colors <- c(
  "IC composite"  = "#333333",
  "Vitality"      = "#E8B33D",
  "Psychological" = "#2E9E8F",
  "Cognitive"     = "#8E6FB0",
  "Locomotion"    = "#6FA86B",
  "Sensory"       = "#B5793F"
)

# ── Domain colour palette ────────────────────────────────────────────────────
domain_colors <- c(
  Vitality      = "#E8B33D",
  Psychological = "#2E9E8F",
  Cognitive     = "#8E6FB0",
  Locomotion    = "#6FA86B",
  Sensory       = "#B5793F"
)

# ── Sex colours ──────────────────────────────────────────────────────────────
sex_colors <- c(Female = "#E36A6A", Male = "#355872")

# ── Input checks ─────────────────────────────────────────────────────────────
require_input <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop(
      label, " was not found at: ", normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }
  invisible(path)
}

require_columns <- function(df, columns, label = deparse(substitute(df))) {
  missing <- setdiff(columns, names(df))
  if (length(missing)) {
    stop(
      label, " is missing required columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(df)
}

read_supp_table <- function(sheet, header_row = 2) {
  require_input(supp_xlsx, "Supplementary Tables workbook")
  readxl::read_excel(
    supp_xlsx,
    sheet = sheet,
    skip = header_row - 1,
    .name_repair = "minimal"
  )
}

read_supp_range <- function(sheet, range) {
  require_input(supp_xlsx, "Supplementary Tables workbook")
  readxl::read_excel(
    supp_xlsx,
    sheet = sheet,
    range = range,
    .name_repair = "minimal"
  )
}

parse_estimate_ci <- function(x) {
  x <- gsub("\u2212", "-", as.character(x), fixed = TRUE)
  x <- gsub("\u2013|\u2014", "-", x)
  parts <- stringr::str_match(
    x,
    "^\\s*([+-]?[0-9.]+)\\s*\\(([+-]?[0-9.]+)\\s*-\\s*([+-]?[0-9.]+)"
  )
  tibble::tibble(
    estimate = suppressWarnings(as.numeric(parts[, 2])),
    lower = suppressWarnings(as.numeric(parts[, 3])),
    upper = suppressWarnings(as.numeric(parts[, 4]))
  )
}

parse_delta_c <- function(x) {
  x <- gsub("\u2212", "-", as.character(x), fixed = TRUE)
  suppressWarnings(as.numeric(stringr::str_extract(x, "[+-]?[0-9.]+")))
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
         width = width, height = height, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300, bg = "white")
  message("Saved: ", name, " (.svg + .png)")
}

# ── GSEA helper ──────────────────────────────────────────────────────────────
clean_pathway_label <- function(x) {
  gsub("_", " ", gsub("^HALLMARK_|^REACTOME_", "", x))
}

message(
  "Helpers loaded. Supplement workbook: ", supp_xlsx,
  "; derived figure inputs: ", derived_dir
)
