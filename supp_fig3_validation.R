# ============================================================================
# supp_fig3_validation.R — Supporting external-validation figure
#
# Inputs:
#   Supplementary Tables workbook, sheets ST10 and ST18
#   <derived_dir>/validation_domain_correlations.csv
#     domain_x, domain_y, correlation
#   <derived_dir>/validation_benchmark_associations.csv
#     domain, benchmark, family, estimate, significant
#   <derived_dir>/validation_attrition.csv
#     exam, age_stratum, retained_percent
#
# Output: output/supp_fig3_validation.svg/.png
# ============================================================================

if (!exists("read_supp_range")) source("00_helpers.R")

corr_path <- file.path(derived_dir, "validation_domain_correlations.csv")
benchmark_path <- file.path(derived_dir, "validation_benchmark_associations.csv")
attrition_path <- file.path(derived_dir, "validation_attrition.csv")

require_input(corr_path, "Validation domain-correlation input")
require_input(benchmark_path, "Validation benchmark-association input")
require_input(attrition_path, "Validation attrition input")

corr_df <- read_csv(corr_path, show_col_types = FALSE)
benchmark_df <- read_csv(benchmark_path, show_col_types = FALSE)
attrition_df <- read_csv(attrition_path, show_col_types = FALSE)

require_columns(
  corr_df, c("domain_x", "domain_y", "correlation"),
  "validation_domain_correlations.csv"
)
require_columns(
  benchmark_df, c("domain", "benchmark", "family", "estimate", "significant"),
  "validation_benchmark_associations.csv"
)
require_columns(
  attrition_df, c("exam", "age_stratum", "retained_percent"),
  "validation_attrition.csv"
)

validation_domains <- c(
  "IC composite", "Vitality", "Locomotion",
  "Cognitive", "Psychological", "Sensory"
)

normalize_domain <- function(x) {
  recode(
    x,
    "IC" = "IC composite",
    "Overall IC" = "IC composite",
    "Cognition" = "Cognitive"
  )
}

# ── A. Domain correlation matrix ─────────────────────────────────────────────
corr_plot_df <- corr_df %>%
  transmute(
    domain_x = factor(normalize_domain(domain_x), levels = validation_domains),
    domain_y = factor(normalize_domain(domain_y), levels = rev(validation_domains)),
    correlation = as.numeric(correlation)
  ) %>%
  filter(!is.na(domain_x), !is.na(domain_y), is.finite(correlation))

p_corr <- ggplot(corr_plot_df, aes(domain_x, domain_y, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.2) +
  scale_fill_gradient2(
    low = "#3B6DA0", mid = "white", high = "#C0392B",
    midpoint = 0, limits = c(-1, 1), name = "Correlation"
  ) +
  labs(x = NULL, y = NULL) +
  theme_heatmap() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# ── B. Current IC level and prior rate of change from ST18 ───────────────────
st18 <- read_supp_range("ST18", "A17:H20")
require_columns(
  st18,
  c("Outcome", "IC level HR (95% CI)", "Prior-slope HR (95% CI)"),
  "ST18"
)
level_ci <- parse_estimate_ci(st18$`IC level HR (95% CI)`)
slope_ci <- parse_estimate_ci(st18$`Prior-slope HR (95% CI)`)

level_slope <- bind_rows(
  bind_cols(
    st18["Outcome"], level_ci,
    tibble::tibble(measure = "Current IC level")
  ),
  bind_cols(
    st18["Outcome"], slope_ci,
    tibble::tibble(measure = "Prior rate of change")
  )
) %>%
  mutate(
    Outcome = factor(Outcome, levels = rev(st18$Outcome)),
    measure = factor(measure, levels = c("Current IC level", "Prior rate of change"))
  )

p_level_slope <- ggplot(level_slope, aes(estimate, Outcome, color = measure)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0, orientation = "y",
    position = position_dodge(width = 0.55)
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.55)) +
  scale_x_log10(breaks = c(0.5, 0.7, 1.0, 1.2, 1.4)) +
  scale_color_manual(
    values = c("Current IC level" = "#355872", "Prior rate of change" = "#C47F34"),
    name = NULL
  ) +
  labs(x = "Hazard ratio (95% CI)", y = NULL) +
  theme_ic() +
  theme(legend.position = "bottom")

# ── C. IC domains and aging benchmarks ───────────────────────────────────────
benchmark_plot_df <- benchmark_df %>%
  transmute(
    domain = factor(normalize_domain(domain), levels = rev(validation_domains)),
    benchmark,
    family,
    estimate = as.numeric(estimate),
    significant = as.logical(significant)
  ) %>%
  filter(!is.na(domain), is.finite(estimate))

p_benchmark <- ggplot(
  benchmark_plot_df,
  aes(benchmark, domain, fill = estimate)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f%s", estimate, if_else(significant, "*", ""))),
    size = 3.1
  ) +
  facet_grid(. ~ family, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#3B6DA0", mid = "white", high = "#C0392B",
    midpoint = 0, name = "Association"
  ) +
  labs(x = NULL, y = NULL) +
  theme_heatmap()

# ── D. Added C-index over DNA-methylation clocks from ST10 ───────────────────
st10_clocks_raw <- read_supp_range("ST10", "A13:D18")
clock_col <- grep("^DNAm clock", names(st10_clocks_raw))
delta_col <- grep("^IC adds over clock", names(st10_clocks_raw))
if (!length(clock_col) || !length(delta_col)) {
  stop("ST10 Panel B is missing the DNAm clock or IC-adds-over-clock column.",
       call. = FALSE)
}

st10_clocks <- tibble::tibble(
  clock = st10_clocks_raw[[clock_col[1]]],
  delta_c = parse_delta_c(st10_clocks_raw[[delta_col[1]]])
) %>%
  filter(!is.na(clock), is.finite(delta_c)) %>%
  arrange(delta_c) %>%
  mutate(clock = factor(clock, levels = clock))

p_clocks <- ggplot(st10_clocks, aes(delta_c, clock)) +
  geom_col(width = 0.62, fill = "#355872") +
  geom_text(
    aes(label = sprintf("+%.3f", delta_c)),
    hjust = -0.15, size = 3.1
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(x = "IC added C-index over clock", y = NULL) +
  theme_ic()

# ── E. Retention across validation visits ────────────────────────────────────
attrition_plot_df <- attrition_df %>%
  mutate(
    exam = factor(exam, levels = unique(exam)),
    age_stratum = factor(age_stratum, levels = unique(age_stratum))
  )

attrition_colors <- setNames(
  c("#4A7C59", "#E0A03A", "#C0392B")[seq_along(levels(attrition_plot_df$age_stratum))],
  levels(attrition_plot_df$age_stratum)
)

p_attrition <- ggplot(
  attrition_plot_df,
  aes(exam, retained_percent, color = age_stratum, group = age_stratum)
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.8) +
  scale_color_manual(values = attrition_colors, name = "Baseline age") +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 100)) +
  labs(x = NULL, y = "Cohort retained") +
  theme_ic() +
  theme(legend.position = "bottom")

supp_fig3 <- (p_corr | p_level_slope) /
  p_benchmark /
  (p_clocks | p_attrition) +
  plot_layout(heights = c(1, 1.15, 0.9))

save_fig(supp_fig3, "supp_fig3_validation", width = 12, height = 13)
