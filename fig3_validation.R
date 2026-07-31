# ============================================================================
# fig3_validation.R — External validation figure
#
# Inputs:
#   Supplementary Tables workbook, sheets ST8–ST10
#   <derived_dir>/validation_longitudinal.csv
#     participant_id, exam, age, sex, domain, score
#   <derived_dir>/validation_participant_outcomes.csv
#     participant_id, baseline_ic, followup_years, death,
#     baseline_condition_count, incident_any
#
# Output: output/fig3_validation.svg/.png
# ============================================================================

if (!exists("read_supp_range")) source("00_helpers.R")

library(mgcv)
library(survival)
library(survminer)

long_path <- file.path(derived_dir, "validation_longitudinal.csv")
outcome_path <- file.path(derived_dir, "validation_participant_outcomes.csv")
require_input(long_path, "Validation longitudinal input")
require_input(outcome_path, "Validation participant-outcome input")

validation_long <- read_csv(long_path, show_col_types = FALSE)
validation_outcomes <- read_csv(outcome_path, show_col_types = FALSE)

require_columns(
  validation_long,
  c("participant_id", "exam", "age", "sex", "domain", "score"),
  "validation_longitudinal.csv"
)
require_columns(
  validation_outcomes,
  c("participant_id", "baseline_ic", "followup_years", "death",
    "baseline_condition_count", "incident_any"),
  "validation_participant_outcomes.csv"
)

validation_domain_order <- c(
  "IC composite", "Vitality", "Cognitive",
  "Locomotion", "Psychological", "Sensory"
)

validation_long <- validation_long %>%
  mutate(
    domain = recode(domain,
                    "Overall IC" = "IC composite",
                    "Cognition" = "Cognitive"),
    domain = factor(domain, levels = validation_domain_order),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>%
  filter(!is.na(domain), is.finite(age), is.finite(score))

# ── A. Age relationships across the validation domains ──────────────────────
repeated_ids <- validation_long %>%
  count(participant_id, domain) %>%
  filter(n >= 2) %>%
  select(participant_id, domain)

validation_lines <- validation_long %>%
  semi_join(repeated_ids, by = c("participant_id", "domain"))

p_age <- ggplot(validation_long, aes(age, score, color = sex)) +
  geom_line(
    data = validation_lines,
    aes(group = participant_id),
    alpha = 0.06, linewidth = 0.18
  ) +
  geom_point(alpha = 0.10, size = 0.35, stroke = 0) +
  geom_smooth(
    method = "gam", formula = y ~ s(x, k = 4),
    se = TRUE, linewidth = 0.85, alpha = 0.10
  ) +
  facet_wrap(~ domain, scales = "free_y", ncol = 3) +
  scale_color_manual(values = sex_colors, name = NULL) +
  scale_fill_manual(values = sex_colors, guide = "none") +
  labs(x = "Age at exam", y = NULL) +
  theme_ic() +
  theme(legend.position = "bottom")

# ── B. IC trajectories by future health ─────────────────────────────────────
health_long <- validation_long %>%
  filter(domain == "IC composite") %>%
  left_join(validation_outcomes, by = "participant_id") %>%
  mutate(
    `Baseline conditions` = cut(
      baseline_condition_count,
      breaks = c(-Inf, 0, 1, 3, Inf),
      labels = c("0", "1", "2-3", "4+")
    ),
    `Incident disease` = factor(
      if_else(incident_any == 1, "Incident", "Disease-free"),
      levels = c("Disease-free", "Incident")
    ),
    Mortality = factor(
      if_else(death == 1, "Deceased", "Alive"),
      levels = c("Alive", "Deceased")
    )
  )

health_specs <- list(
  list(
    column = "Baseline conditions",
    colors = c("0" = "#95D5B2", "1" = "#52B788", "2-3" = "#2D6A4F", "4+" = "#1B4332")
  ),
  list(
    column = "Incident disease",
    colors = c("Disease-free" = "#4A7C59", "Incident" = "#C9973A")
  ),
  list(
    column = "Mortality",
    colors = c("Alive" = "#4A7C59", "Deceased" = "#7A3B2E")
  )
)

health_plots <- lapply(health_specs, function(spec) {
  ggplot(
    health_long %>% filter(!is.na(.data[[spec$column]])),
    aes(age, score, color = .data[[spec$column]], fill = .data[[spec$column]])
  ) +
    geom_smooth(
      method = "gam", formula = y ~ s(x, k = 4),
      linewidth = 0.85, alpha = 0.12
    ) +
    scale_color_manual(values = spec$colors, name = spec$column) +
    scale_fill_manual(values = spec$colors, name = spec$column) +
    labs(x = "Age at exam", y = NULL) +
    theme_ic() +
    theme(legend.position = "bottom")
})
p_health <- wrap_plots(health_plots, nrow = 1)

# ── C. Kaplan–Meier survival by baseline-IC tertile ─────────────────────────
survival_df <- validation_outcomes %>%
  filter(
    is.finite(baseline_ic), is.finite(followup_years),
    !is.na(death), followup_years >= 0
  ) %>%
  mutate(
    ic_tertile = ntile(baseline_ic, 3),
    ic_tertile = factor(
      ic_tertile,
      levels = 1:3,
      labels = c("T1 (lowest IC)", "T2", "T3 (highest IC)")
    )
  )

km_fit <- survfit(Surv(followup_years, death) ~ ic_tertile, data = survival_df)
km <- ggsurvplot(
  km_fit,
  data = survival_df,
  risk.table = TRUE,
  conf.int = FALSE,
  censor = FALSE,
  pval = TRUE,
  palette = c("#A1D99B", "#41AB5D", "#00441B"),
  legend.title = "Baseline IC",
  legend.labs = levels(survival_df$ic_tertile),
  xlab = "Follow-up time (years)",
  ylab = "Overall survival",
  xlim = c(0, floor(max(survival_df$followup_years))),
  break.time.by = 2,
  ggtheme = theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "top")
)
p_km <- wrap_elements(full = km$plot) /
  wrap_elements(full = km$table) +
  plot_layout(heights = c(3, 1))

# ── D. Incident-outcome associations from ST8 ────────────────────────────────
st8 <- read_supp_range("ST8", "A3:F13")
require_columns(st8, c("Outcome", "Model 2 HR (95% CI)"), "ST8")

outcome_order <- c(
  "All-cause mortality", "Cardiovascular death", "Cancer death",
  "Myocardial infarction", "Stroke", "Atrial fibrillation",
  "Hard CVD", "Heart failure", "Incident diabetes"
)
outcome_system <- c(
  "All-cause mortality" = "Mortality",
  "Cardiovascular death" = "Mortality",
  "Cancer death" = "Mortality",
  "Myocardial infarction" = "Cardiovascular",
  "Stroke" = "Cardiovascular",
  "Atrial fibrillation" = "Cardiovascular",
  "Hard CVD" = "Cardiovascular",
  "Heart failure" = "Cardiovascular",
  "Incident diabetes" = "Metabolic"
)

st8_ci <- parse_estimate_ci(st8$`Model 2 HR (95% CI)`)
st8_plot <- bind_cols(st8, st8_ci) %>%
  filter(Outcome %in% outcome_order) %>%
  mutate(
    Outcome = factor(Outcome, levels = rev(outcome_order)),
    system = factor(
      unname(outcome_system[as.character(Outcome)]),
      levels = c("Mortality", "Cardiovascular", "Metabolic")
    )
  )

p_outcomes <- ggplot(st8_plot, aes(estimate, Outcome)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0, orientation = "y", color = "#355872"
  ) +
  geom_point(size = 2, color = "#355872") +
  scale_x_log10(breaks = c(0.3, 0.5, 0.7, 1.0)) +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Hazard ratio per +1 SD IC (95% CI)", y = NULL) +
  theme_ic()

# ── E. Associations across race/ethnic groups from ST9 ──────────────────────
st9 <- read_supp_range("ST9", "A2:H9")
race_levels <- c("White", "Black", "Hispanic", "Chinese")
require_columns(st9, c("Outcome", race_levels), "ST9")
race_colors <- c(
  White = "#4E79A7", Black = "#59A14F",
  Hispanic = "#B07AA1", Chinese = "#F28E2B"
)

st9_long <- st9 %>%
  select(Outcome, all_of(race_levels)) %>%
  pivot_longer(-Outcome, names_to = "group", values_to = "result")
st9_ci <- parse_estimate_ci(st9_long$result)
st9_plot <- bind_cols(st9_long, st9_ci) %>%
  mutate(
    Outcome = factor(Outcome, levels = rev(unique(st9$Outcome))),
    group = factor(group, levels = race_levels)
  )

p_group <- ggplot(st9_plot, aes(estimate, Outcome, color = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0, orientation = "y",
    position = position_dodge(width = 0.65)
  ) +
  geom_point(size = 1.7, position = position_dodge(width = 0.65)) +
  scale_x_log10(breaks = c(0.3, 0.5, 0.7, 1.0, 1.5, 2.0)) +
  scale_color_manual(values = race_colors, name = NULL) +
  labs(x = "Hazard ratio per +1 SD IC (95% CI)", y = NULL) +
  theme_ic() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))

# ── F. Comparison with frailty measures from ST10 ────────────────────────────
st10 <- read_supp_range("ST10", "A3:D10")
measure_levels <- c("IC", "Frailty index", "Fried")
require_columns(st10, c("Outcome", measure_levels), "ST10")

st10_long <- st10 %>%
  pivot_longer(all_of(measure_levels), names_to = "measure", values_to = "result")
st10_ci <- parse_estimate_ci(st10_long$result)
st10_plot <- bind_cols(st10_long, st10_ci) %>%
  mutate(
    Outcome = factor(Outcome, levels = rev(unique(st10$Outcome))),
    measure = factor(measure, levels = measure_levels)
  )

p_frailty <- ggplot(st10_plot, aes(estimate, Outcome, color = measure)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0, orientation = "y",
    position = position_dodge(width = 0.62)
  ) +
  geom_point(size = 1.8, position = position_dodge(width = 0.62)) +
  scale_x_log10(breaks = c(1.0, 1.2, 1.4, 1.6, 1.8)) +
  scale_color_manual(
    values = c(IC = "#333333", `Frailty index` = "#C0392B", Fried = "#3B6DA0"),
    name = NULL
  ) +
  labs(x = "Hazard ratio per 1 SD worse (95% CI)", y = NULL) +
  theme_ic() +
  theme(legend.position = "bottom")

# ── Assemble ─────────────────────────────────────────────────────────────────
fig3_validation <- p_age /
  (p_health | p_km) /
  (p_outcomes | p_group | p_frailty) +
  plot_layout(heights = c(1.05, 1.1, 1.15))

save_fig(fig3_validation, "fig3_validation", width = 15, height = 17)
