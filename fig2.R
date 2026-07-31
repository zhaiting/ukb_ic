# ============================================================================
# fig2.R — Figure 2
#
# Inputs: Supplementary Tables workbook, sheets ST3 and ST6
# Output: output/fig2a.svg/.png
#         output/fig2b.svg/.png
#         output/fig2c.svg/.png
#
# Each panel uses the specified rows and estimates from the workbook.
# ============================================================================

if (!exists("read_supp_table")) source("00_helpers.R")

endpoint_system <- c(
  "Parkinson's" = "Neurological",
  "Dementia" = "Neurological",
  "OCD / stress disorders" = "Psychological",
  "Anxiety" = "Psychological",
  "Depression" = "Psychological",
  "Psychotic disorders" = "Psychological",
  "Bipolar" = "Psychological",
  "Osteoarthritis" = "Locomotion",
  "Dorsopathy" = "Locomotion",
  "Osteoporosis" = "Locomotion",
  "Inflammatory arthritis" = "Locomotion",
  "Glaucoma" = "Sensory",
  "Cataract" = "Sensory",
  "Hearing loss" = "Sensory",
  "Hypertension" = "Cardiovascular",
  "Atrial fibrillation" = "Cardiovascular",
  "Ischaemic heart disease" = "Cardiovascular",
  "Stroke" = "Cardiovascular",
  "Heart failure" = "Cardiovascular",
  "COPD" = "Respiratory",
  "Dyslipidaemia" = "Metabolic",
  "Type 2 diabetes" = "Metabolic",
  "Obesity" = "Metabolic"
)

system_primary_domain <- c(
  Neurological = "Cognitive",
  Psychological = "Psychological",
  Locomotion = "Locomotion",
  Sensory = "Sensory",
  Cardiovascular = "Vitality",
  Respiratory = "Vitality",
  Metabolic = "Vitality"
)

system_order <- names(system_primary_domain)
endpoint_order <- names(endpoint_system)
lane_order <- c("Vitality", "Psychological", "Cognitive", "Locomotion", "Sensory")
lane_offset <- setNames(seq(-0.28, 0.28, length.out = length(lane_order)), lane_order)

st3 <- read_supp_table("ST3", header_row = 2)
require_columns(
  st3,
  c("outcome_class", "endpoint", "exposure", "model", "stratum",
    "hr", "hr_lo", "hr_hi", "p_fdr"),
  "ST3"
)

plot_df <- st3 %>%
  filter(
    outcome_class == "Disease",
    model == "Model 2",
    stratum == "Overall",
    endpoint %in% endpoint_order,
    exposure %in% c("Overall IC", "Vitality", "Psychological",
                    "Cognition", "Locomotion", "Sensory")
  ) %>%
  transmute(
    endpoint,
    system = unname(endpoint_system[endpoint]),
    exposure = recode(exposure,
                      "Overall IC" = "IC composite",
                      "Cognition" = "Cognitive"),
    hr = as.numeric(hr),
    lower = as.numeric(hr_lo),
    upper = as.numeric(hr_hi),
    significant = as.numeric(p_fdr) < 0.05
  ) %>%
  mutate(
    system = factor(system, levels = system_order),
    endpoint = factor(endpoint, levels = rev(endpoint_order))
  )

if (!nrow(plot_df)) stop("No eligible Figure 2a rows were found in ST3.", call. = FALSE)

endpoint_positions <- tibble(
  endpoint = factor(rev(endpoint_order), levels = rev(endpoint_order)),
  y = seq_along(endpoint_order)
) %>%
  mutate(
    endpoint_character = as.character(endpoint),
    system = factor(unname(endpoint_system[endpoint_character]), levels = system_order)
  )

composite_df <- plot_df %>%
  filter(exposure == "IC composite") %>%
  left_join(select(endpoint_positions, endpoint, y), by = "endpoint")

domain_df <- plot_df %>%
  filter(exposure != "IC composite") %>%
  left_join(select(endpoint_positions, endpoint, y), by = "endpoint") %>%
  mutate(
    exposure = factor(exposure, levels = lane_order),
    y_lane = y + lane_offset[as.character(exposure)]
  )

matched_lanes <- endpoint_positions %>%
  mutate(
    primary_domain = unname(system_primary_domain[as.character(system)]),
    y_lane = y + lane_offset[primary_domain]
  )

shared_y <- scale_y_continuous(
  breaks = endpoint_positions$y,
  labels = endpoint_positions$endpoint_character,
  expand = expansion(add = 0.55)
)

p_composite <- ggplot(composite_df, aes(hr, y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper, alpha = significant),
    width = 0, orientation = "y",
    color = score_colors[["IC composite"]], linewidth = 0.55
  ) +
  geom_point(
    aes(alpha = significant),
    color = score_colors[["IC composite"]], size = 2.1
  ) +
  scale_x_log10(breaks = c(0.4, 0.6, 0.8, 1.0, 1.2)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
  shared_y +
  facet_grid(system ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "IC composite HR (95% CI)", y = NULL) +
  theme_ic() +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

p_domains <- ggplot(domain_df, aes(hr, y_lane, color = exposure)) +
  geom_rect(
    data = matched_lanes,
    aes(ymin = y_lane - 0.055, ymax = y_lane + 0.055),
    xmin = -Inf, xmax = Inf, inherit.aes = FALSE,
    fill = "grey91", color = NA
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper, alpha = significant),
    width = 0, orientation = "y", linewidth = 0.45
  ) +
  geom_point(aes(alpha = significant), size = 1.65) +
  scale_x_log10(breaks = c(0.6, 0.8, 1.0, 1.2)) +
  scale_color_manual(values = domain_colors[lane_order], name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = "none") +
  shared_y +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Domain-specific HR (95% CI)", y = NULL) +
  theme_ic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 1))

fig2a <- p_composite + p_domains + plot_layout(widths = c(1.15, 2.6))
save_fig(fig2a, "fig2a", width = 10.5, height = 11.5)

# ── Figure 2b: cancer incidence and cause-specific mortality ────────────────
model_keep <- "Model 2"

# The default display uses endpoints with an Overall stratum. Set this to TRUE
# to add the sex-specific endpoints using their corresponding strata.
include_sex_specific <- FALSE

sex_specific <- c(
  "Breast cancer (female)"  = "Female",
  "Ovary cancer (female)"   = "Female",
  "Uterine cancer (female)" = "Female",
  "Prostate cancer (male)"  = "Male"
)

cancer_keep <- c(
  "Any malignant cancer", "Lung cancer", "Colorectal cancer",
  "Breast cancer (female)", "Prostate cancer (male)", "Pancreatic cancer",
  "Liver cancer", "Kidney cancer", "Bladder cancer", "Oesophageal cancer",
  "Stomach cancer", "Melanoma", "Non-Hodgkin lymphoma", "Leukaemia",
  "Ovary cancer (female)", "Uterine cancer (female)", "Brain cancer",
  "Head & neck cancer", "Myeloma", "Thyroid cancer"
)

if (!include_sex_specific) {
  cancer_keep <- setdiff(cancer_keep, names(sex_specific))
}

mort_order <- c(
  "Suicide", "Infection mortality", "Kidney mortality", "Liver mortality",
  "Diabetes mortality", "Respiratory mortality", "Dementia mortality",
  "Cancer mortality", "CVD mortality", "All-cause mortality"
)

require_columns(
  st3,
  c("outcome_class", "endpoint", "exposure", "model", "stratum",
    "hr", "hr_lo", "hr_hi", "p_fdr"),
  "ST3"
)

composite_rows <- function(df, class_name, endpoints) {
  wanted <- tibble(
    endpoint = endpoints,
    stratum = unname(ifelse(
      endpoints %in% names(sex_specific),
      sex_specific[endpoints],
      "Overall"
    ))
  )

  out <- df %>%
    filter(
      outcome_class == class_name,
      model == model_keep,
      exposure == "Overall IC"
    ) %>%
    inner_join(wanted, by = c("endpoint", "stratum")) %>%
    transmute(
      endpoint,
      stratum,
      hr = as.numeric(hr),
      lower = as.numeric(hr_lo),
      upper = as.numeric(hr_hi),
      significant = as.numeric(p_fdr) < 0.05,
      group = class_name
    ) %>%
    filter(is.finite(hr), is.finite(lower), is.finite(upper))

  missing <- setdiff(endpoints, out$endpoint)
  if (length(missing)) {
    stop(
      class_name, " endpoints absent from ST3 under ", model_keep, ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  out
}

cancer_df <- composite_rows(st3, "Cancer", cancer_keep) %>%
  arrange(hr) %>%
  mutate(endpoint = factor(endpoint, levels = endpoint))

mort_df <- composite_rows(st3, "Mortality", mort_order) %>%
  mutate(endpoint = factor(endpoint, levels = mort_order))

message("Figure 2b rows — cancer: ", nrow(cancer_df),
        "; mortality: ", nrow(mort_df))

forest_ic <- function(df, breaks) {
  ggplot(df, aes(hr, endpoint)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey55", linewidth = 0.4) +
    geom_errorbar(
      aes(xmin = lower, xmax = upper, alpha = significant),
      width = 0, orientation = "y",
      color = score_colors[["IC composite"]], linewidth = 0.5
    ) +
    geom_point(
      aes(alpha = significant),
      color = score_colors[["IC composite"]], size = 2.4
    ) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.28), guide = "none") +
    scale_x_log10(
      breaks = breaks,
      limits = c(min(df$lower) * 0.95, max(max(df$upper), 1) * 1.04),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    facet_grid(group ~ ., switch = "y") +
    labs(title = "IC", x = "Adjusted HR per 1-SD (95% CI)", y = NULL) +
    theme_ic() +
    theme(
      panel.grid.minor = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}

fig2b <- forest_ic(cancer_df, breaks = c(0.5, 0.7, 1.0)) +
  forest_ic(mort_df, breaks = c(0.3, 0.5, 1.0)) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = paste("Prospective associations with cancer incidence",
                  "and cause-specific mortality"),
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    )
  )

save_fig(fig2b, "fig2b", width = 8.6, height = 5.8)

# ── Figure 2c: lifestyle associations ───────────────────────────────────────
factor_map <- tibble::tribble(
  ~st6_name,                                 ~label,                  ~category,
  "Vigorous physical activity (days/week)",  "Vigorous PA (days/wk)", "PA\n(self-report)",
  "Moderate physical activity (days/week)",  "Moderate PA (days/wk)", "PA\n(self-report)",
  "Sleep duration (accel.)",                 "Sleep (accel.)",        "PA\n(accelerometer)",
  "Sedentary time (accel.)",                 "Sedentary",             "PA\n(accelerometer)",
  "Overall acceleration (accel.)",           "Overall acceleration",  "PA\n(accelerometer)",
  "MVPA (accel.)",                           "MVPA",                  "PA\n(accelerometer)",
  "Light activity (accel.)",                 "Light activity",        "PA\n(accelerometer)",
  "Sleep duration (hours)",                  "Sleep duration",        "Sleep",
  "Insomnia",                                "Insomnia",              "Sleep",
  "Tea intake",                              "Tea",                   "Diet &\nalcohol",
  "Raw vegetable intake",                    "Raw vegetables",        "Diet &\nalcohol",
  "Processed meat intake",                   "Processed meat",        "Diet &\nalcohol",
  "Poultry intake",                          "Poultry",               "Diet &\nalcohol",
  "Pork intake",                             "Pork",                  "Diet &\nalcohol",
  "Oily fish intake",                        "Oily fish",             "Diet &\nalcohol",
  "Non-oily fish intake",                    "Non-oily fish",         "Diet &\nalcohol",
  "Fresh fruit intake",                      "Fresh fruit",           "Diet &\nalcohol",
  "Cooked vegetable intake",                 "Cooked vegetables",     "Diet &\nalcohol",
  "Bread intake",                            "Bread",                 "Diet &\nalcohol",
  "Beef intake",                             "Beef",                  "Diet &\nalcohol",
  "Alcohol intake frequency",                "Alcohol frequency",     "Diet &\nalcohol",
  "Smoking status",                          "Smoking status",        "Smoking",
  "Townsend deprivation index",              "Townsend deprivation",  "Socio-\neconomic\n& social",
  "Social engagement composite",             "Social isolation",      "Socio-\neconomic\n& social",
  "Age left full-time education",            "Age at education",      "Socio-\neconomic\n& social",
  "Number of regular medications",           "No. medications",       "Medication",
  "Insulin medication use",                  "Insulin",               "Medication",
  "Cholesterol-lowering medication use",     "Cholesterol meds",      "Medication",
  "Blood pressure medication use",           "BP medication",         "Medication"
)

category_order <- unique(factor_map$category)

domain_labels <- c(
  "Overall IC"    = "IC",
  "Vitality"      = "Vitality",
  "Cognition"     = "Cognition",
  "Locomotion"    = "Locomotion",
  "Psychological" = "Psychological",
  "Sensory"       = "Sensory"
)
domain_order <- unname(domain_labels)

st6 <- read_supp_table("ST6", header_row = 2)
require_columns(
  st6,
  c("Model", "Source", "Lifestyle factor", "IC domain", "beta", "se", "p_fdr"),
  "ST6"
)

lifestyle_df <- st6 %>%
  filter(
    Model == "Model 2",
    `IC domain` %in% names(domain_labels),
    `Lifestyle factor` %in% factor_map$st6_name
  ) %>%
  transmute(
    st6_name = `Lifestyle factor`,
    domain = factor(unname(domain_labels[`IC domain`]), levels = domain_order),
    beta = as.numeric(beta),
    significant = !is.na(as.numeric(p_fdr)) & as.numeric(p_fdr) < 0.05
  ) %>%
  inner_join(factor_map, by = "st6_name") %>%
  mutate(
    label = factor(label, levels = rev(factor_map$label)),
    category = factor(category, levels = category_order)
  ) %>%
  filter(is.finite(beta))

expected_cells <- nrow(factor_map) * length(domain_order)
if (nrow(lifestyle_df) != expected_cells) {
  stop(
    "Expected ", expected_cells, " factor x score cells in ST6, found ",
    nrow(lifestyle_df), ".",
    call. = FALSE
  )
}

message(
  "Figure 2c cells: ", nrow(lifestyle_df),
  " (FDR-significant: ", sum(lifestyle_df$significant),
  "; not significant: ", sum(!lifestyle_df$significant), ")"
)

beta_limit <- 0.20

fig2c <- ggplot(lifestyle_df, aes(domain, label)) +
  geom_tile(fill = NA, colour = "grey90", linewidth = 0.25) +
  geom_point(
    data = filter(lifestyle_df, significant),
    aes(size = abs(beta), fill = beta),
    shape = 21, colour = "grey30", stroke = 0.3
  ) +
  geom_point(
    data = filter(lifestyle_df, !significant),
    aes(size = abs(beta)),
    shape = 21, fill = "white", colour = "grey60", stroke = 0.4
  ) +
  scale_fill_gradient2(
    low = "#8C510A", mid = "white", high = "#01665E", midpoint = 0,
    limits = c(-beta_limit, beta_limit), oob = squish,
    breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
    name = expression(beta ~ "(M2)")
  ) +
  scale_size_area(
    max_size = 5.5, limits = c(0, 0.28), breaks = c(0.05, 0.10, 0.20),
    name = expression(abs(beta))
  ) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Cross-sectional association with lifestyle factors",
    x = NULL, y = NULL,
    caption = paste(
      "Fully adjusted Model 2 standardised beta per 1 SD.",
      "Hollow = not FDR-significant. Model 1 estimates are in ST6."
    )
  ) +
  guides(
    fill = guide_colourbar(
      order = 1, barheight = 0.5, barwidth = 6, title.vjust = 1
    ),
    size = guide_legend(
      order = 2, override.aes = list(fill = "grey50")
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 8.5, colour = "grey15"),
    strip.text.y     = element_text(
      angle = 0, size = 8, face = "bold", lineheight = 0.9
    ),
    strip.background = element_rect(fill = "grey94", colour = "grey80"),
    panel.spacing.y  = unit(3, "pt"),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(size = 9, vjust = 1),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption     = element_text(size = 6.8, colour = "grey40", hjust = 0),
    plot.caption.position = "plot"
  )

save_fig(fig2c, "fig2c", width = 5.0, height = 10.0)
