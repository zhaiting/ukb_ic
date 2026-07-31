# ============================================================================
# fig1.R — Figure 1
#
# Input:  data/ic_trajectories.csv   (long format: participant_id, domain,
#             score, age, sex, visit)
#         data/ic_survival.csv       (ic_score, time_years, event)
#         data/ic_domain_corr.csv    (6×6 correlation matrix of IC domain scores)
#
# Output: output/fig1b.svg/.png
#         output/supp_fig1f.svg/.png
#         output/fig1c.pdf
#         output/fig1e.svg/.png
# ============================================================================

if (!exists("save_fig")) source("00_helpers.R")

library(mgcv)
library(lme4)
library(survival)
library(survminer)
library(corrplot)

# ── Data ───────────────────────────────────────────────────────────────────
trajectory_path <- file.path(data_dir, "ic_trajectories.csv")
survival_path <- file.path(data_dir, "ic_survival.csv")
correlation_path <- file.path(data_dir, "ic_domain_corr.csv")
require_input(trajectory_path, "IC trajectory input")
require_input(survival_path, "IC survival input")
require_input(correlation_path, "IC domain-correlation input")

traj_df <- read_csv(trajectory_path, show_col_types = FALSE)
require_columns(
  traj_df,
  c("participant_id", "domain", "score", "age", "sex", "visit"),
  "ic_trajectories.csv"
)

# ── Plot theme ─────────────────────────────────────────────────────────────
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
  plot.subtitle    = element_text(color = "black", size = 10, face = "italic",
                                  hjust = 0.5),
  legend.background = element_blank(),
  legend.key       = element_blank(),
  legend.text      = element_text(color = "black", size = 9),
  legend.title     = element_text(color = "black", size = 10, face = "bold"),
  legend.position  = "right",
  panel.spacing    = unit(0.1, "lines")
)

sex_traj_colors <- c("Female" = "#E36A6A", "Male" = "#355872")

domain_order <- c("Cognition", "Locomotion", "Psychological",
                  "Sensory", "Vitality", "Overall IC")

traj_df <- traj_df %>%
  mutate(
    sex    = factor(sex, levels = c("Female", "Male")),
    domain = factor(domain, levels = domain_order)
  )

# ── Supplementary Figure 1: cross-sectional pattern ──────────────────────────
set.seed(96)
n_sample <- min(20000, n_distinct(traj_df$participant_id))
ids_sample <- sample(unique(traj_df$participant_id), n_sample)

plot_df_a <- traj_df %>%
  filter(participant_id %in% ids_sample, visit == 0)

supp_fig1f <- ggplot(plot_df_a, aes(x = age, y = score, color = sex)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              se = FALSE, linewidth = 0.7) +
  facet_wrap(~ domain, scales = "free_y", ncol = 3) +
  scale_color_manual(values = sex_traj_colors, name = "Sex") +
  labs(x = "Age", y = NULL) +
  plot.format + theme(legend.position = "none")

# ── Figure 1b: within-person longitudinal estimates ─────────────────────────
ids_repeat <- traj_df %>%
  group_by(participant_id) %>%
  summarise(n_visits = n_distinct(visit), .groups = "drop") %>%
  filter(n_visits >= 3)

set.seed(96)
n_traj <- min(300, nrow(ids_repeat))
ids_traj <- sample(ids_repeat$participant_id, n_traj)

plot_df_b <- traj_df %>% filter(participant_id %in% ids_traj)

plot_line_df <- plot_df_b %>%
  group_by(participant_id, domain) %>%
  filter(n_distinct(age) >= 2) %>%
  ungroup()

fit_withinperson_curve <- function(df) {
  if (nrow(df) < 50 || n_distinct(df$participant_id) < 10 ||
      n_distinct(df$age) < 5) {
    return(tibble())
  }

  model <- suppressWarnings(lme4::lmer(
    score ~ splines::ns(age, df = 4) + (1 | participant_id),
    data = df,
    control = lme4::lmerControl(
      check.conv.singular = "ignore",
      calc.derivs = FALSE
    )
  ))

  grid <- tibble(
    age = seq(
      quantile(df$age, 0.02, na.rm = TRUE),
      quantile(df$age, 0.98, na.rm = TRUE),
      length.out = 80
    )
  )
  grid$fit <- as.numeric(predict(model, newdata = grid, re.form = NA))
  grid
}

curve_df <- traj_df %>%
  filter(is.finite(age), is.finite(score), !is.na(participant_id)) %>%
  group_by(domain, sex) %>%
  group_split() %>%
  lapply(function(df) {
    curve <- fit_withinperson_curve(df)
    if (!nrow(curve)) return(curve)
    mutate(curve, domain = df$domain[[1]], sex = df$sex[[1]])
  }) %>%
  bind_rows() %>%
  mutate(
    domain = factor(domain, levels = domain_order),
    sex = factor(sex, levels = c("Female", "Male"))
  )

fig1b <- ggplot(plot_df_b, aes(x = age, y = score,
                               group = participant_id, color = sex)) +
  geom_line(data = plot_line_df, alpha = 0.1, linewidth = 0.35) +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_line(
    data = curve_df,
    aes(x = age, y = fit, color = sex, group = sex),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  facet_wrap(~ domain, scales = "free_y", ncol = 3) +
  scale_color_manual(values = sex_traj_colors, name = "Sex") +
  labs(x = "Age at repeated visit", y = NULL) +
  plot.format + theme(legend.position = "none")

# ── Figure 1e: Kaplan-Meier survival by IC tertile ───────────────────────────
surv_df <- read_csv(survival_path, show_col_types = FALSE)
require_columns(
  surv_df, c("ic_score", "time_years", "event"), "ic_survival.csv"
)

surv_df <- surv_df %>%
  filter(!is.na(ic_score), !is.na(time_years), !is.na(event), time_years >= 0) %>%
  mutate(
    score_group = ntile(ic_score, 3),
    score_group = factor(score_group, levels = 1:3,
                         labels = c("T1 (lowest IC)", "T2", "T3 (highest IC)"))
  )

km_fit <- survfit(Surv(time_years, event) ~ score_group, data = surv_df)

fig1e <- ggsurvplot(
  fit          = km_fit,
  data         = surv_df,
  risk.table   = TRUE,
  conf.int     = FALSE,
  pval         = TRUE,
  censor       = TRUE,
  ylim         = c(0.75, 1.00),
  palette      = c("#9CD5FF", "#7AAACE", "#355872"),
  xlab         = "Follow-up time (years)",
  ylab         = "Overall survival probability",
  legend.title = "Baseline IC score",
  legend.labs  = levels(surv_df$score_group),
  break.time.by = 2,
  ggtheme      = theme_bw(base_size = 12) +
    theme(
      panel.grid     = element_blank(),
      legend.position = "bottom",
      legend.title   = element_text(face = "bold"),
      plot.title     = element_text(face = "bold", hjust = 0.5)
    )
)

fig1e$table <- fig1e$table +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())

# Save the Kaplan-Meier plot and risk table.
svglite::svglite(file.path(fig_dir, "fig1e.svg"),
                 width = 6.5, height = 5, bg = "white")
print(fig1e)
dev.off()
png(file.path(fig_dir, "fig1e.png"), width = 6.5, height = 5,
    units = "in", res = 300, bg = "white")
print(fig1e)
dev.off()
message("Saved: fig1e (.svg + .png)")

# ── Figure 1c: IC domain correlation matrix ──────────────────────────────────
corr_df <- read_csv(correlation_path, show_col_types = FALSE)
if (
  nrow(corr_df) != ncol(corr_df) ||
  !all(vapply(corr_df, is.numeric, logical(1))) ||
  any(!is.finite(as.matrix(corr_df)))
) {
  stop(
    "ic_domain_corr.csv must be a square, finite numeric matrix.",
    call. = FALSE
  )
}
corr_mat <- as.matrix(corr_df)
rownames(corr_mat) <- colnames(corr_mat)

pdf(file.path(fig_dir, "fig1c.pdf"), width = 5, height = 5)
corrplot(
  corr_mat, type = "lower", diag = FALSE,
  method   = "color",
  col      = colorRampPalette(c("#1B9E77", "white", "#5E3C99"))(200),
  col.lim  = c(0, 1),
  addCoef.col = "black",
  number.cex  = 0.8,
  tl.col      = "black",
  tl.srt      = 45,
  mar         = c(0, 0, 1, 0)
)
dev.off()
message("Saved: fig1c.pdf")

save_fig(supp_fig1f, "supp_fig1f", width = 8, height = 5.6)
save_fig(fig1b, "fig1b", width = 8, height = 5.6)
