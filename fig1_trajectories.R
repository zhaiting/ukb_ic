# ============================================================================
# fig1_trajectories.R — Figure 1
#
# Input:  data/ic_trajectories.csv   (long format: participant_id, domain,
#             score, age, sex, visit, years_since_baseline)
#         data/ic_survival.csv       (participant_id, ic_score, time_years, event)
#         data/ic_domain_corr.csv    (6×6 correlation matrix of IC domain scores)
#
# Output: output/fig1_trajectories.svg/.png   (Panels A + B)
#         output/fig1_km.svg/.png              (Panel C)
#         output/fig1_corrplot.pdf             (Panel D)
#
# Prerequisites: source("00_helpers.R")
#                install.packages(c("survminer", "corrplot"))
# ============================================================================

library(survival)
library(survminer)
library(corrplot)

# ── Data ───────────────────────────────────────────────────────────────────
traj_df <- read_csv(file.path(data_dir, "ic_trajectories.csv"),
                    show_col_types = FALSE)

# ── Theme (matches original plot.format) ───────────────────────────────────
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

# ── Panel A: Cross-sectional ─────────────────────────────────────────────────
set.seed(96)
n_sample <- min(20000, n_distinct(traj_df$participant_id))
ids_sample <- sample(unique(traj_df$participant_id), n_sample)

plot_df_a <- traj_df %>%
  filter(participant_id %in% ids_sample, visit == 0)

fig1a <- ggplot(plot_df_a, aes(x = age, y = score, color = sex)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              se = FALSE, linewidth = 0.7) +
  facet_wrap(~ domain, scales = "free_y", ncol = 3) +
  scale_color_manual(values = sex_traj_colors, name = "Sex") +
  labs(x = "Age", y = NULL) +
  plot.format + theme(legend.position = "none")

# ── Panel B: Longitudinal ───────────────────────────────────────────────────
ids_repeat <- traj_df %>%
  group_by(participant_id) %>%
  summarise(n_visits = n_distinct(visit), .groups = "drop") %>%
  filter(n_visits >= 3)

set.seed(96)
n_traj <- min(300, nrow(ids_repeat))
ids_traj <- sample(ids_repeat$participant_id, n_traj)

plot_df_b <- traj_df %>% filter(participant_id %in% ids_traj)

fig1b <- ggplot(plot_df_b, aes(x = years_since_baseline, y = score,
                                group = participant_id, color = sex)) +
  geom_line(alpha = 0.1, linewidth = 0.35) +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_smooth(aes(group = sex), method = "gam", formula = y ~ s(x, k = 4),
              se = FALSE, linewidth = 0.9) +
  facet_wrap(~ domain, scales = "free_y", ncol = 3) +
  scale_color_manual(values = sex_traj_colors, name = "Sex") +
  labs(x = "Age at repeated visit", y = NULL) +
  plot.format + theme(legend.position = "none")

# ── Panel C: Kaplan-Meier survival by IC tertile ─────────────────────────────
surv_df <- read_csv(file.path(data_dir, "ic_survival.csv"),
                    show_col_types = FALSE)

surv_df <- surv_df %>%
  filter(!is.na(ic_score), !is.na(time_years), !is.na(event), time_years >= 0) %>%
  mutate(
    score_group = ntile(ic_score, 3),
    score_group = factor(score_group, levels = 1:3,
                         labels = c("T1 (lowest IC)", "T2", "T3 (highest IC)"))
  )

km_fit <- survfit(Surv(time_years, event) ~ score_group, data = surv_df)

fig1c <- ggsurvplot(
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

fig1c$table <- fig1c$table +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())

# survminer returns a list; save via print()
svg(file.path(fig_dir, "fig1_km.svg"), width = 6.5, height = 5)
print(fig1c)
dev.off()
png(file.path(fig_dir, "fig1_km.png"), width = 6.5, height = 5,
    units = "in", res = 300)
print(fig1c)
dev.off()
message("Saved: fig1_km (.svg + .png)")

# ── Panel D: IC domain correlation matrix ─────────────────────────────────────
corr_mat <- as.matrix(read_csv(file.path(data_dir, "ic_domain_corr.csv"),
                               show_col_types = FALSE))
rownames(corr_mat) <- colnames(corr_mat)

pdf(file.path(fig_dir, "fig1_corrplot.pdf"), width = 5, height = 5)
corrplot(
  corr_mat, type = "lower", diag = FALSE,
  method   = "color",
  col      = colorRampPalette(c("#1B9E77", "white", "#5E3C99"))(200),
  cl.lim   = c(0, 1),
  addCoef.col = "black",
  number.cex  = 0.8,
  tl.col      = "black",
  tl.srt      = 45,
  mar         = c(0, 0, 1, 0)
)
dev.off()
message("Saved: fig1_corrplot.pdf")

# ── Panels A + B combined ────────────────────────────────────────────────────
fig1ab <- fig1a / fig1b + plot_layout(heights = c(1, 1))

save_fig(fig1ab, "fig1_trajectories", width = 10, height = 14)
