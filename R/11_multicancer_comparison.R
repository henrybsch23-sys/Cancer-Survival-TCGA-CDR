# =============================================================================
# 11_multicancer_comparison.R
# Multi-cancer survival comparison — BRCA vs LUAD
#
# Purpose:
#   Combine the cleaned BRCA and LUAD cohorts to produce cross-cancer
#   comparison outputs: a cohort summary table, an overlaid Kaplan-Meier
#   plot (one curve per cancer type), a side-by-side Cox hazard ratio
#   comparison figure, and a stacked Cox results CSV.
#
# Input:  data/processed/brca_survival_clean.csv  (from 02)
#         data/processed/luad_survival_clean.csv  (from 07)
# Output: data/processed/multicancer_cohort_comparison.csv
#         data/processed/multicancer_cox_results.csv
#         figures/multicancer_km_comparison.png
#         figures/multicancer_cox_comparison.png
# =============================================================================

library(survival)
library(survminer)
library(broom)
library(tidyverse)
library(here)

brca <- read_csv(
  here("data", "processed", "brca_survival_clean.csv"),
  show_col_types = FALSE
)

luad <- read_csv(
  here("data", "processed", "luad_survival_clean.csv"),
  show_col_types = FALSE
)

# ---------------------------------------------------------------------------
# 1. Combined cohort summary table
# ---------------------------------------------------------------------------
make_summary <- function(df, label) {
  tibble(
    cancer_type              = label,
    n_patients               = nrow(df),
    n_os_events              = sum(df$os_event, na.rm = TRUE),
    os_event_rate_pct        = round(100 * mean(df$os_event, na.rm = TRUE), 1),
    median_follow_up_months  = round(median(df$os_time_months, na.rm = TRUE), 1),
    median_age               = round(median(df$age, na.rm = TRUE), 1),
    n_missing_stage          = sum(is.na(df$stage_group))
  )
}

cohort_comparison <- bind_rows(
  make_summary(brca, "TCGA-BRCA"),
  make_summary(luad, "TCGA-LUAD")
)

write_csv(cohort_comparison, here("data", "processed", "multicancer_cohort_comparison.csv"))

cat("\nCohort comparison:\n")
print(cohort_comparison)

# ---------------------------------------------------------------------------
# 2. Combined Cox results (stacked, cancer_type column added)
# ---------------------------------------------------------------------------
fit_cox <- function(df) {
  df_model <- df |>
    filter(
      !is.na(stage_group),
      !is.na(age),
      os_time_months > 0
    ) |>
    mutate(
      stage_group = factor(
        stage_group,
        levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
      )
    )
  coxph(Surv(os_time_months, os_event) ~ age + stage_group, data = df_model)
}

cox_brca <- fit_cox(brca)
cox_luad <- fit_cox(luad)

tidy_cox <- function(model, label) {
  tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    mutate(
      cancer_type = label,
      estimate    = round(estimate,  3),
      conf.low    = round(conf.low,  3),
      conf.high   = round(conf.high, 3),
      p.value     = signif(p.value,  3),
      term_label  = case_when(
        term == "age"                      ~ "Age, per year",
        term == "stage_groupStage II"      ~ "Stage II vs I",
        term == "stage_groupStage III"     ~ "Stage III vs I",
        term == "stage_groupStage IV"      ~ "Stage IV vs I",
        TRUE ~ term
      )
    ) |>
    select(cancer_type, term, term_label, estimate, conf.low, conf.high, p.value)
}

multicancer_cox <- bind_rows(
  tidy_cox(cox_brca, "TCGA-BRCA"),
  tidy_cox(cox_luad, "TCGA-LUAD")
)

write_csv(multicancer_cox, here("data", "processed", "multicancer_cox_results.csv"))

cat("\nMulti-cancer Cox results saved.\n")

# ---------------------------------------------------------------------------
# 3. Overlaid Kaplan-Meier plot (one curve per cancer type)
# ---------------------------------------------------------------------------
# A single curve per cancer type (ignoring stage) gives the clearest visual
# contrast in overall prognosis between BRCA and LUAD.

combined <- bind_rows(
  brca |> select(cancer_type, os_time_months, os_event),
  luad |> select(cancer_type, os_time_months, os_event)
) |>
  mutate(cancer_type = factor(cancer_type, levels = c("BRCA", "LUAD")))

km_fit <- survfit(Surv(os_time_months, os_event) ~ cancer_type, data = combined)

png(here("figures", "multicancer_km_comparison.png"), width = 1100, height = 850)
print(
  ggsurvplot(
    km_fit,
    data         = combined,
    risk.table   = TRUE,
    pval         = TRUE,
    xlab         = "Time since diagnosis (months)",
    ylab         = "Overall survival probability",
    legend.title = "Cancer type",
    legend.labs  = c("TCGA-BRCA", "TCGA-LUAD"),
    title        = "Overall survival: BRCA vs LUAD (TCGA-CDR)",
    palette      = c("#2166ac", "#d7191c"),
    conf.int     = TRUE
  )
)
dev.off()

cat("Multi-cancer KM plot saved to figures/multicancer_km_comparison.png\n")

# ---------------------------------------------------------------------------
# 4. Side-by-side Cox forest plot (BRCA vs LUAD)
# ---------------------------------------------------------------------------
forest_comparison <- multicancer_cox |>
  mutate(
    cancer_type = factor(cancer_type, levels = c("TCGA-BRCA", "TCGA-LUAD")),
    term_label  = factor(
      term_label,
      levels = rev(c("Age, per year", "Stage II vs I", "Stage III vs I", "Stage IV vs I"))
    )
  )

cox_comparison_plot <- ggplot(
  forest_comparison,
  aes(x = estimate, y = term_label, color = cancer_type, shape = cancer_type)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height    = 0.2,
    linewidth = 0.8,
    position  = position_dodge(width = 0.5)
  ) +
  geom_point(
    size     = 3.5,
    position = position_dodge(width = 0.5)
  ) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20),
    limits = c(0.85, 35)
  ) +
  scale_color_manual(values = c("TCGA-BRCA" = "#2166ac", "TCGA-LUAD" = "#d7191c")) +
  scale_shape_manual(values = c("TCGA-BRCA" = 16, "TCGA-LUAD" = 17)) +
  labs(
    title    = "Cox model hazard ratios — BRCA vs LUAD",
    subtitle = "Model: age + pathological stage (Stage I as reference)",
    x        = "Hazard ratio (log scale)",
    y        = NULL,
    color    = "Cancer type",
    shape    = "Cancer type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    text             = element_text(color = "black"),
    axis.text        = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave(
  here("figures", "multicancer_cox_comparison.png"),
  cox_comparison_plot,
  width  = 9,
  height = 5.5,
  dpi    = 300,
  bg     = "white"
)

cat("Multi-cancer Cox comparison plot saved to figures/multicancer_cox_comparison.png\n")
