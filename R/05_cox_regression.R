# =============================================================================
# 05_cox_regression.R
# Cox proportional hazards regression — TCGA-BRCA overall survival
#
# Purpose:
#   Estimate the associations of age at diagnosis and pathological stage with
#   overall survival using a Cox proportional hazards model. Assess the
#   proportional hazards (PH) assumption via Schoenfeld residuals and
#   produce a forest plot of hazard ratios.
#
# Input:  data/processed/tcga_brca_survival_clean.csv
# Output: data/processed/cox_results.csv
#         data/processed/cox_ph_assumption_test.csv
#         figures/cox_forest_plot.png
# =============================================================================

library(survival)
library(broom)
library(tidyverse)
library(here)

# ---------------------------------------------------------------------------
# Load data and prepare model-ready dataset
# ---------------------------------------------------------------------------
df <- read_csv(
  here("data", "processed", "tcga_brca_survival_clean.csv"),
  show_col_types = FALSE
) |>
  filter(
    !is.na(stage_group),     # exclude patients with unknown stage
    !is.na(age),
    !is.na(os_event),
    !is.na(os_time_months),
    os_time_months > 0        # guard against zero follow-up (log undefined)
  ) |>
  mutate(
    # Set Stage I as the reference level so hazard ratios for Stages II–IV
    # are interpreted as relative to the lowest-risk stage group.
    stage_group = factor(
      stage_group,
      levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
    )
  )

# ---------------------------------------------------------------------------
# Fit Cox proportional hazards model
# ---------------------------------------------------------------------------
# Surv() defines the time-to-event outcome: follow-up in months (continuous)
# and the binary event indicator (1 = death, 0 = censored).
# The model adjusts simultaneously for age (continuous) and stage (ordinal),
# so each hazard ratio is conditional on the other predictor.

cox_model <- coxph(
  Surv(os_time_months, os_event) ~ age + stage_group,
  data = df
)

# tidy() extracts model coefficients; exponentiate = TRUE converts log-HR
# to HR scale for direct interpretation (HR = 1 means no association).
cox_results <- tidy(
  cox_model,
  exponentiate = TRUE,
  conf.int = TRUE
) |>
  mutate(
    estimate  = round(estimate,  3),
    conf.low  = round(conf.low,  3),
    conf.high = round(conf.high, 3),
    p.value   = signif(p.value,  3)
  )

# ---------------------------------------------------------------------------
# Proportional hazards assumption test (Schoenfeld residuals)
# ---------------------------------------------------------------------------
# The Cox model assumes the HR is constant over follow-up time (PH assumption).
# cox.zph() tests this by regressing scaled Schoenfeld residuals on time.
# A significant p-value suggests the HR may vary over follow-up,
# in which case results should be interpreted as time-averaged associations.

ph_test <- cox.zph(cox_model)

ph_results <- as.data.frame(ph_test$table) |>
  rownames_to_column("term") |>
  as_tibble() |>
  rename(
    chisq   = chisq,
    p_value = p
  ) |>
  mutate(
    chisq   = round(chisq,   3),
    p_value = signif(p_value, 3)
  )

# ---------------------------------------------------------------------------
# Export model results
# ---------------------------------------------------------------------------
write_csv(cox_results, here("data", "processed", "cox_results.csv"))
write_csv(ph_results,  here("data", "processed", "cox_ph_assumption_test.csv"))

cat("\nCox proportional hazards model results:\n")
print(cox_results)

cat("\nProportional hazards assumption test (Schoenfeld residuals):\n")
print(ph_results)

cat("\nCox model outputs saved to data/processed/\n")

# ---------------------------------------------------------------------------
# Forest plot of hazard ratios
# ---------------------------------------------------------------------------
# Each point is a HR estimate; horizontal bars are 95% CIs.
# The dashed vertical line at HR = 1 marks the null (no association).
# The x-axis is on the log scale so that symmetric CIs appear visually
# symmetric around the point estimate.

forest_data <- cox_results |>
  mutate(
    term_label = case_when(
      term == "age"                      ~ "Age, per year",
      term == "stage_groupStage II"      ~ "Stage II vs Stage I",
      term == "stage_groupStage III"     ~ "Stage III vs Stage I",
      term == "stage_groupStage IV"      ~ "Stage IV vs Stage I",
      TRUE ~ term
    ),
    # Annotation text shown next to each point
    hr_label = paste0(estimate, " (", conf.low, "–", conf.high, ")")
  )

forest_plot <- ggplot(
  forest_data,
  aes(x = estimate, y = reorder(term_label, estimate))
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height   = 0.18,
    linewidth = 0.8,
    color    = "gray30"
  ) +
  geom_point(size = 3, color = "black") +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 20),
    limits = c(0.9, 30)
  ) +
  labs(
    title    = "Cox proportional hazards model",
    subtitle = "Overall survival among TCGA-BRCA patients",
    x        = "Hazard ratio (log scale)",
    y        = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    text             = element_text(color = "black"),
    axis.text        = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(
  here("figures", "cox_forest_plot.png"),
  forest_plot,
  width  = 8,
  height = 5,
  dpi    = 300,
  bg     = "white"
)

cat("\nForest plot saved to figures/cox_forest_plot.png\n")
