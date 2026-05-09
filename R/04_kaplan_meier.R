# =============================================================================
# 04_kaplan_meier.R
# Kaplan-Meier survival estimation and log-rank test — TCGA-BRCA
#
# Purpose:
#   Estimate overall survival curves stratified by pathological tumor stage
#   using the Kaplan-Meier (product-limit) estimator.  A log-rank test is
#   used to evaluate whether survival distributions differ across stage groups.
#   The output figure is used in the Streamlit dashboard.
#
# Input:  data/processed/tcga_brca_survival_clean.csv
# Output: figures/km_stage_plot.png
# =============================================================================

library(survival)
library(survminer)
library(tidyverse)
library(here)

df <- read_csv(here("data", "processed", "tcga_brca_survival_clean.csv"),
               show_col_types = FALSE) |>
  # Remove patients with unknown stage so the KM plot shows only labelled groups.
  # The 24 patients with stage_group = NA (mostly "Stage X" raw codes) would
  # otherwise appear as a separate unlabelled stratum.
  filter(!is.na(stage_group))

# ---------------------------------------------------------------------------
# Kaplan-Meier estimator stratified by pathological stage
# ---------------------------------------------------------------------------
# survfit() computes the product-limit estimate of the survival function
# S(t) = P(T > t) within each stage group separately.
# Patients without an observed event by their last follow-up are censored:
# their exact survival time is unknown but contributes information up to that
# point.

fit <- survfit(Surv(os_time_months, os_event) ~ stage_group, data = df)

# ---------------------------------------------------------------------------
# Plot and export
# ---------------------------------------------------------------------------
# ggsurvplot() overlays all stage-group curves on one panel.
# risk.table = TRUE adds the at-risk count table below the plot, which is
# standard practice in clinical epidemiology publications.
# pval = TRUE displays the log-rank p-value testing the null hypothesis that
# all survival curves are identical.

png(here("figures", "km_stage_plot.png"), width = 1000, height = 800)
print(
  ggsurvplot(
    fit,
    data       = df,
    risk.table = TRUE,   # at-risk counts by time point
    pval       = TRUE,   # log-rank p-value
    xlab       = "Time since diagnosis (months)",
    ylab       = "Overall survival probability",
    legend.title = "Pathological stage",
    palette    = c("#2166ac", "#4dac26", "#d7191c", "#762a83")
  )
)
dev.off()

cat("\nKaplan-Meier plot saved to figures/km_stage_plot.png\n")
# Result: curves show progressively lower survival with advancing stage;
# the log-rank test is statistically significant, motivating the Cox model
# in 05_cox_regression.R.
