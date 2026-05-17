# =============================================================================
# 09_kaplan_meier_luad.R
# Kaplan-Meier survival estimation and log-rank test — TCGA-LUAD
#
# Purpose:
#   Estimate overall survival curves stratified by pathological tumour stage
#   for LUAD patients. Mirrors 04_kaplan_meier.R for direct comparability
#   with BRCA.
#
# Input:  data/processed/luad_survival_clean.csv
# Output: figures/luad_km_stage_plot.png
# =============================================================================

library(survival)
library(survminer)
library(tidyverse)
library(here)

df <- read_csv(
  here("data", "processed", "luad_survival_clean.csv"),
  show_col_types = FALSE
) |>
  filter(!is.na(stage_group))

fit <- survfit(Surv(os_time_months, os_event) ~ stage_group, data = df)

png(here("figures", "luad_km_stage_plot.png"), width = 1000, height = 800)
print(
  ggsurvplot(
    fit,
    data         = df,
    risk.table   = TRUE,
    pval         = TRUE,
    xlab         = "Time since diagnosis (months)",
    ylab         = "Overall survival probability",
    legend.title = "Pathological stage",
    title        = "TCGA-LUAD — Overall survival by pathological stage",
    palette      = c("#2166ac", "#4dac26", "#d7191c", "#762a83")
  )
)
dev.off()

cat("\nKaplan-Meier plot saved to figures/luad_km_stage_plot.png\n")
