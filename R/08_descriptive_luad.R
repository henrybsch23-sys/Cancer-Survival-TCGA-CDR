# =============================================================================
# 08_descriptive_luad.R
# Cohort descriptive statistics — TCGA-LUAD
#
# Purpose:
#   Generate tabular summaries of the cleaned TCGA-LUAD cohort for use in the
#   Streamlit dashboard: overall cohort metrics, stage distribution, gender
#   breakdown, and age-group breakdown. Mirrors 03_descriptive_analysis.R.
#
# Input:  data/processed/luad_survival_clean.csv
# Output: data/processed/luad_cohort_summary.csv
#         data/processed/luad_stage_summary.csv
#         data/processed/luad_gender_summary.csv
#         data/processed/luad_age_summary.csv
# =============================================================================

library(tidyverse)
library(here)

df <- read_csv(
  here("data", "processed", "luad_survival_clean.csv"),
  show_col_types = FALSE
)

cohort_summary <- tibble(
  metric = c(
    "Number of patients",
    "Number of OS events",
    "OS event proportion",
    "Median follow-up, months",
    "Median age at diagnosis",
    "Patients with missing stage"
  ),
  value = c(
    nrow(df),
    sum(df$os_event, na.rm = TRUE),
    round(mean(df$os_event, na.rm = TRUE), 3),
    round(median(df$os_time_months, na.rm = TRUE), 1),
    round(median(df$age, na.rm = TRUE), 1),
    sum(is.na(df$stage_group))
  )
)

stage_summary <- df |>
  count(stage_group, name = "n") |>
  mutate(percent = round(100 * n / sum(n), 1))

gender_summary <- df |>
  count(gender, name = "n") |>
  mutate(percent = round(100 * n / sum(n), 1))

age_summary <- df |>
  count(age_group, name = "n") |>
  mutate(percent = round(100 * n / sum(n), 1))

write_csv(cohort_summary, here("data", "processed", "luad_cohort_summary.csv"))
write_csv(stage_summary,  here("data", "processed", "luad_stage_summary.csv"))
write_csv(gender_summary, here("data", "processed", "luad_gender_summary.csv"))
write_csv(age_summary,    here("data", "processed", "luad_age_summary.csv"))

cat("\nCohort summary (LUAD):\n")
print(cohort_summary)

cat("\nStage summary:\n")
print(stage_summary)

cat("\nDescriptive summary files saved to data/processed/ (luad_ prefix)\n")
