# =============================================================================
# 06_export_for_streamlit.R
# Validate and export all Streamlit app input files — TCGA-BRCA
#
# Purpose:
#   Verify that all required processed files produced by the preceding R
#   scripts (02–05) are present, then write an app_summary.csv containing a
#   high-level cohort overview for display in the Streamlit dashboard.
#
# Input:  data/processed/tcga_brca_survival_clean.csv  (from 02)
#         data/processed/cohort_summary.csv            (from 03)
#         data/processed/stage_summary.csv             (from 03)
#         data/processed/gender_summary.csv            (from 03)
#         data/processed/age_summary.csv               (from 03)
#         data/processed/cox_results.csv               (from 05)
#         data/processed/cox_ph_assumption_test.csv    (from 05)
# Output: data/processed/app_summary.csv
# =============================================================================

library(tidyverse)
library(here)

required_files <- c(
  "tcga_brca_survival_clean.csv",
  "cohort_summary.csv",
  "stage_summary.csv",
  "gender_summary.csv",
  "age_summary.csv",
  "cox_results.csv",
  "cox_ph_assumption_test.csv"
)

processed_dir <- here("data", "processed")

missing_files <- required_files[
  !file.exists(file.path(processed_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing required files: ",
    paste(missing_files, collapse = ", ")
  )
}

df <- read_csv(
  here("data", "processed", "tcga_brca_survival_clean.csv"),
  show_col_types = FALSE
)

app_summary <- tibble(
  metric = c(
    "Cancer type",
    "Patients",
    "Overall survival events",
    "Median follow-up, months",
    "Median age at diagnosis"
  ),
  value = c(
    "TCGA-BRCA",
    as.character(nrow(df)),
    as.character(sum(df$os_event, na.rm = TRUE)),
    as.character(round(median(df$os_time_months, na.rm = TRUE), 1)),
    as.character(round(median(df$age, na.rm = TRUE), 1))
  )
)

write_csv(
  app_summary,
  here("data", "processed", "app_summary.csv")
)

cat("\nAll Streamlit input files are ready.\n")
cat("Saved data/processed/app_summary.csv\n")