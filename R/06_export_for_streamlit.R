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