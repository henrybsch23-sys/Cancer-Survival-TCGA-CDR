# =============================================================================
# 07_clean_luad.R
# Data management and cleaning — TCGA-LUAD cohort
#
# Purpose:
#   Read the raw TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR) Excel file,
#   subset to lung adenocarcinoma patients (TCGA-LUAD), recode and clean the
#   variables needed for survival analysis, and write an analysis-ready CSV.
#   Variable selection mirrors 02_clean_tcga_cdr.R for direct comparability.
#
# Input:  data/raw/TCGA-CDR-SupplementalTableS1.xlsx
# Output: data/processed/luad_survival_clean.csv
# =============================================================================

library(readxl)
library(tidyverse)
library(janitor)
library(here)

file_path <- here("data", "raw", "TCGA-CDR-SupplementalTableS1.xlsx")

raw <- read_excel(file_path, sheet = "TCGA-CDR") |>
  clean_names()

missing_like_values <- c(
  "#N/A",
  "[Not Available]",
  "[Unknown]",
  "Not Available",
  "Unknown",
  "Not Reported",
  ""
)

# ---------------------------------------------------------------------------
# Subset to TCGA-LUAD and select/recode analysis variables
# ---------------------------------------------------------------------------
clean_luad <- raw |>
  filter(type == "LUAD") |>
  transmute(
    patient_id  = bcr_patient_barcode,
    cancer_type = type,

    age         = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),

    gender              = na_if(gender, "#N/A"),
    race                = na_if(race, "#N/A"),
    histological_type   = na_if(histological_type, "#N/A"),

    stage_raw           = na_if(ajcc_pathologic_tumor_stage, "#N/A"),

    os_event        = suppressWarnings(as.numeric(os)),
    os_time_days    = suppressWarnings(as.numeric(os_time)),
    os_time_months  = os_time_days / 30.44
  ) |>
  mutate(
    across(
      where(is.character),
      ~ if_else(.x %in% missing_like_values, NA_character_, .x)
    )
  ) |>
  mutate(
    stage_group = case_when(
      str_detect(stage_raw, "Stage IV")  ~ "Stage IV",
      str_detect(stage_raw, "Stage III") ~ "Stage III",
      str_detect(stage_raw, "Stage II")  ~ "Stage II",
      str_detect(stage_raw, "Stage I")   ~ "Stage I",
      TRUE                               ~ NA_character_
    ),

    age_group = case_when(
      age < 50              ~ "<50",
      age >= 50 & age < 65  ~ "50-64",
      age >= 65             ~ "65+",
      TRUE                  ~ NA_character_
    )
  ) |>
  filter(
    !is.na(os_event),
    !is.na(os_time_days),
    os_time_days > 0
  )

# ---------------------------------------------------------------------------
# Quick data quality checks
# ---------------------------------------------------------------------------
cat("\nNumber of LUAD patients:", nrow(clean_luad), "\n")
cat("Number of OS events:", sum(clean_luad$os_event, na.rm = TRUE), "\n")
cat("Median follow-up in months:", median(clean_luad$os_time_months, na.rm = TRUE), "\n\n")

cat("Stage distribution (including NA):\n")
print(table(clean_luad$stage_group, useNA = "ifany"))

cat("\nGender distribution:\n")
print(table(clean_luad$gender, useNA = "ifany"))

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
write_csv(
  clean_luad,
  here("data", "processed", "luad_survival_clean.csv")
)

cat("\nCleaned LUAD survival dataset saved to data/processed/luad_survival_clean.csv\n")
