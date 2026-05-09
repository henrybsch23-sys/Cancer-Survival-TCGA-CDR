# =============================================================================
# 02_clean_tcga_cdr.R
# Data management and cleaning — TCGA-BRCA cohort
#
# Purpose:
#   Read the raw TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR) Excel file,
#   subset to breast cancer patients (TCGA-BRCA), recode and clean the
#   variables needed for survival analysis, and write an analysis-ready CSV.
#
# Input:  data/raw/TCGA-CDR-SupplementalTableS1.xlsx
# Output: data/processed/tcga_brca_survival_clean.csv
# =============================================================================

library(readxl)
library(tidyverse)
library(janitor)
library(here)

file_path <- here("data", "raw", "TCGA-CDR-SupplementalTableS1.xlsx")

# Read the CDR sheet and standardise all column names to snake_case
raw <- read_excel(file_path, sheet = "TCGA-CDR") |>
  clean_names()

# Inspect available cancer types before subsetting
print(sort(unique(raw$type)))

# Establishing missing values
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
# Subset to TCGA-BRCA and select/recode analysis variables
# ---------------------------------------------------------------------------
# suppressWarnings() is used for as.numeric() coercion because some cells
# contain "#N/A" strings exported from Excel, which produce expected NAs.
# These are also handled explicitly with na_if() for character columns.

clean_brca <- raw |>
  filter(type == "BRCA") |>
  transmute(
    patient_id  = bcr_patient_barcode,
    cancer_type = type,

    # Continuous age; coerce to numeric (Excel sometimes stores as character)
    age         = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),

    # Replace "#N/A" strings (Excel artefact) with proper NA
    gender              = na_if(gender, "#N/A"),
    race                = na_if(race, "#N/A"),
    menopause_status    = na_if(menopause_status, "#N/A"),
    histological_type   = na_if(histological_type, "#N/A"),
    
    stage_raw           = na_if(ajcc_pathologic_tumor_stage, "#N/A"),

    # Overall survival: 1 = death, 0 = censored (alive at last follow-up)
    os_event        = suppressWarnings(as.numeric(os)),
    os_time_days    = suppressWarnings(as.numeric(os_time)),

    # Convert follow-up to months; 30.44 = average days per calendar month
    os_time_months  = os_time_days / 30.44
  ) |>
  mutate(
    # Applying the distinction and missing values
    across(
      where(is.character),
      ~ if_else(.x %in% missing_like_values, NA_character_, .x)
    )
  ) |>
  mutate(
    # Collapse detailed AJCC sub-stages (e.g. "Stage IIA", "Stage IIIC") into
    # four main ordinal groups used in the survival model.
    # str_detect() matches the most specific pattern first (IV before I).
    stage_group = case_when(
      str_detect(stage_raw, "Stage IV")  ~ "Stage IV",
      str_detect(stage_raw, "Stage III") ~ "Stage III",
      str_detect(stage_raw, "Stage II")  ~ "Stage II",
      str_detect(stage_raw, "Stage I")   ~ "Stage I",
      TRUE                               ~ NA_character_
    ),

    # Clinically meaningful age bands for descriptive tables
    age_group = case_when(
      age < 50              ~ "<50",
      age >= 50 & age < 65  ~ "50-64",
      age >= 65             ~ "65+",
      TRUE                  ~ NA_character_
    )
  ) |>
  # Exclude patients with missing survival outcome or zero/negative follow-up.
  # Zero follow-up causes undefined log(time) in survival models.
  filter(
    !is.na(os_event),
    !is.na(os_time_days),
    os_time_days > 0
  )

# ---------------------------------------------------------------------------
# Quick data quality checks
# ---------------------------------------------------------------------------
cat("\nNumber of BRCA patients:", nrow(clean_brca), "\n")
cat("Number of OS events:", sum(clean_brca$os_event, na.rm = TRUE), "\n")
cat("Median follow-up in months:", median(clean_brca$os_time_months, na.rm = TRUE), "\n\n")

cat("Stage distribution (including NA):\n")
print(table(clean_brca$stage_group, useNA = "ifany"))

cat("\nGender distribution:\n")
print(table(clean_brca$gender, useNA = "ifany"))

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
write_csv(
  clean_brca,
  here("data", "processed", "tcga_brca_survival_clean.csv")
)

cat("\nCleaned BRCA survival dataset saved to data/processed/tcga_brca_survival_clean.csv\n")
