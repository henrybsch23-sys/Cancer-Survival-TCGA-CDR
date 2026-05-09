# =============================================================================
# 01_inspect_excel_sheets.R
# Inspect the raw TCGA-CDR Excel workbook
#
# Purpose:
#   List all sheet names in the TCGA Pan-Cancer Clinical Data Resource (CDR)
#   supplemental Excel file and preview the first five rows and column names
#   of each sheet. Used to orient the data before any cleaning is carried out.
#
# Input:  data/raw/TCGA-CDR-SupplementalTableS1.xlsx
# =============================================================================

library(readxl)
library(tidyverse)
library(here)

file_path <- here("data", "raw", "TCGA-CDR-SupplementalTableS1.xlsx")

sheets <- excel_sheets(file_path)
print(sheets)

for (s in sheets) {
  cat("\n\n--- Sheet:", s, "---\n")
  temp <- read_excel(file_path, sheet = s, n_max = 5)
  print(names(temp))
}