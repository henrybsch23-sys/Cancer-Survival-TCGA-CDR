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