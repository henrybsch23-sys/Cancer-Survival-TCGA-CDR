# =============================================================================
# 12_export_for_streamlit_v3.R
# Validate and confirm all v3 Streamlit input files — LUAD + multi-cancer
#
# Purpose:
#   Check that all v3 processed files produced by scripts 07-11 are present
#   before launching the updated Streamlit dashboard.
#
# Input:  All files produced by R/07 through R/11
# =============================================================================

library(here)

v3_required <- c(
  # LUAD-specific files
  "luad_survival_clean.csv",
  "luad_cohort_summary.csv",
  "luad_stage_summary.csv",
  "luad_gender_summary.csv",
  "luad_age_summary.csv",
  "luad_cox_results.csv",
  "luad_cox_ph_assumption_test.csv",
  # Multi-cancer comparison files
  "multicancer_cohort_comparison.csv",
  "multicancer_cox_results.csv"
)

v3_figures <- c(
  "luad_km_stage_plot.png",
  "luad_cox_forest_plot.png",
  "multicancer_km_comparison.png",
  "multicancer_cox_comparison.png"
)

processed_dir <- here("data", "processed")
figures_dir   <- here("figures")

missing_csv <- v3_required[
  !file.exists(file.path(processed_dir, v3_required))
]

missing_fig <- v3_figures[
  !file.exists(file.path(figures_dir, v3_figures))
]

if (length(missing_csv) > 0) {
  stop("Missing v3 CSV files: ", paste(missing_csv, collapse = ", "))
}

if (length(missing_fig) > 0) {
  stop("Missing v3 figure files: ", paste(missing_fig, collapse = ", "))
}

cat("\nAll v3 Streamlit input files are present.\n")
cat("CSV files checked:", length(v3_required), "\n")
cat("Figure files checked:", length(v3_figures), "\n")
cat("\nv3 pipeline complete. You can now launch the Streamlit dashboard.\n")
