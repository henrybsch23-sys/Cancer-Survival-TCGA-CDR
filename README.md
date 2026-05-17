# TCGA Breast Cancer Survival Analysis

## Live dashboard

A deployed version of the Streamlit dashboard is available here:

https://cancer-survival-tcga-cdr-gzxwsq2dbktxmjtlet76rz.streamlit.app/

The dashboard presents cohort summaries, Kaplan–Meier survival curves, Cox proportional hazards regression results, model diagnostics, a machine learning risk prediction demo, and visual interpretations of survival patterns among TCGA-BRCA breast cancer patients.

---

## Overview

This project demonstrates a reproducible survival analysis and machine learning workflow using the open-access **TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR)**. The analysis focuses on breast cancer patients (TCGA-BRCA) and covers the full pipeline from raw data to an interactive dashboard.

The project is split across two languages by design:

- **R** — data cleaning, descriptive statistics, Kaplan–Meier survival estimation, Cox proportional hazards regression, and model diagnostics.
- **Python** — data validation notebooks, machine learning risk prediction, and an interactive Streamlit dashboard.

---

## Motivation

This project was developed as a portfolio project for research roles in cancer epidemiology, biostatistics, precision medicine, and data-driven health research. It demonstrates practical skills in:

- Clinical data cleaning and reproducible pipelines
- Survival analysis (Kaplan–Meier, log-rank, Cox regression)
- Proportional hazards diagnostics (Schoenfeld residuals)
- Censoring-aware machine learning outcome definition
- Supervised classification with imbalanced clinical data
- Model evaluation (AUC, calibration, confusion matrix)
- Interactive communication of epidemiological results

---

## Dataset

The project uses the open-access TCGA Pan-Cancer Clinical Data Resource:

```text
TCGA-CDR-SupplementalTableS1.xlsx
```

Download from the NCI GDC PanCanAtlas publication resources and place in:

```text
data/raw/TCGA-CDR-SupplementalTableS1.xlsx
```

Only open-access clinical and survival outcome data are used. No controlled-access genomic data, raw sequencing data, or identifiable patient data are used.

**Reference:** Liu et al. (2018). An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. *Cell*, 173(2), 400–416.

---

## Research Questions

1. Among TCGA-BRCA patients, are **age at diagnosis** and **pathological tumour stage** associated with overall survival?
2. Can a **censoring-aware** machine learning classifier predict 3-year mortality using clinical variables available at diagnosis?

---

## Methods

### Part 1 — BRCA survival analysis (R, v1)

| Step | Script | Description |
|---|---|---|
| 1 | `R/01_inspect_excel_sheets.R` | Inspect sheet names and preview column headers |
| 2 | `R/02_clean_tcga_cdr.R` | Subset to BRCA, recode variables, handle missing values, export cleaned CSV |
| 3 | `R/03_descriptive_analysis.R` | Generate cohort summary, stage, age, and gender tables |
| 4 | `R/04_kaplan_meier.R` | Kaplan–Meier curves by stage with log-rank test |
| 5 | `R/05_cox_regression.R` | Cox regression (age + stage), Schoenfeld PH test, forest plot |
| 6 | `R/06_export_for_streamlit.R` | Validate all BRCA processed files and export app summary |

### Part 2 — Python data validation and ML (notebooks, v2)

| Notebook | Description |
|---|---|
| `notebooks/01_python_data_check.ipynb` | Load cleaned cohort, inspect missingness, validate distributions |
| `notebooks/02_ml_risk_prediction_demo.ipynb` | Censoring-aware 3-year mortality prediction with Logistic Regression, Random Forest, and Gradient Boosting |

### Part 3 — LUAD pipeline and multi-cancer comparison (R, v3)

| Step | Script | Description |
|---|---|---|
| 7 | `R/07_clean_luad.R` | Subset to LUAD, recode variables, export cleaned CSV |
| 8 | `R/08_descriptive_luad.R` | LUAD cohort summary, stage, age, and gender tables |
| 9 | `R/09_kaplan_meier_luad.R` | Kaplan–Meier curves by stage for LUAD |
| 10 | `R/10_cox_regression_luad.R` | Cox regression (age + stage) for LUAD, forest plot |
| 11 | `R/11_multicancer_comparison.R` | Overlaid KM plot, side-by-side Cox forest plot, combined cohort table |
| 12 | `R/12_export_for_streamlit_v3.R` | Validate all v3 processed files |

The ML outcome definition excludes patients censored before 3 years (53.2 % of the cohort) because their 3-year status is unknown. This avoids outcome misclassification bias and mirrors the landmark analysis approach used in clinical epidemiology.

---

## Analysis Variables

| Original TCGA-CDR column | Cleaned variable | Use |
|---|---|---|
| `bcr_patient_barcode` | `patient_id` | Patient identifier |
| `type` | `cancer_type` | Cancer type filter |
| `age_at_initial_pathologic_diagnosis` | `age` | Continuous predictor |
| `gender` | `gender` | Descriptive variable |
| `race` | `race` | Descriptive / ML feature |
| `menopause_status` | `menopause_status` | Descriptive / ML feature |
| `histological_type` | `histological_type` | Descriptive / ML feature |
| `ajcc_pathologic_tumor_stage` | `stage_group` | Main clinical predictor (grouped I–IV) |
| `OS` | `os_event` | Survival event indicator (1 = death) |
| `OS.time` | `os_time_days`, `os_time_months` | Follow-up time |

Other TCGA-CDR endpoints (`DSS`, `DFI`, `PFI`) are not used in the current version but are candidates for future extensions.

---

## Key Results

### Cohort

| Metric | Value |
|---|---:|
| Patients | 1,083 |
| Overall survival events | 151 |
| Event rate | 13.9 % |
| Median follow-up | 28.3 months |
| Median age at diagnosis | 58 years |
| Patients with missing stage | 24 |

### Stage distribution

| Stage | Patients | Percent |
|---|---:|---:|
| Stage I | 182 | 16.8 % |
| Stage II | 612 | 56.5 % |
| Stage III | 245 | 22.6 % |
| Stage IV | 20 | 1.8 % |
| Missing | 24 | 2.2 % |

### Cox proportional hazards model

Model: `Surv(os_time_months, os_event) ~ age + stage_group` (Stage I as reference)

| Predictor | Hazard ratio | 95 % CI | p-value |
|---|---:|---|---|
| Age, per year | 1.04 | 1.02–1.05 | < 0.001 |
| Stage II vs Stage I | 1.80 | 1.04–3.12 | 0.034 |
| Stage III vs Stage I | 3.61 | 2.03–6.44 | < 0.001 |
| Stage IV vs Stage I | 12.1 | 5.96–24.5 | < 0.001 |

The global Schoenfeld residual test was statistically significant, suggesting the proportional hazards assumption may not hold perfectly for stage. Hazard ratios should be interpreted as average associations over follow-up.

### ML risk prediction (3-year mortality)

After excluding 576 patients censored before 3 years, the ML dataset contained 493 patients (67 deaths, 13.6 % event rate).

| Model | Test AUC | CV AUC (5-fold) |
|---|---:|---:|
| Logistic regression | 0.720 | 0.743 ± 0.032 |
| Random forest | 0.671 | 0.692 ± 0.036 |
| Gradient boosting | — | 0.587 ± 0.040 |

Logistic regression outperformed tree-based models, consistent with the small sample size and the near-linear dose–response relationship between stage and mortality.

---

## Repository Structure

```text
Cancer-Survival-TCGA-CDR/
├── app_1.py                  ← Streamlit dashboard
├── requirements.txt
│
├── data/
│   ├── raw/                  ← Place TCGA-CDR Excel file here (not tracked)
│   └── processed/            ← Generated outputs (brca_*, luad_*, multicancer_*)
│
├── R/                        ← Scripts 01–06 (BRCA), 07–12 (LUAD + comparison)
├── notebooks/                ← 01 data validation, 02 ML risk prediction (BRCA)
├── python/                   ← Utility helpers for the Streamlit app
├── figures/                  ← Plots (brca_*, luad_*, multicancer_*)
└── reports/                  ← Short project report and code sample
```

---

## Setup and Installation

### Prerequisites

- R 4.2 or later
- Python 3.10 or later
- Git

### R package dependencies

```r
install.packages(c(
  "tidyverse",
  "readxl",
  "janitor",
  "survival",
  "survminer",
  "broom",
  "here"
))
```

### Python dependencies

```bash
pip install -r requirements.txt
```

---

## How to Reproduce the Full Analysis

### Step 1 — Download the raw data

Download `TCGA-CDR-SupplementalTableS1.xlsx` from the NCI GDC PanCanAtlas resources and place it in `data/raw/`.

### Step 2 — Run the R pipeline (in order)

**v1 — BRCA analysis:**

```bash
Rscript R/01_inspect_excel_sheets.R
Rscript R/02_clean_tcga_cdr.R
Rscript R/03_descriptive_analysis.R
Rscript R/04_kaplan_meier.R
Rscript R/05_cox_regression.R
Rscript R/06_export_for_streamlit.R
```

**v3 — LUAD and multi-cancer comparison:**

```bash
Rscript R/07_clean_luad.R
Rscript R/08_descriptive_luad.R
Rscript R/09_kaplan_meier_luad.R
Rscript R/10_cox_regression_luad.R
Rscript R/11_multicancer_comparison.R
Rscript R/12_export_for_streamlit_v3.R
```

> **Note:** Each script reads the output of the preceding one. If you update an earlier script, rerun it and all subsequent scripts to keep processed files consistent.

### Step 3 — Run the Python notebooks (in order)

Open and run all cells in:

1. `notebooks/01_python_data_check.ipynb`
2. `notebooks/02_ml_risk_prediction_demo.ipynb`

These generate the ML outputs and figures saved to `data/processed/` and `figures/`.

### Step 4 — Launch the Streamlit app

```bash
streamlit run app_1.py
```

---

## Streamlit Dashboard

The app is organised into six tabs:

| Tab | Content |
|---|---|
| Overview | Analysis summary and preview of the cleaned BRCA dataset |
| Cohort description | Stage, age, and gender distributions with sidebar filters (BRCA) |
| Kaplan–Meier analysis | Survival curves by pathological stage (BRCA) |
| Cox regression | Hazard ratio table, forest plot, and PH assumption test (BRCA) |
| ML risk prediction | 3-year mortality prediction demo with model comparison, ROC curves, calibration, confusion matrix, and model interpretation (BRCA) |
| Methods and limitations | Full methodological description and limitations |
| Multi-cancer comparison | BRCA vs LUAD cohort table, overlaid KM curves, side-by-side Cox forest plot |

---

## Main Outputs

| Output | File |
|---|---|
| Cleaned BRCA cohort | `data/processed/brca_survival_clean.csv` |
| Cohort summary | `data/processed/brca_cohort_summary.csv` |
| Cox model results | `data/processed/brca_cox_results.csv` |
| PH assumption test | `data/processed/brca_cox_ph_assumption_test.csv` |
| ML model comparison | `data/processed/brca_ml_model_comparison.csv` |
| ML test predictions | `data/processed/brca_ml_test_predictions.csv` |
| Kaplan–Meier plot | `figures/brca_km_stage_plot.png` |
| Cox forest plot | `figures/brca_cox_forest_plot.png` |
| ROC curves | `figures/brca_ml_roc_curves.png` |
| Calibration curves | `figures/brca_ml_calibration_curves.png` |
| Confusion matrix | `figures/brca_ml_confusion_matrix_best_model.png` |
| RF feature importance | `figures/brca_ml_feature_importance.png` |
| Logistic odds ratios | `figures/brca_ml_logistic_odds_ratios.png` |
| LUAD cleaned cohort | `data/processed/luad_survival_clean.csv` |
| LUAD Cox model results | `data/processed/luad_cox_results.csv` |
| Multi-cancer cohort table | `data/processed/multicancer_cohort_comparison.csv` |
| Multi-cancer Cox results | `data/processed/multicancer_cox_results.csv` |
| LUAD Kaplan–Meier plot | `figures/luad_km_stage_plot.png` |
| LUAD Cox forest plot | `figures/luad_cox_forest_plot.png` |
| BRCA vs LUAD KM comparison | `figures/multicancer_km_comparison.png` |
| BRCA vs LUAD Cox comparison | `figures/multicancer_cox_comparison.png` |

---

## Limitations

- TCGA is not a population-based cancer registry; results may not generalise to broader populations.
- The analysis is observational and should not be interpreted causally.
- Follow-up time is limited for a substantial fraction of patients (median ~28 months).
- Missingness in race (7.7 %) and menopause status (7.7 %) is handled by complete-case analysis in the survival models and mode imputation in the ML pipeline.
- The Stage IV subgroup is small (n = 20), limiting the precision of Stage IV estimates.
- The proportional hazards assumption is not perfectly satisfied; stage hazard ratios reflect average associations over follow-up.
- ER/PR/HER2 receptor status and tumour grade — the strongest prognostic markers in breast cancer — are not available in the TCGA-CDR clinical table and are absent from all models.
- The ML section is a methodological demonstration. The models have not been externally validated and are not intended for clinical decision support.
- The ML dataset is small (~493 patients, 67 events) after excluding early-censored patients; results should be interpreted as indicative.

---

## Planned Extensions

| Version | Extension | Status |
|---|---|---|
| v1 | BRCA survival analysis, R pipeline, Kaplan–Meier, Cox, dashboard | **Complete** |
| v2 | Python notebooks, ML 3-year mortality prediction | **Complete** |
| v3 | Multi-cancer comparison (BRCA + LUAD) | **Complete** |
| v4 | Alternative endpoints (PFI, DSS) | Planned |
| v5 | Molecular / genomics extension (ER/PR/HER2, PAM50 subtype) | Planned |

---

## Deployment

The deployed Streamlit app uses the processed CSV files and figure outputs committed to the repository. The raw TCGA-CDR Excel file is excluded from version control (see `.gitignore`).
