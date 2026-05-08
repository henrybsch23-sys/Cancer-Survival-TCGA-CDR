# TCGA Breast Cancer Survival Dashboard

## Live dashboard

A deployed version of the Streamlit dashboard is available here:

https://cancer-survival-tcga-cdr-gzxwsq2dbktxmjtlet76rz.streamlit.app/

The app presents cohort summaries, Kaplan–Meier survival curves, Cox proportional hazards model results, proportional hazards diagnostics, and visual interpretation of survival patterns among TCGA-BRCA breast cancer patients.

## Overview

This project demonstrates a reproducible survival analysis workflow using the TCGA Pan-Cancer Clinical Data Resource. The current version focuses on breast cancer patients from TCGA-BRCA and includes data cleaning, descriptive epidemiology, Kaplan–Meier survival estimation, Cox proportional hazards regression, proportional hazards diagnostics, and a Streamlit dashboard.

The statistical analysis is implemented in R, while the dashboard is implemented in Python using Streamlit.

## Motivation

This project was developed as a portfolio project for research roles in cancer epidemiology, biostatistics, precision medicine, and data-driven health research.

The aim is to demonstrate practical skills in:

- Clinical data cleaning
- Survival analysis
- Cox regression
- Statistical interpretation
- Reproducible research workflows
- Interactive communication of epidemiological results

## Dataset

The project uses the open-access TCGA Pan-Cancer Clinical Data Resource:

```text
TCGA-CDR-SupplementalTableS1.xlsx
```

The file should be downloaded from the NCI Genomic Data Commons PanCanAtlas publication resources and placed in:

```text
data/raw/TCGA-CDR-SupplementalTableS1.xlsx
```

Only open-access clinical and survival outcome data are used. No controlled-access genomic data, raw sequencing data, or identifiable patient data are used.

## Research Question

Among TCGA breast cancer patients, are age at diagnosis and pathological tumor stage associated with overall survival?

## Endpoint

The main endpoint is **overall survival**.

| Variable | Meaning |
|---|---|
| `OS` | Overall survival event indicator; 1 = death, 0 = censored/alive at last follow-up |
| `OS.time` | Follow-up time in days |
| `os_time_months` | Follow-up time converted to months |

Patients who did not die during recorded follow-up are treated as censored observations.

## Methods

The analysis includes:

1. Inspection of the TCGA-CDR Excel workbook
2. Filtering to TCGA-BRCA breast cancer patients
3. Cleaning age, stage, survival time, and survival event variables
4. Descriptive cohort analysis
5. Kaplan-Meier survival estimation by pathological stage
6. Log-rank comparison of survival curves
7. Cox proportional hazards regression
8. Proportional hazards assumption testing
9. Export of analysis-ready tables and figures
10. Visualization in a Streamlit dashboard

## Current Analysis Variables

| Original TCGA-CDR column | Cleaned variable | Use |
|---|---|---|
| `bcr_patient_barcode` | `patient_id` | Patient identifier |
| `type` | `cancer_type` | Cancer type filter |
| `age_at_initial_pathologic_diagnosis` | `age` | Predictor |
| `gender` | `gender` | Descriptive variable |
| `race` | `race` | Stored for future descriptive work |
| `ajcc_pathologic_tumor_stage` | `stage_group` | Main clinical predictor |
| `OS` | `os_event` | Survival event indicator |
| `OS.time` | `os_time_days`, `os_time_months` | Follow-up time |

Other TCGA-CDR endpoints such as `DSS`, `DFI`, and `PFI` are not used in the current MVP version, but may be added in future extensions.

## Key Results

The cleaned TCGA-BRCA cohort included:

| Metric | Value |
|---|---:|
| Patients | 1,083 |
| Overall survival events | 151 |
| Overall survival event proportion | 13.9% |
| Median follow-up | 28.3 months |
| Median age at diagnosis | 58 years |
| Patients with missing stage | 24 |

### Stage Distribution

| Stage | Patients | Percent |
|---|---:|---:|
| Stage I | 182 | 16.8% |
| Stage II | 612 | 56.5% |
| Stage III | 245 | 22.6% |
| Stage IV | 20 | 1.8% |
| Missing | 24 | 2.2% |

### Cox Proportional Hazards Model

The primary Cox model was:

```text
Surv(os_time_months, os_event) ~ age + stage_group
```

Stage I was used as the reference group.

| Predictor | Hazard ratio | 95% CI | Interpretation |
|---|---:|---:|---|
| Age, per year | 1.04 | 1.02-1.05 | Each additional year of age was associated with higher mortality hazard |
| Stage II vs Stage I | 1.80 | 1.04-3.12 | Stage II had higher hazard than Stage I |
| Stage III vs Stage I | 3.61 | 2.03-6.44 | Stage III had substantially higher hazard than Stage I |
| Stage IV vs Stage I | 12.1 | 5.96-24.5 | Stage IV had much higher hazard than Stage I |

Kaplan-Meier analysis showed lower overall survival among patients with more advanced pathological stage.

## Model Diagnostics

The proportional hazards assumption was assessed using Schoenfeld residuals.

The global proportional hazards test was statistically significant, suggesting that the Cox proportional hazards assumption may not hold perfectly. The issue appeared mainly related to stage group.

Therefore, stage hazard ratios should be interpreted as average associations over follow-up. Future sensitivity analyses could consider stratified Cox models or time-varying effects.

## Repository Structure

```text
Cancer-Survival-TCGA-CDR/
│
├── app_1.py
├── README.md
├── requirements.txt
├── .gitignore
├── LICENSE
│
├── data/
│   ├── raw/
│   │   └── .gitkeep
│   └── processed/
│       ├── tcga_brca_survival_clean.csv
│       ├── cohort_summary.csv
│       ├── stage_summary.csv
│       ├── gender_summary.csv
│       ├── age_summary.csv
│       ├── cox_results.csv
│       └── cox_ph_assumption_test.csv
│
├── R/
│   ├── 01_inspect_excel_sheets.R
│   ├── 02_clean_tcga_cdr.R
│   ├── 03_descriptive_analysis.R
│   ├── 04_kaplan_meier.R
│   └── 05_cox_regression.R
│
├── notebooks/
│   ├── 01_python_data_check.ipynb
│   └── 02_ml_risk_prediction_demo.ipynb
│
├── python/
│   ├── prepare_app_data.py
│   └── utils.py
│
├── figures/
│   ├── km_stage_plot.png
│   └── cox_forest_plot.png
│
└── reports/
    └── short_project_report.md
```

## Setup and Installation

### Prerequisites

- R 4.5.1 or later
- Python 3.10 or later
- Git
- Streamlit

### R Package Dependencies

Install the required R packages:

```r
install.packages(c(
  "tidyverse",
  "readxl",
  "janitor",
  "survival",
  "survminer",
  "broom",
  "gtsummary",
  "here"
))
```

### Python Dependencies

Install Python dependencies from the project root:

```powershell
pip install -r requirements.txt
```

A minimal `requirements.txt` should include:

```text
streamlit
pandas
numpy
matplotlib
plotly
openpyxl
scikit-learn
jupyter
```

## How to run locally

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

Run the Streamlit dashboard:

```bash
streamlit run app_1.py
```

The app uses the processed CSV files in `data/processed/` and the figures in `figures/`.

## How to reproduce the R analysis

The raw TCGA-CDR Excel file is not included in this repository. To reproduce the full R analysis:

1. Download `TCGA-CDR-SupplementalTableS1.xlsx` from the NCI Genomic Data Commons PanCanAtlas publication resources.
2. Place the file in:

```text
data/raw/TCGA-CDR-SupplementalTableS1.xlsx
```

3. Run the R scripts in order:

```bash
Rscript R/01_inspect_excel_sheets.R
Rscript R/02_clean_tcga_cdr.R
Rscript R/03_descriptive_analysis.R
Rscript R/04_kaplan_meier.R
Rscript R/05_cox_regression.R
```

These scripts generate the cleaned dataset, summary tables, Kaplan–Meier plot, Cox model results, proportional hazards diagnostics, and forest plot.

## Deployment note

The deployed Streamlit app uses the processed CSV files and figure outputs committed to the repository. The raw TCGA-CDR Excel file is not required for the deployed app and is intentionally excluded from version control.

## Main outputs

| Output | File |
|---|---|
| Cleaned BRCA cohort | `data/processed/tcga_brca_survival_clean.csv` |
| Cohort summary | `data/processed/cohort_summary.csv` |
| Stage summary | `data/processed/stage_summary.csv` |
| Age summary | `data/processed/age_summary.csv` |
| Gender summary | `data/processed/gender_summary.csv` |
| Cox model results | `data/processed/cox_results.csv` |
| Proportional hazards diagnostics | `data/processed/cox_ph_assumption_test.csv` |
| Kaplan–Meier plot | `figures/km_stage_plot.png` |
| Cox forest plot | `figures/cox_forest_plot.png` |
| Streamlit dashboard | `app_1.py` |

## Streamlit Dashboard

The Streamlit dashboard includes:

- Cohort overview
- Descriptive summaries
- Kaplan-Meier survival plot by pathological stage
- Cox regression table
- Cox forest plot
- Proportional hazards assumption results
- Methods and limitations

When running locally, the dashboard is available only while the Streamlit process is active. If the terminal or Cursor session is closed, the local dashboard stops. The app can later be deployed through Streamlit Community Cloud.

## Limitations

This project has several important limitations:

- TCGA is not a population-based cancer registry.
- The analysis is observational and should not be interpreted causally.
- Follow-up time is relatively limited for some patients.
- Missingness varies across clinical variables.
- Stage IV includes a small number of BRCA patients in this dataset.
- The proportional hazards assumption may not hold perfectly for stage group.
- Results are intended for educational and portfolio purposes, not clinical decision support.

## Planned Extensions

Potential future extensions include:

1. **Python data validation notebook**
   - Check missingness
   - Validate variable coding
   - Summarize the cleaned BRCA cohort using pandas

2. **Machine learning demo**
   - Predict 3-year mortality status using age and stage
   - Exclude patients censored before 3 years
   - Compare simple logistic regression or random forest performance

3. **Cancer-type comparison**
   - Compare BRCA with LUAD/LUSC or another cancer type
   - Evaluate whether stage-survival patterns differ across cancers

4. **Alternative survival endpoints**
   - Repeat the analysis using PFI or DSS
   - Compare results with OS-based analysis

5. **Molecular extension**
   - Add open-access molecular features later, such as mutation or expression-based summaries
   - This is intentionally left for a later version because it increases project complexity

## Application Relevance

This project demonstrates skills in:

- Epidemiological data management
- R-based statistical analysis
- Survival analysis
- Cox proportional hazards regression
- Model diagnostics
- Python-based dashboard development
- Reproducible research workflows
- Scientific interpretation of clinical data

The project is especially relevant to applications in cancer epidemiology, biostatistics, precision medicine, clinical epidemiology, and data-driven health research.
