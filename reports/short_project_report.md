# Short Project Report

## Background
Cancer survival varies by clinical characteristics such as age and tumor stage.
This project uses open TCGA clinical data to demonstrate survival analysis methods.

## Objective
To evaluate whether age and pathological stage are associated with overall
survival among TCGA-BRCA breast cancer patients.

## Methods
The TCGA-CDR clinical dataset was filtered to BRCA patients. Overall survival was
defined using OS and OS.time. Kaplan-Meier curves were estimated by pathological
stage. Cox proportional hazards regression was used to estimate associations
between age, stage, and mortality hazard.

## Results
The final cohort included 1,083 patients and 151 overall survival events. Median
follow-up was 28.3 months and median age was 58 years. The Kaplan-Meier analysis
showed lower survival among patients with higher stage disease. In the Cox model,
each additional year of age was associated with increased mortality hazard
(HR 1.04, 95% CI 1.02-1.05). Compared with Stage I, hazard ratios were 1.80 for
Stage II, 3.61 for Stage III, and 12.1 for Stage IV.

## Model diagnostics
The proportional hazards assumption test suggested possible non-proportionality
for stage group. Results should therefore be interpreted as average hazard ratios
over follow-up.

## Limitations
TCGA is not population-based, and follow-up may be limited. The analysis is
intended as a reproducible portfolio demonstration rather than a clinical
decision-support model.

## Conclusion
This project demonstrates reproducible data cleaning, survival analysis, Cox
regression, diagnostic checking, and interactive communication of cancer
epidemiology results.