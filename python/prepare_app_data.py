from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data" / "processed"

cohort = pd.read_csv(DATA_DIR / "tcga_brca_survival_clean.csv")

summary = {
    "patients": len(cohort),
    "os_events": int(cohort["os_event"].sum()),
    "median_follow_up_months": round(cohort["os_time_months"].median(), 1),
    "median_age": round(cohort["age"].median(), 1),
}

summary_df = pd.DataFrame(
    [{"metric": key, "value": value} for key, value in summary.items()]
)

summary_df.to_csv(DATA_DIR / "python_app_summary.csv", index=False)

print("Saved python_app_summary.csv")