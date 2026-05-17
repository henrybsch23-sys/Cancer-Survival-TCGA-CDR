from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"


def load_clean_cohort() -> pd.DataFrame:
    """Load the cleaned TCGA-BRCA survival cohort."""
    return pd.read_csv(DATA_PROCESSED / "brca_survival_clean.csv")


def load_cox_results() -> pd.DataFrame:
    """Load Cox proportional hazards model results."""
    return pd.read_csv(DATA_PROCESSED / "brca_cox_results.csv")