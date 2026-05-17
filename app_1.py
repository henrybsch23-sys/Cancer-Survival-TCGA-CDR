from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

st.set_page_config(
    page_title="TCGA Survival Analysis Dashboard",
    page_icon="🧬",
    layout="wide",
)

DATA_DIR = Path("data/processed")
FIG_DIR  = Path("figures")

# ---------------------------------------------------------------------------
# Cancer type registry
# Each entry carries metadata and file stems for every asset.
# ---------------------------------------------------------------------------
CANCER_META = {
    "BRCA": {
        "label":     "TCGA-BRCA",
        "full_name": "Breast Cancer",
        "prefix":    "brca",
        "km_script":  "R/04_kaplan_meier.R",
        "cox_script": "R/05_cox_regression.R",
    },
    "LUAD": {
        "label":     "TCGA-LUAD",
        "full_name": "Lung Adenocarcinoma",
        "prefix":    "luad",
        "km_script":  "R/09_kaplan_meier_luad.R",
        "cox_script": "R/10_cox_regression_luad.R",
    },
}


@st.cache_data
def load_cancer(prefix: str) -> dict | None:
    """Load all processed files for a single cancer type. Returns None if missing."""
    cohort_path = DATA_DIR / f"{prefix}_survival_clean.csv"
    if not cohort_path.exists():
        return None
    return {
        "cohort":         pd.read_csv(cohort_path),
        "cohort_summary": pd.read_csv(DATA_DIR / f"{prefix}_cohort_summary.csv"),
        "stage_summary":  pd.read_csv(DATA_DIR / f"{prefix}_stage_summary.csv"),
        "gender_summary": pd.read_csv(DATA_DIR / f"{prefix}_gender_summary.csv"),
        "age_summary":    pd.read_csv(DATA_DIR / f"{prefix}_age_summary.csv"),
        "cox_results":    pd.read_csv(DATA_DIR / f"{prefix}_cox_results.csv"),
        "ph_test":        pd.read_csv(DATA_DIR / f"{prefix}_cox_ph_assumption_test.csv"),
    }


@st.cache_data
def load_ml() -> dict:
    """Load BRCA ML outputs (all optional)."""
    out: dict = {}
    for key, fname in [
        ("status",        "brca_ml_three_year_status_summary.csv"),
        ("results",       "brca_ml_model_comparison.csv"),
        ("importance",    "brca_ml_feature_importance.csv"),
        ("coefficients",  "brca_ml_logistic_coefficients.csv"),
    ]:
        p = DATA_DIR / fname
        out[key] = pd.read_csv(p) if p.exists() else None
    return out


@st.cache_data
def load_multicancer() -> dict:
    """Load multi-cancer comparison outputs (optional)."""
    out: dict = {}
    for key, fname in [
        ("cohort", "multicancer_cohort_comparison.csv"),
        ("cox",    "multicancer_cox_results.csv"),
    ]:
        p = DATA_DIR / fname
        out[key] = pd.read_csv(p) if p.exists() else None
    return out


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def format_ml_status(ml_status: pd.DataFrame) -> pd.DataFrame:
    label_map = {
        "death_within_3y":           "Death within 3 years",
        "known_no_death_within_3y":  "Alive at 3 years (known)",
        "unknown_censored_before_3y": "Excluded — censored before 3 years",
    }
    d = ml_status.copy()
    d["three_year_status"] = d["three_year_status"].map(label_map).fillna(d["three_year_status"])
    d = d.rename(columns={"three_year_status": "3-year outcome", "n": "Patients (n)", "percent": "Percent (%)"})
    order = ["Death within 3 years", "Alive at 3 years (known)", "Excluded — censored before 3 years"]
    d["_sort"] = d["3-year outcome"].map({v: i for i, v in enumerate(order)}).fillna(99)
    return d.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)


def format_cox_results(cox: pd.DataFrame) -> pd.DataFrame:
    term_labels = {
        "age":                  "Age, per year",
        "stage_groupStage II":  "Stage II vs Stage I",
        "stage_groupStage III": "Stage III vs Stage I",
        "stage_groupStage IV":  "Stage IV vs Stage I",
    }
    d = cox.copy()
    d["Variable"]     = d["term"].map(term_labels).fillna(d["term"])
    d["Hazard ratio"] = d["estimate"].round(2)
    d["95% CI"]       = d["conf.low"].round(2).astype(str) + "–" + d["conf.high"].round(2).astype(str)
    d["p-value"]      = d["p.value"].apply(lambda x: "<0.001" if x < 0.001 else f"{x:.3f}")
    return d[["Variable", "Hazard ratio", "95% CI", "p-value"]]


def format_mc_cox(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["Hazard ratio"] = d["estimate"].round(2)
    d["95% CI"]       = d["conf.low"].round(2).astype(str) + "–" + d["conf.high"].round(2).astype(str)
    d["p-value"]      = d["p.value"].apply(lambda x: "<0.001" if x < 0.001 else f"{x:.3f}")
    return d[["cancer_type", "term_label", "Hazard ratio", "95% CI", "p-value"]].rename(
        columns={"cancer_type": "Cancer type", "term_label": "Predictor"}
    )


def show_image(path: Path, caption: str, script: str = "") -> None:
    if path.exists():
        st.image(str(path), caption=caption, use_container_width=True)
    else:
        msg = f"Figure not found: `{path.name}`."
        if script:
            msg += f" Run `{script}` to generate it."
        st.warning(msg)


def section_divider():
    st.markdown("---")


# ---------------------------------------------------------------------------
# Load all data
# ---------------------------------------------------------------------------
cancer_data: dict[str, dict] = {}
for code, meta in CANCER_META.items():
    result = load_cancer(meta["prefix"])
    if result is not None:
        cancer_data[code] = result

ml   = load_ml()
mc   = load_multicancer()

available_cancers = list(cancer_data.keys())

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.title("TCGA Survival App")

    section_divider()

    if not available_cancers:
        st.error("No processed data found. Run the R pipeline first.")
        st.stop()

    cancer_labels = {k: f"{CANCER_META[k]['label']} — {CANCER_META[k]['full_name']}"
                     for k in available_cancers}

    selected_cancer = st.radio(
        "Cancer type",
        options=available_cancers,
        format_func=lambda k: cancer_labels[k],
    )

    section_divider()

    meta    = CANCER_META[selected_cancer]
    data    = cancer_data[selected_cancer]
    cohort  = data["cohort"]
    prefix  = meta["prefix"]

    st.markdown(
        f"""
        **Dataset:** TCGA-CDR  
        **Selected:** {meta['label']}  
        **Endpoint:** Overall survival  
        **Model:** Cox regression (age + stage)
        """
    )

    if selected_cancer == "BRCA":
        st.markdown("**ML extension:** 3-year mortality prediction")

    section_divider()
    st.caption(
        "Use the radio button above to switch between cancer types. "
        "All tabs update automatically."
    )

# ---------------------------------------------------------------------------
# Page header — dynamic
# ---------------------------------------------------------------------------
st.title(f"TCGA Survival Analysis — {meta['label']} ({meta['full_name']})")

st.markdown(
    f"""
    Reproducible survival analysis using the open-access **TCGA Pan-Cancer Clinical
    Data Resource**. Currently showing **{meta['label']} ({meta['full_name']})**.
    Switch cancer type in the sidebar to explore a different cohort.
    The workflow covers data cleaning, Kaplan-Meier estimation, Cox proportional
    hazards regression, and model diagnostics.
    {"An ML extension for 3-year mortality prediction is available for BRCA." if selected_cancer == "BRCA" else ""}
    """
)

# Top metrics
n_patients   = len(cohort)
n_events     = int(cohort["os_event"].sum())
event_rate   = 100 * cohort["os_event"].mean()
median_fu    = round(cohort["os_time_months"].median(), 1)
median_age   = round(cohort["age"].median(), 1)
n_stage_miss = int(cohort["stage_group"].isna().sum())

m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Patients",         f"{n_patients:,}")
m2.metric("OS events",        f"{n_events:,}", f"{event_rate:.1f}%")
m3.metric("Median follow-up", f"{median_fu} months")
m4.metric("Median age",       f"{median_age} years")
m5.metric("Missing stage",    f"{n_stage_miss}")

section_divider()

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------
tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
    "Overview",
    "Cohort description",
    "Kaplan-Meier analysis",
    "Cox regression",
    "ML risk prediction",
    "Methods and limitations",
    "Multi-cancer comparison",
])

# ── Tab 1: Overview ──────────────────────────────────────────────────────────
with tab1:
    st.header(f"Project overview — {meta['label']}")

    left, right = st.columns([1.3, 1])

    with left:
        st.subheader("Analysis summary")
        st.markdown(
            f"""
            The analysis evaluates whether **age at diagnosis** and
            **pathological tumour stage** are associated with overall survival
            among **{meta['label']} ({meta['full_name']})** patients.

            - **R**: data cleaning, descriptive analysis, Kaplan-Meier curves,
              Cox regression, and model diagnostics.
            - **Python / Streamlit**: interactive dashboard and results communication.
            """
        )
        st.dataframe(data["cohort_summary"], use_container_width=True, hide_index=True)

    with right:
        st.subheader("Analysis-ready dataset")
        st.dataframe(cohort.head(15), use_container_width=True, hide_index=True)
        st.download_button(
            label="Download cleaned cohort CSV",
            data=cohort.to_csv(index=False).encode("utf-8"),
            file_name=f"{prefix}_survival_clean.csv",
            mime="text/csv",
        )

# ── Tab 2: Cohort description ─────────────────────────────────────────────────
with tab2:
    st.header(f"Cohort description — {meta['label']}")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Stage distribution")
        stage_plot_df = (
            cohort["stage_group"]
            .fillna("Missing")
            .value_counts()
            .rename_axis("stage_group")
            .reset_index(name="n")
        )
        fig_stage = px.bar(
            stage_plot_df, x="stage_group", y="n", text="n",
            labels={"stage_group": "Pathological stage", "n": "Patients"},
            title=f"Patients by pathological stage ({meta['label']})",
        )
        fig_stage.update_layout(title_x=0.02, showlegend=False, height=420)
        st.plotly_chart(fig_stage, use_container_width=True)

    with col2:
        st.subheader("Age distribution")
        fig_age = px.histogram(
            cohort, x="age", nbins=25,
            labels={"age": "Age at diagnosis"},
            title=f"Age at diagnosis ({meta['label']})",
        )
        fig_age.update_layout(title_x=0.02, showlegend=False, height=420)
        st.plotly_chart(fig_age, use_container_width=True)

    col3, col4 = st.columns(2)
    with col3:
        st.subheader("Age group summary")
        st.dataframe(data["age_summary"], use_container_width=True, hide_index=True)
    with col4:
        st.subheader("Gender summary")
        st.dataframe(data["gender_summary"], use_container_width=True, hide_index=True)

    st.subheader("Stage summary")
    st.dataframe(data["stage_summary"], use_container_width=True, hide_index=True)

# ── Tab 3: Kaplan-Meier ───────────────────────────────────────────────────────
with tab3:
    st.header(f"Kaplan-Meier survival analysis — {meta['label']}")

    st.markdown(
        """
        Kaplan-Meier curves estimate the probability of remaining alive over time.
        Patients who did not experience death before their last recorded follow-up
        are treated as censored observations.
        """
    )

    show_image(
        FIG_DIR / f"{prefix}_km_stage_plot.png",
        f"Overall survival by pathological stage — {meta['label']}",
        script=meta["km_script"],
    )

    st.info(
        "The Kaplan-Meier curves show lower overall survival among patients with "
        "more advanced pathological stage. The log-rank test evaluates whether "
        "survival distributions differ significantly across stage groups."
    )

# ── Tab 4: Cox regression ─────────────────────────────────────────────────────
with tab4:
    st.header(f"Cox proportional hazards model — {meta['label']}")

    st.markdown(
        f"""
        The Cox model estimates hazard ratios for overall mortality in
        **{meta['label']}** patients. The model includes **age at diagnosis**
        and **pathological stage**, with Stage I as the reference group.
        """
    )

    st.subheader("Hazard ratio table")
    st.dataframe(format_cox_results(data["cox_results"]), use_container_width=True, hide_index=True)

    show_image(
        FIG_DIR / f"{prefix}_cox_forest_plot.png",
        f"Cox model hazard ratios — {meta['label']}",
        script=meta["cox_script"],
    )

    st.subheader("Interpretation")
    st.markdown(
        """
        - A hazard ratio above 1 indicates higher mortality hazard.
        - Each additional year of age is associated with higher estimated hazard.
        - Compared with Stage I, more advanced stages show progressively higher
          mortality hazard.
        - Stage IV shows the highest hazard ratio, but estimates for small subgroups
          should be interpreted with caution.
        """
    )

    st.subheader("Proportional hazards assumption test")
    ph = data["ph_test"]
    st.dataframe(ph, use_container_width=True, hide_index=True)

    global_ph = ph.loc[ph["term"] == "GLOBAL", "p_value"]
    if not global_ph.empty and float(global_ph.iloc[0]) < 0.05:
        st.warning(
            "The global proportional hazards test is statistically significant. "
            "Results should be interpreted as average hazard ratios over follow-up."
        )
    else:
        st.success("No major global proportional hazards issue detected.")

# ── Tab 5: ML risk prediction ─────────────────────────────────────────────────
with tab5:
    st.header("Machine learning risk prediction — TCGA-BRCA")

    if selected_cancer != "BRCA":
        st.info(
            f"The ML risk prediction extension was developed for **TCGA-BRCA** only. "
            f"Switch to BRCA in the sidebar to explore the 3-year mortality prediction demo. "
            f"A LUAD ML extension is planned for a future version."
        )
    else:
        st.markdown(
            """
            This section demonstrates a Python-based ML extension to the survival analysis
            workflow. The prediction task is **death within 3 years after diagnosis**.

            Patients censored before 3 years are excluded because their 3-year outcome
            status is unknown. This avoids incorrectly treating early-censored patients
            as survivors.
            """
        )

        if ml["status"] is None or ml["results"] is None:
            st.warning("ML outputs not found. Run `notebooks/02_ml_risk_prediction_demo.ipynb` to regenerate.")
        else:
            st.subheader("3-year outcome definition")
            st.dataframe(format_ml_status(ml["status"]), use_container_width=True, hide_index=True)

            st.info(
                "Only patients with a known 3-year outcome are used for model training "
                "and evaluation. Early-censored patients are excluded to avoid outcome "
                "misclassification bias."
            )

            st.subheader("Model comparison")
            model_display = ml["results"].copy()
            for col in ["auc","accuracy","balanced_accuracy","sensitivity_recall",
                        "specificity","precision","f1","cv_auc_mean","cv_auc_sd"]:
                if col in model_display.columns:
                    model_display[col] = model_display[col].round(3)

            st.dataframe(model_display, use_container_width=True, hide_index=True)

            selected_model = st.selectbox("Select a model to inspect", model_display["model"].tolist())
            row = model_display.loc[model_display["model"] == selected_model].iloc[0]

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("AUC",         f"{row['auc']:.3f}")
            c2.metric("CV AUC",      f"{row['cv_auc_mean']:.3f}")
            c3.metric("Sensitivity", f"{row['sensitivity_recall']:.3f}")
            c4.metric("Specificity", f"{row['specificity']:.3f}")

            col_left, col_right = st.columns(2)
            with col_left:
                show_image(FIG_DIR / "brca_ml_roc_curves.png", "ROC curves — 3-year mortality")
            with col_right:
                show_image(FIG_DIR / "brca_ml_auc_comparison.png", "AUC comparison across models")

            show_image(FIG_DIR / "brca_ml_calibration_curves.png", "Calibration curves for all models")

            st.subheader("Best model confusion matrix")
            show_image(FIG_DIR / "brca_ml_confusion_matrix_best_model.png", "Confusion matrix — best model")

            st.subheader("Model interpretation")
            col_a, col_b = st.columns(2)
            with col_a:
                show_image(FIG_DIR / "brca_ml_logistic_odds_ratios.png", "Logistic regression odds ratios")
            with col_b:
                show_image(FIG_DIR / "brca_ml_feature_importance.png", "Random forest feature importance")

            if ml["importance"] is not None:
                with st.expander("View random forest feature-importance table"):
                    st.dataframe(ml["importance"], use_container_width=True, hide_index=True)

            if ml["coefficients"] is not None:
                with st.expander("View logistic regression coefficients table"):
                    st.dataframe(ml["coefficients"], use_container_width=True, hide_index=True)

            st.warning(
                "This ML analysis is a methodological demonstration and is not intended "
                "for clinical decision support. Models have not been externally validated."
            )

# ── Tab 6: Methods ────────────────────────────────────────────────────────────
with tab6:
    st.header("Methods and limitations")

    st.subheader("Dataset")
    st.markdown(
        """
        This project uses the open-access **TCGA Pan-Cancer Clinical Data Resource
        (TCGA-CDR)**. The analysis covers **TCGA-BRCA** (breast cancer) and
        **TCGA-LUAD** (lung adenocarcinoma). Only curated clinical and survival
        outcome data are used — no controlled-access genomic data.
        """
    )

    st.subheader("Endpoint")
    st.markdown(
        """
        Overall survival was defined using:

        - `OS`: event indicator — 1 = death, 0 = censored (alive at last follow-up).
        - `OS.time`: follow-up time in days, converted to months for analysis.
        """
    )

    st.subheader("Statistical methods")
    st.markdown(
        """
        - Descriptive cohort summaries
        - Kaplan-Meier survival curves with log-rank tests
        - Cox proportional hazards regression (age + pathological stage)
        - Proportional hazards assumption testing via Schoenfeld residuals
        - Multi-cancer survival comparison (BRCA vs LUAD)
        - Censoring-aware 3-year mortality prediction (BRCA only)
        - Logistic regression, random forest, and gradient boosting
        - Train/test evaluation and stratified cross-validation
        - ROC/AUC, sensitivity, specificity, precision, and F1-score
        """
    )

    st.subheader("Limitations")
    st.markdown(
        """
        - TCGA is not a population-based cancer registry; results may not generalise.
        - The analysis is observational and should not be interpreted causally.
        - Follow-up is limited for a substantial fraction of patients.
        - Missingness varies across clinical variables and is handled by complete-case
          analysis in survival models and mode imputation in the ML pipeline.
        - The Stage IV subgroup is small in both cohorts.
        - The proportional hazards assumption may not hold perfectly for stage.
        - ER/PR/HER2 and tumour grade are unavailable in the TCGA-CDR clinical table.
        - The ML extension is a methodological demonstration and has not been externally
          validated; it is not intended for clinical decision support.
        """
    )

# ── Tab 7: Multi-cancer comparison ───────────────────────────────────────────
with tab7:
    st.header("Multi-cancer comparison: BRCA vs LUAD")

    st.markdown(
        """
        Both cohorts use the same endpoint (overall survival), the same predictor set
        (age + pathological stage), and the same Cox model specification, enabling
        direct comparison of effect estimates across cancer types.
        """
    )

    if mc["cohort"] is None or mc["cox"] is None:
        st.warning(
            "Multi-cancer comparison outputs not found. "
            "Run `R/07` through `R/11` to generate LUAD and comparison files."
        )
    else:
        section_divider()
        st.subheader("Cohort comparison")

        st.dataframe(
            mc["cohort"].rename(columns={
                "cancer_type":             "Cancer type",
                "n_patients":              "Patients (n)",
                "n_os_events":             "OS events (n)",
                "os_event_rate_pct":       "Event rate (%)",
                "median_follow_up_months": "Median follow-up (months)",
                "median_age":              "Median age",
                "n_missing_stage":         "Missing stage (n)",
            }),
            use_container_width=True,
            hide_index=True,
        )

        st.info(
            "LUAD typically has a substantially higher event rate and shorter "
            "median follow-up than BRCA, reflecting the more aggressive natural "
            "history of lung adenocarcinoma."
        )

        section_divider()
        st.subheader("Kaplan-Meier survival comparison")
        st.markdown(
            "Overall survival curves for the full cohort of each cancer type "
            "(not stratified by stage). Shaded bands are 95% confidence intervals."
        )
        show_image(
            FIG_DIR / "multicancer_km_comparison.png",
            "Overall survival: TCGA-BRCA vs TCGA-LUAD",
            script="R/11_multicancer_comparison.R",
        )

        section_divider()
        st.subheader("Cox model hazard ratio comparison")
        st.markdown(
            "Both models use `Surv(os_time_months, os_event) ~ age + stage_group` "
            "with Stage I as the reference. Hazard ratios are overlaid for direct comparison."
        )
        show_image(
            FIG_DIR / "multicancer_cox_comparison.png",
            "Cox hazard ratios: BRCA vs LUAD",
            script="R/11_multicancer_comparison.R",
        )

        st.subheader("Hazard ratio table")
        st.dataframe(format_mc_cox(mc["cox"]), use_container_width=True, hide_index=True)

        section_divider()
        st.subheader("Interpretation")
        st.markdown(
            """
            - **Stage gradient:** Both cancers show a clear dose–response relationship
              between advancing stage and higher mortality hazard.
            - **Magnitude:** Stage III and IV hazard ratios tend to be larger in LUAD,
              reflecting the steeper survival gradient in lung cancer.
            - **Age effect:** The per-year age hazard ratio is broadly similar across
              both cancers.
            - **Baseline survival:** Even within the same stage, LUAD patients have
              substantially lower baseline survival than BRCA patients.
            - **Caveat:** Both models are fitted independently. Differences in HRs do
              not directly imply causal differences — they may reflect variation in
              treatment, follow-up patterns, or unmeasured confounders.
            """
        )
