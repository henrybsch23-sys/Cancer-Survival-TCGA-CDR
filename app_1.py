from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

st.set_page_config(
    page_title="TCGA Breast Cancer Survival Dashboard",
    page_icon="🧬",
    layout="wide",
)

DATA_DIR = Path("data/processed")
FIG_DIR = Path("figures")


@st.cache_data
def load_data():
    cohort = pd.read_csv(DATA_DIR / "brca_survival_clean.csv")
    cohort_summary = pd.read_csv(DATA_DIR / "brca_cohort_summary.csv")
    stage_summary = pd.read_csv(DATA_DIR / "brca_stage_summary.csv")
    gender_summary = pd.read_csv(DATA_DIR / "brca_gender_summary.csv")
    age_summary = pd.read_csv(DATA_DIR / "brca_age_summary.csv")
    cox_results = pd.read_csv(DATA_DIR / "brca_cox_results.csv")
    ph_test = pd.read_csv(DATA_DIR / "brca_cox_ph_assumption_test.csv")

    ml_status = None
    ml_results = None
    ml_feature_importance = None
    ml_logistic_coefficients = None

    if (DATA_DIR / "brca_ml_three_year_status_summary.csv").exists():
        ml_status = pd.read_csv(DATA_DIR / "brca_ml_three_year_status_summary.csv")

    if (DATA_DIR / "brca_ml_model_comparison.csv").exists():
        ml_results = pd.read_csv(DATA_DIR / "brca_ml_model_comparison.csv")

    if (DATA_DIR / "brca_ml_feature_importance.csv").exists():
        ml_feature_importance = pd.read_csv(DATA_DIR / "brca_ml_feature_importance.csv")

    if (DATA_DIR / "brca_ml_logistic_coefficients.csv").exists():
        ml_logistic_coefficients = pd.read_csv(DATA_DIR / "brca_ml_logistic_coefficients.csv")

    # v3 — multi-cancer comparison data
    mc_cohort = None
    mc_cox = None

    if (DATA_DIR / "multicancer_cohort_comparison.csv").exists():
        mc_cohort = pd.read_csv(DATA_DIR / "multicancer_cohort_comparison.csv")

    if (DATA_DIR / "multicancer_cox_results.csv").exists():
        mc_cox = pd.read_csv(DATA_DIR / "multicancer_cox_results.csv")

    return (
        cohort,
        cohort_summary,
        stage_summary,
        gender_summary,
        age_summary,
        cox_results,
        ph_test,
        ml_status,
        ml_results,
        ml_feature_importance,
        ml_logistic_coefficients,
        mc_cohort,
        mc_cox,
    )


def format_ml_status(ml_status: pd.DataFrame) -> pd.DataFrame:
    label_map = {
        "death_within_3y": "Death within 3 years",
        "known_no_death_within_3y": "Alive at 3 years (known)",
        "unknown_censored_before_3y": "Excluded — censored before 3 years",
    }
    display = ml_status.copy()
    display["three_year_status"] = display["three_year_status"].map(label_map).fillna(
        display["three_year_status"]
    )
    display = display.rename(
        columns={
            "three_year_status": "3-year outcome",
            "n": "Patients (n)",
            "percent": "Percent (%)",
        }
    )
    order = ["Death within 3 years", "Alive at 3 years (known)", "Excluded — censored before 3 years"]
    display["_sort"] = display["3-year outcome"].map({v: i for i, v in enumerate(order)}).fillna(99)
    display = display.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
    return display


def format_cox_results(cox_results: pd.DataFrame) -> pd.DataFrame:
    display = cox_results.copy()

    term_labels = {
        "age": "Age, per year",
        "stage_groupStage II": "Stage II vs Stage I",
        "stage_groupStage III": "Stage III vs Stage I",
        "stage_groupStage IV": "Stage IV vs Stage I",
    }

    display["Variable"] = display["term"].map(term_labels).fillna(display["term"])
    display["Hazard ratio"] = display["estimate"].round(2)
    display["95% CI"] = (
        display["conf.low"].round(2).astype(str)
        + "–"
        + display["conf.high"].round(2).astype(str)
    )
    display["p-value"] = display["p.value"].apply(
        lambda x: "<0.001" if x < 0.001 else f"{x:.3f}"
    )

    return display[["Variable", "Hazard ratio", "95% CI", "p-value"]]


def section_divider():
    st.markdown("---")


(
    cohort,
    cohort_summary,
    stage_summary,
    gender_summary,
    age_summary,
    cox_results,
    ph_test,
    ml_status,
    ml_results,
    ml_feature_importance,
    ml_logistic_coefficients,
    mc_cohort,
    mc_cox,
) = load_data()

# Sidebar
with st.sidebar:
    st.title("TCGA Survival App")
    st.markdown(
        """
        **Dataset:** TCGA-CDR  
        **Cancer types:** TCGA-BRCA · TCGA-LUAD  
        **Endpoint:** Overall survival  
        **Main model:** Cox regression  
        **ML extension:** 3-year mortality prediction (BRCA)
        """
    )

    section_divider()

    st.subheader("Cohort filters")

    available_stages = sorted(cohort["stage_group"].dropna().unique())
    selected_stages = st.multiselect(
        "Pathological stage",
        options=available_stages,
        default=available_stages,
    )

    min_age = int(cohort["age"].min())
    max_age = int(cohort["age"].max())

    age_range = st.slider(
        "Age at diagnosis",
        min_value=min_age,
        max_value=max_age,
        value=(min_age, max_age),
    )

    section_divider()

    st.caption(
        "Filters affect descriptive dashboard views only. "
        "The Cox model shown was fitted using the full cleaned analysis cohort."
    )

filtered_cohort = cohort[
    cohort["stage_group"].isin(selected_stages)
    & cohort["age"].between(age_range[0], age_range[1])
].copy()

# Header
st.title("TCGA Breast Cancer Survival Analysis Dashboard")

st.markdown(
    """
    This dashboard presents a reproducible survival analysis using the open-access
    **TCGA Pan-Cancer Clinical Data Resource**. The current version focuses on
    **TCGA-BRCA breast cancer patients** and demonstrates data cleaning,
    Kaplan-Meier survival estimation, Cox proportional hazards regression,
    model diagnostics, and a machine learning extension for short-term mortality
    risk prediction.
    """
)

# Top metrics
n_patients = len(cohort)
n_filtered = len(filtered_cohort)
n_events = int(cohort["os_event"].sum())
event_rate = 100 * cohort["os_event"].mean()
median_follow_up = round(cohort["os_time_months"].median(), 1)
median_age = round(cohort["age"].median(), 1)

metric_col1, metric_col2, metric_col3, metric_col4, metric_col5 = st.columns(5)

metric_col1.metric("Patients", f"{n_patients:,}")
metric_col2.metric("Filtered view", f"{n_filtered:,}")
metric_col3.metric("OS events", f"{n_events:,}", f"{event_rate:.1f}%")
metric_col4.metric("Median follow-up", f"{median_follow_up} months")
metric_col5.metric("Median age", f"{median_age} years")

section_divider()

tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs(
    [
        "Overview",
        "Cohort description",
        "Kaplan-Meier analysis",
        "Cox regression",
        "ML risk prediction",
        "Methods and limitations",
        "Multi-cancer comparison",
    ]
)

with tab1:
    st.header("Project overview")

    left, right = st.columns([1.3, 1])

    with left:
        st.subheader("Analysis summary")

        st.markdown(
            """
            The analysis evaluates whether **age at diagnosis** and
            **pathological tumor stage** are associated with overall survival
            among TCGA-BRCA breast cancer patients.

            The workflow is intentionally split across R and Python:

            - **R**: data cleaning, descriptive analysis, Kaplan-Meier curves,
              Cox regression, and model diagnostics.
            - **Python/Streamlit**: interactive dashboard and communication of results.
            """
        )

        st.info(
            "Main finding: older age and more advanced pathological stage were "
            "associated with higher mortality hazard in the Cox model."
        )

    with right:
        st.subheader("Analysis-ready dataset")

        st.dataframe(
            cohort.head(15),
            use_container_width=True,
            hide_index=True,
        )

        csv = cohort.to_csv(index=False).encode("utf-8")

        st.download_button(
            label="Download cleaned cohort CSV",
            data=csv,
            file_name="brca_survival_clean.csv",
            mime="text/csv",
        )

with tab2:
    st.header("Cohort description")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Stage distribution")

        stage_plot_df = (
            filtered_cohort["stage_group"]
            .value_counts(dropna=False)
            .rename_axis("stage_group")
            .reset_index(name="n")
        )

        stage_plot_df["stage_group"] = stage_plot_df["stage_group"].fillna("Missing")

        fig_stage = px.bar(
            stage_plot_df,
            x="stage_group",
            y="n",
            text="n",
            labels={
                "stage_group": "Pathological stage",
                "n": "Patients",
            },
            title="Patients by pathological stage",
        )

        fig_stage.update_layout(
            title_x=0.02,
            showlegend=False,
            height=420,
        )

        st.plotly_chart(fig_stage, use_container_width=True)

    with col2:
        st.subheader("Age distribution")

        fig_age = px.histogram(
            filtered_cohort,
            x="age",
            nbins=25,
            labels={"age": "Age at diagnosis"},
            title="Age at diagnosis",
        )

        fig_age.update_layout(
            title_x=0.02,
            showlegend=False,
            height=420,
        )

        st.plotly_chart(fig_age, use_container_width=True)

    col3, col4 = st.columns(2)

    with col3:
        st.subheader("Age group summary")
        st.dataframe(age_summary, use_container_width=True, hide_index=True)

    with col4:
        st.subheader("Gender summary")
        st.dataframe(gender_summary, use_container_width=True, hide_index=True)

    st.subheader("Stage summary")
    st.dataframe(stage_summary, use_container_width=True, hide_index=True)

with tab3:
    st.header("Kaplan-Meier survival analysis")

    st.markdown(
        """
        Kaplan-Meier curves estimate the probability of remaining alive over time.
        Patients who did not experience death before their last recorded follow-up
        are treated as censored observations.
        """
    )

    km_path = FIG_DIR / "brca_km_stage_plot.png"

    if km_path.exists():
        st.image(
            str(km_path),
            caption="Overall survival by pathological stage",
            use_container_width=True,
        )
    else:
        st.warning("Kaplan-Meier plot not found. Run R/04_kaplan_meier.R to regenerate.")

    st.info(
        "The Kaplan-Meier curves show lower overall survival among patients with "
        "more advanced pathological stage. The log-rank test indicated a strong "
        "difference between stage groups."
    )

with tab4:
    st.header("Cox proportional hazards model")

    st.markdown(
        """
        The Cox model estimates hazard ratios for overall mortality.
        The model included **age at diagnosis** and **pathological stage**.
        Stage I was used as the reference group.
        """
    )

    formatted_cox = format_cox_results(cox_results)

    st.subheader("Hazard ratio table")

    st.dataframe(
        formatted_cox,
        use_container_width=True,
        hide_index=True,
    )

    forest_path = FIG_DIR / "brca_cox_forest_plot.png"

    if forest_path.exists():
        st.subheader("Hazard ratio forest plot")
        st.image(
            str(forest_path),
            caption="Cox model hazard ratios",
            use_container_width=True,
        )
    else:
        st.warning("Cox forest plot not found. Run R/05_cox_regression.R to regenerate.")

    st.subheader("Interpretation")

    st.markdown(
        """
        - A hazard ratio above 1 indicates higher mortality hazard.
        - Each additional year of age was associated with higher estimated hazard.
        - Compared with Stage I, more advanced stages showed progressively higher
          mortality hazard.
        - Stage IV showed the highest hazard ratio, but this estimate should be
          interpreted carefully because the Stage IV subgroup is small.
        """
    )

    st.subheader("Proportional hazards assumption test")

    st.dataframe(ph_test, use_container_width=True, hide_index=True)

    global_ph = ph_test.loc[ph_test["term"] == "GLOBAL", "p_value"]

    if not global_ph.empty and float(global_ph.iloc[0]) < 0.05:
        st.warning(
            "The global proportional hazards test is statistically significant. "
            "This suggests that the Cox proportional hazards assumption may not "
            "hold perfectly. Results should be interpreted as average hazard "
            "ratios over follow-up."
        )
    else:
        st.success("No major global proportional hazards issue detected.")

with tab5:
    st.header("Machine learning risk prediction demo")

    st.markdown(
        """
        This section demonstrates a Python-based machine learning extension to the
        survival analysis workflow. The prediction task is **death within 3 years
        after diagnosis**.

        Patients censored before 3 years are excluded from this prediction target
        because their 3-year outcome status is unknown. This avoids incorrectly
        treating early-censored patients as survivors.
        """
    )

    if ml_status is None or ml_results is None:
        st.warning(
            "ML outputs not found. Run notebooks/02_ml_risk_prediction_demo.ipynb to regenerate."
        )
    else:
        st.subheader("3-year outcome definition")

        st.dataframe(
            format_ml_status(ml_status),
            use_container_width=True,
            hide_index=True,
        )

        st.info(
            "Only patients with a known 3-year outcome status are used for model training "
            "and evaluation. Patients censored before 3 years are excluded because their "
            "3-year status is genuinely unknown — retaining them as 'survivors' would "
            "introduce outcome misclassification bias."
        )

        st.subheader("Model comparison")

        model_display = ml_results.copy()

        numeric_cols = [
            "auc",
            "accuracy",
            "balanced_accuracy",
            "sensitivity_recall",
            "specificity",
            "precision",
            "f1",
            "cv_auc_mean",
            "cv_auc_sd",
        ]

        for col in numeric_cols:
            if col in model_display.columns:
                model_display[col] = model_display[col].round(3)

        st.dataframe(
            model_display,
            use_container_width=True,
            hide_index=True,
        )

        selected_model = st.selectbox(
            "Select a model to inspect",
            model_display["model"].tolist(),
        )

        selected_row = model_display.loc[
            model_display["model"] == selected_model
        ].iloc[0]

        col1, col2, col3, col4 = st.columns(4)

        col1.metric("AUC", f"{selected_row['auc']:.3f}")
        col2.metric("CV AUC", f"{selected_row['cv_auc_mean']:.3f}")
        col3.metric("Sensitivity", f"{selected_row['sensitivity_recall']:.3f}")
        col4.metric("Specificity", f"{selected_row['specificity']:.3f}")

        st.markdown(
            """
            **Interpretation.** AUC measures discrimination, while sensitivity and
            specificity show the trade-off between detecting patients who died within
            3 years and correctly identifying those who did not. Because the outcome
            is imbalanced, precision is expected to be limited.
            """
        )

        roc_path = FIG_DIR / "brca_ml_roc_curves.png"
        auc_path = FIG_DIR / "brca_ml_auc_comparison.png"

        col_left, col_right = st.columns(2)

        with col_left:
            if roc_path.exists():
                st.image(
                    str(roc_path),
                    caption="ROC curves for 3-year mortality prediction",
                    use_container_width=True,
                )
            else:
                st.warning("ROC curve plot not found.")

        with col_right:
            if auc_path.exists():
                st.image(
                    str(auc_path),
                    caption="AUC comparison across ML models",
                    use_container_width=True,
                )
            else:
                st.warning("AUC comparison plot not found.")

        calibration_path = FIG_DIR / "brca_ml_calibration_curves.png"

        if calibration_path.exists():
            st.subheader("Calibration curves")
            st.markdown(
                "Calibration measures whether predicted probabilities match observed event "
                "rates. A well-calibrated model follows the diagonal. The Brier score "
                "summarises calibration error (lower = better)."
            )
            col_cal, _ = st.columns([1, 1])
            with col_cal:
                st.image(
                    str(calibration_path),
                    caption="Calibration curves for all models",
                    use_container_width=True,
                )

        st.subheader("Best model confusion matrix")

        confusion_path = FIG_DIR / "brca_ml_confusion_matrix_best_model.png"

        if confusion_path.exists():
            col_cm, _ = st.columns([1, 1])
            with col_cm:
                st.image(
                    str(confusion_path),
                    caption="Confusion matrix for the best-performing model",
                    use_container_width=True,
                )
        else:
            st.warning("Confusion matrix plot not found.")

        st.subheader("Model interpretation")

        logistic_path = FIG_DIR / "brca_ml_logistic_odds_ratios.png"
        importance_path = FIG_DIR / "brca_ml_feature_importance.png"

        col_a, col_b = st.columns(2)

        with col_a:
            if logistic_path.exists():
                st.image(
                    str(logistic_path),
                    caption="Logistic regression age and stage associations",
                    use_container_width=True,
                )

        with col_b:
            if importance_path.exists():
                st.image(
                    str(importance_path),
                    caption="Random forest feature importance",
                    use_container_width=True,
                )

        if ml_feature_importance is not None:
            with st.expander("View random forest feature-importance table"):
                st.dataframe(
                    ml_feature_importance,
                    use_container_width=True,
                    hide_index=True,
                )

        if ml_logistic_coefficients is not None:
            with st.expander("View logistic regression coefficients table"):
                st.dataframe(
                    ml_logistic_coefficients,
                    use_container_width=True,
                    hide_index=True,
                )

        st.warning(
            "This ML analysis is a methodological demonstration. It is not intended "
            "for clinical decision support. The models use a limited set of TCGA "
            "clinical variables and would require external validation before any "
            "clinical interpretation."
        )

with tab6:
    st.header("Methods and limitations")

    st.subheader("Dataset")

    st.markdown(
        """
        This project uses the open-access TCGA Pan-Cancer Clinical Data Resource,
        focusing on TCGA-BRCA breast cancer patients. The analysis uses curated
        clinical and survival outcome data only.
        """
    )

    st.subheader("Endpoint")

    st.markdown(
        """
        Overall survival was defined using:

        - `OS`: event indicator, where 1 indicates death and 0 indicates censoring.
        - `OS.time`: follow-up time in days, converted to months for analysis.
        """
    )

    st.subheader("Statistical methods")

    st.markdown(
        """
        - Descriptive cohort summaries
        - Kaplan-Meier survival curves
        - Log-rank comparison of survival curves
        - Cox proportional hazards regression
        - Proportional hazards assumption testing using Schoenfeld residuals
        - Censoring-aware 3-year mortality prediction
        - Logistic regression, random forest, and gradient boosting
        - Train/test evaluation and stratified cross-validation
        - ROC/AUC, sensitivity, specificity, precision, and F1-score
        """
    )

    st.subheader("Limitations")

        st.markdown(
        """
        - TCGA is not a population-based cancer registry.
        - The analysis is observational and should not be interpreted causally.
        - Follow-up may be limited for some patients.
        - Missingness varies across clinical variables.
        - The Stage IV subgroup is small in the TCGA-BRCA subset.
        - The proportional hazards assumption may not hold perfectly for stage.
        - Results are for educational and portfolio purposes, not clinical decision support.
        - The ML extension is a methodological risk-prediction demo and has not been externally validated.
        """
    )

with tab7:
    st.header("Multi-cancer comparison: BRCA vs LUAD")

    st.markdown(
        """
        This section extends the analysis to **TCGA-LUAD** (lung adenocarcinoma) and
        compares survival patterns with **TCGA-BRCA** (breast cancer). Both cohorts use
        the same endpoint (overall survival), the same predictor set (age + pathological
        stage), and the same Cox model specification, enabling direct comparison of
        effect estimates across cancer types.
        """
    )

    if mc_cohort is None or mc_cox is None:
        st.warning(
            "Multi-cancer comparison outputs not found. "
            "Run R/07 through R/11 to generate LUAD and comparison files."
        )
    else:
        section_divider()
        st.subheader("Cohort comparison")

        display_cohort = mc_cohort.rename(columns={
            "cancer_type":             "Cancer type",
            "n_patients":              "Patients (n)",
            "n_os_events":             "OS events (n)",
            "os_event_rate_pct":       "Event rate (%)",
            "median_follow_up_months": "Median follow-up (months)",
            "median_age":              "Median age",
            "n_missing_stage":         "Missing stage (n)",
        })

        st.dataframe(display_cohort, use_container_width=True, hide_index=True)

        st.info(
            "LUAD typically has a substantially higher event rate and shorter "
            "median follow-up than BRCA, reflecting the more aggressive natural "
            "history of lung adenocarcinoma."
        )

        section_divider()
        st.subheader("Kaplan-Meier survival comparison")

        st.markdown(
            "Overall survival curves are shown for the full cohort of each cancer "
            "type (not stratified by stage). The shaded bands are 95% confidence intervals."
        )

        km_mc_path = FIG_DIR / "multicancer_km_comparison.png"

        if km_mc_path.exists():
            st.image(
                str(km_mc_path),
                caption="Overall survival: TCGA-BRCA vs TCGA-LUAD",
                use_container_width=True,
            )
        else:
            st.warning("Multi-cancer KM plot not found. Run R/11_multicancer_comparison.R.")

        section_divider()
        st.subheader("Cox model hazard ratio comparison")

        st.markdown(
            """
            Both models use the same specification: `Surv(os_time_months, os_event) ~ age + stage_group`
            with Stage I as the reference. Hazard ratios from each cancer type are overlaid
            for direct comparison.
            """
        )

        cox_mc_path = FIG_DIR / "multicancer_cox_comparison.png"

        if cox_mc_path.exists():
            st.image(
                str(cox_mc_path),
                caption="Cox hazard ratios: BRCA vs LUAD",
                use_container_width=True,
            )
        else:
            st.warning("Multi-cancer Cox plot not found. Run R/11_multicancer_comparison.R.")

        st.subheader("Hazard ratio table")

        def format_mc_cox(df: pd.DataFrame) -> pd.DataFrame:
            display = df.copy()
            display["Hazard ratio"] = display["estimate"].round(2)
            display["95% CI"] = (
                display["conf.low"].round(2).astype(str)
                + "\u2013"
                + display["conf.high"].round(2).astype(str)
            )
            display["p-value"] = display["p.value"].apply(
                lambda x: "<0.001" if x < 0.001 else f"{x:.3f}"
            )
            return display[["cancer_type", "term_label", "Hazard ratio", "95% CI", "p-value"]].rename(
                columns={"cancer_type": "Cancer type", "term_label": "Predictor"}
            )

        st.dataframe(
            format_mc_cox(mc_cox),
            use_container_width=True,
            hide_index=True,
        )

        section_divider()
        st.subheader("Interpretation")

        st.markdown(
            """
            - **Stage gradient:** Both cancers show a clear dose–response relationship
              between advancing stage and higher mortality hazard, validating stage as a
              strong prognostic marker across cancer types.
            - **Magnitude of stage effects:** Stage III and IV hazard ratios tend to be
              larger in LUAD than BRCA, reflecting the steeper survival gradient in lung cancer.
            - **Age effect:** The per-year age hazard ratio is typically similar across
              both cancers, suggesting age-related risk is relatively consistent.
            - **Baseline survival:** Even within the same stage, LUAD patients have
              substantially lower baseline survival than BRCA patients.
            - **Interpretation caveat:** Both models are fitted independently within each
              cancer type and share no patients. Differences in HRs do not directly imply
              causal differences — they may reflect differences in treatment, follow-up
              patterns, or unmeasured confounders.
            """
        )