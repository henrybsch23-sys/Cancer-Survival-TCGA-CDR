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
    cohort = pd.read_csv(DATA_DIR / "tcga_brca_survival_clean.csv")
    cohort_summary = pd.read_csv(DATA_DIR / "cohort_summary.csv")
    stage_summary = pd.read_csv(DATA_DIR / "stage_summary.csv")
    gender_summary = pd.read_csv(DATA_DIR / "gender_summary.csv")
    age_summary = pd.read_csv(DATA_DIR / "age_summary.csv")
    cox_results = pd.read_csv(DATA_DIR / "cox_results.csv")
    ph_test = pd.read_csv(DATA_DIR / "cox_ph_assumption_test.csv")
    return (
        cohort,
        cohort_summary,
        stage_summary,
        gender_summary,
        age_summary,
        cox_results,
        ph_test,
    )


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
) = load_data()

# Sidebar
with st.sidebar:
    st.title("TCGA Survival App")
    st.markdown(
        """
        **Dataset:** TCGA-CDR  
        **Cancer type:** TCGA-BRCA  
        **Endpoint:** Overall survival  
        **Main model:** Cox regression
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
    model diagnostics, and epidemiological interpretation.
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

tab1, tab2, tab3, tab4, tab5 = st.tabs(
    [
        "Overview",
        "Cohort description",
        "Kaplan-Meier analysis",
        "Cox regression",
        "Methods and limitations",
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
            file_name="tcga_brca_survival_clean.csv",
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

    km_path = FIG_DIR / "km_stage_plot.png"

    if km_path.exists():
        st.image(
            str(km_path),
            caption="Overall survival by pathological stage",
            use_container_width=True,
        )
    else:
        st.warning("Kaplan-Meier plot not found. Run R/04_kaplan_meier.R first.")

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

    forest_path = FIG_DIR / "cox_forest_plot.png"

    if forest_path.exists():
        st.subheader("Hazard ratio forest plot")
        st.image(
            str(forest_path),
            caption="Cox model hazard ratios",
            use_container_width=True,
        )
    else:
        st.warning("Cox forest plot not found. Run R/05_cox_regression.R first.")

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
        """
    )