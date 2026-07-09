import sqlite3
import pandas as pd
import plotly.express as px
from scipy import stats
import streamlit as st

DB_PATH = "teiko.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

st.set_page_config(
    page_title="Immune Cell Dashboard",
    page_icon="🧬",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Data loaders (cached)
# ---------------------------------------------------------------------------

@st.cache_data
def load_frequencies():
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql(
        """
        SELECT
            s.sample_id  AS Sample,
            su.project_id AS project,
            su.condition,
            su.treatment,
            su.response,
            su.sex,
            s.sample_type,
            s.time_from_treatment_start,
            cc.population AS Population,
            cc.count      AS Count
        FROM samples s
        JOIN subjects su ON s.subject_id = su.subject_id
        JOIN cell_counts cc ON s.sample_id = cc.sample_id
        """,
        conn,
    )
    conn.close()
    totals = df.groupby("Sample")["Count"].sum().rename("Total_count")
    df = df.join(totals, on="Sample")
    df["Percentage"] = (df["Count"] / df["Total_count"] * 100).round(2)
    return df


@st.cache_data
def compute_stats():
    df = load_frequencies()
    subset = df[
        (df["condition"] == "melanoma")
        & (df["treatment"] == "miraclib")
        & (df["sample_type"] == "PBMC")
    ].copy()

    n_tests = len(POPULATIONS)
    rows = []
    for pop in POPULATIONS:
        pop_data = subset[subset["Population"] == pop]
        resp = pop_data[pop_data["response"] == "yes"]["Percentage"].values
        non_resp = pop_data[pop_data["response"] == "no"]["Percentage"].values
        stat, p = stats.mannwhitneyu(resp, non_resp, alternative="two-sided")
        p_bonf = min(p * n_tests, 1.0)
        rows.append(
            {
                "Population": pop,
                "n (responders)": len(resp),
                "n (non-responders)": len(non_resp),
                "Median — responders": round(float(pd.Series(resp).median()), 2),
                "Median — non-responders": round(float(pd.Series(non_resp).median()), 2),
                "Mann-Whitney U": round(stat, 2),
                "p-value": round(p, 4),
                "Significant (raw, p<0.05)": "Yes" if p < 0.05 else "No",
                "p-value (Bonferroni)": round(p_bonf, 4),
                "Significant (Bonferroni)": "Yes" if p_bonf < 0.05 else "No",
            }
        )
    return pd.DataFrame(rows), subset


@st.cache_data
def load_baseline():
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql(
        """
        SELECT
            su.subject_id,
            su.project_id AS project,
            su.sex,
            su.response,
            s.sample_id,
            cc.population,
            cc.count
        FROM samples s
        JOIN subjects su ON s.subject_id = su.subject_id
        JOIN cell_counts cc ON s.sample_id = cc.sample_id
        WHERE su.condition              = 'melanoma'
          AND s.sample_type             = 'PBMC'
          AND s.time_from_treatment_start = 0
          AND su.treatment              = 'miraclib'
        """,
        conn,
    )
    conn.close()
    return df


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

def main():
    st.title("Immune Cell Population Analysis")
    st.caption("Clinical trial data — Loblaw Bio")

    tab1, tab2, tab3 = st.tabs([
        "Part 2 — Cell Frequencies",
        "Part 3 — Statistical Analysis",
        "Part 4 — Subset Analysis",
    ])

    # -----------------------------------------------------------------------
    # Tab 1 — Frequencies
    # -----------------------------------------------------------------------
    with tab1:
        st.header("Relative frequency of each cell population per sample")
        df = load_frequencies()

        col1, col2, col3, col4 = st.columns(4)
        conditions  = ["All"] + sorted(df["condition"].dropna().unique().tolist())
        treatments  = ["All"] + sorted(df["treatment"].dropna().unique().tolist())
        stypes      = ["All"] + sorted(df["sample_type"].dropna().unique().tolist())
        pop_opts    = ["All"] + POPULATIONS

        sel_cond   = col1.selectbox("Condition",    conditions)
        sel_treat  = col2.selectbox("Treatment",    treatments)
        sel_stype  = col3.selectbox("Sample type",  stypes)
        sel_pop    = col4.selectbox("Population",   pop_opts)

        filt = df.copy()
        if sel_cond  != "All": filt = filt[filt["condition"]   == sel_cond]
        if sel_treat != "All": filt = filt[filt["treatment"]   == sel_treat]
        if sel_stype != "All": filt = filt[filt["sample_type"] == sel_stype]
        if sel_pop   != "All": filt = filt[filt["Population"]  == sel_pop]

        st.dataframe(
            filt[["Sample", "Total_count", "Population", "Count", "Percentage"]].reset_index(drop=True),
            use_container_width=True,
            height=380,
        )

        if not filt.empty:
            fig = px.box(
                filt,
                x="Population",
                y="Percentage",
                color="condition",
                title="Relative frequency by population",
                labels={"Percentage": "Relative frequency (%)"},
            )
            st.plotly_chart(fig, use_container_width=True)

    # -----------------------------------------------------------------------
    # Tab 2 — Statistical analysis
    # -----------------------------------------------------------------------
    with tab2:
        st.header("Responders vs Non-responders")
        st.markdown(
            "Filter applied: **melanoma** patients on **miraclib**, **PBMC** samples only.  \n"
            "Test: Mann-Whitney U (two-sided), Bonferroni-corrected for "
            f"{len(POPULATIONS)} comparisons (adjusted α = 0.05 / {len(POPULATIONS)} = "
            f"{0.05 / len(POPULATIONS):.3f})."
        )

        results_df, subset = compute_stats()

        sig_pops_raw = results_df[results_df["Significant (raw, p<0.05)"] == "Yes"]["Population"].tolist()
        sig_pops_bonf = results_df[results_df["Significant (Bonferroni)"] == "Yes"]["Population"].tolist()

        if sig_pops_bonf:
            st.success(
                f"Significant after Bonferroni correction: **{', '.join(sig_pops_bonf)}**"
            )
        elif sig_pops_raw:
            st.warning(
                f"**{', '.join(sig_pops_raw)}** significant at raw p < 0.05, but does **not** "
                f"survive Bonferroni correction for {len(POPULATIONS)} comparisons — "
                "treat as a weak/unconfirmed signal, not a robust finding."
            )
        else:
            st.info("No populations showed a statistically significant difference at p < 0.05.")

        st.subheader("Summary table")
        st.dataframe(results_df, use_container_width=True)

        st.subheader("Boxplots — all populations")
        fig_box = px.box(
            subset,
            x="Population",
            y="Percentage",
            color="response",
            color_discrete_map={"yes": "#4C9BE8", "no": "#E85C5C"},
            points="all",
            labels={"Percentage": "Relative frequency (%)", "response": "Response"},
            title="Cell population frequencies: Responders vs Non-responders",
            category_orders={"response": ["yes", "no"]},
        )
        fig_box.update_layout(legend_title_text="Responder")
        st.plotly_chart(fig_box, use_container_width=True)

        st.subheader("Single population view")
        sel_pop2 = st.selectbox("Select population", POPULATIONS, key="stat_pop")
        pop_sub = subset[subset["Population"] == sel_pop2]
        pop_row = results_df[results_df["Population"] == sel_pop2].iloc[0]
        p_val, p_bonf_val = pop_row["p-value"], pop_row["p-value (Bonferroni)"]
        fig2 = px.box(
            pop_sub,
            x="response",
            y="Percentage",
            color="response",
            color_discrete_map={"yes": "#4C9BE8", "no": "#E85C5C"},
            points="all",
            labels={"Percentage": "Relative frequency (%)", "response": "Response"},
            title=f"{sel_pop2}  —  p = {p_val}  (Bonferroni p = {p_bonf_val})",
            category_orders={"response": ["yes", "no"]},
        )
        fig2.update_layout(showlegend=False)
        st.plotly_chart(fig2, use_container_width=True)

    # -----------------------------------------------------------------------
    # Tab 3 — Subset analysis
    # -----------------------------------------------------------------------
    with tab3:
        st.header("Baseline melanoma PBMC samples on miraclib")
        st.markdown(
            "Filter: **melanoma** · **PBMC** · **time_from_treatment_start = 0** · **miraclib**"
        )

        baseline = load_baseline()
        samples_df  = baseline[["sample_id", "subject_id", "project", "sex", "response"]].drop_duplicates("sample_id")
        subjects_df = samples_df.drop_duplicates("subject_id")

        m1, m2 = st.columns(2)
        m1.metric("Total samples",  len(samples_df))
        m2.metric("Unique subjects", len(subjects_df))

        c1, c2, c3 = st.columns(3)

        with c1:
            st.subheader("Samples per project")
            by_proj = (
                samples_df.groupby("project")["sample_id"]
                .count()
                .reset_index()
                .rename(columns={"sample_id": "Samples", "project": "Project"})
            )
            fig_proj = px.bar(by_proj, x="Project", y="Samples", text="Samples")
            fig_proj.update_traces(textposition="outside")
            fig_proj.update_layout(showlegend=False, margin=dict(t=30))
            st.plotly_chart(fig_proj, use_container_width=True)

        with c2:
            st.subheader("Subjects by response")
            by_resp = (
                subjects_df.groupby("response")["subject_id"]
                .count()
                .reset_index()
                .rename(columns={"subject_id": "Subjects", "response": "Response"})
            )
            fig_resp = px.pie(
                by_resp,
                values="Subjects",
                names="Response",
                color="Response",
                color_discrete_map={"yes": "#4C9BE8", "no": "#E85C5C"},
            )
            st.plotly_chart(fig_resp, use_container_width=True)

        with c3:
            st.subheader("Subjects by sex")
            by_sex = (
                subjects_df.groupby("sex")["subject_id"]
                .count()
                .reset_index()
                .rename(columns={"subject_id": "Subjects", "sex": "Sex"})
            )
            fig_sex = px.pie(
                by_sex,
                values="Subjects",
                names="Sex",
                color="Sex",
                color_discrete_map={"M": "#6CBDE8", "F": "#E8A06C"},
            )
            st.plotly_chart(fig_sex, use_container_width=True)

        st.divider()
        st.subheader("Average B cell count — melanoma males, responders, baseline")
        bcell_avg = baseline[
            (baseline["population"] == "b_cell")
            & (baseline["sex"] == "M")
            & (baseline["response"] == "yes")
        ]["count"].mean()
        st.metric("Average B cells", f"{bcell_avg:.2f}")


if __name__ == "__main__":
    main()
