import os
import sqlite3
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

DB_PATH = "teiko.db"
OUTPUT_DIR = "outputs"
PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

os.makedirs(PLOTS_DIR, exist_ok=True)


def get_connection():
    return sqlite3.connect(DB_PATH)


# ---------------------------------------------------------------------------
# Part 2 — relative frequencies
# ---------------------------------------------------------------------------

def part2_cell_frequencies(conn):
    df = pd.read_sql(
        """
        SELECT
            s.sample_id  AS Sample,
            cc.population AS Population,
            cc.count      AS Count
        FROM samples s
        JOIN cell_counts cc ON s.sample_id = cc.sample_id
        ORDER BY s.sample_id, cc.population
        """,
        conn,
    )

    totals = df.groupby("Sample")["Count"].sum().rename("Total_count")
    df = df.join(totals, on="Sample")
    df["Percentage"] = (df["Count"] / df["Total_count"] * 100).round(2)
    df = df[["Sample", "Total_count", "Population", "Count", "Percentage"]]

    out = os.path.join(OUTPUT_DIR, "cell_frequencies.csv")
    df.to_csv(out, index=False)
    print(f"[Part 2] {len(df)} rows saved to {out}")
    return df


# ---------------------------------------------------------------------------
# Part 3 — statistical analysis
# ---------------------------------------------------------------------------

def part3_statistical_analysis(conn, freq_df):
    meta = pd.read_sql(
        """
        SELECT s.sample_id, su.condition, su.treatment, su.response, s.sample_type
        FROM samples s
        JOIN subjects su ON s.subject_id = su.subject_id
        """,
        conn,
    )

    df = freq_df.merge(meta, left_on="Sample", right_on="sample_id")
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
                "n_responders": len(resp),
                "n_non_responders": len(non_resp),
                "median_responders": round(float(pd.Series(resp).median()), 2),
                "median_non_responders": round(float(pd.Series(non_resp).median()), 2),
                "MannWhitneyU": round(stat, 4),
                "p_value": round(p, 4),
                "significant": p < 0.05,
                "p_value_bonferroni": round(p_bonf, 4),
                "significant_bonferroni": p_bonf < 0.05,
            }
        )

    results_df = pd.DataFrame(rows)
    out = os.path.join(OUTPUT_DIR, "statistical_results.csv")
    results_df.to_csv(out, index=False)
    print(f"\n[Part 3] Statistical results saved to {out}")
    print(results_df.to_string(index=False))

    sig_raw = results_df[results_df["significant"]]["Population"].tolist()
    sig_bonf = results_df[results_df["significant_bonferroni"]]["Population"].tolist()
    print(f"\n[Part 3] Significant at raw p<0.05: {sig_raw or 'none'}")
    print(
        f"[Part 3] Significant after Bonferroni correction "
        f"(n={n_tests} tests, adjusted alpha={0.05 / n_tests:.3f}): {sig_bonf or 'none'}"
    )

    _save_boxplots(subset, results_df)
    return results_df


def _save_boxplots(subset, results_df):
    fig, axes = plt.subplots(1, 5, figsize=(20, 6))
    colors = {"yes": "#4C9BE8", "no": "#E85C5C"}
    labels = ["Responders", "Non-resp."]

    for ax, pop in zip(axes, POPULATIONS):
        pop_data = subset[subset["Population"] == pop]
        resp_vals = pop_data[pop_data["response"] == "yes"]["Percentage"].values
        non_resp_vals = pop_data[pop_data["response"] == "no"]["Percentage"].values

        bp = ax.boxplot(
            [resp_vals, non_resp_vals],
            tick_labels=labels,
            patch_artist=True,
            medianprops=dict(color="black", linewidth=2),
            flierprops=dict(marker="o", markersize=3, alpha=0.5),
        )
        bp["boxes"][0].set_facecolor(colors["yes"])
        bp["boxes"][1].set_facecolor(colors["no"])

        row = results_df[results_df["Population"] == pop].iloc[0]
        p, p_bonf = row["p_value"], row["p_value_bonferroni"]
        if p_bonf < 0.05:
            sig_label = "** sig. (Bonferroni)"
        elif p < 0.05:
            sig_label = "* sig. (uncorrected)"
        else:
            sig_label = "ns"
        ax.set_title(f"{pop}\np={p:.4f}  {sig_label}", fontsize=9)
        ax.set_ylabel("Relative frequency (%)")
        ax.tick_params(axis="x", labelsize=8)

    n_tests = len(POPULATIONS)
    fig.suptitle(
        "Immune cell frequencies — Responders vs Non-responders\n"
        "(Melanoma · Miraclib · PBMC)  ·  Mann-Whitney U, "
        f"Bonferroni-corrected for {n_tests} tests (adjusted α={0.05 / n_tests:.3f})",
        fontsize=12,
        fontweight="bold",
    )
    plt.tight_layout()

    out = os.path.join(PLOTS_DIR, "boxplots_responders_vs_nonresponders.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[Part 3] Boxplot saved to {out}")


# ---------------------------------------------------------------------------
# Part 4 — subset analysis
# ---------------------------------------------------------------------------

def part4_subset_analysis(conn):
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
        WHERE su.condition            = 'melanoma'
          AND s.sample_type           = 'PBMC'
          AND s.time_from_treatment_start = 0
          AND su.treatment            = 'miraclib'
        """,
        conn,
    )

    samples_df = df[["sample_id", "subject_id", "project", "sex", "response"]].drop_duplicates("sample_id")
    subjects_df = samples_df.drop_duplicates("subject_id")

    print(f"\n[Part 4] Baseline melanoma PBMC miraclib samples: {len(samples_df)}")

    # 4.2a — samples per project
    by_project = (
        samples_df.groupby("project")["sample_id"]
        .count()
        .reset_index()
        .rename(columns={"sample_id": "sample_count"})
    )
    print("\n  Samples per project:")
    print(by_project.to_string(index=False))
    by_project.to_csv(os.path.join(OUTPUT_DIR, "part4_by_project.csv"), index=False)

    # 4.2b — subjects by response
    by_response = (
        subjects_df.groupby("response")["subject_id"]
        .count()
        .reset_index()
        .rename(columns={"subject_id": "subject_count"})
    )
    print("\n  Subjects by response:")
    print(by_response.to_string(index=False))
    by_response.to_csv(os.path.join(OUTPUT_DIR, "part4_by_response.csv"), index=False)

    # 4.2c — subjects by sex
    by_sex = (
        subjects_df.groupby("sex")["subject_id"]
        .count()
        .reset_index()
        .rename(columns={"subject_id": "subject_count"})
    )
    print("\n  Subjects by sex:")
    print(by_sex.to_string(index=False))
    by_sex.to_csv(os.path.join(OUTPUT_DIR, "part4_by_sex.csv"), index=False)

    # Answer: avg B cells for melanoma male responders at time=0
    bcell_avg = df[
        (df["population"] == "b_cell")
        & (df["sex"] == "M")
        & (df["response"] == "yes")
    ]["count"].mean()

    answer = f"Average B cell count (melanoma, male, responder, time=0): {bcell_avg:.2f}"
    print(f"\n  {answer}")

    with open(os.path.join(OUTPUT_DIR, "part4_bcell_answer.txt"), "w") as f:
        f.write(answer + "\n")

    return by_project, by_response, by_sex, bcell_avg


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    conn = get_connection()
    try:
        print("=" * 60)
        print("PART 2 — Cell Frequency Table")
        print("=" * 60)
        freq_df = part2_cell_frequencies(conn)

        print("\n" + "=" * 60)
        print("PART 3 — Statistical Analysis")
        print("=" * 60)
        part3_statistical_analysis(conn, freq_df)

        print("\n" + "=" * 60)
        print("PART 4 — Subset Analysis")
        print("=" * 60)
        part4_subset_analysis(conn)
    finally:
        conn.close()
