# Immune Cell Population Analysis

Clinical trial analysis tool for exploring how a drug candidate affects immune cell populations across patient samples.

---

## Setup and running

### Requirements

Python 3.9+ is required. All dependencies are listed in `requirements.txt`.

### Install dependencies

```bash
make setup
# or: pip install -r requirements.txt
```

### Run the full pipeline

```bash
make pipeline
```

This runs `load_data.py` (initializes the SQLite database and loads `cell-count.csv`) followed by `analysis.py` (generates all output tables and plots). Outputs are written to the `outputs/` directory.

### Start the dashboard

```bash
make dashboard
# or: streamlit run dashboard.py
```

The dashboard will open at `http://localhost:8501`.

> **Dashboard link:** https://teiko-technical-sedveykyy9rbdsfwrfyogn.streamlit.app/

---

## Database schema

The schema is split into four normalized tables:

```
projects    (project_id PK)
subjects    (subject_id PK, project_id FK, condition, age, sex, treatment, response)
samples     (sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
cell_counts (id PK, sample_id FK, population TEXT, count INTEGER)
```

Cell counts are stored in a long/tidy format — one row per population per sample — rather than as five separate columns on the sample row.

**Rationale and scalability:**

- **Separation of entities** — projects, subjects, samples, and measurements are distinct concepts with different cardinalities and update frequencies. Keeping them separate avoids redundancy (e.g., a subject's age and sex don't need to repeat across every sample row) and makes updates cleaner.

- **Unpivoted cell counts** — storing counts as `(sample_id, population, count)` rows rather than wide columns means adding a new cell population (e.g., `treg_cell`) requires inserting new rows, not altering the schema. At scale — thousands of assays with dozens of marker panels — this matters a lot. It also makes aggregation queries (sum counts per sample, filter by population) straightforward with standard SQL.

- **Indexed foreign keys** — indexes on `cell_counts(sample_id)`, `cell_counts(population)`, `samples(subject_id)`, and `subjects(project_id)` keep joins fast as rows grow into the tens of millions.

- **Scalability path** — at hundreds of projects and thousands of samples, the core schema holds without changes. You would add: a `panels` table to group populations by assay type, a `timepoints` table if time metadata becomes richer, and potentially partition `cell_counts` by project or time range. The current foreign key structure makes those extensions additive.

---

## Code structure

```
.
├── load_data.py      # Part 1 — initializes SQLite DB, loads CSV
├── analysis.py       # Parts 2-4 — generates tables, stats, and plots
├── dashboard.py      # Interactive Streamlit dashboard
├── requirements.txt
├── Makefile
├── cell-count.csv
└── outputs/
    ├── cell_frequencies.csv          # Part 2
    ├── statistical_results.csv       # Part 3
    ├── part4_by_project.csv          # Part 4
    ├── part4_by_response.csv         # Part 4
    ├── part4_by_sex.csv              # Part 4
    ├── part4_bcell_answer.txt        # Part 4 answer
    └── plots/
        └── boxplots_responders_vs_nonresponders.png
```

**Design decisions:**

- `load_data.py` is kept minimal — its only job is schema creation and data ingestion. It drops and recreates the database on each run so the pipeline is fully reproducible.

- `analysis.py` contains three clearly separated functions (`part2_*`, `part3_*`, `part4_*`), each responsible for one section. This makes it easy to re-run or modify a single part without touching the others.

- `dashboard.py` queries the database directly (using the same SQL as `analysis.py`) rather than reading the output CSVs. This keeps the dashboard live and filterable, and means it stays correct even if the underlying data changes.

- Streamlit was chosen for the dashboard because it requires no frontend code, is trivially deployable to Streamlit Cloud, and handles the interactive filtering and charting with minimal boilerplate.

---

## Part 3 findings

**Do relative frequencies of any cell population differ between responders and non-responders?**

Filtered to: condition=melanoma, treatment=miraclib, sample_type=PBMC. For each of the 5 populations, responders vs. non-responders were compared with a Mann-Whitney U test (two-sided, non-parametric — no normality assumption on relative frequency).

Since 5 populations are tested simultaneously, raw p-values were also Bonferroni-corrected (adjusted α = 0.05 / 5 = 0.010) to control the false-positive rate across the family of tests.

| Population | p-value | p-value (Bonferroni) | Significant (raw) | Significant (Bonferroni) |
|---|---|---|---|---|
| cd4_t_cell | 0.0134 | 0.0670 | Yes | **No** |
| b_cell | 0.0557 | 0.2787 | No | No |
| nk_cell | 0.1211 | 0.6053 | No | No |
| monocyte | 0.1635 | 0.8175 | No | No |
| cd8_t_cell | 0.6392 | 1.0000 | No | No |

**Conclusion:** `cd4_t_cell` is the only population with a raw p < 0.05, but it does not survive Bonferroni correction. At this sample size, none of the 5 populations show a statistically robust difference in relative frequency between responders and non-responders — the `cd4_t_cell` signal should be treated as a weak, unconfirmed lead worth following up with more data rather than a finding to act on.

Full results: `outputs/statistical_results.csv`. Boxplots: `outputs/plots/boxplots_responders_vs_nonresponders.png`.

---

## Part 4 answer

**Considering melanoma males, what is the average number of B cells for responders at time=0?**

**10401.28**

(Filtered to: condition=melanoma, sex=M, response=yes, time_from_treatment_start=0, sample_type=PBMC, treatment=miraclib)
