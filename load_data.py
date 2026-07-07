import sqlite3
import csv
import os

DB_PATH = "teiko.db"
CSV_PATH = "cell-count.csv"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def init_db(conn):
    conn.executescript("""
        PRAGMA foreign_keys = ON;

        CREATE TABLE IF NOT EXISTS projects (
            project_id  TEXT PRIMARY KEY
        );

        CREATE TABLE IF NOT EXISTS subjects (
            subject_id  TEXT PRIMARY KEY,
            project_id  TEXT NOT NULL REFERENCES projects(project_id),
            condition   TEXT,
            age         INTEGER,
            sex         TEXT,
            treatment   TEXT,
            response    TEXT
        );

        CREATE TABLE IF NOT EXISTS samples (
            sample_id                   TEXT PRIMARY KEY,
            subject_id                  TEXT NOT NULL REFERENCES subjects(subject_id),
            sample_type                 TEXT,
            time_from_treatment_start   INTEGER
        );

        CREATE TABLE IF NOT EXISTS cell_counts (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id   TEXT NOT NULL REFERENCES samples(sample_id),
            population  TEXT NOT NULL,
            count       INTEGER NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_cc_sample     ON cell_counts(sample_id);
        CREATE INDEX IF NOT EXISTS idx_cc_population ON cell_counts(population);
        CREATE INDEX IF NOT EXISTS idx_samples_subj  ON samples(subject_id);
        CREATE INDEX IF NOT EXISTS idx_subj_project  ON subjects(project_id);
    """)
    conn.commit()


def load_csv(conn, csv_path):
    cur = conn.cursor()
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cur.execute(
                "INSERT OR IGNORE INTO projects VALUES (?)",
                (row["project"],),
            )
            cur.execute(
                "INSERT OR IGNORE INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    row["subject"],
                    row["project"],
                    row["condition"] or None,
                    int(row["age"]) if row["age"] else None,
                    row["sex"] or None,
                    row["treatment"] or None,
                    row["response"] or None,
                ),
            )
            cur.execute(
                "INSERT OR IGNORE INTO samples VALUES (?, ?, ?, ?)",
                (
                    row["sample"],
                    row["subject"],
                    row["sample_type"] or None,
                    int(row["time_from_treatment_start"]) if row["time_from_treatment_start"] else None,
                ),
            )
            for pop in POPULATIONS:
                cur.execute(
                    "INSERT INTO cell_counts (sample_id, population, count) VALUES (?, ?, ?)",
                    (row["sample"], pop, int(row[pop])),
                )
    conn.commit()


if __name__ == "__main__":
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    conn = sqlite3.connect(DB_PATH)
    try:
        init_db(conn)
        load_csv(conn, CSV_PATH)
        print(f"Database created: {DB_PATH}")
    finally:
        conn.close()
