## Setup & Running the Code

### Requirements

Python 3.9+ is required. Install dependencies with:

```
pip install pandas numpy scipy dash plotly
```

### Steps to Run

1. Ensure `cell-count.csv` and `load_data.py` are in the same directory
2. Run the script:

```bash
python load_data.py
```

This will:
- Connect to `cell_counts.db` (SQLite database)
- Load and compute all analysis results for Parts 2, 3, and 4
- Launch an interactive dashboard at `http://127.0.0.1:8050/`


**Note:** `load_data.py` reads from an existing database. If `cell_counts.db` does not exist yet, run `load_data.py` first to create and populate it:


## Database Schema

### Tables

**`subjects`** — one row per patient

| Column | Type | Description |
|---|---|---|
| `subject_id` | TEXT (PK) | Unique patient identifier |
| `project` | TEXT | Clinical trial project |
| `condition` | TEXT | Diagnosis (e.g. melanoma, carcinoma) |
| `age` | INTEGER | Patient age |
| `sex` | TEXT | Patient sex |
| `treatment` | TEXT | Drug administered |
| `response` | TEXT | Treatment response (yes/no) |

**`samples`** — one row per biological sample

| Column | Type | Description |
|---|---|---|
| `sample_id` | TEXT (PK) | Unique sample identifier |
| `subject_id` | TEXT (FK) | Links to `subjects` |
| `sample_type` | TEXT | Sample type (e.g. PBMC) |
| `time_from_treatment_start` | INTEGER | Days since treatment began |

**`cell_counts`** — one row per sample × cell type (long format)

| Column | Type | Description |
|---|---|---|
| `auto_id` | INTEGER (PK) | Auto-incremented surrogate key |
| `sample_id` | TEXT (FK) | Links to `samples` |
| `cell_type` | TEXT | Immune cell population name |
| `count` | INTEGER | Raw cell count |

### Design Rationale

**Separation of subjects and samples.** Subject-level facts (age, sex, condition, treatment, response) are constant across all of a patient's samples. Storing them once in `subjects` and referencing via foreign key avoids repeating the same values across every timepoint row, which reduces redundancy and prevents update anomalies so if a subject's metadata needs correcting, only one row changes.

**Long format for cell counts.** Rather than storing each cell type as its own column (`b_cell`, `cd8_t_cell`, etc.), counts are unpivoted into long format with a `cell_type` column. This has several advantages such as adding a new cell population requires no schema change (just new rows), queries that aggregate or filter by cell type are simpler (`WHERE cell_type = 'cd8_t_cell'` vs selecting a specific column), and the table remains well-normalised regardless of how many populations are measured.

**Foreign keys with referential integrity.** `PRAGMA foreign_keys = ON` is set at connection time so SQLite enforces that every `sample_id` in `cell_counts` exists in `samples`, and every `subject_id` in `samples` exists in `subjects`.

**`UNIQUE(sample_id, cell_type)` constraint on `cell_counts`.** Ensures a given cell type can only be recorded once per sample, and combined with `INSERT OR IGNORE` makes the loading script safe to re-run without creating duplicate rows.

### Scalability Considerations

With hundreds of projects, thousands of samples, and diverse analytics needs, the schema holds up well with some targeted additions:

**Indexing.** At scale the most expensive queries will be joins and filtered aggregations. Adding indexes on the foreign key columns and commonly filtered fields keeps query times fast:


**A `projects` lookup table.** Currently `project` is a plain text column on `subjects`. At hundreds of projects this should become its own table with a foreign key reference, so project-level metadata (start date, sponsor, therapeutic area) can be stored without repeating it on every subject row.

**Materialised views for derived analytics.** The `cell_counts` table in long format is ideal for flexible querying but requires aggregation on every read. For high-frequency dashboard queries, a pre-computed `sample_frequencies` table (sample × cell type × percentage) maintained by a scheduled job would significantly reduce query cost.

**Migration to a production database.** SQLite is appropriate for this analysis but would be replaced by PostgreSQL or similar at scale. Because the schema uses standard SQL with no SQLite-specific features beyond `PRAGMA`, migration is straightforward — the same `CREATE TABLE` statements and queries work with minor syntax adjustments.

---

## Code Structure

```
.
├── load_data.py      
├── cell-count.csv     
└── cell_counts.db     
```



### `load_data.py` 

Initialises the database schema and loads all rows from `cell-count.csv`. Each insert function handles one table (`insert_subjects`, `insert_samples`, `insert_cellCounts`). The CSV's five wide cell-type columns are unpivoted into long format using `pd.melt()` before inserting into `cell_counts`. `INSERT OR IGNORE` is used throughout so the script is safe to re-run against an existing database without creating duplicates.

Structured in three layers:

**Data functions** (`get_freq`, `melanoma_pbmc`, `cal_statistics`, `baseline_query`, etc.) each own a single SQL query and return a clean DataFrame. All cohort filtering (melanoma, miraclib, PBMC) happens at the database level so only relevant rows are loaded into memory. `load_all()` calls each of these once at startup and stores the results in a dict, so the database is never queried again during the session.

**Chart functions** (`make_boxplot`, `make_freq_bar`, `make_pie`) each take a DataFrame and return a Plotly figure. Keeping chart logic separate from layout logic means visuals can be updated without touching the Dash component tree.

**Layout and app** (`build_app`) assembles the three tab views from the pre-built DataFrames and figures, registers two callbacks (tab switching and the sample dropdown), and returns the configured app. Only the per-sample bar chart uses a callback since it changes based on user input — everything else is built once at startup and served as static HTML.

