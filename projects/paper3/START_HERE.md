# Paper 3 Intern Start Guide (Detailed)

Read this file top to bottom. Copy and paste commands exactly.

## Project goal

Identify failure of ARG detection under sequence divergence.

You will produce:

- recall-vs-identity figure (BLAST + embedding lines)
- identity-stratified evaluation table (30-50, 50-70, 70-90)
- benchmark summary table (overall, low-ID, high-ID)
- identity-clustered split files

## Where to run commands

Always run commands from the repository root folder (`OPALS_Group19`) on your own computer.

Activate your local virtual environment first.

Windows PowerShell:

```powershell
./.venv/Scripts/Activate.ps1
```

Linux/WSL/macOS:

```bash
source .venv/bin/activate
```

## Before you start

Make sure these files exist from Paper 2:

- `projects/paper2/results/low_identity_per_query.csv`
- `projects/paper2/results/query_labels_card.csv`
- `projects/paper2/results/db_labels_card.csv`
- `projects/paper2/proteins/card_query_homolog.faa`
- `projects/paper2/proteins/card_db_homolog.faa`

## Step 1: Create identity-clustered split (required for robust evaluation)

Run:

```powershell
python projects/paper3/analysis/identity_cluster_split.py --pairwise-identity-csv projects/paper3/results/query_query_identity.csv --labels-csv projects/paper2/results/query_labels_card.csv --identity-threshold 70 --test-fraction 0.2 --out-assignments projects/paper3/results/split_assignments.csv --out-summary projects/paper3/results/split_summary.csv
```

Note:

- `query_query_identity.csv` must include `query_id,subject_id,pident` from a query-vs-query BLAST run.

## Step 2: Run divergence benchmark pipeline

```powershell
python projects/paper3/analysis/build_divergence_results.py --blast-per-query projects/paper2/results/low_identity_per_query.csv --query-labels projects/paper2/results/query_labels_card.csv --db-labels projects/paper2/results/db_labels_card.csv --query-fasta projects/paper2/proteins/card_query_homolog.faa --db-fasta projects/paper2/proteins/card_db_homolog.faa --embedding-model esm2 --hf-model-id facebook/esm2_t6_8M_UR50D --batch-size 16 --max-length 192 --identity-bin-out projects/paper3/results/identity_bin_recall.csv --benchmark-table-out projects/paper3/results/benchmark_table.csv --figure-out projects/paper3/analysis/figures/recall_vs_identity.png --embedding-per-query-out projects/paper3/results/embedding_per_query.csv
```

Success check:


- `projects/paper3/results/identity_bin_recall.csv` exists.
- `projects/paper3/results/benchmark_table.csv` exists.
- `projects/paper3/analysis/figures/recall_vs_identity.png` exists.

## Step 3: Optional real embedding (ProtBERT or ESM)

If `torch` and `transformers` are installed, run:

```powershell
python projects/paper3/analysis/build_divergence_results.py --blast-per-query projects/paper2/results/low_identity_per_query.csv --query-labels projects/paper2/results/query_labels_card.csv --db-labels projects/paper2/results/db_labels_card.csv --query-fasta projects/paper2/proteins/card_query_homolog.faa --db-fasta projects/paper2/proteins/card_db_homolog.faa --embedding-model protbert --identity-bin-out projects/paper3/results/identity_bin_recall_protbert.csv --benchmark-table-out projects/paper3/results/benchmark_table_protbert.csv --figure-out projects/paper3/analysis/figures/recall_vs_identity_protbert.png --embedding-per-query-out projects/paper3/results/embedding_per_query_protbert.csv
```

## Files to submit

- `projects/paper3/results/identity_bin_recall.csv`
- `projects/paper3/results/benchmark_table.csv`
- `projects/paper3/analysis/figures/recall_vs_identity.png`
- `projects/paper3/results/split_assignments.csv`
- `projects/paper3/results/split_summary.csv`

## Common mistakes

- Using random row-level split instead of identity-clustered split.
- Reporting only overall recall and not identity-stratified recall.
- Claiming "performed better" without quantifying recall decrease under divergence.
- Running ProtBERT/ESM mode without `torch` and `transformers` installed.
