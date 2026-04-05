# Paper 3: ARG Detection Failure Under Sequence Divergence

For interns with limited coding experience: start with `START_HERE.md` in this folder.

## Objective

Quantify where ARG detection fails as sequence divergence increases.

Core contribution:

- identify divergence-driven failure regimes, not just average method ranking
- enforce identity-aware evaluation to reduce train/test leakage
- report stratified recall (30-50, 50-70, 70-90 identity bins)

Hypotheses:

- H1: BLAST performance drops with divergence
- H2: embeddings are more robust
- H3: standard evaluation overestimates performance

## Inputs

- ../paper2/results/low_identity_per_query.csv
- ../paper2/results/query_labels_card.csv
- ../paper2/results/db_labels_card.csv
- ../paper2/proteins/card_query_homolog.faa
- ../paper2/proteins/card_db_homolog.faa

## Outputs

- results/identity_bin_recall.csv
- results/benchmark_table.csv
- results/embedding_per_query.csv
- analysis/figures/recall_vs_identity.png
- results/split_assignments.csv (from identity-clustered split)
- results/split_summary.csv (from identity-clustered split)

## Workflow

1. Build identity-clustered split assignments (to prevent leakage):

```powershell
python projects/paper3/analysis/identity_cluster_split.py --pairwise-identity-csv projects/paper3/results/query_query_identity.csv --labels-csv projects/paper2/results/query_labels_card.csv --identity-threshold 70 --test-fraction 0.2 --out-assignments projects/paper3/results/split_assignments.csv --out-summary projects/paper3/results/split_summary.csv
```

2. Generate divergence benchmark outputs and main figure (ESM2):

```powershell
python projects/paper3/analysis/build_divergence_results.py --blast-per-query projects/paper2/results/low_identity_per_query.csv --query-labels projects/paper2/results/query_labels_card.csv --db-labels projects/paper2/results/db_labels_card.csv --query-fasta projects/paper2/proteins/card_query_homolog.faa --db-fasta projects/paper2/proteins/card_db_homolog.faa --embedding-model esm2 --hf-model-id facebook/esm2_t6_8M_UR50D --batch-size 16 --max-length 192 --identity-bin-out projects/paper3/results/identity_bin_recall.csv --benchmark-table-out projects/paper3/results/benchmark_table.csv --figure-out projects/paper3/analysis/figures/recall_vs_identity.png --embedding-per-query-out projects/paper3/results/embedding_per_query.csv
```

3. Optional real embedding run (requires torch + transformers):

```powershell
python projects/paper3/analysis/build_divergence_results.py --blast-per-query projects/paper2/results/low_identity_per_query.csv --query-labels projects/paper2/results/query_labels_card.csv --db-labels projects/paper2/results/db_labels_card.csv --query-fasta projects/paper2/proteins/card_query_homolog.faa --db-fasta projects/paper2/proteins/card_db_homolog.faa --embedding-model protbert --identity-bin-out projects/paper3/results/identity_bin_recall_protbert.csv --benchmark-table-out projects/paper3/results/benchmark_table_protbert.csv --figure-out projects/paper3/analysis/figures/recall_vs_identity_protbert.png --embedding-per-query-out projects/paper3/results/embedding_per_query_protbert.csv
```

## Figures

1. Recall vs sequence identity (BLAST line + embedding line).

## Expected signal

If divergence failure exists, BLAST recall decreases by a measurable percentage as identity drops. Aggregate-only reporting can mask these low-identity failure modes.

