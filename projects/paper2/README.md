# Paper 2: ARG Detection Under Sequence Divergence

For new contributors, start with `START_HERE.md` in this folder.

## Objective

Evaluate ARG detection tradeoffs under sequence divergence: BLAST has higher absolute accuracy, while embedding recall is more stable across identity bins (including `<40%` identity).

## Main Artifacts

- `projects/paper2/results/blast_vs_ml_metrics.csv`
- `projects/paper2/results/low_identity_comparison.csv`
- `projects/paper2/results/identity_bin_recall.csv`
- `projects/paper2/results/recall_gap_by_identity.csv`
- `projects/paper2/results/low_identity_per_query.csv`
- `projects/paper2/analysis/figures/embedding_pca.png`
- `projects/paper2/analysis/figures/method_comparison_bar.png`
- `projects/paper2/analysis/figures/identity_bin_recall.png`
- `projects/paper2/analysis/figures/recall_gap_vs_identity.png`

## CARD-Backed Reproducible Workflow

1. Download CARD archive payload to local file:

```powershell
python -c "import urllib.request, pathlib; p=pathlib.Path('tmp_card_latest_data.bin'); p.write_bytes(urllib.request.urlopen('https://card.mcmaster.ca/latest/data', timeout=120).read()); print(p, p.stat().st_size)"
```

2. Prepare query/db FASTA files, labels, and embeddings from CARD:

```powershell
python projects/paper2/analysis/prepare_card_benchmark_data.py `
	--card-archive tmp_card_latest_data.bin `
	--query-member ./protein_fasta_protein_homolog_model.fasta `
	--db-member ./protein_fasta_protein_homolog_model.fasta `
	--label-category "Resistance Mechanism" `
	--max-query 220 `
	--max-db 2200 `
	--min-class-count 10 `
	--out-query-fasta projects/paper2/proteins/card_query_homolog.faa `
	--out-db-fasta projects/paper2/proteins/card_db_homolog_train.faa `
	--out-query-labels projects/paper2/results/query_labels_card.csv `
	--out-db-labels projects/paper2/results/db_labels_card.csv `
	--out-embeddings projects/paper2/results/embeddings/protein_embeddings.npy `
	--out-training-labels projects/paper2/results/training_labels_template.csv
```

3. Run full benchmark pipeline (training + BLAST + figures + low-identity):

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
./projects/paper2/analysis/run_paper2_full_pipeline.ps1 `
	-PythonExe ./.venv-1/Scripts/python.exe `
	-RunLowIdentity `
	-QueryFasta projects/paper2/proteins/card_query_homolog.faa `
	-DbFasta projects/paper2/proteins/card_db_homolog_train.faa `
	-QueryLabels projects/paper2/results/query_labels_card.csv `
	-DbLabels projects/paper2/results/db_labels_card.csv
```

## Demo / Toy Run

Use demo mode only for pipeline smoke tests:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
./projects/paper2/analysis/run_paper2_full_pipeline.ps1 -UseDemoBaselineData -RunLowIdentity
```

## Notes

- The current embedding artifact in this repository is a 128-d proxy feature matrix generated from sequence k-mer hashing in `prepare_card_benchmark_data.py`.
- This run is feature-based model vs alignment, not a full protein foundation-model vs alignment comparison.
- The intended future upgrade is full protein language-model embeddings (target 1280 dimensions).

