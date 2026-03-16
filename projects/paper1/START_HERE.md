# Paper 1 Intern Start Guide (Detailed)

Read this file top to bottom. Copy and paste commands exactly.

## Project goal

Measure how much ARG variation is explained by:

- pipeline
- environment
- pipeline x environment interaction

## Before you start

0. Open a terminal in the repository root folder (`OPALS_Group19`).
1. Activate your local virtual environment.

Windows PowerShell:

```powershell
./.venv/Scripts/Activate.ps1
```

Linux/WSL/macOS:

```bash
source .venv/bin/activate
```

2. Confirm these files exist:
	- `projects/paper1/metadata/metadata_final.csv`
	- `projects/paper1/metadata/download_list.txt`
3. Confirm required tools are installed for this paper:
	- SRA toolkit (`prefetch`, `fastq-dump`)
	- `fastqc`
	- `megahit`
	- `prodigal`
	- `diamond`

## Step 1: Download reads

```bash
bash scripts/download_from_list.sh projects/paper1/metadata/download_list.txt projects/paper1/raw_reads
```

Success check:

- `projects/paper1/raw_reads/` contains `.fastq` files.

## Step 2: Run FastQC

```bash
bash scripts/run_fastqc_all.sh projects/paper1/raw_reads projects/paper1/results/fastqc
```

Success check:

- `projects/paper1/results/fastqc/` contains HTML reports.

## Step 3: Run Pipeline A for each sample

Use one sample first (test run):

```bash
bash scripts/run_arg_pipeline.sh SRR13853495 projects/paper1/raw_reads/SRR13853495_1.fastq projects/paper1/raw_reads/SRR13853495_2.fastq db/CARD.fasta db/CARD projects/paper1/results
```

Then run all samples by repeating the same pattern.

Success check per sample:

- `projects/paper1/results/<sample_id>/arg_hits.tsv` exists.

## Step 4: Build Pipeline A matrix

```bash
python scripts/build_arg_matrix.py --hits-glob "projects/paper1/results/*/arg_hits.tsv" --output projects/paper1/results/ARG_matrix_pipelineA.csv
```

Success check:

- `projects/paper1/results/ARG_matrix_pipelineA.csv` has one row per sample.

## Step 5: Add Pipeline B and Pipeline C matrices

Populate these files:

- `projects/paper1/results/ARG_matrix_pipelineB.csv`
- `projects/paper1/results/ARG_matrix_pipelineC.csv`

Required format for all pipeline matrices:

- First column: `sample_id`
- Remaining columns: ARG names
- Values: normalized counts (`ARG_hits / total_reads`)

## Step 6: Build long table for model

Fill `projects/paper1/results/pipeline_long_table.csv` with columns:

- `sample_id`
- `environment`
- `pipeline`
- `ARG_total`

Example rows:

| sample_id | environment | pipeline | ARG_total |
|---|---|---|---:|
| SRR13853495 | wastewater | pipelineA | 0.0042 |
| SRR13853495 | wastewater | pipelineB | 0.0059 |
| SRR13853495 | wastewater | pipelineC | 0.0038 |

## Step 7: Run variance decomposition

```bash
python projects/paper1/analysis/run_variance_decomposition.py --input projects/paper1/results/pipeline_long_table.csv --output projects/paper1/results/variance_decomposition.csv
```

Success check:

- `projects/paper1/results/variance_decomposition.csv` exists.
- It includes factors for pipeline and environment.

## Files to submit

- `projects/paper1/results/ARG_matrix_pipelineA.csv`
- `projects/paper1/results/ARG_matrix_pipelineB.csv`
- `projects/paper1/results/ARG_matrix_pipelineC.csv`
- `projects/paper1/results/pipeline_long_table.csv`
- `projects/paper1/results/variance_decomposition.csv`

## Common mistakes

- Running commands from the wrong folder.
- Missing `_1.fastq` and `_2.fastq` file names.
- Mixing raw counts and normalized counts in pipeline matrices.
- Misspelling column names in `pipeline_long_table.csv`.
