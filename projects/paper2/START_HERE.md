# Paper 2 Intern Start Guide (Detailed)

Read this file top to bottom. Copy and paste commands exactly.

## Project goal

Test whether embedding-based machine learning detects ARGs better than alignment.

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

1. Understand what a protein FASTA file is.

Protein FASTA files usually end with `.faa`. They contain protein sequences in text format.

Each protein entry looks like this:

```text
>protein_name_or_id
MKTIIALSYIFCLVFADYKDDDDK
```

Important:

- Each protein starts with a line beginning with `>`.
- The next lines are amino-acid letters (`A, C, D, E, ...`).

2. Put protein FASTA files in `projects/paper2/proteins/`.

Option A (fastest): copy existing `.faa` files from another project.

```powershell
Copy-Item projects/paper1/proteins/*.faa projects/paper2/proteins/ -Force
```

Option B (if you do not already have `.faa` files): generate one from an assembly with Prodigal.

```bash
prodigal -i projects/paper1/assemblies/SRR13853495/final.contigs.fa -a projects/paper2/proteins/SRR13853495.faa -p meta
```

3. Check that at least one `.faa` file exists.

```powershell
Get-ChildItem projects/paper2/proteins/*.faa
```

4. Make sure `projects/paper2/results/embeddings/` exists.

```powershell
New-Item -ItemType Directory -Force -Path projects/paper2/results/embeddings | Out-Null
```

5. Install required Python packages (one time):

```powershell
python -m pip install fair-esm numpy pandas scikit-learn
```

## Step 1: Generate embeddings

Create or run your embedding extraction workflow and save output to:

- `projects/paper2/results/embeddings/protein_embeddings.npy`

Success check:

- File exists and is not empty.

Quick check command:

```powershell
Get-Item projects/paper2/results/embeddings/protein_embeddings.npy
```

## Step 2: Label training rows

Edit:

- `projects/paper2/results/training_labels_template.csv`

Required columns:

- `row_index` (integer index in embeddings array)
- `label` (ARG class)

Add at least 10 labeled rows before training.

## Step 3: Train baseline classifier

```powershell
python projects/paper2/analysis/train_arg_classifier.py --embeddings projects/paper2/results/embeddings/protein_embeddings.npy --labels projects/paper2/results/training_labels_template.csv --metrics-out projects/paper2/results/blast_vs_ml_metrics.csv
```

Success check:

- `projects/paper2/results/blast_vs_ml_metrics.csv` is created.

## Step 4: Compare against alignment baseline

Run BLAST/DIAMOND workflow and record baseline metrics in:

- `projects/paper2/results/blast_vs_ml_metrics.csv`

Columns:

- `method`
- `precision`
- `recall`
- `f1`
- `roc_auc`

## Step 5: Make figures

Produce and save under `projects/paper2/analysis/figures/`:

1. Embedding clustering plot (UMAP).
2. ROC curve.
3. BLAST vs ML comparison.
4. Novel predicted ARG clusters.

## Files to submit

- `projects/paper2/results/embeddings/protein_embeddings.npy`
- `projects/paper2/results/training_labels_template.csv`
- `projects/paper2/results/blast_vs_ml_metrics.csv`

## Common mistakes

- Using row indexes that do not exist in the `.npy` array.
- Training with too few labeled rows.
- Mixing class names (for example `beta_lactam` and `beta-lactam`) in the same label column.
