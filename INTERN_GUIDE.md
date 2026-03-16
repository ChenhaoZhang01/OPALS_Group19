# Intern Quick Start (No Build Needed)

Goal: start your assigned paper safely on a new computer.

## Default path (recommended)

1. Pull or clone the repository.
2. Use `metadata/master_metadata.csv` for metadata.
3. Open your assigned paper folder under `projects/`.
4. Follow that folder's `START_HERE.md`.

## Start a New Paper Project (example: Paper 1)

Run this once from the repository root:

```powershell
./scripts/init_paper_project.ps1 -ProjectName paper1 -MetadataPath metadata/master_metadata.csv -MinReads 200000
```

Then:

1. Open `projects/paper1/`.
2. Use `projects/paper1/metadata/download_list.txt` for SRA download.
3. Save all outputs only under `projects/paper1/`.
4. Follow `projects/paper1/README.md` checklist.

## Paper 1 on a New Computer

1. Clone and enter the repository.

```powershell
git clone <repo-url>
cd OPALS_Group19
```

2. Create local Python environment and metadata bootstrap.

Windows:

```powershell
./START_PHASE1.bat
```

Linux/WSL/macOS:

```bash
bash scripts/quickstart_phase1.sh
```

3. Activate local environment.

Windows PowerShell:

```powershell
./.venv/Scripts/Activate.ps1
```

Linux/WSL/macOS:

```bash
source .venv/bin/activate
```

4. Verify Paper 1 starter files exist.

```powershell
Test-Path projects/paper1/metadata/metadata_final.csv
Test-Path projects/paper1/metadata/download_list.txt
```

Both should return `True`.

5. Start Paper 1 workflow.

```text
Open projects/paper1/START_HERE.md and run steps in order.
```

6. Run one sample first (from repository root) before any full batch run.

```bash
bash scripts/run_arg_pipeline.sh SRR13853495 projects/paper1/raw_reads/SRR13853495_1.fastq projects/paper1/raw_reads/SRR13853495_2.fastq db/CARD.fasta db/CARD projects/paper1/results
```

7. Build Paper 1 matrix after arg_hits files exist.

```bash
python scripts/build_arg_matrix.py --hits-glob "projects/paper1/results/*/arg_hits.tsv" --output projects/paper1/results/ARG_matrix_pipelineA.csv
```

If your machine is missing `megahit`, `prodigal`, `diamond`, `prefetch`, or `fastq-dump`, install the pipeline tools yourself in WSL/Linux:

1. Install WSL (Windows only, one time):

```powershell
wsl --install
```

2. Open your Linux shell (Ubuntu), go to the repository, and create the conda environment from this repo file:

```bash
conda env create -f environment-linux.yml
conda activate resistome
```

3. Verify required tools:

```bash
which prefetch fastq-dump fastqc megahit prodigal diamond
```

4. If all tools print paths, re-run Step 6 and Step 7 from this guide in the Linux shell.

Note: if `db/CARD.fasta` is missing, download the CARD protein FASTA and place it at `db/CARD.fasta`, then re-run the pipeline command.


## Rebuild option A (Windows, only if asked)

1. Open the project folder.
2. Double-click `START_PHASE1.bat`.
3. Wait for completion. You should see: `Success. Metadata is ready`.

## What this phase does

- creates `.venv` if missing
- builds `metadata/master_metadata.csv`
- does not run heavy assembly or ARG detection tools

## Common issues

- Python not found:
  - install Python 3.10+
  - rerun `START_PHASE1.bat`
- ENA network hiccup:
  - rerun the same command
- Not enough irrigation samples in one run:
  - rerun with a higher query limit:

```powershell
./scripts/quickstart_phase1.ps1 -PerQueryLimit 10000
```
