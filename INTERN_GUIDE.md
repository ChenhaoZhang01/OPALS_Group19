# Intern Quick Start (No Build Needed)

Goal: start analysis immediately with included files.

## Default path (recommended)

1. Pull or clone the repository.
2. Use `metadata/master_metadata.csv` for metadata.
3. Use `results/ARG_matrix.csv` for ARG analysis.

You do not need to run setup for this step.

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


## Rebuild option A (Windows, only if asked)

1. Open the project folder.
2. Double-click `START_PHASE1.bat`.
3. Wait for completion. You should see: `Success. Metadata is ready`.

## What this phase does

- creates `.venv` if missing
- builds `metadata/master_metadata.csv`
- does not run heavy assembly or ARG detection tools

To rebuild `results/ARG_matrix.csv`, your team must first generate sample ARG hit files (`results/*/arg_hits.tsv`) via the full pipeline.

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
