# Resistome Metagenome ARG Pipeline

This repository builds a master ENA metadata library and a core ARG abundance matrix for ~100 environmental metagenomes.

## Intern Quick Start (5 Minutes, No Build Needed)

These instructions are portable to any computer after cloning this repository.

1. Pull or clone this repository.
2. Use `metadata/master_metadata.csv` for sample metadata (already included).
3. Start inside your assigned paper folder under `projects/` and follow `START_HERE.md`.

No setup is required for this default path.

If `.venv` is not present on a new machine, run one of these first:

- Windows: `START_PHASE1.bat`
- Linux/WSL/macOS: `bash scripts/quickstart_phase1.sh`

Current repository already includes:

- `metadata/master_metadata.csv`
- `projects/paper1/` through `projects/paper4/` with initialized subfolders and metadata/download files

Refresh option (Windows, only if your mentor asks):

1. Double-click `START_PHASE1.bat`.
2. Wait for `Success. Metadata is ready`.

If you prefer terminal commands:

PowerShell (Windows):

```powershell
./scripts/quickstart_phase1.ps1
```

Bash (Linux/WSL/macOS):

```bash
bash scripts/quickstart_phase1.sh
```

See `INTERN_GUIDE.md` for plain-language instructions and troubleshooting.

Phase 1 quickstart will:

- create a local `.venv` if needed
- generate `metadata/master_metadata.csv`
- skip heavy metagenomics tools (`diamond`, `megahit`, `prodigal`) until Phase 2

## Project Folder Per Paper

Use one folder per paper to keep intern work isolated and organized.

Pre-created paper folders:

- `projects/paper1/` (Pipeline Uncertainty Decomposition)
- `projects/paper2/` (Foundation Model ARG Detection)
- `projects/paper3/` (Early Warning Signals)
- `projects/paper4/` (Wastewater Treatment Natural Experiment)

Each folder has a paper-specific `README.md` and a beginner-friendly `START_HERE.md`.

Initialize a new paper workspace (only for future/new papers):

```powershell
./scripts/init_paper_project.ps1 -ProjectName paper1 -MetadataPath metadata/master_metadata.csv -MinReads 200000
```

This creates:

- `projects/paper1/metadata/metadata_final.csv`
- `projects/paper1/metadata/download_list.txt`
- `projects/paper1/raw_reads/`
- `projects/paper1/assemblies/`
- `projects/paper1/proteins/`
- `projects/paper1/results/arg_hits/`
- `projects/paper1/analysis/figures/`

For interns starting Paper 1 now:

1. Work only inside `projects/paper1/`.
2. Start with `projects/paper1/START_HERE.md`.
3. Follow `projects/paper1/README.md` for full project context.

## Project Targets

| Environment | Target Samples |
|---|---:|
| wastewater | 40 |
| soil | 25 |
| river | 20 |
| irrigation | 15 |

## Outputs

- `metadata/master_metadata.csv`
- `results/ARG_matrix.csv`

## Research Build Status

Heavy build steps (environment/tool setup, read download, assembly, ARG calling, and matrix aggregation) are handled by maintainers, not interns.

In this workspace run:

- metadata filter step is complete (`metadata_final.csv` and `download_list.txt` were generated)
- core matrix rebuild requires ARG hit files under `results/*/arg_hits.tsv` and Linux-first tools/data not present in this repo snapshot

`master_metadata.csv` columns:

- `sample_id`
- `environment`
- `country`
- `year`
- `study`
- `timepoint`
- `treatment`
- `read_count`

## Optional Analysis Commands (After ARG Matrix Exists)

These are optional downstream commands once a real ARG matrix is available.

1. Build normalized ARG matrix and merged dataset:

```bash
python scripts/build_arg_dataset.py --metadata metadata_final.csv --arg-matrix results/ARG_matrix.csv --normalized-out results/ARG_matrix_normalized.csv --dataset-out results/ARG_dataset.csv
```

2. Generate first summary plots:

```bash
python scripts/plot_first_summary.py --dataset results/ARG_dataset.csv --outdir analysis/figures
```

3. Check dataset quality targets:

```bash
python scripts/check_dataset_quality.py --dataset results/ARG_dataset.csv --arg-matrix results/ARG_matrix.csv --report results/dataset_quality_report.txt
```

## Collaboration Guide

### Branching and commits

- Create feature branches from `main`.
- Keep commits focused: metadata curation, pipeline code, analysis notebooks, or manuscript artifacts.
- Open PRs with:
  - what changed
  - why
  - sample count impact
  - validation notes

### Data conventions

- Do not edit `metadata/master_metadata.csv` manually without documenting the reason in your PR.
- Keep raw FASTQ files out of git. Use local storage or shared data storage.
- Store intermediate outputs under `results/<sample_id>/`.

### Reproducibility checklist

- Confirm tool versions (`diamond`, `megahit`, `prodigal`, `blastn/blastp`, `fastqc`, `sra-tools`).
- Save full command lines used for production runs.
- Record CARD database version and date in PR notes.

### Suggested folder layout

```text
config/
metadata/
scripts/
results/
  <sample_id>/
    assembly/
    proteins.faa
    genes.fna
    arg_hits.tsv
```

## Notes

- On native Windows, prefer using WSL2 for heavy metagenomics tools.
- `scripts/build_master_metadata.py` is designed to bootstrap the library quickly, then you can manually curate final inclusion/exclusion before downloading all reads.
