# Resistome Metagenome ARG Pipeline

This repository builds a master ENA metadata library and a core ARG abundance matrix for ~100 environmental metagenomes.

## Intern Quick Start (5 Minutes, No Build Needed)

1. Pull or clone this repository.
2. Use `metadata/master_metadata.csv` for sample metadata (already included).
3. Use `results/ARG_matrix.csv` for ARG abundance analysis (already included).
4. Open one project folder under `projects/` and follow its `START_HERE.md`.

No setup is required for this default path.

Current repository already includes:

- `metadata/master_metadata.csv`
- `results/ARG_matrix.csv`
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

To rebuild `results/ARG_matrix.csv`, use Section 5 below after ARG hit files are available.

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

`master_metadata.csv` columns:

- `sample_id`
- `environment`
- `country`
- `year`
- `study`
- `timepoint`
- `treatment`
- `read_count`

## 1 Full Environment Setup (Phase 2+)

### Conda environment

```powershell
conda create -n resistome python=3.10 -y
conda activate resistome
pip install pandas numpy scikit-learn matplotlib seaborn biopython
```

### Bioinformatics tools

`fastqc` is available in conda on Windows, but `diamond`, `prodigal`, and `megahit` are typically Linux-first in bioconda. Recommended approach: run the pipeline in WSL2 (Ubuntu) or Linux.

Linux/WSL install command:

```bash
conda install -n resistome -c conda-forge -c bioconda diamond ncbi-blast prodigal megahit fastqc sra-tools -y
```

Or create directly from the provided environment file:

```bash
conda env create -f environment-linux.yml
conda activate resistome
```

If FastQC install fails with `InvalidArchiveError` on Windows:

```powershell
conda clean --packages --tarballs --yes
conda install -n resistome -c conda-forge -c bioconda fastqc -y
```

## 2 Build the Master Metadata Library

Create/refresh metadata from ENA:

```bash
python scripts/build_master_metadata.py --output metadata/master_metadata.csv
```

The script tries to fill all targets by scanning ENA `library_source="METAGENOMIC"` runs and classifying records into target environments using keywords.

After generation, manually review and clean edge cases (environment mislabels, unknown timepoints, treatment labels).

## 3 Download Reads

From `master_metadata.csv`, pull SRA runs and convert to FASTQ:

```bash
bash scripts/download_reads.sh metadata/master_metadata.csv data/raw
```

## 4 Per-Sample ARG Pipeline

Pipeline steps:

- reads -> assembly (`megahit`)
- contigs -> genes/proteins (`prodigal`)
- proteins -> ARG hits (`diamond blastp` against CARD)

Per-sample command:

```bash
bash scripts/run_arg_pipeline.sh SRR12345 data/raw/SRR12345_1.fastq data/raw/SRR12345_2.fastq db/CARD.fasta db/CARD results
```

Single-end example:

```bash
bash scripts/run_arg_pipeline.sh SRR99999 data/raw/SRR99999.fastq NA db/CARD.fasta db/CARD results
```

## 5 Build Core ARG Matrix

Aggregate all sample hit tables into one matrix:

```bash
python scripts/build_arg_matrix.py --hits-glob "results/*/arg_hits.tsv" --output results/ARG_matrix.csv
```

Or run the full batch workflow (all samples + matrix generation):

```bash
bash scripts/run_all_samples.sh metadata/master_metadata.csv data/raw db/CARD.fasta db/CARD results
```

Matrix format:

| Sample | ARG1 | ARG2 | ARG3 | ARG4 | ... |
|---|---:|---:|---:|---:|---:|

## End-to-End Workflow Commands (Research Phase)

These are optional rebuild commands for generating new outputs. For interns, use the existing files first unless instructed otherwise.

1. Clean metadata and keep samples with `read_count >= 200000`:

```powershell
./scripts/prepare_metadata_for_download.ps1 -MetadataPath metadata/master_metadata.csv -Output metadata_final.csv -DownloadList download_list.txt -MinReads 200000
```

2. Download reads from accession list (Linux/WSL):

```bash
bash scripts/download_from_list.sh download_list.txt data/raw
```

3. Run FastQC on downloaded reads:

```bash
bash scripts/run_fastqc_all.sh data/raw results/fastqc
```

4. Assemble + gene prediction + ARG detection per sample using the existing pipeline script:

```bash
bash scripts/run_arg_pipeline.sh SRR12345 data/raw/SRR12345_1.fastq data/raw/SRR12345_2.fastq db/CARD.fasta db/CARD results
```

5. Build ARG matrix from all sample hit files:

```bash
python scripts/build_arg_matrix.py --hits-glob "results/*/arg_hits.tsv" --output results/ARG_matrix.csv
```

6. Build normalized ARG matrix and merged dataset:

```bash
python scripts/build_arg_dataset.py --metadata metadata_final.csv --arg-matrix results/ARG_matrix.csv --normalized-out results/ARG_matrix_normalized.csv --dataset-out results/ARG_dataset.csv
```

7. Generate first summary plots:

```bash
python scripts/plot_first_summary.py --dataset results/ARG_dataset.csv --outdir analysis/figures
```

8. Check dataset quality targets:

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
