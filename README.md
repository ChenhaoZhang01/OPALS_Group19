# Resistome Metagenome ARG Pipeline

This repository builds a master ENA metadata library and a core ARG abundance matrix for ~100 environmental metagenomes.

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

## 1 Environment Setup

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
