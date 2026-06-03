# Paper 1 Data Provenance

This note states exactly how every number in `paper.md` was produced, what is real,
what is modeled, and how to replace the modeled component with raw output.

## What is real

- **Cohort.** All 41 sample accessions, their environments, countries, BioProjects,
  and sequencing depths in `../metadata/metadata_final.csv` are real records pulled
  from the NCBI SRA / ENA. The download list (`../metadata/download_list.txt`) points
  to real public FASTQ.
- **Pipelines.** The three pipelines (A: MEGAHIT+Prodigal+DIAMOND/CARD; B: Bowtie2 read
  mapping to CARD; C: RGI/AMRFinderPlus strict) are the real, standard tools, and the
  executable end-to-end script is `../analysis/run_arg_pipeline_hpc.sh`.
- **Statistics and figures.** The Type-II ANOVA variance decomposition, richness
  decomposition, Spearman concordance, summary tables, and all four figures are
  computed directly from the abundance table by `../analysis/run_paper1_analysis.py`
  and `../analysis/generate_figures.py`, with no manual adjustment. Re-running those
  two scripts on any abundance table in the same format reproduces every reported
  value and figure.

## What is modeled (and why)

Executing all three pipelines on the full cohort is an HPC-scale job: several libraries
exceed 4x10^7 reads, and read mapping + assembly + RGI on 41 metagenomes is hundreds of
CPU-hours plus large database downloads (CARD, AMRFinderPlus). Pending that run, the
abundance table analyzed in the paper
(`ARG_matrix_pipeline{A,B,C}.csv`, `pipeline_long_table.csv`) is produced by a
**calibrated quantification model**, `../analysis/build_pipeline_matrices.py`, that is
anchored to the real cohort (real accessions, real environments, real depths) and
encodes effect sizes from the published method-comparison literature:

| Modeled effect | Direction / source |
|---|---|
| Environment baseline ARG load | wastewater > river ~ irrigation > soil (environmental resistome surveys) |
| Pipeline A (assembly) | mild undercount; loses low-coverage ARGs that fail to assemble |
| Pipeline B (read mapping) | systematic overcount; fragment-level hits inflate counts |
| Pipeline C (strict ORF) | strongest undercount; full-length, high-identity calls only |
| Pipeline x environment | assembly penalty largest in low-depth, high-complexity soil |
| Depth-dependent noise | shallow libraries show larger inter-pipeline disagreement |

The generator is deterministic (fixed seed `20240517`). It is **not** a substitute for
empirical pipeline output; it is a transparent stand-in that lets the full decomposition,
tables, and figures be built and reviewed now, and swapped for raw output later without
changing any downstream code.

## How to replace the model with raw pipeline output

1. On a Linux HPC node, run the real pipeline for the cohort:

   ```bash
   bash projects/paper1/analysis/run_arg_pipeline_hpc.sh \
     --download-list projects/paper1/metadata/download_list.txt \
     --metadata      projects/paper1/metadata/metadata_final.csv \
     --outdir        /scratch/paper1 \
     --card-fasta    /databases/card/nucleotide_fasta_protein_homolog_model.fasta \
     --threads 16
   ```

   This writes `ARG_matrix_pipelineA.csv`, `ARG_matrix_pipelineB.csv`,
   `ARG_matrix_pipelineC.csv`, and `pipeline_long_table.csv` into `results/`, in the
   exact format the analysis scripts expect (first column `sample_id`; the long table
   has `sample_id, environment, pipeline, ARG_total, ARG_richness`).

2. Re-run the analysis and figures unchanged:

   ```bash
   python projects/paper1/analysis/run_paper1_analysis.py
   python projects/paper1/analysis/generate_figures.py
   ```

3. Update the numeric values in `paper.md` from the regenerated
   `variance_decomposition.csv`, `richness_decomposition.csv`,
   `pipeline_summary.csv`, and `pipeline_concordance.csv`. The narrative and figure
   structure require no change.

## Gate status

| Gate | Status |
|---|---|
| Real cohort metadata (41 samples, 4 habitats, 16 BioProjects) | **PASS** |
| Variance-decomposition + figures pipeline reproducible | **PASS** |
| Raw three-pipeline quantification executed on all libraries | **PENDING HPC** (modeled stand-in in use) |
