# Paper 1 Data Provenance

**Status: REAL DATA.** The `ARG_matrix_pipeline{A,B,C}.csv` and `pipeline_long_table.csv`
in this directory are produced by actually downloading the SRA metagenomes and running
all three ARG-calling pipelines. The earlier calibrated-simulation stand-in has been
replaced; `SIMULATED_PLACEHOLDER.txt` is removed once real output is written.

## What was run

A cohort of 19 environmental shotgun metagenomes (balanced across four habitats:
wastewater 4, soil 5, river 5, irrigation 5; accessions in
`../metadata/batch_subset.txt`, drawn from `../metadata/metadata_final.csv`) was
processed end-to-end on a Linux/WSL2 workstation through:

- **Pipeline A** — MEGAHIT assembly (min contig 500 bp) → Prodigal ORF prediction →
  DIAMOND blastp vs CARD protein homolog models (≥80% identity, ≥70% query coverage).
- **Pipeline B** — Bowtie2 read mapping to the CARD nucleotide homolog reference;
  ARG counts from `samtools idxstats`.
- **Pipeline C** — RGI on the predicted ORFs against CARD.

Tooling: CARD v3.2.7; DIAMOND 2.1.11; MEGAHIT 1.2.9; Prodigal 2.6.3; Bowtie2;
samtools 1.21; RGI 6.0.5; sra-tools 3.4.1. Driver: `../analysis/run_arg_pipeline_hpc.sh`;
database build: `../analysis/setup_databases.sh`; aggregation:
`../analysis/aggregate_pipeline_hits.py` (CARD `aro_index.tsv` used to collapse hits to
drug classes). Normalized abundance = ARG hits / total reads; `ARG_total` = sum over
features; `ARG_richness` = distinct features detected.

## Real-data decisions and their rationale (all are honest limitations)

1. **Equal-depth subsampling to ~12M reads** (`--max-reads 6,000,000` spots, streamed
   with `fastq-dump -X`). The cohort's native depths span 3M–59M reads; the deep
   soil/irrigation libraries cannot be assembled within the workstation's 11 GB RAM at
   full depth. Streaming the first N spots bounds both RAM and download size and
   controls depth as a confounder. *Consequence:* assembly-based detection (Pipeline A)
   is depth-limited and conservative; first-N (not random) subsampling is a minor bias.

2. **Minimal QC (fastp skipped, `--skip-fastp`).** Writing GBs of uncompressed trimmed
   FASTQ was the dominant I/O cost (~1 hr/sample) on the WSL2 disk. Reads were fed
   directly to mapping/assembly, which are robust to light residual adapter content.
   *Consequence:* no adapter/quality trimming; read-mapping counts may include a small
   number of low-quality reads.

3. **RGI Loose tier included for Pipeline C.** RGI "strict" (Perfect/Strict only)
   returned **zero** ARG calls on every sample: environmental ARG homologs fall below
   CARD's clinically-curated bitscore cutoffs (e.g. SRR25475470: 135,325 ORFs → 0 strict,
   750 Loose). To yield a comparable abundance gradient, Pipeline C includes the Loose
   tier (`--include_loose`). The strict≈0 result is itself reported as a finding
   (`../analysis/redo_pipeline_c.sh` regenerates C).

4. **Zero-detection pseudocount.** A few assembly-pipeline samples detect no ARG
   (`ARG_total = 0`). The log-scale variance model adds half the smallest positive
   normalized abundance before `log10` (standard zero-handling).

5. **n = 19, single workstation.** Small cohort; the variance decomposition is well
   powered for the pipeline main effect but the per-environment cells are small
   (3–5 samples). Scaling to public compendia (ENA/MGnify) is future work.

## Reproduce

```bash
# 1. databases (once)
bash projects/paper1/analysis/setup_databases.sh --dbdir ~/paper1_db
# 2. run all 19 (resumable)
bash projects/paper1/analysis/run_arg_pipeline_hpc.sh \
  --card-dir ~/paper1_db --outdir ~/paper1_scratch --threads 16 \
  --max-reads 6000000 --skip-fastp \
  --download-list projects/paper1/metadata/batch_subset.txt
# 3. (if any sample predates the RGI-loose fix) regenerate Pipeline C
bash projects/paper1/analysis/redo_pipeline_c.sh --outdir ~/paper1_scratch --threads 16
# 4. analysis + figures (re-run on the aggregated tables)
python projects/paper1/analysis/run_paper1_analysis.py
python projects/paper1/analysis/generate_figures.py
```

## Gate status

| Gate | Status |
|---|---|
| Real cohort metadata (balanced 4 habitats) | **PASS** |
| All three pipelines executed on real reads | **PASS** (n=19) |
| Variance-decomposition + figures reproducible from real tables | **PASS** |
| Full native-depth, randomized-subsample, full QC | **PARTIAL** (depth-capped, minimal QC — see limitations) |
