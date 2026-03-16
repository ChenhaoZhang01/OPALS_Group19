# Paper 1: Pipeline Uncertainty Decomposition in ARG Ecology

For interns with limited coding experience: start with `START_HERE.md` in this folder.

## Objective

Quantify variation in ARG abundance explained by:

- pipeline
- environment
- pipeline x environment interaction

Target model:

`ARG_total ~ pipeline + environment + pipeline:environment`

## Inputs

- metadata/metadata_final.csv
- metadata/download_list.txt
- raw_reads/

## Outputs

- results/ARG_matrix_pipelineA.csv
- results/ARG_matrix_pipelineB.csv
- results/ARG_matrix_pipelineC.csv
- results/pipeline_long_table.csv
- results/variance_decomposition.csv
- analysis/figures/

## Workflow

1. Run Pipeline A (assembly + prodigal + DIAMOND).

Assembly example:

```bash
megahit -1 raw_reads/SRRXXXX_1.fastq -2 raw_reads/SRRXXXX_2.fastq -o assemblies/SRRXXXX
```

Gene prediction example:

```bash
prodigal -i assemblies/SRRXXXX/final.contigs.fa -a proteins/SRRXXXX.faa -p meta
```

ARG search example:

```bash
diamond blastp -d card -q proteins/SRRXXXX.faa -o results/arg_hits/SRRXXXX_pipelineA.tsv
```

2. Run Pipeline B (read mapping to CARD).

```bash
bowtie2-build card.fasta card_index
bowtie2 -x card_index -1 raw_reads/SRRXXXX_1.fastq -2 raw_reads/SRRXXXX_2.fastq -S results/SRRXXXX_pipelineB.sam
```

Convert each SAM to ARG counts and build `results/ARG_matrix_pipelineB.csv`.

3. Run Pipeline C (AMRFinderPlus or RGI).

RGI example:

```bash
rgi main -i proteins/SRRXXXX.faa -o results/SRRXXXX_pipelineC
```

4. Build three normalized ARG matrices.

Normalization:

`normalized_ARG = ARG_hits / total_reads`

5. Combine to long format table:

| sample_id | environment | pipeline | ARG_total |
|---|---|---|---:|

6. Fit variance model in Python.

```python
import statsmodels.formula.api as smf

model = smf.ols(
	"ARG_total ~ pipeline + environment + pipeline:environment",
	data=df,
).fit()
```

7. Export variance explained per factor to `results/variance_decomposition.csv`.

## Figures

1. ARG abundance across pipelines.
2. Variance partition bar chart.
3. Pipeline x environment interaction plot.
4. ARG richness comparison.

## Expected signal

- pipeline explains about 20 to 40 percent.
- environment explains about 40 to 60 percent.

