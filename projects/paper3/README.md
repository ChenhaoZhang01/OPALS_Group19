# Paper 3: Early Warning Signals for ARG Expansion

For interns with limited coding experience: start with `START_HERE.md` in this folder.

## Objective

Test whether genomic features, especially mobile genetic elements (MGEs), rise before ARG expansion.

Hypothesis:

`MGE(t)` predicts `ARG(t+1)`

## Inputs

- results/ARG_dataset.csv
- assemblies/
- proteins/
- metadata/metadata_final.csv

## Outputs

- results/features/mge_abundance.csv
- results/features/kmer_entropy.csv
- results/features/time_ordered_table.csv
- results/models/lag_regression_results.csv
- results/models/cross_correlation.csv
- analysis/figures/

## Workflow

1. Detect mobile element markers.

Target genes include:

- transposase
- integrase
- plasmid replication proteins

Use mobileOG reference with DIAMOND:

```bash
diamond blastp -d mobileog -q proteins/SRRXXXX.faa -o results/features/SRRXXXX_mobileog.tsv
```

2. Compute sequence complexity metrics.

Python entropy example:

```python
from scipy.stats import entropy
```

Write entropy and diversity features per sample.

3. Align time series by study and sampling order.

Priority example cohort:

- PRJNA1149575

4. Run lag analysis.

Regression example:

```python
ARG_t1 ~ MGE_t
```

Also compute cross-correlation across lag values.

## Figures

1. ARG time series.
2. MGE time series.
3. Lag correlation plot.
4. Prediction performance summary.

## Expected signal

MGE abundance may increase before ARG expansion, supporting early warning use in surveillance.

