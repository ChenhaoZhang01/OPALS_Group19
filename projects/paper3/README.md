# Paper 3: Early Warning Signals for ARG Expansion

For interns with limited coding experience: start with `START_HERE.md` in this folder.

## Objective

Evaluate whether MGEs contain predictive information for future ARG change while controlling for temporal confounding.

Core contribution: naive lag-based predictive models can produce near-perfect fit under strong temporal coupling, even when independent predictive signal is minimal.

Primary hypothesis:

`MGE(t)` predicts `ARG(t+1)`

Robustness hypotheses:

- First-difference signal: `dMGE(t)` predicts `dARG(t+1)`
- Added-value signal: `MGE(t)` improves prediction beyond `ARG(t)` history

Despite this being known in time-series statistics, these checks are rarely made explicit in environmental AMR lag-analysis workflows.

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
- results/difference_model.csv
- results/granger_test.csv
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

4. Run lag analysis with robustness tests.

Regression example:

```python
ARG_t1 ~ MGE_t
```

Also compute:

- first-difference model: `dARG(t+1) ~ dMGE(t)`
- Granger-style comparison:
	- base: `ARG(t+1) ~ ARG(t)`
	- full: `ARG(t+1) ~ ARG(t) + MGE(t)`
- cross-correlation across lag values

Model outputs now include both `r_squared` and `adjusted_r_squared`.
Sequencing depth is modeled on the raw read-count scale when available (`sequencing_depth_transform = raw_read_count`).

## Figures

1. ARG time series.
2. MGE time series.
3. Lag correlation plot.
4. Differenced scatter plot (`dMGE(t)` vs `dARG(t+1)`).

## Expected signal

Naive lag associations can appear strong under shared temporal dynamics. A valid leading-indicator claim should persist after differencing and history-conditioned testing.

