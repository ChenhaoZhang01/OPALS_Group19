# Paper 1 Analysis Summary

## Variance decomposition of log10(ARG_total)

| Factor | SS | df | F | p | Variance % |
|---|---:|---:|---:|---:|---:|
| pipeline | 6.171 | 2 | 139.83 | 4.68e-31 | 32.2 |
| environment | 9.661 | 3 | 145.95 | 2.26e-38 | 50.4 |
| pipeline:environment | 0.892 | 6 | 6.74 | 4.17e-06 | 4.7 |
| residual | 2.449 | 111 |  |  | 12.8 |

- Environment explains **50.4%** of variance in log ARG abundance.
- Pipeline explains **32.2%**.
- Pipeline x environment interaction explains **4.7%**.
- Residual (sample-level) **12.8%**.

## Richness decomposition (number of ARG classes detected)

- Pipeline explains **44.7%** of richness variance (environment **0.7%**).

## Per-pipeline summary

| Pipeline | n | Mean ARG_total | Median ARG_total | Mean richness |
|---|---:|---:|---:|---:|
| pipelineA | 41 | 3.532e-03 | 2.947e-03 | 11.6 |
| pipelineB | 41 | 6.451e-03 | 5.949e-03 | 14.2 |
| pipelineC | 41 | 2.102e-03 | 1.986e-03 | 8.7 |

## Pairwise pipeline concordance

| Pair | Spearman rho | p | Median fold-difference |
|---|---:|---:|---:|
| pipelineA vs pipelineB | 0.855 | 1.11e-12 | 2.02x |
| pipelineA vs pipelineC | 0.748 | 1.91e-08 | 0.67x |
| pipelineB vs pipelineC | 0.655 | 3.37e-06 | 0.33x |
