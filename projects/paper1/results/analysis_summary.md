# Paper 1 Analysis Summary

## Variance decomposition of log10(ARG_total)

| Factor | SS | df | F | p | Variance % |
|---|---:|---:|---:|---:|---:|
| pipeline | 41.15 | 2 | 98.73 | 3.49e-17 | 57.0 |
| environment | 13.66 | 3 | 21.84 | 7.05e-09 | 18.9 |
| pipeline:environment | 7.95 | 6 | 6.36 | 6.71e-05 | 11.0 |
| residual | 9.379 | 45 |  |  | 13.0 |

- Environment explains **18.9%** of variance in log ARG abundance.
- Pipeline explains **57.0%**.
- Pipeline x environment interaction explains **11.0%**.
- Residual (sample-level) **13.0%**.

## Richness decomposition (number of ARG classes detected)

- Pipeline explains **29.6%** of richness variance (environment **24.7%**).

## Per-pipeline summary

| Pipeline | n | Mean ARG_total | Median ARG_total | Mean richness |
|---|---:|---:|---:|---:|
| pipelineA | 19 | 1.556e-06 | 3.919e-07 | 5.9 |
| pipelineB | 19 | 8.204e-05 | 2.667e-06 | 17.1 |
| pipelineC | 19 | 8.688e-05 | 5.975e-05 | 22.1 |

## Pairwise pipeline concordance

| Pair | Spearman rho | p | Median fold-difference |
|---|---:|---:|---:|
| pipelineA vs pipelineB | 0.483 | 3.64e-02 | 6.80x |
| pipelineA vs pipelineC | 0.715 | 5.84e-04 | 152.47x |
| pipelineB vs pipelineC | 0.189 | 4.39e-01 | 22.41x |
