# Paper 4 (corrected) — ARG & class-1 integron removal across quaternary lines

Data: projects/paper4/data/external/oa_scan/PMC12735332/PMC12735332/supp_unzipped/Supplementary-tables-CEBEDEAU.xlsx (Table S10); campaigns ['C10-I', 'C12-I', 'C2-I', 'C5-I']; genes ['blaAmpC', 'ermB', 'sul1', 'sul2', 'tetW', 'intI1'].

Positive log10 values = reduction (line has fewer copies than its OUT WWTP input).

## Per-line mean removal across ARGs (log10)

| Line | Mean ARG removal (abs) | Mean ARG removal (16S-norm) | intI1 removal (abs) | intI1 removal (norm) |
|---|---:|---:|---:|---:|
| Constructed wetland | 2.00 | 0.55 | 2.26 | 0.81 |
| Ozonation + GAC | 1.52 | 1.11 | 0.46 | 0.05 |
| GAC | 0.95 | 0.30 | 0.94 | 0.29 |

## Removal by line × gene (absolute log10)

| gene    |   CW |   GAC |   AOP |
|:--------|-----:|------:|------:|
| blaAmpC | 1.83 |  1.02 |  1.11 |
| ermB    | 2.28 |  1.44 |  2.6  |
| intI1   | 2.26 |  0.94 |  0.46 |
| sul1    | 2.06 |  0.79 |  1.33 |
| sul2    | 1.63 |  0.75 |  0.31 |
| tetW    | 2.21 |  0.73 |  2.23 |

## 16S-normalized removal by line × gene (log10)

| gene    |   CW |   GAC |   AOP |
|:--------|-----:|------:|------:|
| blaAmpC | 0.38 |  0.38 |  0.7  |
| ermB    | 0.83 |  0.79 |  2.19 |
| intI1   | 0.81 |  0.29 |  0.05 |
| sul1    | 0.62 |  0.14 |  0.92 |
| sul2    | 0.18 |  0.11 | -0.1  |
| tetW    | 0.76 |  0.08 |  1.82 |

## Conventional plant removal (IN WWTP -> OUT WWTP, absolute log10)

| gene    |   in_out_log_removal_mean |   n |
|:--------|--------------------------:|----:|
| blaAmpC |                      2.93 |   2 |
| ermB    |                      3.53 |   2 |
| sul1    |                      2.5  |   2 |
| sul2    |                      2.68 |   2 |
| tetW    |                      3.11 |   2 |
| intI1   |                      2.13 |   2 |

## Mechanism (biomass vs selective)

biomass_vs_selective_gap = abs removal - 16S-normalized removal; large positive = removal driven mostly by biomass loss, not selective ARG removal.

| line   |   arg_abs_removal |   intI1_abs_removal |   arg_norm_removal |   intI1_norm_removal |   biomass_vs_selective_gap |
|:-------|------------------:|--------------------:|-------------------:|---------------------:|---------------------------:|
| CW     |              2    |                2.26 |               0.55 |                 0.81 |                       1.45 |
| GAC    |              0.95 |                0.94 |               0.3  |                 0.29 |                       0.65 |
| AOP    |              1.52 |                0.46 |               1.11 |                 0.05 |                       0.41 |