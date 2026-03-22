# Evaluating Predictive Signal Versus Temporal Coupling in Resistome-Mobilome Dynamics

## Abstract
Environmental ARG surveillance often seeks leading indicators that could warn of future resistome expansion. We evaluated whether mobile genetic element (MGE) abundance at time t predicts ARG burden at time t+1, while explicitly testing whether apparent prediction is confounded by shared temporal structure. In this repository run (33 lag pairs), a naive forward lag model showed strong association (ARG(t+1) ~ MGE(t) + C(study) + sequencing_depth; coefficient 18.554, p < 1e-8, R-squared = 0.998). However, lag correlations were uniformly high across offsets (t to t: 0.997; t to t+1: 0.995; t to t+2: 0.993), and reverse-direction modeling was similarly strong (ARG(t) to MGE(t+1), R-squared = 0.999). After trend-removing differencing, the effect weakened (dARG(t+1) ~ dMGE(t): coefficient 4.545, p = 0.0616, R-squared = 0.119). A Granger-style added-value test further showed that MGE does not significantly improve prediction beyond ARG history in this dataset (delta R-squared = 0.000027; p = 0.372). These findings support a methodological contribution: naive lag models can overstate predictive interpretation under strong temporal coupling, and robust surveillance inference requires differencing and history-conditioned tests.

Naive lag-based predictive models can produce near-perfect fit (R-squared approximately 1) even when no independent predictive signal exists, due to strong temporal coupling and autocorrelation in resistome dynamics.

## Introduction
Mobilome-resistome coupling is biologically plausible because MGEs mediate horizontal transfer and can co-vary with resistance burden. This motivates a key surveillance question: does MGE variation provide predictive information about future ARG change, beyond contemporaneous co-movement?

Simple lag models are attractive but vulnerable to confounding by shared drift, autocorrelation, and synchronized study-level trajectories. Strong fit alone does not establish directional influence.

This study therefore evaluates predictive claims under increasing rigor: (1) lag regression, (2) reverse-direction check, (3) first-difference regression, and (4) Granger-style added-value testing.

While temporal autocorrelation is well understood in statistical theory, its consequences for predictive interpretation are rarely operationalized in resistome studies, where high model fit is frequently interpreted as evidence of leading indicators without formal diagnostic testing.

## Materials and Methods
### Data and preprocessing
Input was projects/paper3/results/features_table.csv with columns study, sample_id, order, mge_abundance, entropy, arg_total, and sequencing_depth. Data were sorted by study and order. Lag features were built within study:

- ARG(t+1) from arg_total shifted by -1
- ARG(t+2) from arg_total shifted by -2
- MGE(t+1) from mge_abundance shifted by -1

First-difference variables were constructed as:

- dMGE(t) = MGE(t) - MGE(t-1)
- dARG(t+1) = ARG(t+1) - ARG(t)

### Models
Primary lag model:

ARG(t+1) ~ MGE(t) + C(study) + sequencing_depth

Sequencing depth was modeled on the raw read-count scale (no log transform) when available.

Future work will assess robustness to log-transformed depth and alternative scaling, although the magnitude of the observed collapse suggests that conclusions are unlikely to depend on this choice.

Entropy-augmented model:

ARG(t+1) ~ MGE(t) + entropy(t) + C(study) + sequencing_depth

Reverse-direction model:

MGE(t+1) ~ ARG(t) + C(study) + sequencing_depth

First-difference model:

dARG(t+1) ~ dMGE(t)

Granger-style comparison:

- Base: ARG(t+1) ~ ARG(t) + C(study) + sequencing_depth
- Full: ARG(t+1) ~ ARG(t) + MGE(t) + C(study) + sequencing_depth

Nested model improvement was assessed with F-test and delta R-squared.

Pearson correlations were computed for t to t, t to t+1, and t to t+2.

### Outputs
The pipeline writes:

- projects/paper3/results/lag_regression_results.csv
- projects/paper3/results/model_comparison.csv
- projects/paper3/results/directionality_test.csv
- projects/paper3/results/cross_correlation.csv
- projects/paper3/results/difference_model.csv
- projects/paper3/results/granger_test.csv

## Results
### Step 1: Naive lag model shows extremely strong signal

Forward model:

- coefficient: 18.553611
- standard error: 0.589673
- p-value: < 1e-8
- 95% CI: [17.345719, 19.761502]
- R-squared: 0.998303
- adjusted R-squared: 0.998060
- n: 33

At face value, this result would suggest an exceptionally strong and potentially actionable predictive relationship.

### Step 2: Uniformly high lag correlations raise suspicion

| lag | correlation | n |
|---|---:|---:|
| t to t | 0.996761 | 36 |
| t to t+1 | 0.995173 | 33 |
| t to t+2 | 0.992825 | 30 |

Despite high model fit, lag correlations remain high across multiple time offsets, suggesting strong temporal coupling.

These results demonstrate that lagged association, even when statistically significant and near-perfect in fit, is insufficient to establish predictive or causal direction.

### Supporting check: Directionality remains symmetric

| direction | coefficient | p-value | R-squared | n |
|---|---:|---:|---:|---:|
| MGE(t) to ARG(t+1) | 18.554 | < 1e-8 | 0.998 | 33 |
| ARG(t) to MGE(t+1) | 0.054 | < 1e-8 | 0.999 | 33 |

The symmetry of forward and reverse models indicates that lagged association alone is insufficient to establish directional influence.

### Step 3: After removing temporal trends, predictive signal weakens

| model | coefficient | p-value | R-squared | adjusted R-squared | n |
|---|---:|---:|---:|---:|---:|
| dARG(t+1) ~ dMGE(t) | 4.545 | 0.061649 | 0.119225 | 0.087769 | 30 |

After removing shared temporal trends, the apparent predictive signal collapses, with R-squared decreasing from 0.998 to 0.119 and statistical significance disappearing, indicating that the original signal was largely driven by shared temporal structure.

### Step 4: Final test shows no added predictive value beyond ARG history

| metric | value |
|---|---:|
| Base R-squared | 0.999097 |
| Base adjusted R-squared | 0.998720 |
| Full R-squared | 0.999124 |
| Full adjusted R-squared | 0.998731 |
| Delta R-squared | 0.000027 |
| Added-value F | 0.823097 |
| Added-value p-value | 0.372303 |

In this run, MGE does not add significant predictive value beyond ARG history.

### Figures
Generated figures in projects/paper3/analysis/figures/:

1. scatter_mge_t_vs_arg_t1.png
2. timeseries_mge_arg_by_study.png
3. lag_comparison_correlation.png
4. differenced_scatter_dmge_t_vs_darg_t1.png

## Discussion
The key contribution is not a claim of confirmed early warning signal, but a framework for evaluating predictive claims under confounding-prone longitudinal structure.

Naive lag models in this dataset produce extremely high fit, but differencing and history-conditioned tests indicate that much of the apparent predictability is explained by shared temporal dynamics. This is a stronger and more generalizable conclusion for surveillance-method development.

These findings highlight that predictive surveillance in resistome systems requires disentangling true leading indicators from shared temporal dynamics, which can otherwise produce misleadingly strong lag correlations.

These findings indicate that predictive surveillance in resistome systems may require fundamentally different signals, such as mechanistic markers or intervention-driven perturbations, rather than relying on endogenous temporal dynamics alone.

Residual diagnostics were not exhaustively profiled in this run, and formal heteroskedasticity checks are a recommended next extension before applying non-linear or multi-horizon predictive models.

## Position in Program
Within the four-paper sequence, this study contributes the predictive-inference layer:

- Paper 1: measurement reliability
- Paper 2: detection under domain shift
- Paper 3: limits and validation of predictive inference
- Paper 4: causal intervention modeling

## Limitations
The current analysis remains constrained by modest lag-pair count, curated feature trajectories, and lack of external validation cohorts. Granger-style testing here is linear and one-step; broader causal inference would require longer series, richer covariates, and intervention-aware designs.

## Conclusion
We introduce a general evaluation framework for predictive claims in longitudinal resistome analysis, applicable beyond ARG-MGE systems, and demonstrate that naive lag models can be confounded by strong temporal coupling. In this run, differencing and Granger-style testing do not support a strong independent leading-indicator effect of MGE beyond ARG history. This reframes the result from over-claimed prediction to rigorous evaluation of predictive limits, which is the more robust methodological contribution.

This study provides a general diagnostic framework for distinguishing true predictive signals from temporal coupling in longitudinal omics data. This failure mode is likely to affect any longitudinal omics or ecological system in which multiple variables co-evolve over time, suggesting broad implications for predictive modeling beyond resistome surveillance.
