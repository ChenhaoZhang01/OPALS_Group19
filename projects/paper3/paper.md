<p align="center"><b>Evaluating Predictive Signal Versus Temporal Coupling in Resistome-Mobilome Dynamics</b></p>
<p align="center"><b>Chenhao Zhang¹, Eliana Wong¹*, Ashley Fang¹*, Wanze Tang¹*, Linda Shi², Yujie Men²</b></p>
<p align="center">¹Institute of Engineering in Medicine, University of California, San Diego, San Diego, La Jolla, CA 92095<br>²Department of Chemical and Environmental Engineering, University of California, Riverside, California 92521<br>
*High school students participating in IEM OPALS program</p>

**Abstract -** *Environmental ARG surveillance often seeks leading indicators that could warn of future resistome expansion. We evaluate whether mobile genetic element (MGE) abundance at time t predicts ARG burden at time t+1, while explicitly testing whether apparent predictive performance is confounded by shared temporal structure. Using a curated feature table with 33 lag pairs across multiple studies, a naive forward lag model (ARG(t+1) ~ MGE(t) + study + sequencing depth) yields an extremely strong association (coefficient 18.554, p < 1e-8, R-squared = 0.998). However, lag correlations are uniformly high across offsets (t→t: 0.997; t→t+1: 0.995; t→t+2: 0.993), and reverse-direction modeling is similarly strong (ARG(t) → MGE(t+1), R-squared = 0.999). After trend-removing differencing, the association collapses (dARG(t+1) ~ dMGE(t): p = 0.0616, R-squared = 0.119). A Granger-style added-value test shows no significant improvement beyond ARG history (delta R-squared = 0.000027; p = 0.372). These results demonstrate that naive lag models can overstate predictive interpretation under strong temporal coupling. We provide a rigorous, reproducible framework for assessing predictive claims in longitudinal resistome studies.*

**Keywords:** ARG surveillance, mobilome, lag regression, temporal confounding, Granger-style test, time-series diagnostics

## 1. Introduction
Mobile genetic elements (MGEs) mediate horizontal gene transfer and can co-vary with antibiotic resistance gene (ARG) abundance. This biological coupling motivates a surveillance question: does MGE abundance contain predictive information about future ARG change, beyond contemporaneous co-movement?

Simple lag models are attractive because they are easy to interpret and implement. However, longitudinal biological systems often exhibit shared drift, autocorrelation, and study-level synchronization. Under these conditions, a high R-squared in a lag model does not imply predictive directionality or independent predictive signal.

This study introduces a reproducible workflow for evaluating predictive claims under confounding-prone temporal structure. The broader microbiome time-series literature has recently emphasized standardized evaluation and the risk of over-interpreting interaction or causality claims from limited or strongly autocorrelated series [4], and recent methodological work has re-examined when Granger-style inference is valid for ecological count time series [5]. We operationalize a set of robustness checks that are common in time-series statistics but rarely applied in resistome studies: (1) forward lag regression, (2) reverse-direction testing, (3) first-difference regression, and (4) Granger-style added-value testing. We also note that new forecasting-focused resistome work highlights the applied demand for reliable longitudinal modeling of ARG abundance [6].

## 2. Methods

### Data and Preprocessing
The input is a features table with columns: `study`, `sample_id`, `order`, `mge_abundance`, `entropy`, `arg_total`, and `sequencing_depth`. Observations are ordered within study by sampling time. Lagged variables are computed within each study:
- ARG(t+1): arg_total shifted by -1
- ARG(t+2): arg_total shifted by -2
- MGE(t+1): mge_abundance shifted by -1

First-difference variables are defined as:
- dMGE(t) = MGE(t) - MGE(t-1)
- dARG(t+1) = ARG(t+1) - ARG(t)

Sequencing depth is modeled on the raw read-count scale when available. Missing sequencing depth is left as null and excluded from models that require it.

### Models
We evaluate the following models:

1. Forward lag model:
	- ARG(t+1) ~ MGE(t) + C(study) + sequencing_depth

2. Entropy-augmented model:
	- ARG(t+1) ~ MGE(t) + entropy(t) + C(study) + sequencing_depth

3. Reverse-direction model:
	- MGE(t+1) ~ ARG(t) + C(study) + sequencing_depth

4. First-difference model:
	- dARG(t+1) ~ dMGE(t)

5. Granger-style added-value test:
	- Base: ARG(t+1) ~ ARG(t) + C(study) + sequencing_depth
	- Full: ARG(t+1) ~ ARG(t) + MGE(t) + C(study) + sequencing_depth

Nested model improvement is assessed using F-tests and delta R-squared. Pearson correlations are computed for t→t, t→t+1, and t→t+2 to assess the stability of correlation across lags.

## 3. Results

### Forward Lag Model
The forward lag model shows a strong association:
- coefficient = 18.553611
- p-value < 1e-8
- R-squared = 0.998303 (n=33)

Taken alone, this suggests a strong predictive relationship between MGE(t) and ARG(t+1).

<p align="center"><img src="analysis/figures/scatter_mge_t_vs_arg_t1.png" width="700"></p>
<p align="center"><b>Fig. 1.</b> Lagged scatter showing MGE(t) vs. ARG(t+1).</p>

<p align="center"><img src="analysis/figures/timeseries_mge_arg_by_study.png" width="700"></p>
<p align="center"><b>Fig. 2.</b> Study-specific time-series trajectories for MGE abundance and ARG burden.</p>

### Lag Correlation Diagnostics
Lag correlations remain high across offsets:

**Table 1. Lag correlations across offsets.**

| Lag | Correlation | n |
|---|---:|---:|
| t→t | 0.996761 | 36 |
| t→t+1 | 0.995173 | 33 |
| t→t+2 | 0.992825 | 30 |

This uniformity indicates strong temporal coupling rather than directional prediction.

<p align="center"><img src="analysis/figures/lag_comparison_correlation.png" width="700"></p>
<p align="center"><b>Fig. 3.</b> Lag-correlation comparison across offsets (t→t, t→t+1, t→t+2).</p>

### Directionality Check
Reverse-direction modeling is similarly strong:

**Table 2. Forward vs. reverse directionality check.**

| Direction | Coefficient | R-squared | n |
|---|---:|---:|---:|
| MGE(t) → ARG(t+1) | 18.554 | 0.998 | 33 |
| ARG(t) → MGE(t+1) | 0.054 | 0.999 | 33 |

The symmetry of forward and reverse fits suggests that the lag signal is not directional.

### Differenced Model
After removing shared trends, the signal collapses:

**Table 3. First-difference regression results.**

| Model | Coefficient | p-value | R-squared | n |
|---|---:|---:|---:|---:|
| dARG(t+1) ~ dMGE(t) | 4.545 | 0.0616 | 0.119 | 30 |

This indicates the original predictive signal is largely driven by shared temporal structure.

<p align="center"><img src="analysis/figures/differenced_scatter_dmge_t_vs_darg_t1.png" width="700"></p>
<p align="center"><b>Fig. 4.</b> Differenced scatter showing dMGE(t) vs. dARG(t+1).</p>

### Granger-Style Added Value
Adding MGE(t) to a history-based model does not yield significant improvement:

**Table 4. Granger-style added-value test summary.**

| Metric | Value |
|---|---:|
| Delta R-squared | 0.000027 |
| Added-value F | 0.823097 |
| Added-value p-value | 0.372303 |

MGE does not provide independent predictive gain beyond ARG history in this run.

## 4. Discussion
The main contribution is methodological: naive lag models can appear nearly perfect in longitudinal resistome data even when predictive signal is not independent. The combination of lag correlation diagnostics, reverse-direction checks, differencing, and Granger-style testing provides a practical framework to avoid over-interpretation.

In this dataset, the strong forward association is explained by shared temporal dynamics rather than leading-indicator effects. This reframes the result from a predictive claim to a cautionary demonstration: in the presence of strong temporal coupling, high fit is not evidence of predictive direction.

## 5. Future Work and Limitations
Limitations include modest lag-pair counts, curated feature trajectories, and the absence of external validation cohorts. The Granger-style test is linear and one-step; longer series and additional covariates are needed to assess multi-horizon or nonlinear relationships.

Future work should include:
1. Longer time series with more independent studies.
2. Incorporation of mechanistic or intervention-linked covariates.
3. Nonlinear or multi-horizon predictive models with robust diagnostic checks.
4. Formal residual diagnostics, including heteroskedasticity and autocorrelation tests.

## 6. Conclusion
We present a reproducible framework for evaluating predictive claims in longitudinal resistome–mobilome data. In this run, naive lag models show near-perfect fit, but differencing, lag-correlation diagnostics, reverse-direction testing, and Granger-style evaluation demonstrate that the apparent predictive signal is largely due to temporal coupling. These results emphasize that surveillance claims about leading indicators must be validated with explicit robustness tests, not inferred from high R-squared alone.

## 7. References
1. statsmodels documentation for linear models and ANOVA. https://www.statsmodels.org
2. pandas documentation for tabular data processing. https://pandas.pydata.org
3. General guidance on temporal autocorrelation and differencing in time-series analysis. https://otexts.com/fpp3
4. Schluter, J., Hussey, G., Valeriano, J., Zhang, C., Sullivan, A., & Fenyö, D. (2024). The MTIST platform: a microbiome time series inference standardized test. Research Square (preprint). https://doi.org/10.21203/rs.3.rs-4343683/v1
5. Papaspyropoulos, K. G., & Kugiumtzis, D. (2024). On the Validity of Granger Causality for Ecological Count Time Series. Econometrics. https://doi.org/10.3390/econometrics12020013
6. Choi, J. M., Rumi, M. A., Brown, C. L., Vikesland, P. J., Pruden, A., & Zhang, L. (2026). ARGfore: A Multivariate Framework for Forecasting Antibiotic Resistance Gene Abundances Using Time-Series Metagenomic Datasets. IEEE Access. https://doi.org/10.1109/access.2026.3667074
