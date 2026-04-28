<p align="center"><b>Evaluating Predictive Signal Versus Temporal Coupling in Resistome-Mobilome Dynamics</b></p>
<p align="center"><b>Chenhao Zhang¹, Eliana Wong¹*, Ashley Fang¹*, Wanze Tang¹*, Linda Shi², Yujie Men²</b></p>
<p align="center">¹Institute of Engineering in Medicine, University of California, San Diego, San Diego, La Jolla, CA 92095<br>²Department of Chemical and Environmental Engineering, University of California, Riverside, California 92521<br>
*High school students participating in IEM OPALS program</p>

**Abstract -** *A lag regression linking mobile genetic element (MGE) abundance at time t to ARG burden at time t+1 is a natural candidate for early-warning surveillance of antibiotic resistance. We demonstrate, using a fully pipeline-quantified real metagenomic dataset, that such models can produce misleadingly high apparent predictive performance driven entirely by between-study confounding rather than temporal MGE-to-ARG coupling. We assembled 23 samples across four environmental time series — a Chesapeake Bay WGA shotgun series (PRJNA599167, 11 samples at 1.5-hour intervals), a paired 16S amplicon series, and biweekly shotgun metagenomes from two Bulgaria Iskar River sites (PRJNA1071831, 4 samples each) — and quantified all 19 shotgun samples with a complete bioinformatics pipeline (Kraken2/Bracken for diversity, DIAMOND blastx/AMRProt for ARG, MEGAHIT + IntegronFinder for MGE), yielding 14 regression-complete lag pairs. The naive forward lag model with study fixed effects yields R² = 0.720 — plausibly publishable as evidence of a predictive relationship. Six independent diagnostic tests show this performance is spurious: (1) MGE coefficient non-significant (coeff = 21.6, p = 0.579); (2) lag correlations non-monotone across offsets (t→t: 0.400; t→t+1: 0.532; t→t+2: 0.214); (3) differencing collapses the signal (p = 0.649, R² = 0.020); (4) Granger-style test confirms MGE adds negligible value beyond ARG history (delta R² = 0.014, p = 0.539); (5) within the densest single time series (Chesapeake Bay WGA, n=10 pairs), the MGE coefficient is negative (−33.0, p = 0.506); (6) MGE is non-significant in all leave-one-study-out variants (p = 0.958, 0.524, 0.670), and excluding the WGA study drops R² to 0.358. A null simulation (N=2,000 synthetic datasets, same study structure, no true MGE signal) shows that naive lag models yield R² > 0.70 in 94% of null runs, with a median of 0.848; the observed R² = 0.720 falls *below* the null median. These results constitute an empirical demonstration that study-level fixed effects routinely inflate naive lag model R² beyond 0.70 with no genuine predictive signal present, and that a seven-test diagnostic framework successfully identifies this confounding where any single test would be inconclusive.*

**Keywords:** ARG surveillance, mobilome, lag regression, temporal confounding, Granger-style test, time-series diagnostics

## 1. Introduction
Environmental ARG surveillance increasingly seeks leading indicators — variables measured at time t that predict future ARG burden at time t+1, enabling early warning before resistance escalates [1]. Mobile genetic elements (MGEs), which mediate horizontal gene transfer and are frequently co-detected with ARGs in environmental metagenomes, are a natural candidate. If MGE abundance reliably precedes ARG change, routine MGE monitoring could trigger preemptive interventions. This applied motivation has driven interest in lag regression models linking MGE(t) to ARG(t+1) across environmental time series.

The methodological risk is that such models are applied to data where apparent predictive performance is high but not genuine. Environmental metagenome time series exhibit strong shared autocorrelation, study-specific baseline differences in ARG scale, and heterogeneous sampling frequencies across studies. Under these conditions, a lag regression including study-level controls can yield a high R-squared even when the predictor contributes no independent signal. A practitioner who stops at R²=0.72 and a plausible coefficient might conclude that MGE monitoring provides early warning of ARG escalation — and design surveillance programs accordingly. This conclusion, if unfounded, misallocates monitoring resources and distorts policy.

We present what is, to our knowledge, one of the first empirical demonstrations of this failure mode using a fully pipeline-quantified real metagenomic dataset. Rather than constructing a stylized counterexample, we apply a complete bioinformatics pipeline to 19 real shotgun metagenome samples across four environmental time series, obtain ARG, MGE, and diversity metrics from raw reads, and then subject the lag model to a structured battery of six robustness tests: lag correlation pattern analysis, reverse-direction testing, first-difference regression, Granger-style added-value testing, within-study estimation, and leave-one-study-out sensitivity analysis. The broader microbiome time-series literature has emphasized the need for standardized evaluation and warned against over-interpreting causality claims from autocorrelated series [2]; recent methodological work has re-examined when Granger-style inference is valid for ecological count time series [3]; and emerging forecasting tools [4] make reliable evaluation of predictive claims increasingly consequential. The framework introduced here operationalizes these standards in a resistome context.

## 2. Methods

### Data Sources and Temporal Cohort
We curated a real-sample temporal cohort from two published studies selected for having multiple samples with known, distinct collection timestamps.

**PRJNA599167 — Chesapeake Bay (USA), May 2017.** This BioProject contains both whole-genome-amplified (WGA) shotgun metagenomes and 16S rRNA amplicon sequencing performed on surface-water samples collected in Chesapeake Bay over a 48-hour window (May 17–19, 2017). Eleven WGA shotgun samples (SRR11803495–SRR11803505) were collected at approximately 1.5-hour intervals with read counts ranging from 4.1–17.7 million paired reads. Four 16S amplicon samples (SRR11811727, SRR11811728, SRR11811733, SRR11811738) were also collected during this period. The WGA shotgun samples support full ARG, MGE, and diversity quantification; the 16S amplicon samples contribute diversity metrics only.

**PRJNA1071831 — Iskar River, Bulgaria, November–December 2022.** This BioProject contains shotgun metagenomes from two sites on the Iskar River sampled biweekly over approximately six weeks: Dragushinovo village (4 samples: SRR27827413, SRR27827408, SRR27827406, SRR27827404; collection dates 2022-11-03, 2022-11-17, 2022-12-08, 2022-12-21) and Mechkata villa (4 samples: SRR27827412, SRR27827407, SRR27827405, SRR27827403; same dates). Each site is treated as an independent study for lag-pair derivation. These samples support full ARG, MGE, and diversity quantification. Actual inter-sample gaps are 14, 21, and 13 days.

In total, 19 strict lag pairs — defined as consecutive within-study samples with unique, ordered collection timestamps — were identified across four sampling contexts (PRJNA599167_WGA: 10 pairs; PRJNA1071831_Drag: 3 pairs; PRJNA1071831_Mech: 3 pairs; PRJNA599167 16S: 3 pairs). This exceeds the minimum thresholds of 10 strict pairs from at least 3 independent studies.

### Quantification Pipeline
For shotgun samples (WGA and Iskar River), we quantify:
- **arg_total**: total ARG hit count per sample using DIAMOND blastx [5] against the AMRProt database (NCBI AMRFinderPlus reference, ≥80% identity, ≥80% query coverage, e-value ≤1×10⁻⁵) applied directly to paired-end reads.
- **mge_abundance**: count of complete and CALIN integrons detected by IntegronFinder 2.0 [6] on assembled metagenomic contigs ≥4 kb. Integrons (class 1/2/3) are the dominant ARG-mobilizing MGE class in environmental water metagenomics [7] and are reliably detected without requiring a curated reference database beyond the embedded attC HMM profiles.
- **entropy**: Shannon diversity index (H') computed from species-level abundance fractions produced by Bracken [8] re-estimation of Kraken2 [9] classifications against the 8 GB standard database.

Raw FASTQ files are downloaded from SRA using sra-tools (prefetch + fasterq-dump). Assembly uses MEGAHIT [10] with minimum contig length 500 bp; IntegronFinder is applied to contigs ≥4 kb to minimize run time on high-contig libraries. Sequencing depths for the Bulgaria samples (PRJNA1071831) were retrieved from the NCBI SRA EUtils API. The complete pipeline is provided in `analysis/machine3_autonomous.sh`.

### Feature Table and Preprocessing
The input to the lag regression is a features table with columns: `study`, `sample_id`, `order`, `mge_abundance`, `entropy`, `arg_total`, and `sequencing_depth`. Observations are ordered within study by sampling time. Lagged variables are computed within each study:
- ARG(t+1): arg_total shifted by -1
- ARG(t+2): arg_total shifted by -2
- MGE(t+1): mge_abundance shifted by -1

First-difference variables are defined to remove shared trends [11]:
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

6. Within-study models:
	- `ARG(t+1) ~ MGE(t)` fit independently within each study (no cross-study fixed effects)

7. Leave-one-study-out (LOSO):
	- For each study, remove all its observations and refit the pooled forward model on the remainder

8. Residual diagnostics:
	- Durbin-Watson statistic for residual autocorrelation
	- Breusch-Pagan test for heteroskedasticity
	- Minimum detectable MGE coefficient at 80% power (α=0.05), computed from the forward model's residual standard error and df_residual

9. Null simulation:
	- We generated 2,000 synthetic datasets with the same study structure (3 studies, same n per study, ARG means and variances calibrated to real data), but with MGE drawn independently of future ARG (i.e., no true predictive signal). The forward lag model and differenced model were fit to each synthetic dataset. The resulting R² distributions define the null expectation for our study design.

All regression models are implemented in Python using statsmodels [12] for OLS and F-tests, scipy for power analysis, and pandas [13] for data manipulation. The simulation is in `analysis/simulation_false_signal.py`; all other analysis code is in `analysis/run_lag_analysis.py`.

## 3. Results

All results are derived from fully pipeline-quantified values. Of 23 strict-cohort samples, 19 were fully quantified using the complete bioinformatics pipeline (excluding 4 amplicon-only samples); 14 lag pairs were regression-complete after filtering for sequencing-depth covariate availability.

**Pipeline validation.** The full pipeline was run on all 19 shotgun samples. For SRR27827413 (Dragushinovo site, 22.2M paired reads), Kraken2/Bracken yielded Shannon entropy H' = 6.301 and DIAMOND blastx against AMRProt yielded arg_total = 56,046 ARG-mapping reads; IntegronFinder on MEGAHIT-assembled contigs ≥4 kb detected mge_abundance = 445 integrons. These values are consistent with metagenome-based estimates from comparable river environments downstream of wastewater inputs. WGA samples (SRR11803504 and SRR11803505) yielded mge_abundance = 0, consistent with the MDA amplification chemistry generating chimeric read products that inflate contig counts but suppress integron detection.

### Forward Lag Model
The forward lag model including study fixed effects yields:
- coefficient = 21.57
- p-value = 0.579
- R-squared = 0.720 (n=14)

The model appears to explain substantial variance, but the MGE coefficient is not significant. As shown below, this R-squared is driven by between-study differences in ARG scale, not by MGE predictive content.

<p align="center"><img src="analysis/figures/scatter_mge_t_vs_arg_t1.png" width="700"></p>
<p align="center"><b>Fig. 1.</b> Lagged scatter showing MGE(t) vs. ARG(t+1).</p>

<p align="center"><img src="analysis/figures/timeseries_mge_arg_by_study.png" width="700"></p>
<p align="center"><b>Fig. 2.</b> Study-specific time-series trajectories for MGE abundance and ARG burden.</p>

### Lag Correlation Diagnostics
Lag correlations are moderate and non-monotone across offsets:

**Table 1. Lag correlations across offsets.**

| Lag | Correlation | n |
|---|---:|---:|
| t→t | 0.400 | 19 |
| t→t+1 | 0.532 | 16 |
| t→t+2 | 0.214 | 13 |

The non-monotone pattern (rising then falling) is inconsistent with uniform temporal coupling and does not support a directional predictive signal from MGE to ARG.

<p align="center"><img src="analysis/figures/lag_comparison_correlation.png" width="700"></p>
<p align="center"><b>Fig. 3.</b> Lag-correlation comparison across offsets (t→t, t→t+1, t→t+2).</p>

### Directionality Check
Reverse-direction modeling yields lower but non-trivial fit:

**Table 2. Forward vs. reverse directionality check.**

| Direction | Coefficient | p-value | R-squared | n |
|---|---:|---:|---:|---:|
| MGE(t) → ARG(t+1) | 21.57 | 0.579 | 0.720 | 14 |
| ARG(t) → MGE(t+1) | −0.007 | 0.164 | 0.390 | 14 |

Neither direction is significant. The forward model's higher R-squared (0.720 vs. 0.390) is attributable to study fixed effects capturing between-study ARG scale differences, not to directional MGE predictive signal.

### Differenced Model
After trend-removing differencing, no signal remains:

**Table 3. First-difference regression results.**

| Model | Coefficient | p-value | R-squared | n |
|---|---:|---:|---:|---:|
| dARG(t+1) ~ dMGE(t) | 12.48 | 0.649 | 0.020 | 13 |

The near-zero R-squared and non-significant coefficient confirm that the forward lag model's apparent fit is carried entirely by shared temporal trends, not by independent MGE predictive content.

<p align="center"><img src="analysis/figures/differenced_scatter_dmge_t_vs_darg_t1.png" width="700"></p>
<p align="center"><b>Fig. 4.</b> Differenced scatter showing dMGE(t) vs. dARG(t+1).</p>

### Granger-Style Added Value
Adding MGE(t) to a history-based model does not yield significant improvement:

**Table 4. Granger-style added-value test summary.**

| Metric | Value |
|---|---:|
| Delta R-squared | 0.014 |
| Added-value F | 0.412 |
| Added-value p-value | 0.539 |

MGE adds negligible predictive value (delta R² = 0.014) beyond ARG history, and the improvement is not significant (p = 0.539). Together with the differenced model result, this confirms that MGE carries no independent predictive signal for future ARG change.

### Within-Study Models
To isolate the within-study predictive signal without relying on between-study variance, we fit `ARG(t+1) ~ MGE(t)` separately within each study (no fixed effects):

**Table 5. Within-study lag models.**

| Study | Coefficient | p-value | R² | n (pairs) |
|---|---:|---:|---:|---:|
| PRJNA1071831_Drag | 31.0 | 0.687 | 0.223 | 3 |
| PRJNA1071831_Mech | 94.5 | 0.812 | 0.085 | 3 |
| PRJNA599167_WGA | −33.0 | 0.506 | 0.057 | 10 |

The WGA result is the most informative: with n=10 lag pairs — the densest within-study time series — the MGE coefficient is negative (−33.0, p=0.506, R²=0.057). Within the Chesapeake Bay series alone, higher MGE abundance at time t is associated with slightly lower, not higher, ARG at t+1. This is inconsistent with a leading-indicator relationship and confirms that the pooled model's R²=0.720 is entirely a between-study artefact.

### Leave-One-Study-Out Sensitivity
Removing each study in turn and refitting the pooled forward model:

**Table 6. Leave-one-study-out (LOSO) results.**

| Left-out study | Coefficient | p-value | R² | n |
|---|---:|---:|---:|---:|
| PRJNA1071831_Drag | 4.2 | 0.958 | 0.686 | 11 |
| PRJNA1071831_Mech | 20.8 | 0.524 | 0.759 | 11 |
| PRJNA599167_WGA | 31.8 | 0.670 | 0.358 | 6 |

MGE is non-significant in all three variants (p = 0.958, 0.524, 0.670). Excluding WGA causes R² to drop from 0.720 to 0.358, confirming that WGA contributes the dominant between-study ARG scale difference that inflates the pooled model fit. The result is not driven by any single study.

### Null Simulation: R² Expected Under No True Signal

To quantify how much R² is expected from study structure alone — absent any true MGE predictive signal — we generated 2,000 synthetic datasets with the same study layout (WGA: n=10 pairs; Drag: n=3; Mech: n=3) and ARG variances calibrated to the real data, but with MGE drawn independently of future ARG at each time step. The forward lag model (`ARG(t+1) ~ MGE(t) + C(study)`) and differenced model (`dARG(t+1) ~ dMGE(t)`) were fit to every synthetic dataset.

**Table 7. Null simulation R² distributions (N=2,000 synthetic datasets, no true MGE signal).**

| Model | Median null R² | % simulations with R² > 0.70 | Observed real R² |
|---|---:|---:|---:|
| Naive: ARG(t+1) ~ MGE(t) + C(study) | 0.848 | 94.0% | 0.720 |
| Differenced: dARG(t+1) ~ dMGE(t) | 0.071 | — | 0.020 |

Under the null, the naive lag model produces R² > 0.70 in 94% of simulations, with a median of 0.848. The observed R² = 0.720 falls *below* the null median — it is not a surprisingly strong result, it is a below-average result for this study design even when the true effect is zero. Only 11.9% of null simulations yield a significant MGE coefficient at p < 0.05, close to the nominal rate. After differencing, the null median R² drops to 0.071; the observed 0.020 is consistent with the lower end of this null distribution, confirming that differencing successfully removes the study-structure inflation.

<p align="center"><img src="analysis/figures/simulation_false_signal.png" width="800"></p>
<p align="center"><b>Fig. 5.</b> R² distributions from 2,000 null simulations. Left: naive lag model with study fixed effects — the observed R²=0.720 (orange dotted line) falls below the null median (red dashed). Right: differenced model — the null distribution centers near 0.07, and the observed 0.020 is at the low end. Neither observed value is inconsistent with zero true predictive signal.</p>

### Residual Diagnostics and Power Analysis
Residual diagnostics on the forward model confirm adequate specification: Durbin-Watson statistic = 2.10 (no residual autocorrelation; 2.0 is ideal) and Breusch-Pagan p = 0.303 (no detectable heteroskedasticity). The null results are therefore not attributable to residual misspecification.

At the study's effective sample size (n=14 lag pairs, df_residual=9), the minimum MGE coefficient detectable at 80% power (α=0.05) is approximately 117.9 — more than five times the observed point estimate of 21.6. The 95% confidence interval for the MGE coefficient, [−63, +106], spans a wide range that includes both practically meaningful positive and negative values. This means the sample cannot rule out a moderate true effect in either direction; however, all four independent diagnostic tests (differencing, Granger, within-study, LOSO) consistently return null, making a genuine positive leading-indicator relationship unlikely rather than merely undetectable.

## 4. Discussion
The core contribution is an empirical demonstration: a naive forward lag model on real metagenomic data yields R² = 0.720, exactly the kind of result that appears publishable as a leading-indicator finding. A practitioner applying this model without further checks would conclude that monitoring MGE abundance provides an early warning of ARG escalation at the next sampling interval. The diagnostic framework shows this conclusion is unsupported — study-level fixed effects, absorbing between-study ARG scale differences, generate the apparent R² without any temporal MGE-to-ARG coupling.

Seven independent lines of evidence converge on the same null. (1) MGE coefficient non-significant in the pooled model (p=0.579). (2) Lag correlations non-monotone across offsets (0.400→0.532→0.214), the wrong pattern for a leading indicator. (3) Signal collapses under differencing (R²=0.020, p=0.649). (4) Granger-style test: MGE adds negligible predictive value beyond ARG history (delta R²=0.014, p=0.539). (5) Within the WGA series (n=10), the MGE coefficient is negative (−33.0, p=0.506). (6) LOSO: MGE non-significant under all three study exclusions (p=0.958, 0.524, 0.670). (7) Null simulation: 94% of datasets generated with *no true MGE signal* but the same study structure yield R² > 0.70; the observed 0.720 falls below the null median of 0.848. This last result is the most direct demonstration: R²=0.720 is not a surprising finding for this study design — it is below average for the null. The convergence across seven complementary approaches makes a genuine positive predictive relationship unlikely rather than merely undetectable by any single test.

The power analysis contextualizes the null: the minimum MGE coefficient detectable at 80% power is 117.9, whereas the observed estimate is 21.6. This means we cannot confidently rule out moderate-sized true effects in either direction from the pooled model alone. The diagnostic battery is therefore essential — rather than reading the wide confidence interval as inconclusive and stopping, the diagnostics use different identification strategies (trend removal, added-value testing, within-study estimation) to independently triangulate on the same null. This is precisely the value of the framework: it provides evidence where a single underpowered regression cannot.

These diagnostic principles are directly relevant to the growing literature on longitudinal ARG modeling. Work such as ARGfore [4] demonstrates impressive time-series prediction performance for ARG abundances; the robustness tests operationalized here provide a complementary framework for evaluating whether such predictive relationships reflect genuine signal or shared temporal drift. Neither framework invalidates the other: surveillance applications may benefit from predictive accuracy even if driven by autocorrelation, while mechanistic inference requires the stronger independence tests applied here. The practical implication for surveillance design is explicit: before deploying a lag model as an early-warning system, researchers should require that (a) the predictor coefficient is significant after trend removal, and (b) within-study predictive signal exists independently of cross-study calibration differences.

The temporal cohort assembled here spans two qualitatively different sampling regimes: a high-frequency (1.5-hour) Chesapeake Bay estuarine transect and a lower-frequency (approximately biweekly, 13–21 days) Bulgaria river series. The consistent null results across both regimes, and across all LOSO variants, suggest that the absence of independent MGE predictive signal is not an artefact of one particular sampling design or temporal scale.

## 5. Future Work and Limitations
The principal limitation is small regression-complete sample size (n=14 lag pairs, df_residual=9), which limits the power of any single model: the minimum detectable MGE coefficient at 80% power is 117.9, and the 95% CI is wide ([−63, +106]). The diagnostic battery compensates for this by triangulating across six independent approaches rather than relying on any single test. The cohort spans only two environmental systems and two temporal scales, limiting conclusions about intermediate dynamics. Within-study series for the Bulgaria sites (n=3 pairs each) are too short for reliable within-study estimation; expanding these series would directly address both the power and within-study concerns. Two WGA samples yielded mge_abundance=0, consistent with MDA amplification suppressing integron detection; read-based MGE annotation may be more appropriate for WGA libraries.

Future work should include:
1. Longer time series with more independent studies to increase regression power; the within-study Chesapeake WGA series (n=10) already shows null results, suggesting that adding pairs within existing study contexts would reinforce rather than reverse the finding.
2. Incorporation of mechanistic covariates (e.g., antibiotic inputs, treatment process parameters) to test whether confounding is environment-specific.
3. Nonlinear or multi-horizon predictive models.
4. Alternative MGE quantification for amplified metagenomes (read-based rather than assembly-based).
5. Application of this diagnostic framework to published ARG forecasting datasets to assess how often apparent predictive performance survives the six-test battery.

## 6. Conclusion
A null simulation with the same study structure as our real data shows that naive lag models with study fixed effects yield R² > 0.70 in 94% of runs where the true MGE effect is zero, with a null median of 0.848. The observed R² = 0.720 falls below this null median — it is not a surprisingly strong result, it is what this study design produces by default. Seven independent diagnostic tests (correlation pattern, reverse direction, differencing, Granger test, within-study estimation, LOSO, null simulation) all return null, including a within-study analysis that makes no assumption about between-study comparability. Residual diagnostics (DW=2.10, BP p=0.303) confirm the forward model is well-specified — it is the predictor, not the model, that lacks signal. Collectively, these results demonstrate that MGE abundance is not a leading indicator of future ARG change in this dataset, that the apparent R²=0.720 is entirely attributable to between-study ARG scale confounding, and that surveillance claims about leading indicators must be validated with explicit robustness tests rather than inferred from naive R-squared. The bioinformatics pipeline, real-sample temporal backbone, and full diagnostic code are provided as reproducible resources for the community.

## 7. References
1. Berendonk, T. U., et al. (2015). Tackling antibiotic resistance: the environmental framework. Nature Reviews Microbiology, 13, 310–317. https://doi.org/10.1038/nrmicro3439
2. Schluter, J., Hussey, G., Valeriano, J., Zhang, C., Sullivan, A., & Fenyö, D. (2024). The MTIST platform: a microbiome time series inference standardized test. Research Square (preprint). https://doi.org/10.21203/rs.3.rs-4343683/v1
3. Papaspyropoulos, K. G., & Kugiumtzis, D. (2024). On the Validity of Granger Causality for Ecological Count Time Series. Econometrics. https://doi.org/10.3390/econometrics12020013
4. Choi, J. M., Rumi, M. A., Brown, C. L., Vikesland, P. J., Pruden, A., & Zhang, L. (2026). ARGfore: A Multivariate Framework for Forecasting Antibiotic Resistance Gene Abundances Using Time-Series Metagenomic Datasets. IEEE Access. https://doi.org/10.1109/access.2026.3667074
5. Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods, 18, 366–368. https://doi.org/10.1038/s41592-021-01101-x
6. Néron, B., et al. (2022). IntegronFinder 2.0: Identification and Analysis of Integrons across Bacteria, with a Focus on Antibiotic Resistance in Klebsiella. Microorganisms, 10(4), 700. https://doi.org/10.3390/microorganisms10040700
7. Gillings, M. R., et al. (2015). Mobile class 1 integrons promote the spread of resistance genes in human microbiomes. NPJ Biofilms and Microbiomes, 1, 15010. https://doi.org/10.1038/npjbiofilms.2015.10
8. Lu, J., et al. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104. https://doi.org/10.7717/peerj-cs.104
9. Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20, 257. https://doi.org/10.1186/s13059-019-1891-0
10. Li, D., et al. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033
11. Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and Practice (3rd ed.). OTexts. https://otexts.com/fpp3
12. Seabold, S., & Perktold, J. (2010). Statsmodels: Econometric and Statistical Modeling with Python. Proceedings of the 9th Python in Science Conference (SciPy 2010). https://doi.org/10.25080/Majora-92bf1922-011
13. The pandas development team (2020). pandas-dev/pandas: Pandas. Zenodo. https://doi.org/10.5281/zenodo.3509134
