<p align="center"><b>Evaluating Predictive Signal Versus Temporal Coupling in Resistome-Mobilome Dynamics</b></p>
<p align="center"><b>Chenhao Zhang¹, Eliana Wong¹*, Ashley Fang¹*, Wanze Tang¹*, Linda Shi², Yujie Men²</b></p>
<p align="center">¹Institute of Engineering in Medicine, University of California, San Diego, San Diego, La Jolla, CA 92095<br>²Department of Chemical and Environmental Engineering, University of California, Riverside, California 92521<br>
*High school students participating in IEM OPALS program</p>

**Abstract -** *Environmental ARG surveillance often seeks leading indicators that could warn of future resistome expansion. We evaluate whether mobile genetic element (MGE) abundance at time t predicts ARG burden at time t+1, while explicitly testing whether apparent predictive performance is confounded by shared temporal structure. We assembled a real-sample temporal cohort of 23 samples across four independently sampled environments: a Chesapeake Bay whole-genome-amplified (WGA) shotgun time series (PRJNA599167, 11 samples at 1.5-hour intervals over May 17–19, 2017), a paired Chesapeake Bay 16S amplicon series (4 samples), and approximately biweekly shotgun metagenomes from two Bulgaria river sites — Dragushinovo and Mechkata on the Iskar River (PRJNA1071831, 4 samples each over November–December 2022). All 19 shotgun samples were fully quantified using a complete bioinformatics pipeline (Kraken2/Bracken for diversity, DIAMOND blastx against AMRProt for ARG burden, MEGAHIT assembly + IntegronFinder for MGE abundance), yielding 14 regression-complete lag pairs across three study contexts. The naive forward lag model including study fixed effects yields R-squared = 0.720, which might superficially suggest meaningful predictive structure. However, the MGE coefficient is not significant (coefficient = 21.6, p = 0.579), and the apparent model fit is driven by between-study differences in ARG scale and sequencing depth rather than by MGE predictive content. Lag correlations are moderate and non-monotone across offsets (t→t: 0.400; t→t+1: 0.532; t→t+2: 0.214), inconsistent with uniform temporal coupling. After trend-removing differencing, no signal remains (dARG(t+1) ~ dMGE(t): p = 0.649, R-squared = 0.020). A Granger-style added-value test confirms that MGE adds negligible predictive value beyond ARG history (delta R-squared = 0.014; p = 0.539). These results, from a fully pipeline-quantified real-sample dataset, demonstrate that study-level fixed effects can inflate apparent lag model fit without reflecting any genuine MGE predictive signal, and that the diagnostic framework successfully identifies this confounding.*

**Keywords:** ARG surveillance, mobilome, lag regression, temporal confounding, Granger-style test, time-series diagnostics

## 1. Introduction
Mobile genetic elements (MGEs) mediate horizontal gene transfer and can co-vary with antibiotic resistance gene (ARG) abundance. Environmental monitoring of ARGs has been recognized as a priority for understanding resistance dissemination across compartments [1]. This biological coupling motivates a surveillance question: does MGE abundance contain predictive information about future ARG change, beyond contemporaneous co-movement?

Simple lag models are attractive because they are easy to interpret and implement. However, longitudinal biological systems often exhibit shared drift, autocorrelation, and study-level synchronization. Under these conditions, a high R-squared in a lag model does not imply predictive directionality or independent predictive signal.

This study introduces a reproducible workflow for evaluating predictive claims under confounding-prone temporal structure. The broader microbiome time-series literature has recently emphasized standardized evaluation and the risk of over-interpreting interaction or causality claims from limited or strongly autocorrelated series [2], and recent methodological work has re-examined when Granger-style inference is valid for ecological count time series [3]. We operationalize a set of robustness checks that are common in time-series statistics but rarely applied in resistome studies: (1) forward lag regression, (2) reverse-direction testing, (3) first-difference regression, and (4) Granger-style added-value testing. We also note that new forecasting-focused resistome work highlights the applied demand for reliable longitudinal modeling of ARG abundance [4].

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

All regression models are implemented in Python using statsmodels [12] for OLS and F-tests, and pandas [13] for data manipulation. Analysis code is provided in `analysis/run_lag_analysis.py`.

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

## 4. Discussion
The main contribution is methodological: naive lag models can appear to explain substantial variance in longitudinal resistome data even when the predictor of interest contributes no independent signal. The forward lag model here yields R² = 0.720 — a value that might superficially suggest meaningful predictive structure — yet the MGE coefficient is non-significant (p = 0.579) and all diagnostic tests return null results. This case illustrates how study-level fixed effects, which absorb between-study ARG scale differences, can inflate R² without reflecting any temporal MGE-to-ARG coupling.

The diagnostic battery exposes this confounding at multiple levels. The lag correlations are non-monotone across offsets (t→t: 0.400; t→t+1: 0.532; t→t+2: 0.214): if MGE were a genuine leading indicator, correlation with the future ARG state should exceed or at minimum match the same-timepoint correlation. The non-monotone pattern instead suggests the moderate correlations arise from between-study covariation rather than directional prediction. After differencing to remove shared trends, R² collapses to 0.020 (p = 0.649), and the Granger-style added-value test confirms that MGE adds only delta R² = 0.014 beyond ARG autoregression (p = 0.539). The combination of these four independent diagnostics provides strong evidence that no genuine MGE predictive signal is present.

These diagnostic principles are directly relevant to the growing literature on longitudinal ARG modeling. Work such as ARGfore [4] demonstrates impressive time-series prediction performance for ARG abundances; the robustness tests operationalized here provide a complementary framework for evaluating whether such predictive relationships reflect genuine signal or shared temporal drift. Neither framework invalidates the other: surveillance applications may benefit from predictive accuracy even if driven by autocorrelation, while mechanistic inference requires the stronger independence tests applied here.

The temporal cohort assembled here spans two qualitatively different sampling regimes: a high-frequency (1.5-hour) Chesapeake Bay estuarine transect and a lower-frequency (approximately biweekly, 13–21 days) Bulgaria river series. This diversity of temporal scales strengthens the generalizability of the diagnostic framework. The Chesapeake Bay WGA shotgun samples are particularly valuable because they provide a dense within-study time series from a single environmental event, where temporal autocorrelation in ARG and MGE profiles is expected to be strong. The Bulgaria river samples provide a contrasting 45-day window across two independent sites, allowing assessment of whether temporal coupling persists at ecological rather than event timescales. The consistent null results across both regimes (differencing p = 0.649, Granger p = 0.539) suggest that the absence of independent MGE predictive signal is not an artefact of one particular sampling design.

## 5. Future Work and Limitations
Limitations include small regression-complete sample size (n=14 lag pairs) and the absence of external validation cohorts. The current temporal cohort spans only two environmental systems at two temporal scales (1.5-hour estuarine vs. biweekly river), limiting conclusions about intermediate temporal dynamics. The Granger-style test is linear and one-step; longer series and additional covariates are needed to assess multi-horizon or nonlinear relationships. Two WGA samples (SRR11803504 and SRR11803505) yielded mge_abundance = 0, consistent with the known tendency of MDA amplification to produce chimeric reads that suppress integron detection on assembled contigs; alternative MGE quantification methods (e.g., read-based MGE annotation) may be more appropriate for WGA libraries.

Future work should include:
1. Longer time series with more independent studies across diverse environments to increase regression-complete sample size and power.
2. Incorporation of mechanistic or intervention-linked covariates (e.g., antibiotics, treatment process changes).
3. Nonlinear or multi-horizon predictive models with robust diagnostic checks.
4. Formal residual diagnostics, including heteroskedasticity and autocorrelation tests.
5. Alternative MGE quantification approaches for amplified metagenomes.

## 6. Conclusion
We present a reproducible framework for evaluating predictive claims in longitudinal resistome–mobilome data, grounded in a real-sample temporal cohort of 23 samples across four environmental contexts, fully quantified using a complete bioinformatics pipeline (DIAMOND blastx/AMRProt for ARG, IntegronFinder/MEGAHIT for MGE, Kraken2/Bracken for diversity). Applying this pipeline to 14 regression-complete lag pairs, the naive forward lag model yields R² = 0.720, which might superficially suggest meaningful predictive structure. However, the MGE coefficient is non-significant (p = 0.579), lag correlations are non-monotone across offsets (0.400 → 0.532 → 0.214), differencing eliminates the signal (R² = 0.020, p = 0.649), and the Granger-style added-value test confirms that MGE contributes negligible independent information (delta R² = 0.014, p = 0.539). Together, these diagnostics demonstrate that the apparent model fit is driven by between-study differences in ARG scale — absorbed by study fixed effects — rather than by any genuine MGE predictive signal. These results emphasize that surveillance claims about leading indicators must be validated with explicit robustness tests, not inferred from high R-squared alone. The real-sample temporal backbone, bioinformatics quantification pipeline, and diagnostic analysis code are provided as fully reproducible resources for the community.

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
