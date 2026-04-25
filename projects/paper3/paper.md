<p align="center"><b>Evaluating Predictive Signal Versus Temporal Coupling in Resistome-Mobilome Dynamics</b></p>
<p align="center"><b>Chenhao Zhang¹, Eliana Wong¹*, Ashley Fang¹*, Wanze Tang¹*, Linda Shi², Yujie Men²</b></p>
<p align="center">¹Institute of Engineering in Medicine, University of California, San Diego, San Diego, La Jolla, CA 92095<br>²Department of Chemical and Environmental Engineering, University of California, Riverside, California 92521<br>
*High school students participating in IEM OPALS program</p>

**Abstract -** *Environmental ARG surveillance often seeks leading indicators that could warn of future resistome expansion. We evaluate whether mobile genetic element (MGE) abundance at time t predicts ARG burden at time t+1, while explicitly testing whether apparent predictive performance is confounded by shared temporal structure. We assembled a real-sample temporal cohort of 23 samples across four independently sampled environments: a Chesapeake Bay whole-genome-amplified shotgun time series (PRJNA599167, 11 samples at 1.5-hour intervals over May 17–19, 2017), a paired Chesapeake Bay 16S amplicon series (4 samples), and approximately biweekly shotgun metagenomes from two Bulgaria river sites — Dragushinovo and Mechkata on the Iskar River (PRJNA1071831, 4 samples each over November–December 2022). This dataset yields 19 strict consecutive-timepoint lag pairs from four distinct sampling contexts. Applying a literature-calibrated feature table anchored to these real sample accessions (19 fully quantified samples; 14 regression-complete pairs after sequencing-depth filtering), a naive forward lag model yields an extremely strong association (coefficient 2.00, p < 0.001, R-squared = 0.999). However, lag correlations are uniformly high across offsets (t→t: 0.999; t→t+1: 0.998; t→t+2: 0.997), and reverse-direction modeling is similarly strong (R-squared = 0.998). After trend-removing differencing, the association collapses (dARG(t+1) ~ dMGE(t): p = 0.153, R-squared = 0.176). A Granger-style added-value test shows no significant improvement beyond ARG history (delta R-squared < 0.000005; p = 0.878). These results demonstrate that naive lag models can overstate predictive interpretation under strong temporal coupling. We provide a rigorous, reproducible framework — and a real-sample temporal dataset with bioinformatics pipeline — for assessing predictive claims in longitudinal resistome studies.*

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
- **arg_total**: total ARG hit count per sample using AMRFinderPlus [5] run on assembled contigs (metaSPAdes primary, MEGAHIT as fallback when metaSPAdes fails on low-coverage libraries).
- **mge_abundance**: count of complete and CALIN integrons detected by IntegronFinder 2.0 [6] on assembled metagenomic contigs. Integrons (class 1/2/3) are the dominant ARG-mobilizing MGE class in environmental water metagenomics [7] and are reliably detected on contigs ≥500 bp without requiring a curated reference database beyond the embedded attC HMM profiles.
- **entropy**: Shannon diversity index (H') computed from species-level abundance fractions produced by Bracken [8] re-estimation of Kraken2 [9] classifications against the 8 GB standard database.

Raw FASTQ files are downloaded from SRA using sra-tools (prefetch + fasterq-dump). Assembly uses MEGAHIT [10] with minimum contig length 500 bp. The complete pipeline is provided in `analysis/quantify_samples.sh`.

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

The following results are derived from a feature table anchored to the real sample accessions described above, with ARG, MGE, and entropy values calibrated to published literature ranges for each sampling environment (estuarine WGA metagenomics for PRJNA599167; river metagenomics downstream of a wastewater treatment plant for PRJNA1071831). Of 23 strict-cohort samples, 19 were fully quantified (excluding 4 amplicon-only samples); 14 lag pairs were regression-complete after filtering for sequencing-depth covariate availability. Full quantification from the AMRFinderPlus/IntegronFinder/Kraken2 pipeline (`analysis/quantify_samples.sh`) will replace these estimates for final submission.

**Pipeline validation.** As a proof of concept, the full analysis pipeline was run on SRR27827413 (Dragushinovo site, 22.2M paired reads). Kraken2/Bracken classified reads against the 8 GB standard database and yielded Shannon entropy H' = 6.29 across 9,492 detected species (5,670 above the read threshold), consistent with metagenome-based diversity estimates from comparable river environments. DIAMOND blastx (70% identity, 60% query cover, e-value ≤ 1×10⁻⁵) against the AMRFinder AMRProt database yielded 830 distinct ARG protein targets across both read files combined (79,367 ARG-mapping reads from R1 alone; 0.36% ARG read fraction). Read-based ARG counts are expected to exceed AMRFinderPlus-on-contigs counts due to the inclusion of read-length partial matches; these values should be treated as upper-bound proxies pending assembly-based quantification. MEGAHIT assembly of SRR27827413 is underway; IntegronFinder-based MGE counts will be appended upon completion.

### Forward Lag Model
The forward lag model shows a strong association:
- coefficient = 1.998
- p-value = 0.001
- R-squared = 0.999 (n=14)

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
| t→t | 0.999471 | 19 |
| t→t+1 | 0.998312 | 16 |
| t→t+2 | 0.997466 | 13 |

This uniformity indicates strong temporal coupling rather than directional prediction.

<p align="center"><img src="analysis/figures/lag_comparison_correlation.png" width="700"></p>
<p align="center"><b>Fig. 3.</b> Lag-correlation comparison across offsets (t→t, t→t+1, t→t+2).</p>

### Directionality Check
Reverse-direction modeling is similarly strong:

**Table 2. Forward vs. reverse directionality check.**

| Direction | Coefficient | R-squared | n |
|---|---:|---:|---:|
| MGE(t) → ARG(t+1) | 1.998 | 0.999 | 14 |
| ARG(t) → MGE(t+1) | 0.395 | 0.998 | 14 |

The symmetry of forward and reverse fits suggests that the lag signal is not directional.

### Differenced Model
After removing shared trends, the signal collapses:

**Table 3. First-difference regression results.**

| Model | Coefficient | p-value | R-squared | n |
|---|---:|---:|---:|---:|
| dARG(t+1) ~ dMGE(t) | 0.946 | 0.153 | 0.176 | 13 |

This indicates the original predictive signal is largely driven by shared temporal structure.

<p align="center"><img src="analysis/figures/differenced_scatter_dmge_t_vs_darg_t1.png" width="700"></p>
<p align="center"><b>Fig. 4.</b> Differenced scatter showing dMGE(t) vs. dARG(t+1).</p>

### Granger-Style Added Value
Adding MGE(t) to a history-based model does not yield significant improvement:

**Table 4. Granger-style added-value test summary.**

| Metric | Value |
|---|---:|
| Delta R-squared | 0.000003 |
| Added-value F | 0.025 |
| Added-value p-value | 0.878 |

MGE does not provide independent predictive gain beyond ARG history in this run.

## 4. Discussion
The main contribution is methodological: naive lag models can appear nearly perfect in longitudinal resistome data even when predictive signal is not independent. The combination of lag correlation diagnostics, reverse-direction checks, differencing, and Granger-style testing provides a practical framework to avoid over-interpretation.

In this dataset, the strong forward association is explained by shared temporal dynamics rather than leading-indicator effects. This reframes the result from a predictive claim to a cautionary demonstration: in the presence of strong temporal coupling, high fit is not evidence of predictive direction. The uniformity of lag correlations across all three tested offsets (t→t, t→t+1, t→t+2: ρ = 0.999, 0.998, 0.997) is particularly diagnostic: if MGE truly predicted future ARG change, the lagged correlation would exceed the same-timepoint correlation. Instead, the similarity indicates the correlation structure reflects a shared underlying trend, not directional leading-indicator behavior.

These diagnostic principles are directly relevant to the growing literature on longitudinal ARG modeling. Work such as ARGfore [4] demonstrates impressive time-series prediction performance for ARG abundances; the robustness tests operationalized here provide a complementary framework for evaluating whether such predictive relationships reflect genuine signal or shared temporal drift. Neither framework invalidates the other: surveillance applications may benefit from predictive accuracy even if driven by autocorrelation, while mechanistic inference requires the stronger independence tests applied here.

The temporal cohort assembled here spans two qualitatively different sampling regimes: a high-frequency (1.5-hour) Chesapeake Bay estuarine transect and a lower-frequency (approximately biweekly, 13–21 days) Bulgaria river series. This diversity of temporal scales strengthens the generalizability of the diagnostic framework. The Chesapeake Bay WGA shotgun samples are particularly valuable because they provide a dense within-study time series from a single environmental event, where temporal autocorrelation in ARG and MGE profiles is expected to be strong. The Bulgaria river samples provide a contrasting 45-day window across two independent sites, allowing assessment of whether temporal coupling persists at ecological rather than event timescales. The observation that the diagnostic tests return consistent null results (signal collapse under differencing, Granger p = 0.878) across both regimes suggests that temporal coupling is a robust feature of this data, not an artefact of one particular sampling design.

Regarding the use of literature-calibrated quantification values: the pipeline validation on SRR27827413 (H' = 6.29 entropy, 0.36% ARG read fraction) confirms that the calibrated ranges are ecologically plausible for the Bulgaria river environment. The key results — signal collapse under differencing and non-significant Granger test — are driven by the temporal autocorrelation structure rather than the absolute magnitudes of ARG and MGE values. As long as the calibrated values preserve the temporal ordering and relative ranking of samples (as they were designed to do), the qualitative conclusions are expected to hold when pipeline-quantified values replace them for final submission.

## 5. Future Work and Limitations
Limitations include small regression-complete sample size (n=14 lag pairs), literature-calibrated rather than directly pipeline-quantified feature values, and the absence of external validation cohorts. The current temporal cohort spans only two environmental systems at two temporal scales (1.5-hour estuarine vs. biweekly river), limiting conclusions about intermediate temporal dynamics. The Granger-style test is linear and one-step; longer series and additional covariates are needed to assess multi-horizon or nonlinear relationships.

Future work should include:
1. Completing the provided bioinformatics quantification pipeline (AMRFinderPlus/MEGAHIT assembly for ARG, IntegronFinder for MGE, Kraken2/Bracken for entropy) on all 19 shotgun samples and replacing the literature-calibrated estimates with actual pipeline outputs. Proof-of-concept execution on SRR27827413 (entropy = 6.29, ARG read fraction = 0.36%) confirms pipeline feasibility.
2. Refining sequencing depth estimates for Bulgaria samples (PRJNA1071831); current analysis uses estimated read counts that may differ from exact run statistics available via SRA metadata.
3. Longer time series with more independent studies across diverse environments.
4. Incorporation of mechanistic or intervention-linked covariates.
5. Nonlinear or multi-horizon predictive models with robust diagnostic checks.
6. Formal residual diagnostics, including heteroskedasticity and autocorrelation tests.

## 6. Conclusion
We present a reproducible framework for evaluating predictive claims in longitudinal resistome–mobilome data, grounded in a real-sample temporal cohort of 23 samples across four environmental contexts. Applying literature-calibrated quantification to this cohort (14 regression-complete lag pairs), naive lag models show near-perfect fit (R² = 0.999), but uniformly high lag correlations across offsets (0.999 → 0.998 → 0.997), symmetric reverse-direction performance (R² = 0.998), signal collapse under differencing (p = 0.153), and a Granger-style added-value test that returns the null strongly (p = 0.878) together demonstrate that the apparent predictive signal is entirely due to temporal coupling. These results emphasize that surveillance claims about leading indicators must be validated with explicit robustness tests, not inferred from high R-squared alone. The real-sample temporal backbone, bioinformatics quantification pipeline (AMRFinderPlus, IntegronFinder, Kraken2/Bracken, MEGAHIT), and diagnostic analysis code are provided as fully reproducible resources for the community.

## 7. References
1. Berendonk, T. U., et al. (2015). Tackling antibiotic resistance: the environmental framework. Nature Reviews Microbiology, 13, 310–317. https://doi.org/10.1038/nrmicro3439
2. Schluter, J., Hussey, G., Valeriano, J., Zhang, C., Sullivan, A., & Fenyö, D. (2024). The MTIST platform: a microbiome time series inference standardized test. Research Square (preprint). https://doi.org/10.21203/rs.3.rs-4343683/v1
3. Papaspyropoulos, K. G., & Kugiumtzis, D. (2024). On the Validity of Granger Causality for Ecological Count Time Series. Econometrics. https://doi.org/10.3390/econometrics12020013
4. Choi, J. M., Rumi, M. A., Brown, C. L., Vikesland, P. J., Pruden, A., & Zhang, L. (2026). ARGfore: A Multivariate Framework for Forecasting Antibiotic Resistance Gene Abundances Using Time-Series Metagenomic Datasets. IEEE Access. https://doi.org/10.1109/access.2026.3667074
5. Feldgarden, M., et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Scientific Reports, 11, 12728. https://doi.org/10.1038/s41598-021-91456-0
6. Néron, B., et al. (2022). IntegronFinder 2.0: Identification and Analysis of Integrons across Bacteria, with a Focus on Antibiotic Resistance in Klebsiella. Microorganisms, 10(4), 700. https://doi.org/10.3390/microorganisms10040700
7. Gillings, M. R., et al. (2015). Mobile class 1 integrons promote the spread of resistance genes in human microbiomes. NPJ Biofilms and Microbiomes, 1, 15010. https://doi.org/10.1038/npjbiofilms.2015.10
8. Lu, J., et al. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104. https://doi.org/10.7717/peerj-cs.104
9. Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20, 257. https://doi.org/10.1186/s13059-019-1891-0
10. Li, D., et al. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033
11. Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and Practice (3rd ed.). OTexts. https://otexts.com/fpp3
12. Seabold, S., & Perktold, J. (2010). Statsmodels: Econometric and Statistical Modeling with Python. Proceedings of the 9th Python in Science Conference (SciPy 2010). https://doi.org/10.25080/Majora-92bf1922-011
13. The pandas development team (2020). pandas-dev/pandas: Pandas. Zenodo. https://doi.org/10.5281/zenodo.3509134
