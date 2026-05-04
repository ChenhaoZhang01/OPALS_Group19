<p align="center"><b>When High R² Misleads: Diagnosing False Predictive Signals in Resistome–Mobilome Time Series</b></p>
<p align="center"><b>Chenhao Zhang¹, Eliana Wong¹*, Ashley Fang¹*, Wanze Tang¹*, Linda Shi², Yujie Men²</b></p>
<p align="center">¹Institute of Engineering in Medicine, University of California San Diego, La Jolla, CA 92093<br>²Department of Chemical and Environmental Engineering, University of California, Riverside, CA 92521<br>
*High school students participating in IEM OPALS program</p>

**Abstract -** *A lag regression linking mobile genetic element (MGE) abundance at time t to antibiotic resistance gene (ARG) burden at t+1 is a natural early-warning model for resistance surveillance. We demonstrate that such models can produce misleadingly high apparent predictive performance driven by between-study confounding rather than genuine temporal coupling. From 19 shotgun metagenomes across four environmental time series (Chesapeake Bay WGA and two Bulgaria Iskar River sites), we derive 14 regression-complete lag pairs using a complete bioinformatics pipeline. The naive forward lag model yields R² = 0.720, yet six independent diagnostic tests all return null: differencing collapses the signal (R²=0.020), Granger added-value is negligible (delta R²=0.014, p=0.539), and the within-study coefficient is negative in the densest series (−33.0, n=10). A calibration simulation (N=2,000) shows naive lag models yield R²>0.70 in 94% of null runs (median 0.848); the observed R²=0.720 falls below the null median. The central finding is that standard lag model evidence cannot distinguish null from true-signal conditions in multi-study metagenome designs. The diagnostic battery can, and code is provided as a reproducible template.*

**Keywords:** ARG surveillance, mobilome, lag regression, observational equivalence, temporal confounding, Granger-style test, calibration simulation, time-series diagnostics

## 1. Introduction

Environmental ARG surveillance increasingly seeks leading indicators — variables at time t predicting ARG burden at t+1 — and mobile genetic elements (MGEs) are a natural candidate, since they mediate horizontal gene transfer and are frequently co-detected with ARGs [1]. The methodological risk is that lag regressions linking MGE(t) to ARG(t+1) can yield high R² and plausible coefficients even when the predictor contributes no independent signal, due to shared autocorrelation and between-study ARG scale differences. A practitioner stopping at R²=0.72 would rationally conclude that MGE monitoring provides early warning and design surveillance infrastructure accordingly — but if the apparent signal is entirely confounded, those decisions rest on a false foundation.

We present what is, to our knowledge, the first empirical demonstration that standard lag model evidence is observationally equivalent under null and true-signal conditions in a real metagenomic resistome dataset, paired with a calibration simulation that proves the diagnostic battery can distinguish the two. We apply a complete bioinformatics pipeline to 19 real shotgun metagenome samples and subject the lag model to a nine-test diagnostic battery. The broader microbiome time-series literature has emphasized the need for standardized evaluation [2]; recent work has re-examined when Granger-style inference is valid for ecological count data [3]; and emerging forecasting tools [4] make reliable evaluation increasingly consequential. The framework introduced here operationalizes these standards in a resistome context.

## 2. Methods

### 2.1. Data Sources and Temporal Cohort

We assembled a temporal cohort from two published studies with known, ordered collection timestamps. **PRJNA599167 (Chesapeake Bay, USA, May 2017)** provides 11 WGA shotgun metagenomes at ~1.5-hour intervals over a 48-hour window (4.1–17.7 M paired reads each; SRR11803495–SRR11803505) and 4 paired 16S amplicon samples; shotgun samples support full ARG/MGE/diversity quantification while amplicon samples contribute diversity only. **PRJNA1071831 (Iskar River, Bulgaria, Nov–Dec 2022)** provides 4 biweekly shotgun metagenomes each from Dragushinovo and Mechkata sites (inter-sample gaps: 14, 21, and 13 days); each site is treated as an independent study. In total, 19 strict lag pairs are identified across four sampling contexts (WGA: 10; Drag: 3; Mech: 3; 16S: 3), exceeding the minimum threshold of 10 pairs from ≥3 independent studies.

### 2.2. Quantification Pipeline

For shotgun samples (WGA and Iskar River), we quantify:
- **arg_total**: total ARG hit count per sample using DIAMOND blastx [5] against the AMRProt database (NCBI AMRFinderPlus reference, ≥80% identity, ≥80% query coverage, e-value ≤1×10⁻⁵) applied directly to paired-end reads.
- **mge_abundance**: count of complete and CALIN integrons detected by IntegronFinder 2.0 [6] on assembled metagenomic contigs ≥4 kb. Class 1 integrons are the dominant ARG-mobilizing MGE in environmental water metagenomics [7].
- **entropy**: Shannon diversity index (H') computed from species-level abundance fractions produced by Bracken [8] re-estimation of Kraken2 [9] classifications against the 8 GB standard database.

Raw FASTQ files are downloaded from SRA using sra-tools (prefetch + fasterq-dump). Assembly uses MEGAHIT [10] with minimum contig length 500 bp; IntegronFinder is applied to contigs ≥4 kb to minimize run time on high-contig libraries. Sequencing depths for the Bulgaria samples were retrieved from the NCBI SRA EUtils API. The complete pipeline is provided in `analysis/machine3_autonomous.sh`.

### 2.3. Feature Table and Preprocessing

Observations are ordered within study by sampling time. Lagged predictors (ARG(t+1), MGE(t+1)) and first-difference variables (dMGE(t) = MGE(t)−MGE(t−1); dARG(t+1) = ARG(t+1)−ARG(t)) are computed within each study to remove shared trends [11]. A key structural heterogeneity is the sampling interval: WGA samples are at ~1.5-hour intervals while Bulgaria samples are biweekly (13–21 days apart), so all lag regressions reflect "one sampling interval ahead" rather than a fixed time horizon.

### 2.4. Models

We evaluate: (1) forward lag model ARG(t+1) ~ MGE(t) + C(study) + sequencing_depth; (2) reverse-direction model MGE(t+1) ~ ARG(t) + C(study); (3) first-difference model dARG(t+1) ~ dMGE(t); (4) Granger-style added-value test comparing base model ARG(t+1) ~ ARG(t) + C(study) vs. full model adding MGE(t); (5) within-study models ARG(t+1) ~ MGE(t) fit independently per study; (6) leave-one-study-out (LOSO) refits; (7) residual diagnostics (Durbin-Watson, Breusch-Pagan, power analysis); and (8) a null simulation (N=2,000 synthetic datasets, same study structure, no true MGE effect). Models use statsmodels [12] and pandas [13].

## 3. Results

All results are derived from fully pipeline-quantified values. Of 23 strict-cohort samples, 19 were fully quantified using the complete bioinformatics pipeline (excluding 4 amplicon-only samples); 14 lag pairs were regression-complete after filtering for sequencing-depth covariate availability.

**Pipeline validation.** For SRR27827413 (Dragushinovo, 22.2M reads): H' = 6.301, arg_total = 56,046, mge_abundance = 445 integrons — consistent with comparable wastewater-impacted river environments.

### 3.1. Forward Lag Model

The forward lag model yields coefficient = 21.57, p-value = 0.579, R² = 0.720 (n=14). R² is driven by between-study ARG scale differences, not MGE predictive content.

<p align="center"><img src="analysis/figures/scatter_mge_t_vs_arg_t1.png" width="750"></p>
<p align="center">Fig. 1: MGE(t) vs. ARG(t+1) colored by study. The pooled OLS fit (black dashed) produces apparent R²=0.720, but this reflects between-study ARG scale differences (three distinct clusters) rather than within-study predictive signal. Dotted lines show within-study fits, which are flat or negative.</p>

<p align="center"><img src="analysis/figures/timeseries_mge_arg_by_study.png" width="750"></p>
<p align="center">Fig. 2: Study-specific ARG (solid) and MGE (dashed) time series on dual y-axes. ARG baselines differ ~3× between the WGA series and Bulgaria sites, generating the between-study variance that inflates pooled R².</p>

### 3.2. Lag Correlation Diagnostics

<p align="center">Table 1: Lag correlations across offsets.</p>

| Lag | Correlation | n |
|---|---:|---:|
| t→t | 0.400 | 19 |
| t→t+1 | 0.532 | 16 |
| t→t+2 | 0.214 | 13 |

The non-monotone pattern (rising then falling) is inconsistent with uniform temporal coupling and does not support a directional predictive signal from MGE to ARG.

<p align="center"><img src="analysis/figures/lag_comparison_correlation.png" width="650"></p>
<p align="center">Fig. 3: Lag correlation across offsets. A genuine leading indicator would show monotonically increasing correlation with lag offset; the non-monotone pattern (0.400 → 0.532 → 0.214) is inconsistent with directional prediction.</p>

The reverse-direction model MGE(t+1) ~ ARG(t) yields R²=0.390, p=0.164 — also non-significant — confirming the forward model's higher R² reflects between-study ARG scale differences, not directionality.

### 3.3. Differenced Model

<p align="center">Table 3: First-difference and Granger-style results.</p>

| Test | Coefficient / Metric | p-value | R-squared | n |
|---|---:|---:|---:|---:|
| dARG(t+1) ~ dMGE(t) | 12.48 | 0.649 | 0.020 | 13 |
| Granger delta R² | 0.014 | 0.539 | — | 14 |

After trend removal, R² collapses to 0.020; adding MGE to a model that already includes ARG history adds delta R²=0.014 (p=0.539). Both confirm the forward model's apparent fit carries no independent MGE predictive content.

<p align="center"><img src="analysis/figures/differenced_scatter_dmge_t_vs_darg_t1.png" width="800"></p>
<p align="center">Fig. 4: Effect of differencing. Left: raw scatter (R²=0.720) showing study-cluster structure. Right: after differencing (R²=0.020), no signal remains.</p>

### 3.4. Within-Study Models

<p align="center">Table 5: Within-study lag models.</p>

| Study | Coefficient | p-value | R² | n (pairs) |
|---|---:|---:|---:|---:|
| PRJNA1071831_Drag | 31.0 | 0.687 | 0.223 | 3 |
| PRJNA1071831_Mech | 94.5 | 0.812 | 0.085 | 3 |
| PRJNA599167_WGA | −33.0 | 0.506 | 0.057 | 10 |

With n=10 lag pairs in the WGA series, the MGE coefficient is negative (−33.0, p=0.506, R²=0.057). Within the Chesapeake Bay series, higher MGE abundance at time t is associated with slightly lower, not higher, ARG at t+1 — inconsistent with a leading-indicator relationship.

### 3.5. Leave-One-Study-Out Sensitivity

<p align="center">Table 6: Leave-one-study-out (LOSO) results.</p>

| Left-out study | Coefficient | p-value | R² | n |
|---|---:|---:|---:|---:|
| PRJNA1071831_Drag | 4.2 | 0.958 | 0.686 | 11 |
| PRJNA1071831_Mech | 20.8 | 0.524 | 0.759 | 11 |
| PRJNA599167_WGA | 31.8 | 0.670 | 0.358 | 6 |

MGE is non-significant in all three variants (p = 0.958, 0.524, 0.670), and excluding WGA drops R² from 0.720 to 0.358. Residual diagnostics on the forward model confirm adequate specification: Durbin-Watson = 2.10 and Breusch-Pagan p = 0.303 (no autocorrelation or heteroskedasticity). The minimum detectable MGE coefficient at 80% power is 117.9 — five times the observed estimate of 21.6 — so the diagnostic battery triangulates across multiple approaches rather than relying on any single test.

### 3.6. Calibration Simulation: Null vs. True Signal Scenarios

To test whether the framework can distinguish genuine signal from confounding, we ran 2,000 simulations under two scenarios (WGA: n=10 pairs; Drag: n=3; Mech: n=3; ARG baselines calibrated to real data):

- **Scenario A (null):** MGE drawn independently of future ARG (beta=0). No true predictive signal.
- **Scenario B (true signal):** ARG(t+1) = study_baseline + 150·MGE(t) + noise. A true causal effect of realistic magnitude.

<p align="center">Table 7: Two-scenario simulation results (N=2,000 synthetic datasets each).</p>

| Metric | Scenario A (null, β=0) | Scenario B (true signal, β=150) | Observed real data |
|---|---:|---:|---:|
| Naive model R² median | 0.856 | 0.940 | 0.720 |
| % naive runs with R² > 0.70 | 94% | 99% | — |
| % naive runs with MGE p < 0.05 | 11% | 90% | p=0.579 |
| Differenced model R² median | 0.070 | 0.692 | 0.020 |
| Within-WGA coeff median | −0.1 | +149 | −33.0 |
| % within-WGA coeff positive | 50% | 98% | — |

The naive model yields R²>0.70 in 94% of null runs and 99% of true-signal runs — it cannot distinguish the two scenarios. The differenced model and within-study coefficient decisively separate them. The observed real data — differenced R²=0.020, within-WGA coeff=−33.0 — matches Scenario A closely and is entirely inconsistent with Scenario B.

<p align="center"><img src="analysis/figures/simulation_false_signal.png" width="850"></p>
<p align="center">Fig. 5: Calibration simulation. Top row: naive lag model R² distributions — both scenarios produce R²>0.70 most of the time; the observed 0.720 (orange dotted line) falls in the Scenario A range. Bottom row: differenced model R² distributions — Scenario A clusters near 0 (observed 0.020 matches) while Scenario B spreads toward 1 (median 0.69). The naive model is uninformative; the differenced model correctly separates null from signal.</p>

## 4. Discussion

The core contribution is an empirical demonstration of observational equivalence: in multi-study environmental metagenome designs, naive lag regression models produce results that are statistically indistinguishable under null and true-signal conditions. The real dataset yields a naive lag model R²=0.720 with a plausible MGE coefficient; the calibration simulation shows that 94% of null datasets and 99% of true-signal datasets generate R²>0.70 under the same study structure. Standard empirical evidence cannot place this result in either category.

Seven independent diagnostics converge on null (Tables 1–7): non-significant pooled coefficient (p=0.579), non-monotone lag correlations, collapsed differenced R² (0.020), negligible Granger added-value (delta R²=0.014), negative WGA within-study coefficient (−33.0), non-significant LOSO across all variants, and 94% of null simulations yielding R²>0.70 — placing the observed 0.720 below the null median. Convergence across seven approaches makes a genuine positive relationship unlikely rather than merely undetectable by any single test.

The calibration simulation proves the framework has discriminatory power: in Scenario B, differenced R² rises to a median of 0.69 and the within-WGA coefficient recovers the true beta (median=149, 98% positive). The real data — differenced R²=0.020, within-WGA coeff=−33 — matches Scenario A on both criteria. A researcher who observed only R²=0.720 would have no basis for this classification; one applying the full battery would.

The practical consequence extends to existing forecasting frameworks such as ARGfore [4]: applying differencing and within-location estimation to published results would directly test whether high R² reflects genuine temporal coupling or between-study scale differences. Measurement conditions here — integron-only MGE proxy, MDA amplification artifacts in WGA samples, heterogeneous sequencing depths — attenuate true coefficients toward zero, so consistent null results hold despite conditions that work against the null. Null results across both the 1.5-hour WGA and biweekly Bulgaria regimes strengthen this conclusion.

The principal limitation is small sample size (n=14 lag pairs), which the diagnostic battery compensates for by triangulating across six independent approaches. The cohort spans only two environmental systems; Bulgaria within-study series (n=3 pairs each) are too short for reliable estimation. Future work should extend to longer time series, mechanistic covariates, and application of this battery to published forecasting datasets.

## 5. Conclusion

**High R² is not evidence of predictive signal in multi-study resistome surveillance data.** Naive lag models yield R²>0.70 in 94% of null runs (median 0.856); the observed R²=0.720 falls below the null median. Seven diagnostics return null while the calibration simulation confirms the framework has power to detect real signals. The practical recommendation is explicit: before deploying a lag regression as an early-warning system, require that (1) the predictor coefficient is significant after trend removal and (2) within-study predictive signal exists independently of between-study scale differences. Neither condition is met here. The full pipeline and diagnostic code are provided as a reproducible template.

## Acknowledgements

The authors thank the OPALS program at the Institute of Engineering in Medicine, UC San Diego, for supporting this research. C.Z., L.S., and Y.M. designed the study and developed the analysis framework. E.W., A.F., and W.T. contributed to data analysis and interpretation. Pipeline development and computational work were conducted using publicly available SRA data from NCBI BioProjects PRJNA599167 and PRJNA1071831.

## References

[1] T. U. Berendonk et al., "Tackling antibiotic resistance: the environmental framework," Nature Reviews Microbiology, vol. 13, pp. 310–317, 2015. doi: 10.1038/nrmicro3439

[2] J. Schluter, G. Hussey, J. Valeriano, C. Zhang, A. Sullivan, and D. Fenyö, "The MTIST platform: a microbiome time series inference standardized test," Research Square, preprint, 2024. doi: 10.21203/rs.3.rs-4343683/v1

[3] K. G. Papaspyropoulos and D. Kugiumtzis, "On the validity of Granger causality for ecological count time series," Econometrics, vol. 12, no. 2, p. 13, 2024. doi: 10.3390/econometrics12020013

[4] J. M. Choi, M. A. Rumi, C. L. Brown, P. J. Vikesland, A. Pruden, and L. Zhang, "ARGfore: A multivariate framework for forecasting antibiotic resistance gene abundances using time-series metagenomic datasets," IEEE Access, 2026. doi: 10.1109/access.2026.3667074

[5] B. Buchfink, K. Reuter, and H.-G. Drost, "Sensitive protein alignments at tree-of-life scale using DIAMOND," Nature Methods, vol. 18, pp. 366–368, 2021. doi: 10.1038/s41592-021-01101-x

[6] B. Néron et al., "IntegronFinder 2.0: identification and analysis of integrons across bacteria, with a focus on antibiotic resistance in Klebsiella," Microorganisms, vol. 10, no. 4, p. 700, 2022. doi: 10.3390/microorganisms10040700

[7] M. R. Gillings et al., "Mobile class 1 integrons promote the spread of resistance genes in human microbiomes," NPJ Biofilms and Microbiomes, vol. 1, p. 15010, 2015. doi: 10.1038/npjbiofilms.2015.10

[8] J. Lu et al., "Bracken: estimating species abundance in metagenomics data," PeerJ Computer Science, vol. 3, p. e104, 2017. doi: 10.7717/peerj-cs.104

[9] D. E. Wood, J. Lu, and B. Langmead, "Improved metagenomic analysis with Kraken 2," Genome Biology, vol. 20, p. 257, 2019. doi: 10.1186/s13059-019-1891-0

[10] D. Li et al., "MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph," Bioinformatics, vol. 31, no. 10, pp. 1674–1676, 2015. doi: 10.1093/bioinformatics/btv033

[11] R. J. Hyndman and G. Athanasopoulos, Forecasting: Principles and Practice, 3rd ed. OTexts, 2021. [Online]. Available: https://otexts.com/fpp3

[12] S. Seabold and J. Perktold, "Statsmodels: econometric and statistical modeling with Python," in Proc. 9th Python in Science Conf. (SciPy 2010), 2010. doi: 10.25080/Majora-92bf1922-011

[13] The pandas development team, "pandas-dev/pandas: Pandas," Zenodo, 2020. doi: 10.5281/zenodo.3509134
