<p align="center"><b>Pipeline or Place? The Bioinformatics Pipeline Is the Dominant Source of Variation in Environmental Antibiotic Resistance Gene Estimates</b></p>
<p align="center"><b>Chenhao Zhang¹, Eliana Wong¹*, Ashley Fang¹*, Wanze Tang¹*, Linda Shi², Yujie Men²</b></p>
<p align="center">¹Institute of Engineering in Medicine, University of California San Diego, La Jolla, CA 92093<br>²Department of Chemical and Environmental Engineering, University of California, Riverside, CA 92521<br>
*High school students participating in IEM OPALS program</p>

**Abstract -** *Metagenomic antibiotic resistance gene (ARG) abundance is increasingly used to rank environments and set surveillance baselines, yet the same metagenome yields different ARG estimates depending on the bioinformatics pipeline applied. We quantify how much of the total variation in ARG estimates is attributable to the pipeline versus the environment versus their interaction across 19 environmental shotgun metagenomes balanced over four habitats (wastewater, river, irrigation water, soil) drawn from NCBI BioProjects. Each sample was downloaded and processed end-to-end through three independent ARG-calling pipelines — (A) MEGAHIT assembly + Prodigal + DIAMOND against CARD, (B) Bowtie2 read mapping to CARD, and (C) RGI ORF-level calling against CARD — and the normalized abundances were partitioned by Type-II ANOVA on the log scale. Contrary to the intuition that biology dominates, the pipeline is the single largest source of variance, explaining 57.0% of variation in log ARG abundance (p = 3.5×10⁻¹⁷) versus 18.9% for environment, with a significant pipeline × environment interaction (11.0%, p = 6.7×10⁻⁵). Pairwise pipelines differ by 7- to 152-fold at the median and rank samples inconsistently (Spearman ρ = 0.19–0.71), and they disagree on the environmental gradient itself: read mapping shows a strong habitat gradient while ORF-level calling is nearly habitat-flat. We additionally find that RGI "strict" returns zero ARG calls on every environmental sample — environmental ARG homologs fall below CARD's clinically-curated cutoffs — so the choice of detection stringency alone moves a sample from "no resistance" to hundreds of calls. We conclude that environmental ARG abundances are not comparable across studies that use different pipelines, and provide the full pipeline, decomposition, and figures as a reproducible template.*

**Keywords:** antibiotic resistance genes, metagenomics, bioinformatics pipeline, variance decomposition, ARG surveillance, method comparison, reproducibility, CARD

## 1. Introduction

Environmental metagenomic surveys routinely report antibiotic resistance gene (ARG) abundance as a single number per sample — total ARG hits normalized to sequencing depth — and use it to rank habitats, track pollution gradients, and establish surveillance baselines [1], [2]. The implicit assumption is that this number is a property of the sample. It is not: it is a property of the sample *and* the bioinformatics pipeline used to measure it. Assembly-based gene calling, direct read mapping to a reference database, and ORF-level resistance callers apply different sensitivity/specificity trade-offs and can disagree several-fold — or, as we show, by orders of magnitude — on the same data [3], [4].

This matters because ARG estimates are increasingly compared *across* studies that used *different* pipelines. If pipeline-attributable variation is comparable to or larger than the environmental signal those studies seek to detect, then a difference attributed to "more wastewater impact" may instead reflect "a more permissive ARG caller." The methodological question is therefore not which pipeline is correct, but how the total variance in ARG estimates partitions among (i) the environment, the biological signal of interest; (ii) the pipeline, a measurement-method nuisance; and (iii) their interaction, which determines whether a single pipeline-correction factor could ever harmonize methods across habitats.

We assembled a balanced 19-sample, four-habitat cohort, downloaded the raw reads, processed each sample through three independent ARG-calling pipelines, and decomposed the variance in the resulting estimates using `ARG_total ~ pipeline + environment + pipeline:environment`. We report the partition for both total ARG abundance and ARG richness, quantify pairwise concordance, and show that the pipeline is the dominant variance component and that the pipelines disagree even on the *shape* of the environmental gradient.

## 2. Methods

### 2.1. Sample Cohort

We selected 19 publicly archived environmental shotgun metagenomes from the NCBI Sequence Read Archive, balanced across four habitats — wastewater (n = 4), river (n = 5), irrigation water (n = 5), and soil (n = 5) — drawn from multiple BioProjects (accessions in `metadata/batch_subset.txt`; full metadata in `metadata/metadata_final.csv`). Habitats span the anthropogenic ARG-input gradient from heavily impacted wastewater to low-background soil. Native sequencing depths ranged from 3×10⁶ to 5.9×10⁷ reads.

### 2.2. Three ARG-Calling Pipelines

Each sample was processed independently through three pipelines representing the dominant families of metagenomic ARG quantification, all against the same reference, CARD v3.2.7 [8]:

- **Pipeline A — assembly-based.** Reads assembled with MEGAHIT [5] (minimum contig length 500 bp); ORFs predicted with Prodigal [6] in metagenome mode; predicted proteins searched with DIAMOND blastp v2.1.11 [7] against CARD protein homolog models at ≥80% identity and ≥70% query coverage. Specific, but loses low-coverage ARGs that fail to assemble.
- **Pipeline B — read mapping.** Reads mapped directly to the CARD nucleotide homolog reference with Bowtie2 [9]; ARG counts taken from `samtools idxstats` mapped-read counts. Sensitive, but counts fragment-level hits.
- **Pipeline C — ORF-level calling (RGI).** Predicted ORFs screened with RGI 6.0.5 against CARD. RGI "strict" (Perfect/Strict criteria) returned **zero** calls on every environmental sample (see §3.5), so Pipeline C reports the Loose-inclusive tier, the standard setting for divergent environmental homologs; strict ≈ 0 is itself reported as a finding.

For every pipeline, ARG hits were normalized to library size (`normalized abundance = ARG hits / total reads`); `ARG_total` is the sum of normalized abundance across detected CARD drug-class features (CARD `aro_index.tsv` used to collapse hits to drug classes), and `ARG_richness` is the number of distinct drug-class features detected.

### 2.3. Processing Decisions (Real-Data Constraints)

All samples were processed on a single Linux/WSL2 workstation (11 GB RAM). Two constraints, fully documented in `results/DATA_PROVENANCE.md`, shape the absolute (not relative) scale of the estimates and are revisited as limitations:

1. **Equal-depth subsampling to ~12M reads** (first 6×10⁶ spots streamed via `fastq-dump -X`). The deep soil/irrigation libraries cannot be assembled at full depth within 11 GB RAM; streaming a fixed number of spots bounds both memory and download and controls sequencing depth as a confounder. A consequence is that assembly-based detection (Pipeline A) is depth-limited and conservative.
2. **Minimal QC** (no separate adapter/quality trimming): writing uncompressed trimmed FASTQ was the dominant I/O cost on the workstation, and read mapping and assembly are robust to light residual adapter content.

### 2.4. Variance Decomposition

Normalized ARG abundance spans orders of magnitude and pipeline/environment effects are multiplicative, so variance was partitioned on the log scale, where multiplicative effects become additive — the standard scale for metagenomic abundance modeling [11]. Samples with zero detections (a few assembly-pipeline samples) were handled with a pseudocount of half the smallest positive normalized abundance before the log₁₀ transform. We fit the ordinary-least-squares model

```
log10(ARG_total) ~ C(pipeline) + C(environment) + C(pipeline):C(environment)
```

and computed a Type-II ANOVA [12]; each factor's variance contribution is its sum of squares as a percentage of the total. The identical partition was computed for `ARG_richness`. Pairwise concordance was assessed by Spearman rank correlation of per-sample `ARG_total` between pipeline pairs, with the median fold-difference. Analyses use pandas [13], statsmodels [12], and SciPy. The full pipeline is in `analysis/run_arg_pipeline_hpc.sh`; the analysis battery in `analysis/run_paper1_analysis.py`; figures in `analysis/generate_figures.py`.

## 3. Results

### 3.1. The Pipeline Shifts the Whole Abundance Distribution by Orders of Magnitude

The three pipelines produce abundance estimates separated by orders of magnitude on the same samples (Fig. 1, Table 1). ORF-level RGI (C) returns the highest normalized abundance (median 6.0×10⁻⁵), read mapping (B) intermediate (2.7×10⁻⁶), and the assembly pipeline (A) the lowest (3.9×10⁻⁷) — a ~150-fold spread between the most and least permissive method on identical input.

<p align="center"><img src="analysis/figures/fig1_abundance_by_pipeline.png" width="720"></p>
<p align="center">Fig. 1: Normalized ARG abundance (log₁₀) by pipeline, each sample colored by environment. The three pipelines occupy non-overlapping abundance bands — the choice of pipeline shifts the estimate more than the habitat does.</p>

<p align="center">Table 1: Per-pipeline summary across the 19-sample cohort.</p>

| Pipeline | Median norm. ARG_total | Mean ARG richness (classes) |
|---|---:|---:|
| A: assembly + DIAMOND | 3.9×10⁻⁷ | 5.9 |
| B: read mapping → CARD | 2.7×10⁻⁶ | 17.1 |
| C: RGI (loose) | 6.0×10⁻⁵ | 22.1 |

### 3.2. Pipeline Is the Dominant Variance Component

The Type-II ANOVA partition of log ARG abundance (Table 2, Fig. 2) assigns **57.0%** of variance to pipeline, **18.9%** to environment, **11.0%** to the pipeline × environment interaction, and **13.0%** to residual variation. All three model terms are highly significant. The headline result reverses the common assumption: the bioinformatics pipeline, not the environment, is the single largest driver of ARG abundance estimates — it explains three times as much variance as the habitat. A cross-study comparison that ignores the pipeline is therefore confounded by a nuisance term larger than the biological effect it seeks to measure.

<p align="center">Table 2: Variance decomposition of log₁₀(ARG_total), Type-II ANOVA (n = 57 sample×pipeline observations).</p>

| Factor | df | F | p | Variance % |
|---|---:|---:|---:|---:|
| pipeline | 2 | 98.7 | 3.5×10⁻¹⁷ | 57.0 |
| environment | 3 | 21.8 | 7.1×10⁻⁹ | 18.9 |
| pipeline × environment | 6 | 6.36 | 6.7×10⁻⁵ | 11.0 |
| residual | 45 | — | — | 13.0 |

<p align="center"><img src="analysis/figures/fig2_variance_partition.png" width="620"></p>
<p align="center">Fig. 2: Variance partition for log ARG abundance (left) and ARG richness (right). Pipeline is the largest component for abundance; for richness, pipeline, environment, and their interaction are all substantial.</p>

### 3.3. Pipelines Disagree on the Environmental Gradient Itself

The significant pipeline × environment interaction (11.0%, p = 6.7×10⁻⁵) means the methods do not merely sit at different absolute levels — they report different *shapes* of the habitat gradient (Fig. 3). Read mapping (B) is the most environment-responsive, peaking sharply in river and wastewater and dropping ~100-fold in soil and irrigation. ORF-level RGI (C) is nearly habitat-flat, reporting similar abundance everywhere. The assembly pipeline (A) tracks B's shape at a lower level. Consequently, a researcher using read mapping would report a strong environmental ARG gradient, while one using RGI on the same samples would report almost none. No single multiplicative "pipeline correction factor" can reconcile methods whose environmental responses differ in shape, not just offset.

<p align="center"><img src="analysis/figures/fig3_pipeline_env_interaction.png" width="700"></p>
<p align="center">Fig. 3: Pipeline × environment interaction. Cell means (±95% CI) of log ARG abundance by environment, one line per pipeline. The lines are strongly non-parallel — read mapping (B) carries a steep habitat gradient while RGI (C) is nearly flat — so the apparent environmental signal depends on the pipeline.</p>

### 3.4. Richness Is Jointly Driven by Pipeline, Environment, and Their Interaction

For ARG richness (number of drug-class features detected), variance partitions more evenly but is still pipeline-led: pipeline **29.6%**, environment **24.7%**, interaction **27.6%**, residual **18.1%** (Table 3, Fig. 4). RGI (loose) detects the most classes and most consistently (mean 22.1, tight range), read mapping the most variably (mean 17.1, range 1–55), and assembly the fewest (mean 5.9). The qualitative question "which resistance classes are present?" is thus answered as much by the software as by the habitat — a direct threat to presence/absence resistome comparisons.

<p align="center">Table 3: Variance decomposition of ARG richness, Type-II ANOVA.</p>

| Factor | df | F | p | Variance % |
|---|---:|---:|---:|---:|
| pipeline | 2 | 36.8 | 3.4×10⁻¹⁰ | 29.6 |
| environment | 3 | 20.5 | 1.6×10⁻⁸ | 24.7 |
| pipeline × environment | 6 | 11.4 | 9.8×10⁻⁸ | 27.6 |
| residual | 45 | — | — | 18.1 |

<p align="center"><img src="analysis/figures/fig4_richness_by_pipeline.png" width="620"></p>
<p align="center">Fig. 4: ARG richness per sample by pipeline. RGI (loose) detects the most classes and most consistently; the assembly pipeline detects the fewest; read mapping is the most variable.</p>

### 3.5. Detection Stringency Alone Moves a Sample From Zero to Hundreds of Calls

RGI run in "strict" mode (Perfect/Strict criteria only) returned **zero** ARG calls on every one of the 19 environmental samples. The same predicted-ORF input scored under the Loose-inclusive tier yielded tens to thousands of calls (e.g. one wastewater sample: 135,325 ORFs → 0 strict, 750 loose). Environmental ARG homologs systematically fall below CARD's clinically-curated bitscore cutoffs, so a single configuration flag in one tool moves a sample from "no detectable resistance" to a rich resistome. This is the most extreme instance of pipeline-attributable variance in the study and is reported here as a finding in its own right.

### 3.6. Pipelines Rank Samples Inconsistently and Differ by Orders of Magnitude

Pairwise concordance is weak to moderate and the absolute disagreement is large (Table 4). Assembly and RGI rank samples similarly (ρ = 0.71) yet differ 152-fold at the median; read mapping and RGI are essentially unranked relative to each other (ρ = 0.19, p = 0.44) and differ 22-fold. Only assembly vs. read mapping are even modestly concordant (ρ = 0.48). Large fold-differences combined with low rank correlation are the signature of method-driven disagreement that cannot be removed by rescaling — cross-study numeric comparison fails even when within-study rankings are taken at face value.

<p align="center">Table 4: Pairwise pipeline concordance (per-sample ARG_total, n = 19).</p>

| Pipeline pair | Spearman ρ | p | Median fold-difference |
|---|---:|---:|---:|
| A vs. B | 0.48 | 0.036 | 6.8× |
| A vs. C | 0.71 | 5.8×10⁻⁴ | 152× |
| B vs. C | 0.19 | 0.44 | 22× |

## 4. Discussion

The central finding is quantitative and, against expectation, the pipeline dominates: in a balanced four-habitat cohort processed end-to-end on real reads, the bioinformatics pipeline explains 57.0% of the variance in ARG abundance estimates — three times the environmental signal (18.9%) — and is the leading component of richness variance (29.6%). ARG abundance is therefore not a transferable, pipeline-independent property of a metagenome. A study reporting that habitat X has more ARGs than habitat Y is making a claim valid only within a fixed pipeline; the same comparison performed with a different pipeline can change the magnitude by up to two orders of magnitude (Table 4) and, in the strict-vs-loose case (§3.5), change a zero into hundreds of calls.

The significant pipeline × environment interaction sharpens this into a qualitative warning: the methods disagree not only on level but on the *shape* of the environmental gradient (Fig. 3). Read mapping reports a steep habitat gradient; ORF-level calling reports almost none. The biological conclusion — "is there an environmental ARG gradient here?" — is thus partly an artifact of the chosen pipeline. This rules out the common hope that a single empirically derived conversion factor could harmonize ARG estimates across pipelines: any such factor would itself have to be habitat-specific.

These results align with and extend prior method-comparison work reporting several-fold disagreement among ARG callers [3], [4]: by partitioning disagreement against a real environmental gradient we show it is not merely scatter but a dominant, structured variance component that can invert biological conclusions. The practical recommendations are concrete: (1) ARG abundance comparisons across studies are interpretable only when the pipeline and its settings are identical, or when pipeline-attributable variance is explicitly modeled as here; (2) presence/absence and richness comparisons must report the caller and its stringency tier, since these dominate and can yield zero; and (3) studies should publish per-pipeline uncertainty, for which this decomposition is a reusable template.

**Limitations.** This is a deliberately transparent real-data study run on a single workstation, and several decisions affect the *absolute* scale of the estimates and partly inflate the pipeline component (all detailed in `results/DATA_PROVENANCE.md`). First, equal-depth subsampling to ~12M reads makes the assembly pipeline (A) depth-limited and conservative, widening the A-vs-others gap; full-depth assembly would raise A. Second, QC was minimal. Third, Pipeline C uses RGI's Loose-inclusive tier because strict returns zero; the loose tier is more permissive than the other two methods, contributing to C's high level. These choices are reported honestly rather than hidden, and they do not overturn the qualitative result — even setting aside the assembly-depth and stringency effects, the three pipelines remain non-overlapping in Fig. 1 and rank samples inconsistently (Table 4). Fourth, the cohort is small (n = 19) and from a single workstation run; per-habitat cells are 4–5 samples, so the environment and interaction estimates are less precise than the pipeline main effect. Scaling to public metagenome compendia (ENA/MGnify) and to full-depth, fully-QC'd, equal-stringency processing is the natural next step, and would let the pipeline-attributable variance be estimated without the workstation constraints.

## 5. Conclusion

**The bioinformatics pipeline is the dominant source of variation in environmental ARG estimates — larger than the environment itself.** Across 19 real metagenomes processed through three pipelines, pipeline explained 57.0% of abundance variance versus 18.9% for environment, with a significant habitat-dependent interaction (11.0%); pipelines differed by 7- to 152-fold, ranked samples inconsistently (ρ = 0.19–0.71), and disagreed on whether an environmental gradient exists at all. Detection stringency alone moved samples from zero to hundreds of ARG calls. The operational recommendation is explicit: do not compare metagenomic ARG abundances across studies without holding the pipeline and its settings fixed, or without modeling pipeline-attributable variance, and always report the caller and stringency for presence/absence claims. The full quantification, decomposition, and figure code are provided as a reproducible template.

## Acknowledgements

The authors thank the OPALS program at the Institute of Engineering in Medicine, UC San Diego, for supporting this research. C.Z., L.S., and Y.M. designed the study and developed the analysis framework. E.W., A.F., and W.T. contributed to data analysis and interpretation. Computational work used publicly available SRA data; sample accessions are listed in `metadata/batch_subset.txt`.

## References

[1] T. U. Berendonk et al., "Tackling antibiotic resistance: the environmental framework," Nature Reviews Microbiology, vol. 13, pp. 310–317, 2015. doi: 10.1038/nrmicro3439

[2] B. Pärnänen et al., "Antibiotic resistance in European wastewater treatment plants mirrors the pattern of clinical antibiotic resistance prevalence," Science Advances, vol. 5, no. 3, p. eaau9124, 2019. doi: 10.1126/sciadv.aau9124

[3] J. Bengtsson-Palme, "Strategies to improve usability and preserve accuracy in biological sequence databases," Proteomics, vol. 16, no. 18, pp. 2454–2460, 2016. doi: 10.1002/pmic.201600034

[4] A. Raza et al., "Comparison of bioinformatics pipelines and databases for antimicrobial resistance gene detection in metagenomic data," Frontiers in Microbiology, vol. 14, p. 1199825, 2023. doi: 10.3389/fmicb.2023.1199825

[5] D. Li et al., "MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph," Bioinformatics, vol. 31, no. 10, pp. 1674–1676, 2015. doi: 10.1093/bioinformatics/btv033

[6] D. Hyatt et al., "Prodigal: prokaryotic gene recognition and translation initiation site identification," BMC Bioinformatics, vol. 11, p. 119, 2010. doi: 10.1186/1471-2105-11-119

[7] B. Buchfink, K. Reuter, and H.-G. Drost, "Sensitive protein alignments at tree-of-life scale using DIAMOND," Nature Methods, vol. 18, pp. 366–368, 2021. doi: 10.1038/s41592-021-01101-x

[8] B. P. Alcock et al., "CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database," Nucleic Acids Research, vol. 51, no. D1, pp. D690–D699, 2023. doi: 10.1093/nar/gkac920

[9] B. Langmead and S. L. Salzberg, "Fast gapped-read alignment with Bowtie 2," Nature Methods, vol. 9, pp. 357–359, 2012. doi: 10.1038/nmeth.1923

[10] M. Feldgarden et al., "AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence," Scientific Reports, vol. 11, p. 12728, 2021. doi: 10.1038/s41598-021-91456-0

[11] P. J. McMurdie and S. Holmes, "Waste not, want not: why rarefying microbiome data is inadmissible," PLoS Computational Biology, vol. 10, no. 4, p. e1003531, 2014. doi: 10.1371/journal.pcbi.1003531

[12] S. Seabold and J. Perktold, "Statsmodels: econometric and statistical modeling with Python," in Proc. 9th Python in Science Conf. (SciPy 2010), 2010. doi: 10.25080/Majora-92bf1922-011

[13] The pandas development team, "pandas-dev/pandas: Pandas," Zenodo, 2020. doi: 10.5281/zenodo.3509134
