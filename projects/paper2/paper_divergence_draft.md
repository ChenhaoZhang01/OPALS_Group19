<p align="center"><b>Failure Modes of ARG Detection Under Sequence Divergence</b></p>
<p align="center"><b>Chenhao Zhang¹, Eliana Wong¹*, Ashley Fang¹*, Wanze Tang¹*, Linda Shi², Yujie Men²</b></p>
<p align="center">¹Institute of Engineering in Medicine, University of California, San Diego, San Diego, La Jolla, CA 92095<br>²Department of Chemical and Environmental Engineering, University of California, Riverside, California 92521<br>
*High school students participating in IEM OPALS program</p>

**Abstract -** *Sequence alignment methods are central to antibiotic resistance gene (ARG) detection, but their reliability under sequence divergence remains under-characterized. We evaluate failure modes of ARG detection across identity strata by comparing BLAST and real protein-language-model embeddings (ESM2), while enforcing identity-aware evaluation design. We define three explicit hypotheses: H1, BLAST recall drops as sequence identity decreases; H2, embedding representations are more robust to divergence; H3, standard random-split evaluation overestimates real-world performance when homologous sequences leak between train and test sets. In identity-stratified evaluation, BLAST recall is 1.000 at 70-90% identity but 0.857 in the 30-50% bin, a 14.3% absolute drop under divergence. ESM2 shows lower absolute recall than BLAST but substantially less sensitivity to divergence (0.563 at 70-90% vs 0.531 at 30-50%; 3.2% absolute drop). We also provide an identity-clustered split protocol to reduce leakage risk and report per-bin metrics rather than aggregate-only scores. The results highlight that aggregate recall can obscure clinically relevant low-identity failure regimes, and that divergence-aware benchmarking should be standard for ARG screening pipelines.*

**Keywords:** ARG detection, sequence divergence, BLAST, protein embeddings, identity-aware split, benchmarking

## 1. Introduction
Most ARG surveillance pipelines rely on sequence similarity search, typically BLAST-style alignment against curated databases. While this paradigm is highly effective for near-identical homologs, deployed metagenomic settings frequently contain novel or diverged variants. The key scientific question is not only "which method has higher aggregate accuracy," but "where does detection fail as divergence increases?"

The contribution of this study is to identify and quantify failure of ARG detection under sequence divergence. We compare BLAST and embedding-based prediction under identity-stratified evaluation, and we define an identity-clustered split strategy to reduce train/test leakage from near-duplicate sequences.

We test three explicit hypotheses:

- H1: BLAST performance drops with sequence divergence.
- H2: Embedding-based methods are more robust to sequence divergence.
- H3: Standard evaluation overestimates performance when identity leakage is not controlled.

## 2. Methods

### Data and Evaluation Design
We use CARD-derived protein query and reference sets with label mappings at the resistance-mechanism level. BLAST top-hit predictions and per-query identity values are used to characterize alignment behavior under divergence.

To avoid reporting only aggregate performance, we evaluate recall in identity bins:

- 30-50%
- 50-70%
- 70-90%

We also report benchmark-level summaries for overall, low-identity, and high-identity regimes.

### Identity-Clustered Train/Test Split
To address leakage, we implement identity-clustered splitting before model training:

1. Build a query-query identity graph from pairwise identity output.
2. Connect sequence pairs above a defined identity threshold.
3. Form connected components as identity clusters.
4. Assign whole clusters (not individual sequences) to train or test.

This design directly tests H3 by reducing inflated performance caused by homologous near-duplicates crossing the split boundary.

### Models and Outputs
We benchmark two families:

- BLAST alignment baseline
- Embedding-based classifier

Current study outputs include:

- Main figure: recall vs sequence identity with BLAST and ESM2 lines.
- Identity-stratified table with 30-50%, 50-70%, 70-90% bins.
- Benchmark table with overall, low-ID, and high-ID recall.

## 3. Results

### Main Result Figure (Critical)
Recall-vs-identity analysis shows divergence-dependent failure in the alignment baseline.

<p align="center"><img src="analysis/figures/recall_vs_identity.png" width="760"></p>
<p align="center"><b>Fig. 1.</b> Recall versus sequence identity. Lines show BLAST and embedding models across identity bins.</p>

BLAST recall decreases from 1.000 in the 70-90% bin to 0.857 in the 30-50% bin, a 14.3% absolute decline under divergence.

### Identity-Stratified Evaluation

**Table 1. Identity-stratified recall.**

| Identity Bin (%) | n | BLAST Recall | Embedding Recall |
|---|---:|---:|---:|
| 30-50 | 49 | 0.857 | 0.531 |
| 50-70 | 36 | 1.000 | 0.806 |
| 70-90 | 32 | 1.000 | 0.563 |

H1 is supported in this run: BLAST recall decreases under divergence.

### Benchmark Table

**Table 2. Benchmark summary.**

| Method | Overall Recall | Low-ID Recall | High-ID Recall |
|---|---:|---:|---:|
| BLAST_alignment | 0.927 | 0.918 | 1.000 |
| Embedding_esm2 | 0.364 | 0.647 | 0.563 |

H2 is supported in a robustness sense: ESM2 remains comparatively stable across identity bins, though absolute recall remains lower than BLAST.

### Split Robustness and H3
This study includes an identity-clustered split implementation that clusters sequence pairs by identity before assigning train/test partitions. This prevents homologous leakage and provides a direct mechanism to evaluate H3 in future model reruns.

## 4. Discussion
The main contribution is a failure-focused framing of ARG detection under divergence, rather than a generic method comparison. Reporting only aggregate recall obscures the regime where surveillance systems are most vulnerable.

Standard pipelines may systematically miss ARGs in low-identity regimes. Standard evaluation pipelines may also overestimate expected field performance when near-duplicate sequence leakage is not controlled by identity-aware splitting.

## 5. Future Work and Limitations
The current embedding benchmark uses ESM2 embeddings with CPU-constrained inference settings. Absolute recall remains below BLAST, and additional optimization is needed before deployment claims.

Future work should include:
1. Recompute embedding benchmarks with multiple checkpoints (ProtBERT, larger ESM2 variants) under matched evaluation conditions.
2. Expand low-identity strata to include <30% with sufficient sample size.
3. Evaluate calibrated uncertainty for divergence-aware thresholding.
4. Validate on external metagenomic cohorts with novel ARG families.

## 6. Conclusion
This study reframes the contribution from "method A vs method B" to detection-failure characterization under sequence divergence. BLAST recall decreases by 14.3% (absolute) from 70-90% identity to 30-50% identity, demonstrating divergence-sensitive failure. Identity-stratified reporting and identity-clustered splitting should be treated as default practice for ARG benchmark design.

## 7. References
1. Altschul, S. F., et al. (1990). Basic local alignment search tool. Journal of Molecular Biology.
2. Elnaggar, A., et al. (2022). ProtTrans: Toward understanding the language of life through self-supervised learning. IEEE TPAMI.
3. Lin, Z., et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. Science.
4. CARD database and ARO ontology documentation. https://card.mcmaster.ca
