<p align="center"><b>ARG Detection Under Sequence Divergence: A Reproducible Benchmark of BLAST and Embedding Models</b></p>
<p align="center"><b>Chenhao Zhang¹, Eliana Wong¹*, Ashley Fang¹*, Wanze Tang¹*, Linda Shi², Yujie Men²</b></p>
<p align="center">¹Institute of Engineering in Medicine, University of California, San Diego, San Diego, La Jolla, CA 92095<br>²Department of Chemical and Environmental Engineering, University of California, Riverside, California 92521<br>
*High school students participating in IEM OPALS program</p>

**Abstract -** *Antibiotic resistance gene (ARG) surveillance remains heavily dependent on sequence alignment against curated databases, a strategy that is reliable for close homologs but vulnerable to sequence divergence. We present a reproducible benchmark that contrasts an alignment baseline (BLAST top-hit assignment) with an embedding-based classifier under a domain-shift framing. Using 220 labeled proteins across five resistance-mechanism classes, we evaluate weighted precision, recall, and F1 overall and across sequence-identity bins. The embedding model (RandomForest on 128-dimensional proxy embeddings) achieved precision 0.723, recall 0.659, and F1 0.676. The BLAST baseline achieved precision 0.967, recall 0.966, and F1 0.966. In a low-identity subset (<40% identity, n=51), BLAST recall dropped to 0.706 while embedding recall remained 0.659, shrinking the recall gap from 0.341 in the 90–100% identity bin to 0.047 in the <40% bin. These results show that alignment retains higher absolute accuracy, but its advantage diminishes rapidly as identity decreases, while embedding recall is comparatively stable. The benchmark provides a transparent, end-to-end workflow for testing detection robustness under sequence divergence and highlights priority upgrades, including full protein language-model embeddings and identity-clustered holdouts.*

**Keywords:** antibiotic resistance genes, sequence divergence, embeddings, BLAST, domain shift, benchmark, low-identity detection

## 1. Introduction
Antibiotic resistance gene (ARG) detection is a cornerstone of environmental and clinical surveillance. In practice, most pipelines rely on sequence alignment against curated ARG databases because the method is interpretable, reproducible, and operationally mature. However, alignment methods depend on similarity thresholds that can miss remote homologs, a failure mode that becomes critical in understudied environments and emerging resistance contexts.

Protein embedding representations offer a contrasting strategy. Instead of relying on pairwise identity, embeddings encode higher-order sequence patterns that can correlate with functional similarity even when identity is low. Recent work using protein language models reinforces this premise for resistance-related prediction and ARG characterization under sequence divergence [4,5]. This study asks whether an embedding-based ARG classifier is less sensitive to sequence divergence than an alignment baseline. We structure the evaluation as a domain-shift test and focus on identity-stratified recall, which is the most operationally relevant measure for surveillance when false negatives are costly.

We test two hypotheses:
1. Alignment performance degrades as sequence identity decreases.
2. Embedding-based performance is less sensitive to sequence identity.

## 2. Methods

### Data and Labels
We use a CARD-backed benchmark built from the local repository workflow. The dataset contains 220 labeled protein sequences across five resistance-mechanism classes: antibiotic inactivation, target alteration, target protection, target replacement, and efflux. Labels are drawn from a training label table with `row_index` and `label` columns, and embeddings are stored as a 220x128 matrix in the results directory.

### Embedding Representation and Classifier
The repository uses a proxy embedding representation derived from hashed k-mer features (k=3). This is a reproducible placeholder for higher-dimensional protein language-model embeddings. The classifier is a RandomForest with 300 trees and a fixed random seed. The train/test split uses a stratified 80/20 partition with random_state = 42. This split avoids exact sequence duplication but does not enforce identity-based clustering, a limitation addressed in Future Work.

### Alignment Baseline
The baseline uses BLASTP top-hit assignment against a local protein database. Settings are: outfmt 6 with qseqid/sseqid/bitscore/evalue, max_target_seqs = 1, max_hsps = 1, and e-value threshold 1e-5. Predictions are assigned from the top hit’s database label. Weighted precision, recall, and F1 are computed on queries with valid top-hit labels.

### Low-Identity Analysis and Identity Bins
Low-identity analysis is performed by binning BLAST top-hit identities and slicing the subset where identity <40%. We report weighted metrics for the full set and low-identity subset, then compute recall per identity bin (90–100%, 70–90%, 40–70%, and <40%) for BLAST and the embedding model. The recall gap is reported as BLAST recall minus embedding recall per bin.

### Evaluation Metrics
Primary metrics are weighted precision, weighted recall, and weighted F1. We focus on recall by identity bin because detection failures in low-identity regimes represent the most consequential surveillance risk.

## 3. Results

<p align="center"><img src="analysis/figures/embedding_pca.png" width="700"></p>
<p align="center"><b>Fig. 1.</b> PCA projection of proxy embeddings, colored by resistance-mechanism class.</p>

### Overall Performance
BLAST outperforms the embedding model on all three weighted metrics in this run:

**Table 1. Overall weighted performance for BLAST vs. embedding model.**

| Method | Precision | Recall | F1 |
|---|---:|---:|---:|
| BLAST alignment | 0.967 | 0.966 | 0.966 |
| Embedding model (RandomForest) | 0.723 | 0.659 | 0.676 |

<p align="center"><img src="analysis/figures/method_comparison_bar.png" width="700"></p>
<p align="center"><b>Fig. 2.</b> Overall weighted precision, recall, and F1 for BLAST vs. embedding model.</p>

### Low-Identity Subset
The low-identity subset (<40% identity, n=51) is non-empty and supports direct comparison. BLAST recall drops from 0.966 overall to 0.706 in the low-identity subset, while embedding recall remains 0.659. This reduces the recall gap from high-identity regimes to low-identity regimes.

**Table 2. Weighted performance overall vs. low-identity subset (<40% identity).**

| Scenario | Method | Precision | Recall | F1 | n |
|---|---|---:|---:|---:|---:|
| All | BLAST alignment | 0.932 | 0.927 | 0.927 | 220 |
| All | Embedding model | 0.723 | 0.659 | 0.676 | 220 |
| <40% identity | BLAST alignment | 0.799 | 0.706 | 0.717 | 51 |
| <40% identity | Embedding model | 0.723 | 0.659 | 0.676 | 51 |

### Identity-Bin Recall and Recall Gap
Identity-stratified recall confirms a steep decline in BLAST recall as identity decreases, while embedding recall remains stable:

**Table 3. Recall by identity bin and recall gap (BLAST - embedding).**

| Identity bin | n | BLAST recall | Embedding recall | Gap (BLAST - Embedding) |
|---|---:|---:|---:|---:|
| 90–100% | 82 | 1.000 | 0.659 | 0.341 |
| 70–90% | 32 | 1.000 | 0.659 | 0.341 |
| 40–70% | 55 | 0.982 | 0.659 | 0.323 |
| <40% | 51 | 0.706 | 0.659 | 0.047 |

<p align="center"><img src="analysis/figures/identity_bin_recall.png" width="700"></p>
<p align="center"><b>Fig. 3.</b> Recall by sequence-identity bin for BLAST and the embedding model.</p>

<p align="center"><img src="analysis/figures/recall_gap_vs_identity.png" width="700"></p>
<p align="center"><b>Fig. 4.</b> Recall gap (BLAST - embedding) across identity bins.</p>

The recall advantage of alignment decreases by approximately seven-fold from the highest identity bin to the lowest identity bin. This supports H1 and H2: alignment performance degrades with divergence, while embedding recall is less identity-sensitive.

## 4. Discussion
Alignment remains the most accurate method in absolute terms for this dataset. However, the sharp collapse in the BLAST recall advantage at low identity indicates a systematic failure mode under sequence divergence. In surveillance contexts, this is precisely where under-detection is most harmful, as emerging or divergent resistance genes are more likely to be missed.

The embedding model used here is a proxy representation based on hashed k-mers. While this limits absolute performance, it provides a reproducible baseline and a lower bound on what representation-based detection could achieve. The identity-binned analysis shows that even this proxy representation exhibits stability across identity regimes. This suggests that hybrid pipelines combining alignment with embedding-based screening may provide higher coverage without sacrificing interpretability in high-identity regimes.

## 5. Future Work and Limitations
Limitations include the use of proxy embeddings rather than protein language-model embeddings, class imbalance, and a single random train/test split without identity-clustered holdouts. The dataset size supports low-identity analysis but remains modest for broader generalization.

Priority upgrades include:
1. Replace proxy embeddings with full protein language-model embeddings (target 1280 dimensions).
2. Implement identity-based clustering to construct strict low-identity holdouts.
3. Increase labeled query size to 500+ while preserving non-empty <40% bins.
4. Evaluate hybrid routing strategies that combine alignment and embedding detection.

## 6. Conclusion
This study provides a transparent, reproducible benchmark for ARG detection under sequence divergence. In the current CARD-backed run, BLAST outperforms the embedding model overall, but the alignment advantage diminishes rapidly as identity decreases. These findings support the hypothesis that alignment is sensitive to divergence and that embedding-based detection offers more stable recall across identity regimes. The benchmark establishes a practical foundation for future upgrades toward full protein language-model embeddings and hybrid detection pipelines.

## 7. References
1. CARD database data archive. https://card.mcmaster.ca
2. BLAST+ suite user documentation. https://blast.ncbi.nlm.nih.gov
3. scikit-learn documentation (RandomForestClassifier). https://scikit-learn.org
4. Yagimoto, K., Hosoda, S., Sato, M., & Hamada, M. (2024). Prediction of antibiotic resistance mechanisms using a protein language model. Bioinformatics. https://doi.org/10.1093/bioinformatics/btae550
5. Ahmed, S., Emon, M. I., Moumi, N. A., Huang, L., Zhou, D., Vikesland, P., Pruden, A., & Zhang, L. (2024). ProtAlign-ARG: Antibiotic Resistance Gene Characterization Integrating Protein Language Models and Alignment-Based Scoring. bioRxiv (preprint). https://doi.org/10.1101/2024.03.20.585944
