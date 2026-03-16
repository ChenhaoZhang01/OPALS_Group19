# Paper 2: Foundation Model ARG Detection

For interns with limited coding experience: start with `START_HERE.md` in this folder.

## Objective

Test whether protein language model embeddings detect ARGs better than sequence alignment.

## Inputs

- proteins/*.faa (environmental proteins)
- CARD protein sequences

## Outputs

- results/embeddings/protein_embeddings.npy
- results/training/training_table.csv
- results/models/
- results/predictions/predicted_ARG_classes.csv
- results/benchmark/blast_vs_ml_metrics.csv
- analysis/figures/

## Workflow

1. Generate protein embeddings.

Install package:

```bash
pip install fair-esm
```

Python setup example:

```python
import esm

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
```

Save one 1280-dimensional embedding per protein in `results/embeddings/protein_embeddings.npy`.

2. Build labeled training table from CARD.

| protein | embedding | label |
|---|---|---|

3. Train classifiers.

- Random Forest
- SVM
- Neural Network

4. Predict ARG classes for environmental proteins.

5. Run alignment benchmark (BLAST or DIAMOND) on same proteins.

6. Compare precision, recall, and ROC across methods.

## Figures

1. Embedding clustering plot (UMAP).
2. Classifier ROC curve.
3. BLAST versus ML comparison chart.
4. Novel predicted ARG cluster overview.

## Expected signal

Embedding models may recover remote ARG homologs that alignment-only methods miss.

