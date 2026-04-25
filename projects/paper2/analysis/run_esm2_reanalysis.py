#!/usr/bin/env python3
"""Full ESM2 reanalysis for Paper 2.

Steps:
  1. Generate ESM2 embeddings for all query sequences (tmp_q3.faa, 1827 seqs)
  2. Train RandomForest on labeled subset (card_query_homolog.faa, 564 seqs)
  3. Predict labels for all 1827 query sequences
  4. Match predictions with existing BLAST per-query data (pident, blast_correct)
  5. Recompute identity-binned recall
  6. Update blast_vs_ml_metrics.csv and identity_bin_recall_every10_moredata.csv
  7. Regenerate figures
"""

from __future__ import annotations

import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch

MAX_SEQ_LEN = 1022
BATCH_SIZE = 8
MIN_N_RELIABLE = 50
BIN_EDGES = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

REPO = Path(__file__).resolve().parents[3]
PAPER2 = REPO / "projects/paper2"
RESULTS = PAPER2 / "results"
PROTEINS = PAPER2 / "proteins"
FIGURES = PAPER2 / "analysis/figures"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> list[tuple[str, str]]:
    seqs: list[tuple[str, str]] = []
    header: str | None = None
    buf: list[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(buf)))
                header = line[1:]
                buf = []
            elif line:
                buf.append(line)
    if header is not None:
        seqs.append((header, "".join(buf)))
    return seqs


def generate_esm2_embeddings(
    sequences: list[tuple[str, str]], model_name: str = "esm2_t6_8M_UR50D"
) -> np.ndarray:
    import esm  # noqa: PLC0415

    loader = getattr(esm.pretrained, model_name)
    model, alphabet = loader()
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    num_layers = model.num_layers

    truncated = sum(1 for _, s in sequences if len(s) > MAX_SEQ_LEN)
    if truncated:
        print(f"  Truncating {truncated} sequences to {MAX_SEQ_LEN} aa")

    all_embeddings: list[np.ndarray] = []
    for start in range(0, len(sequences), BATCH_SIZE):
        batch = [(h, s[:MAX_SEQ_LEN]) for h, s in sequences[start : start + BATCH_SIZE]]
        _, _, tokens = batch_converter(batch)
        with torch.no_grad():
            results = model(tokens, repr_layers=[num_layers], return_contacts=False)
        reps = results["representations"][num_layers]
        for i, (_, seq) in enumerate(batch):
            seq_len = min(len(seq), tokens.shape[1] - 2)
            emb = reps[i, 1 : seq_len + 1].mean(dim=0).float().numpy()
            all_embeddings.append(emb)
        done = min(start + BATCH_SIZE, len(sequences))
        print(f"  {done}/{len(sequences)}", end="\r", flush=True)
    print()
    return np.stack(all_embeddings, axis=0).astype(np.float32)


def weighted_metrics(
    y_true: list[str], y_pred: list[str]
) -> tuple[float, float, float]:
    from collections import Counter  # noqa: PLC0415
    classes = sorted(set(y_true) | set(y_pred))
    support = Counter(y_true)
    n = len(y_true)
    wp = wr = wf = 0.0
    for c in classes:
        tp = sum(1 for t, p in zip(y_true, y_pred) if t == c and p == c)
        fp = sum(1 for t, p in zip(y_true, y_pred) if t != c and p == c)
        fn = sum(1 for t, p in zip(y_true, y_pred) if t == c and p != c)
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) else 0.0
        w = support[c]
        wp += prec * w
        wr += rec * w
        wf += f1 * w
    return wp / n, wr / n, wf / n


def bin_label(pident: float) -> str | None:
    for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
        if lo <= pident < hi:
            return f"{lo}_to_{hi}"
    return None


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main() -> int:
    from sklearn.ensemble import RandomForestClassifier  # noqa: PLC0415
    from sklearn.metrics import f1_score, precision_score, recall_score  # noqa: PLC0415
    from sklearn.model_selection import train_test_split  # noqa: PLC0415

    # ------------------------------------------------------------------
    # Step 1: ESM2 embeddings for all query sequences (tmp_q3.faa)
    # ------------------------------------------------------------------
    q3_fasta = PROTEINS / "tmp_q3.faa"
    q3_emb_path = RESULTS / "embeddings" / "protein_embeddings_q3.npy"

    q3_sequences = read_fasta(q3_fasta)
    print(f"Step 1: Generating ESM2 embeddings for {len(q3_sequences)} sequences in tmp_q3.faa ...")
    q3_embeddings = generate_esm2_embeddings(q3_sequences)
    q3_emb_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(q3_emb_path, q3_embeddings)
    print(f"  Saved: {q3_emb_path}  shape={q3_embeddings.shape}")

    # header → row index in q3
    q3_header_to_idx: dict[str, int] = {h: i for i, (h, _) in enumerate(q3_sequences)}

    # ------------------------------------------------------------------
    # Step 2: Train RandomForest on labeled 564-sequence set
    # ------------------------------------------------------------------
    print("\nStep 2: Training RandomForest on labeled set (card_query_homolog.faa) ...")
    train_emb_path = RESULTS / "embeddings" / "protein_embeddings.npy"
    train_embeddings = np.load(train_emb_path)  # 564 x 320, ESM2

    train_labels_path = RESULTS / "training_labels_template.csv"
    train_labels: list[str] = []
    train_indices: list[int] = []
    with train_labels_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            idx = int(row["row_index"].strip())
            lbl = row["label"].strip()
            if 0 <= idx < len(train_embeddings) and lbl:
                train_indices.append(idx)
                train_labels.append(lbl)

    X = train_embeddings[train_indices]
    y = train_labels

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    clf = RandomForestClassifier(n_estimators=300, random_state=42)
    clf.fit(X_train, y_train)
    y_pred_test = clf.predict(X_test)

    precision = precision_score(y_test, y_pred_test, average="weighted", zero_division=0)
    recall = recall_score(y_test, y_pred_test, average="weighted", zero_division=0)
    f1 = f1_score(y_test, y_pred_test, average="weighted", zero_division=0)
    print(f"  ESM2 RF — precision={precision:.4f}  recall={recall:.4f}  f1={f1:.4f}")

    # Update blast_vs_ml_metrics.csv (keep BLAST row, replace RF row)
    metrics_path = RESULTS / "blast_vs_ml_metrics.csv"
    fieldnames = ["method", "precision", "recall", "f1", "roc_auc"]
    existing_rows: list[dict[str, str]] = []
    if metrics_path.exists():
        with metrics_path.open("r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                if (row.get("method") or "").strip() != "RandomForest_embeddings":
                    existing_rows.append(dict(row))
    existing_rows.append({
        "method": "RandomForest_ESM2",
        "precision": f"{precision:.6f}",
        "recall": f"{recall:.6f}",
        "f1": f"{f1:.6f}",
        "roc_auc": "NA",
    })
    with metrics_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(existing_rows)
    print(f"  Updated: {metrics_path}")

    # ------------------------------------------------------------------
    # Step 3: Predict labels for all 1,827 q3 sequences
    # ------------------------------------------------------------------
    print("\nStep 3: Predicting ESM2 labels for all query sequences ...")
    q3_preds = clf.predict(q3_embeddings)  # 1827 predictions
    q3_pred_map: dict[str, str] = {
        q3_sequences[i][0]: q3_preds[i] for i in range(len(q3_sequences))
    }

    # ------------------------------------------------------------------
    # Step 4: Rebuild per-query file with ESM2 predictions
    # ------------------------------------------------------------------
    print("\nStep 4: Rebuilding per-query ESM2 analysis ...")
    base_perquery_path = RESULTS / "embedding_per_query_every10.csv"
    with base_perquery_path.open("r", encoding="utf-8", newline="") as f:
        base_rows = list(csv.DictReader(f))

    new_perquery_rows: list[dict[str, str]] = []
    missing = 0
    for row in base_rows:
        qid = row["query_id"]
        if qid not in q3_pred_map:
            missing += 1
            continue
        esm2_pred = q3_pred_map[qid]
        true_lbl = row["true_label"]
        new_perquery_rows.append({
            "query_id": qid,
            "pident": row["pident"],
            "true_label": true_lbl,
            "pred_label_embedding": esm2_pred,
            "embedding_correct": str(int(esm2_pred == true_lbl)),
            "blast_correct": row["blast_correct"],
        })
    print(f"  Rows matched: {len(new_perquery_rows)}  missing ESM2 embedding: {missing}")

    perquery_out = RESULTS / "embedding_per_query_every10_moredata.csv"
    fieldnames_pq = ["query_id", "pident", "true_label", "pred_label_embedding",
                     "embedding_correct", "blast_correct"]
    with perquery_out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames_pq)
        writer.writeheader()
        writer.writerows(new_perquery_rows)
    print(f"  Saved: {perquery_out}")

    # ------------------------------------------------------------------
    # Step 5: Compute identity-binned recall
    # ------------------------------------------------------------------
    print("\nStep 5: Computing identity-binned recall ...")
    bin_counts: dict[str, dict[str, int]] = defaultdict(lambda: {"n": 0, "blast": 0, "esm2": 0})
    for row in new_perquery_rows:
        pident = float(row["pident"])
        b = bin_label(pident)
        if b is None:
            continue
        bin_counts[b]["n"] += 1
        bin_counts[b]["blast"] += int(row["blast_correct"])
        bin_counts[b]["esm2"] += int(row["embedding_correct"])

    bin_recall_path = RESULTS / "identity_bin_recall_every10_moredata.csv"
    with bin_recall_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["identity_bin", "n", "blast_recall", "embedding_recall"])
        for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
            key = f"{lo}_to_{hi}"
            d = bin_counts.get(key, {"n": 0, "blast": 0, "esm2": 0})
            n = d["n"]
            blast_r = f"{d['blast'] / n:.6f}" if n > 0 else "0.000000"
            esm2_r = f"{d['esm2'] / n:.6f}" if n > 0 else "0.000000"
            writer.writerow([key, n, blast_r, esm2_r])
            if n > 0:
                print(f"  {key}: n={n}  BLAST={float(blast_r):.3f}  ESM2={float(esm2_r):.3f}")
    print(f"  Saved: {bin_recall_path}")

    # ------------------------------------------------------------------
    # Step 6: Update benchmark table
    # ------------------------------------------------------------------
    print("\nStep 6: Updating benchmark table ...")
    bm_path = RESULTS / "benchmark_table_every10_moredata.csv"
    all_blast = sum(d["blast"] for d in bin_counts.values())
    all_esm2 = sum(d["esm2"] for d in bin_counts.values())
    all_n = sum(d["n"] for d in bin_counts.values())
    low_bins = [f"{lo}_to_{hi}" for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]) if hi <= 60]
    high_bins = [f"{lo}_to_{hi}" for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]) if lo >= 60]
    low_blast = sum(bin_counts[b]["blast"] for b in low_bins if b in bin_counts)
    low_esm2 = sum(bin_counts[b]["esm2"] for b in low_bins if b in bin_counts)
    low_n = sum(bin_counts[b]["n"] for b in low_bins if b in bin_counts)
    high_blast = sum(bin_counts[b]["blast"] for b in high_bins if b in bin_counts)
    high_esm2 = sum(bin_counts[b]["esm2"] for b in high_bins if b in bin_counts)
    high_n = sum(bin_counts[b]["n"] for b in high_bins if b in bin_counts)

    with bm_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["method", "overall_recall", "low_id_recall", "high_id_recall"])
        writer.writerow(["BLAST_alignment",
                         f"{all_blast/all_n:.6f}" if all_n else "NA",
                         f"{low_blast/low_n:.6f}" if low_n else "NA",
                         f"{high_blast/high_n:.6f}" if high_n else "NA"])
        writer.writerow(["ESM2_RandomForest",
                         f"{all_esm2/all_n:.6f}" if all_n else "NA",
                         f"{low_esm2/low_n:.6f}" if low_n else "NA",
                         f"{high_esm2/high_n:.6f}" if high_n else "NA"])
    print(f"  Saved: {bm_path}")

    # ------------------------------------------------------------------
    # Step 7: Regenerate figures
    # ------------------------------------------------------------------
    print("\nStep 7: Regenerating figures ...")
    fig_script = Path(__file__).parent / "generate_final_figures.py"
    subprocess.run([sys.executable, str(fig_script)], check=True)

    print("\nDone. ESM2 reanalysis complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
