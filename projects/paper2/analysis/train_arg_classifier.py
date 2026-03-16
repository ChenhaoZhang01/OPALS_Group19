#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train a starter ARG classifier for Paper 2")
    parser.add_argument("--embeddings", required=True, help="Path to protein_embeddings.npy")
    parser.add_argument("--labels", required=True, help="CSV with columns row_index,label")
    parser.add_argument("--metrics-out", required=True, help="CSV output for model metrics")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        import numpy as np
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.metrics import f1_score, precision_score, recall_score
        from sklearn.model_selection import train_test_split
    except ImportError as exc:
        raise SystemExit(
            "This script needs numpy and scikit-learn. Install with: pip install numpy scikit-learn"
        ) from exc

    embeddings = np.load(args.embeddings)

    labels = []
    indices = []
    with open(args.labels, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            idx_text = (row.get("row_index") or "").strip()
            label = (row.get("label") or "").strip()
            if not idx_text or not label:
                continue
            idx = int(idx_text)
            if idx < 0 or idx >= len(embeddings):
                continue
            indices.append(idx)
            labels.append(label)

    if len(indices) < 10:
        raise ValueError("Need at least 10 labeled rows for a starter train/test split")

    X = embeddings[indices]
    y = labels

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    model = RandomForestClassifier(n_estimators=300, random_state=42)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    precision = precision_score(y_test, y_pred, average="weighted", zero_division=0)
    recall = recall_score(y_test, y_pred, average="weighted", zero_division=0)
    f1 = f1_score(y_test, y_pred, average="weighted", zero_division=0)

    os.makedirs(os.path.dirname(args.metrics_out), exist_ok=True)
    with open(args.metrics_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "precision", "recall", "f1", "roc_auc"])
        writer.writerow(["RandomForest_embeddings", f"{precision:.6f}", f"{recall:.6f}", f"{f1:.6f}", "NA"])

    print(f"Labeled proteins used: {len(indices)}")
    print(f"Wrote metrics: {args.metrics_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
