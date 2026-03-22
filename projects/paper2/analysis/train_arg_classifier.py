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

    fieldnames = ["method", "precision", "recall", "f1", "roc_auc"]
    rows = []
    if os.path.exists(args.metrics_out):
        with open(args.metrics_out, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                method = (row.get("method") or "").strip()
                if not method:
                    continue
                rows.append(
                    {
                        "method": method,
                        "precision": (row.get("precision") or "").strip(),
                        "recall": (row.get("recall") or "").strip(),
                        "f1": (row.get("f1") or "").strip(),
                        "roc_auc": (row.get("roc_auc") or "").strip(),
                    }
                )

    rf_row = {
        "method": "RandomForest_embeddings",
        "precision": f"{precision:.6f}",
        "recall": f"{recall:.6f}",
        "f1": f"{f1:.6f}",
        "roc_auc": "NA",
    }

    updated = False
    for idx, row in enumerate(rows):
        if row["method"] == rf_row["method"]:
            rows[idx] = rf_row
            updated = True
            break
    if not updated:
        rows.append(rf_row)

    with open(args.metrics_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Labeled proteins used: {len(indices)}")
    print(f"Wrote metrics: {args.metrics_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
