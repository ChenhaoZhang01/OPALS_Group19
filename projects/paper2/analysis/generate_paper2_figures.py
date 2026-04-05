#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_labels(labels_csv: Path) -> tuple[np.ndarray, np.ndarray]:
    indices: list[int] = []
    labels: list[str] = []
    with labels_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            idx_text = (row.get("row_index") or "").strip()
            label = (row.get("label") or "").strip()
            if not idx_text or not label:
                continue
            indices.append(int(idx_text))
            labels.append(label)
    return np.array(indices, dtype=int), np.array(labels)


def pca_2d(x: np.ndarray) -> np.ndarray:
    x_centered = x - x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(x_centered, full_matrices=False)
    return u[:, :2] * s[:2]


def save_embedding_pca_plot(embeddings: np.ndarray, labels_idx: np.ndarray, labels: np.ndarray, out_png: Path) -> None:
    x = embeddings[labels_idx]
    x2 = pca_2d(x)

    classes = sorted(set(labels.tolist()))
    palette = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]
    color_map = {cls: palette[i % len(palette)] for i, cls in enumerate(classes)}

    plt.figure(figsize=(7.5, 5.5), dpi=150)
    for cls in classes:
        mask = labels == cls
        plt.scatter(x2[mask, 0], x2[mask, 1], s=45, alpha=0.9, c=color_map[cls], label=cls)

    plt.title("Embedding PCA by ARG Class")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(frameon=False)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png)
    plt.close()


def save_method_comparison(metrics_csv: Path, out_png: Path) -> None:
    methods: list[str] = []
    precision: list[float] = []
    recall: list[float] = []
    f1: list[float] = []

    with metrics_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            p = (row.get("precision") or "").strip()
            r = (row.get("recall") or "").strip()
            f = (row.get("f1") or "").strip()
            if p in {"", "pending", "NA"} or r in {"", "pending", "NA"} or f in {"", "pending", "NA"}:
                continue
            methods.append((row.get("method") or "").strip())
            precision.append(float(p))
            recall.append(float(r))
            f1.append(float(f))

    if not methods:
        return

    x = np.arange(len(methods))
    width = 0.24

    plt.figure(figsize=(8.2, 5.5), dpi=150)
    plt.bar(x - width, precision, width=width, label="Precision", color="#1f77b4")
    plt.bar(x, recall, width=width, label="Recall", color="#ff7f0e")
    plt.bar(x + width, f1, width=width, label="F1", color="#2ca02c")

    plt.ylim(0, 1.0)
    plt.xticks(x, methods, rotation=15, ha="right")
    plt.ylabel("Score")
    plt.title("Method Comparison")
    plt.legend(frameon=False)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png)
    plt.close()


def save_identity_bin_recall(identity_bin_csv: Path, out_png: Path) -> None:
    if not identity_bin_csv.exists():
        return

    bins: list[str] = []
    blast_recall: list[float] = []
    embedding_recall: list[float] = []

    with identity_bin_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            n_text = (row.get("n") or "").strip()
            if not n_text or int(n_text) <= 0:
                continue
            b = (row.get("identity_bin") or "").strip()
            br = (row.get("blast_recall") or "").strip()
            er = (row.get("embedding_recall") or "").strip()
            if b and br not in {"", "NA"} and er not in {"", "NA"}:
                bins.append(b)
                blast_recall.append(float(br))
                embedding_recall.append(float(er))

    if not bins:
        return

    x = np.arange(len(bins))
    width = 0.34
    plt.figure(figsize=(8.4, 5.2), dpi=150)
    plt.bar(x - width / 2, blast_recall, width=width, label="BLAST recall", color="#3b6fb6")
    plt.bar(x + width / 2, embedding_recall, width=width, label="Embedding recall", color="#c95a49")
    plt.ylim(0, 1.0)
    plt.xticks(x, bins)
    plt.ylabel("Recall")
    plt.xlabel("Sequence identity bin")
    plt.title("Recall by Sequence Identity Bin")
    plt.legend(frameon=False)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png)
    plt.close()


def save_recall_gap_by_identity(identity_bin_csv: Path, out_csv: Path, out_png: Path) -> None:
    if not identity_bin_csv.exists():
        return

    bins: list[str] = []
    blast_recall: list[float] = []
    embedding_recall: list[float] = []

    with identity_bin_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            n_text = (row.get("n") or "").strip()
            if not n_text or int(n_text) <= 0:
                continue
            b = (row.get("identity_bin") or "").strip()
            br = (row.get("blast_recall") or "").strip()
            er = (row.get("embedding_recall") or "").strip()
            if b and br not in {"", "NA"} and er not in {"", "NA"}:
                bins.append(b)
                blast_recall.append(float(br))
                embedding_recall.append(float(er))

    if not bins:
        return

    gap = [b - e for b, e in zip(blast_recall, embedding_recall)]

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["identity_bin", "blast_recall", "embedding_recall", "delta_recall_blast_minus_embedding"])
        for b, br, er, dr in zip(bins, blast_recall, embedding_recall, gap):
            writer.writerow([b, f"{br:.6f}", f"{er:.6f}", f"{dr:.6f}"])

    x = np.arange(len(bins))
    plt.figure(figsize=(8.2, 5.0), dpi=150)
    plt.plot(x, gap, marker="o", linewidth=2.0, color="#2c7fb8")
    plt.axhline(0.0, color="#666666", linewidth=1.0, linestyle="--")
    plt.xticks(x, bins)
    plt.ylabel("Recall gap (BLAST - embedding)")
    plt.xlabel("Sequence identity bin")
    plt.title("BLAST-Embedding Recall Gap vs Identity")
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png)
    plt.close()


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    embeddings_path = repo / "projects/paper2/results/embeddings/protein_embeddings.npy"
    labels_path = repo / "projects/paper2/results/training_labels_template.csv"
    metrics_path = repo / "projects/paper2/results/blast_vs_ml_metrics.csv"
    identity_bin_path = repo / "projects/paper2/results/identity_bin_recall.csv"
    recall_gap_path = repo / "projects/paper2/results/recall_gap_by_identity.csv"
    figures_dir = repo / "projects/paper2/analysis/figures"

    embeddings = np.load(embeddings_path)
    label_idx, labels = load_labels(labels_path)

    save_embedding_pca_plot(
        embeddings=embeddings,
        labels_idx=label_idx,
        labels=labels,
        out_png=figures_dir / "embedding_pca.png",
    )
    save_method_comparison(
        metrics_csv=metrics_path,
        out_png=figures_dir / "method_comparison_bar.png",
    )
    save_identity_bin_recall(
        identity_bin_csv=identity_bin_path,
        out_png=figures_dir / "identity_bin_recall.png",
    )
    save_recall_gap_by_identity(
        identity_bin_csv=identity_bin_path,
        out_csv=recall_gap_path,
        out_png=figures_dir / "recall_gap_vs_identity.png",
    )

    print("Wrote figures to", figures_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
