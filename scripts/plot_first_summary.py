#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
from collections import Counter, defaultdict

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate first summary plots from ARG dataset")
    parser.add_argument(
        "--dataset",
        default="results/ARG_dataset.csv",
        help="Path to merged ARG dataset (default: results/ARG_dataset.csv)",
    )
    parser.add_argument(
        "--outdir",
        default="analysis/figures",
        help="Directory for output figures (default: analysis/figures)",
    )
    return parser.parse_args()


def parse_float(value: str) -> float:
    if value is None:
        return 0.0
    text = value.strip()
    if not text:
        return 0.0
    try:
        return float(text)
    except ValueError:
        return 0.0


def safe_figure_with_message(path: str, title: str, message: str) -> None:
    plt.figure(figsize=(8, 5))
    plt.title(title)
    plt.text(0.5, 0.5, message, ha="center", va="center")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def main() -> int:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.dataset, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [row for row in reader]

    fig1 = os.path.join(args.outdir, "figure1_dataset_composition.png")
    fig2 = os.path.join(args.outdir, "figure2_arg_abundance_by_environment.png")
    fig3 = os.path.join(args.outdir, "figure3_arg_richness_histogram.png")

    if not rows:
        safe_figure_with_message(fig1, "Dataset Composition", "No rows in dataset")
        safe_figure_with_message(fig2, "ARG Abundance by Environment", "No rows in dataset")
        safe_figure_with_message(fig3, "ARG Richness Distribution", "No rows in dataset")
        print("WARNING: dataset is empty. Wrote placeholder figures.")
        return 0

    env_counts = Counter((row.get("environment") or "NA").strip() or "NA" for row in rows)

    plt.figure(figsize=(8, 5))
    envs = list(env_counts.keys())
    counts = [env_counts[e] for e in envs]
    plt.bar(envs, counts)
    plt.title("Figure 1: Dataset Composition")
    plt.ylabel("Number of Samples")
    plt.xlabel("Environment")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(fig1, dpi=200)
    plt.close()

    by_env = defaultdict(list)
    for row in rows:
        env = (row.get("environment") or "NA").strip() or "NA"
        by_env[env].append(parse_float(row.get("ARG_total", "0")))

    plt.figure(figsize=(8, 5))
    labels = list(by_env.keys())
    values = [by_env[l] for l in labels]
    plt.boxplot(values, tick_labels=labels)
    plt.title("Figure 2: ARG Abundance by Environment")
    plt.ylabel("ARG_total")
    plt.xlabel("Environment")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(fig2, dpi=200)
    plt.close()

    richness = [parse_float(row.get("ARG_richness", "0")) for row in rows]
    plt.figure(figsize=(8, 5))
    plt.hist(richness, bins=20)
    plt.title("Figure 3: ARG Richness Distribution")
    plt.xlabel("ARG_richness")
    plt.ylabel("Number of Samples")
    plt.tight_layout()
    plt.savefig(fig3, dpi=200)
    plt.close()

    print(f"Wrote: {fig1}")
    print(f"Wrote: {fig2}")
    print(f"Wrote: {fig3}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
