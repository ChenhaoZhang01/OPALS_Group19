#!/usr/bin/env python3
"""Generate publication-quality figures for Paper 2."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


def load_identity_bins(
    csv_path: Path,
) -> tuple[list[str], list[int], list[float], list[float]]:
    bins: list[str] = []
    ns: list[int] = []
    blast: list[float] = []
    embed: list[float] = []
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            n = int(row["n"])
            if n == 0:
                continue
            raw = row["identity_bin"]
            label = raw.replace("_to_", "–").replace("_", "")
            bins.append(label)
            ns.append(n)
            blast.append(float(row["blast_recall"]))
            embed.append(float(row["embedding_recall"]))
    return bins, ns, blast, embed


def load_benchmark_metrics(
    csv_path: Path,
) -> tuple[list[str], list[float], list[float], list[float]]:
    methods: list[str] = []
    precision: list[float] = []
    recall: list[float] = []
    f1: list[float] = []
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            p = (row.get("precision") or "").strip()
            r = (row.get("recall") or "").strip()
            ff = (row.get("f1") or "").strip()
            if p in {"", "NA", "pending"} or r in {"", "NA", "pending"}:
                continue
            methods.append((row.get("method") or "").strip())
            precision.append(float(p))
            recall.append(float(r))
            f1.append(float(ff))
    return methods, precision, recall, f1


def make_identity_recall_figure(
    bins: list[str],
    ns: list[int],
    blast: list[float],
    embed: list[float],
    out_png: Path,
    min_n_reliable: int = 50,
) -> None:
    """Fig. 1 — recall vs. identity bin, only bins with n >= min_n_reliable shown."""
    filtered = [
        (b, n, bv, ev)
        for b, n, bv, ev in zip(bins, ns, blast, embed)
        if n >= min_n_reliable
    ]
    if not filtered:
        print(f"Warning: no bins with n >= {min_n_reliable}, skipping figure.")
        return
    bins_f, ns_f, blast_f, embed_f = zip(*filtered)

    fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=180)

    x = np.arange(len(bins_f))
    width = 0.36

    blast_color = "#2166ac"
    embed_color = "#d6604d"

    for i, (bv, ev, n) in enumerate(zip(blast_f, embed_f, ns_f)):
        ax.bar(x[i] - width / 2, bv, width, color=blast_color,
               edgecolor="white", linewidth=0.8, zorder=3)
        ax.bar(x[i] + width / 2, ev, width, color=embed_color,
               edgecolor="white", linewidth=0.8, zorder=3)

        ax.text(x[i] - width / 2, bv + 0.008, f"{bv:.3f}",
                ha="center", va="bottom", fontsize=8.5, color="#222222")
        ax.text(x[i] + width / 2, ev + 0.008, f"{ev:.3f}",
                ha="center", va="bottom", fontsize=8.5, color="#222222")

        ax.text(x[i], -0.09, f"n={n}", ha="center", va="top",
                fontsize=8, color="#555555",
                transform=ax.get_xaxis_transform())

    ax.set_ylim(0, 0.40)
    ax.set_xlim(-0.6, len(bins_f) - 0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{b}%" for b in bins_f], fontsize=10)
    ax.set_ylabel("Recall", fontsize=11)
    ax.set_xlabel("Sequence identity bin (query vs. top BLAST hit)",
                  fontsize=11, labelpad=26)
    ax.set_title(
        "ARG Detection Recall Across Sequence Identity Strata",
        fontsize=12, pad=10,
    )

    ax.yaxis.grid(True, linestyle="--", linewidth=0.6, color="#cccccc", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        mpatches.Patch(color=blast_color, label="BLAST alignment"),
        mpatches.Patch(color=embed_color, label="RF embedding classifier"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=10, loc="upper right")

    fig.tight_layout(rect=[0, 0.06, 1, 1])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")


def make_benchmark_figure(
    methods: list[str],
    precision: list[float],
    recall: list[float],
    f1: list[float],
    out_png: Path,
) -> None:
    """Fig. 2 — overall benchmark bar chart (precision / recall / F1)."""
    clean_names = []
    for m in methods:
        if "BLAST" in m:
            clean_names.append("BLAST alignment")
        elif "RandomForest" in m or "RF" in m:
            clean_names.append("RF embedding classifier")
        else:
            clean_names.append(m)

    x = np.arange(len(clean_names))
    width = 0.24

    fig, ax = plt.subplots(figsize=(7.5, 5.0), dpi=180)
    ax.bar(x - width, precision, width, label="Precision", color="#4dac26", zorder=3)
    ax.bar(x,         recall,    width, label="Recall",    color="#2166ac", zorder=3)
    ax.bar(x + width, f1,        width, label="F1",        color="#d6604d", zorder=3)

    for xi, (p, r, f) in enumerate(zip(precision, recall, f1)):
        for offset, val in [(-width, p), (0, r), (width, f)]:
            ax.text(xi + offset, val + 0.01, f"{val:.3f}", ha="center",
                    va="bottom", fontsize=8, color="#333333")

    ax.set_ylim(0, 1.12)
    ax.set_xticks(x)
    ax.set_xticklabels(clean_names, fontsize=10)
    ax.set_ylabel("Score", fontsize=11)
    ax.set_title("Overall Benchmark: Weighted Precision, Recall, and F1",
                 fontsize=12, pad=10)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.6, color="#cccccc", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=9.5)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    results = repo / "projects/paper2/results"
    figures = repo / "projects/paper2/analysis/figures"

    # Fig. 1: identity-stratified recall (bins with n >= 50 only)
    bins, ns, blast_r, embed_r = load_identity_bins(
        results / "identity_bin_recall_every10_moredata.csv"
    )
    make_identity_recall_figure(
        bins, ns, blast_r, embed_r,
        out_png=figures / "recall_vs_identity_every10_moredata.png",
    )

    # Fig. 2: overall benchmark comparison
    methods, prec, rec, f1 = load_benchmark_metrics(results / "blast_vs_ml_metrics.csv")
    if methods:
        make_benchmark_figure(
            methods, prec, rec, f1,
            out_png=figures / "method_comparison_bar.png",
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
