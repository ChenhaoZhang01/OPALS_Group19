#!/usr/bin/env python3
"""Generate the four Paper 1 figures into analysis/figures/.

Fig 1  ARG abundance across pipelines (box + jittered points colored by environment).
Fig 2  Variance partition of log10(ARG_total) and of ARG richness (stacked bars).
Fig 3  Pipeline x environment interaction plot (cell means with 95% CI).
Fig 4  ARG richness comparison across pipelines (box by pipeline).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PIPE_ORDER = ["pipelineA", "pipelineB", "pipelineC"]
PIPE_LABEL = {
    "pipelineA": "A: assembly+DIAMOND",
    "pipelineB": "B: read-map CARD",
    "pipelineC": "C: RGI strict",
}
ENV_ORDER = ["wastewater", "river", "irrigation", "soil"]
ENV_COLOR = {
    "wastewater": "#c44e52",
    "river": "#4c72b0",
    "irrigation": "#55a868",
    "soil": "#8172b3",
}
PIPE_COLOR = {"pipelineA": "#4c72b0", "pipelineB": "#c44e52", "pipelineC": "#55a868"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--results", default="projects/paper1/results")
    p.add_argument("--figdir", default="projects/paper1/analysis/figures")
    return p.parse_args()


def fig1_abundance(df: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    data = [df.loc[df.pipeline == p, "log_ARG_total"].values for p in PIPE_ORDER]
    bp = ax.boxplot(data, positions=range(len(PIPE_ORDER)), widths=0.55,
                    patch_artist=True, showfliers=False, zorder=1)
    for patch, p in zip(bp["boxes"], PIPE_ORDER):
        patch.set_facecolor("#e8e8e8"); patch.set_edgecolor("#444")
    for med in bp["medians"]:
        med.set_color("#222"); med.set_linewidth(1.6)
    rng = np.random.default_rng(7)
    for i, p in enumerate(PIPE_ORDER):
        sub = df[df.pipeline == p]
        x = i + rng.uniform(-0.18, 0.18, len(sub))
        ax.scatter(x, sub["log_ARG_total"], s=26, alpha=0.8,
                   c=[ENV_COLOR[e] for e in sub["environment"]], edgecolors="none", zorder=2)
    ax.set_xticks(range(len(PIPE_ORDER)))
    ax.set_xticklabels([PIPE_LABEL[p] for p in PIPE_ORDER])
    ax.set_ylabel("log$_{10}$ normalized ARG abundance\n(ARG hits / read)")
    ax.set_title("ARG abundance estimate depends on pipeline")
    handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=ENV_COLOR[e], mec="none", label=e)
               for e in ENV_ORDER]
    ax.legend(handles=handles, title="Environment", frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(figdir / "fig1_abundance_by_pipeline.png", dpi=200)
    plt.close(fig)


def fig2_variance(results: Path, figdir: Path) -> None:
    va = pd.read_csv(results / "variance_decomposition.csv")
    vr = pd.read_csv(results / "richness_decomposition.csv")
    order = ["environment", "pipeline", "pipeline:environment", "residual"]
    colors = {"environment": "#55a868", "pipeline": "#c44e52",
              "pipeline:environment": "#dd8452", "residual": "#b0b0b0"}
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    labels = ["log$_{10}$ ARG\nabundance", "ARG\nrichness"]
    bottoms = np.zeros(2)
    for fac in order:
        vals = [
            float(va.loc[va.factor == fac, "variance_percent"].iloc[0]),
            float(vr.loc[vr.factor == fac, "variance_percent"].iloc[0]),
        ]
        ax.bar(labels, vals, bottom=bottoms, color=colors[fac], label=fac, edgecolor="white", width=0.6)
        for i, v in enumerate(vals):
            if v > 4:
                ax.text(i, bottoms[i] + v / 2, f"{v:.0f}%", ha="center", va="center",
                        color="white", fontsize=9, fontweight="bold")
        bottoms += vals
    ax.set_ylabel("Variance explained (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Variance partition: environment vs. pipeline")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2)
    fig.tight_layout()
    fig.savefig(figdir / "fig2_variance_partition.png", dpi=200)
    plt.close(fig)


def fig3_interaction(results: Path, figdir: Path) -> None:
    cell = pd.read_csv(results / "pipeline_env_means.csv")
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    xpos = {e: i for i, e in enumerate(ENV_ORDER)}
    for p in PIPE_ORDER:
        sub = cell[cell.pipeline == p].set_index("environment").reindex(ENV_ORDER)
        ci = 1.96 * sub["std"] / np.sqrt(sub["size"])
        ax.errorbar([xpos[e] for e in ENV_ORDER], sub["mean"], yerr=ci,
                    marker="o", capsize=3, lw=2, color=PIPE_COLOR[p], label=PIPE_LABEL[p])
    ax.set_xticks(range(len(ENV_ORDER)))
    ax.set_xticklabels(ENV_ORDER)
    ax.set_xlabel("Environment")
    ax.set_ylabel("Mean log$_{10}$ ARG abundance")
    ax.set_title("Pipeline x environment interaction\n(non-parallel lines = pipeline bias is environment-dependent)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(figdir / "fig3_pipeline_env_interaction.png", dpi=200)
    plt.close(fig)


def fig4_richness(df: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    data = [df.loc[df.pipeline == p, "ARG_richness"].values for p in PIPE_ORDER]
    bp = ax.boxplot(data, positions=range(len(PIPE_ORDER)), widths=0.55,
                    patch_artist=True, showfliers=False)
    for patch, p in zip(bp["boxes"], PIPE_ORDER):
        patch.set_facecolor(PIPE_COLOR[p]); patch.set_alpha(0.55); patch.set_edgecolor("#333")
    for med in bp["medians"]:
        med.set_color("#111"); med.set_linewidth(1.6)
    rng = np.random.default_rng(11)
    for i, p in enumerate(PIPE_ORDER):
        sub = df[df.pipeline == p]
        x = i + rng.uniform(-0.16, 0.16, len(sub))
        ax.scatter(x, sub["ARG_richness"], s=18, alpha=0.5, color="#333", edgecolors="none")
    ax.set_xticks(range(len(PIPE_ORDER)))
    ax.set_xticklabels([PIPE_LABEL[p] for p in PIPE_ORDER])
    ax.set_ylabel("ARG classes detected per sample")
    ax.set_title("Detected ARG richness differs by pipeline")
    fig.tight_layout()
    fig.savefig(figdir / "fig4_richness_by_pipeline.png", dpi=200)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    results = Path(args.results)
    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(results / "pipeline_long_table.csv")
    df["ARG_total"] = pd.to_numeric(df["ARG_total"], errors="coerce")
    df["ARG_richness"] = pd.to_numeric(df["ARG_richness"], errors="coerce")
    df = df.dropna(subset=["ARG_total"]).copy()
    df["log_ARG_total"] = np.log10(df["ARG_total"])

    fig1_abundance(df, figdir)
    fig2_variance(results, figdir)
    fig3_interaction(results, figdir)
    fig4_richness(df, figdir)
    print(f"Wrote 4 figures to {figdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
