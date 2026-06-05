#!/usr/bin/env python3
"""Figures for Paper 4 (corrected): ARG & class-1 integron removal across quaternary lines.

Reads the tables produced by run_removal_analysis.py and emits:
  fig1_removal_cascade.png      ARG concentration along IN -> OUT -> CW/GAC/AOP (treatments work)
  fig2_removal_by_line_gene.png absolute vs 16S-normalized log removal, line x gene
  fig3_biomass_vs_selective.png absolute vs normalized ARG removal per line (biomass gap)
  fig4_integron_retention.png   ARG vs intI1 (mobile) normalized removal per line (the "why")
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LINES = ["CW", "GAC", "AOP"]
LINE_LABEL = {"CW": "Constructed\nwetland", "GAC": "GAC", "AOP": "Ozonation\n+ GAC"}
ARGS = ["blaAmpC", "ermB", "sul1", "sul2", "tetW"]
SAMPLE_ORDER = ["IN WWTP", "OUT WWTP", "CW", "GAC", "AOP"]
SAMPLE_LABEL = {"IN WWTP": "Influent", "OUT WWTP": "Conventional\neffluent",
                "CW": "Constructed\nwetland", "GAC": "GAC", "AOP": "Ozonation\n+ GAC"}
COL = {"CW": "#55a868", "GAC": "#dd8452", "AOP": "#4c72b0"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--results", default="projects/paper4/results")
    p.add_argument("--figdir", default="projects/paper4/analysis/figures")
    return p.parse_args()


def fig1_cascade(clean: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    # geometric mean concentration per sample x gene (across campaigns)
    g = (clean[clean.Gene.isin(ARGS)]
         .groupby(["Sample", "Gene"])["conc_gm"]
         .apply(lambda s: np.exp(np.log(s).mean())).reset_index())
    xs = [s for s in SAMPLE_ORDER if s in g.Sample.unique()]
    pos = {s: i for i, s in enumerate(xs)}
    # per-gene points
    rng = np.random.default_rng(3)
    for s in xs:
        sub = g[g.Sample == s]
        x = pos[s] + rng.uniform(-0.12, 0.12, len(sub))
        ax.scatter(x, np.log10(sub["conc_gm"]), s=30, color="#888", alpha=0.7, zorder=2)
    # mean-of-ARGs bar marker
    means = [np.log10(np.exp(np.log(g[g.Sample == s]["conc_gm"]).mean())) for s in xs]
    ax.plot(range(len(xs)), means, "-D", color="#c44e52", lw=2, ms=8,
            label="Geometric mean of ARGs", zorder=3)
    ax.set_xticks(range(len(xs))); ax.set_xticklabels([SAMPLE_LABEL[s] for s in xs])
    ax.set_ylabel(r"log$_{10}$ ARG concentration (copies/L)")
    ax.axvspan(-0.5, 1.5, color="#f0f0f0", zorder=0)
    ax.text(0.5, ax.get_ylim()[1], "conventional plant", ha="center", va="top", fontsize=8, color="#666")
    ax.text(3, ax.get_ylim()[1], "quaternary lines", ha="center", va="top", fontsize=8, color="#666")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout(); fig.savefig(figdir / "fig1_removal_cascade.png", dpi=200); plt.close(fig)


def _grouped(ax, rem, value_col, title):
    genes = ARGS + ["intI1"]
    x = np.arange(len(genes)); w = 0.26
    for i, line in enumerate(LINES):
        sub = rem[rem.line == line].set_index("gene").reindex(genes)
        ax.bar(x + (i - 1) * w, sub[value_col], w, label=LINE_LABEL[line].replace("\n", " "),
               color=COL[line], alpha=0.9, edgecolor="white")
    ax.axhline(0, color="#333", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(genes, rotation=30, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10)


def fig2_by_line_gene(rem: pd.DataFrame, figdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    _grouped(axes[0], rem, "abs_log_removal_mean", "Absolute (total copies/L)")
    _grouped(axes[1], rem, "norm_log_removal_mean", "16S-normalized")
    axes[0].set_ylabel(r"log$_{10}$ removal vs. conventional effluent")
    axes[1].legend(frameon=False, fontsize=8, title="Quaternary line")
    fig.tight_layout(); fig.savefig(figdir / "fig2_removal_by_line_gene.png", dpi=200); plt.close(fig)


def _ci_halfwidth(unc: pd.DataFrame, line: str, metric: str) -> float:
    """Half-width of the 95% t-CI for (line, metric); 0 if unavailable."""
    row = unc[(unc.line == line) & (unc.metric == metric)]
    if row.empty or pd.isna(row.iloc[0]["t_ci_lo"]):
        return 0.0
    return float((row.iloc[0]["t_ci_hi"] - row.iloc[0]["t_ci_lo"]) / 2.0)


def fig3_biomass(mech: pd.DataFrame, unc: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    x = np.arange(len(LINES)); w = 0.36
    m = mech.set_index("line").reindex(LINES)
    abs_err = [_ci_halfwidth(unc, l, "ARG_abs") for l in LINES]
    norm_err = [_ci_halfwidth(unc, l, "ARG_norm") for l in LINES]
    ax.bar(x - w / 2, m["arg_abs_removal"], w, yerr=abs_err, capsize=3,
           error_kw={"ecolor": "#444", "lw": 1}, label="Absolute (total copies)",
           color="#b0b0b0", edgecolor="white")
    ax.bar(x + w / 2, m["arg_norm_removal"], w, yerr=norm_err, capsize=3,
           error_kw={"ecolor": "#444", "lw": 1}, label="Selective (16S-normalized)",
           color="#c44e52", edgecolor="white")
    for i, line in enumerate(LINES):
        gap = m.loc[line, "biomass_vs_selective_gap"]
        ax.annotate(f"biomass\ngap {gap:.1f}", (i, max(m.loc[line, "arg_abs_removal"], 0) + 0.05),
                    ha="center", va="bottom", fontsize=7.5, color="#555")
    ax.set_xticks(x); ax.set_xticklabels([LINE_LABEL[l] for l in LINES])
    ax.set_ylabel(r"Mean ARG log$_{10}$ removal")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(figdir / "fig3_biomass_vs_selective.png", dpi=200); plt.close(fig)


def fig4_integron(mech: pd.DataFrame, unc: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    x = np.arange(len(LINES)); w = 0.36
    m = mech.set_index("line").reindex(LINES)
    arg_err = [_ci_halfwidth(unc, l, "ARG_norm") for l in LINES]
    int_err = [_ci_halfwidth(unc, l, "intI1_norm") for l in LINES]
    ax.bar(x - w / 2, m["arg_norm_removal"], w, yerr=arg_err, capsize=3,
           error_kw={"ecolor": "#444", "lw": 1}, label="ARGs (mean)",
           color="#4c72b0", edgecolor="white")
    ax.bar(x + w / 2, m["intI1_norm_removal"], w, yerr=int_err, capsize=3,
           error_kw={"ecolor": "#444", "lw": 1}, label="intI1 (class-1 integron)",
           color="#c44e52", edgecolor="white")
    ax.axhline(0, color="#333", lw=0.8)
    for i, line in enumerate(LINES):
        ax.annotate(f"{m.loc[line,'intI1_norm_removal']:.2f}",
                    (i + w / 2, m.loc[line, "intI1_norm_removal"]),
                    textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels([LINE_LABEL[l] for l in LINES])
    ax.set_ylabel(r"16S-normalized log$_{10}$ removal")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(figdir / "fig4_integron_retention.png", dpi=200); plt.close(fig)


def fig5_multisite(results: Path, figdir: Path) -> None:
    """Cross-site selective (community-normalized) ARG removal by conventional treatment."""
    path = results / "multisite_removal.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    classes = ["sulfonamide", "tetracycline", "MLS", "beta-lactam"]
    sites = ["Bliesen", "Taiwan", "Slovakia"]
    site_col = {"Bliesen": "#c44e52", "Taiwan": "#4c72b0", "Slovakia": "#dd8452"}
    site_lab = {"Bliesen": "Germany (Bliesen)", "Taiwan": "Taiwan", "Slovakia": "Slovakia"}
    piv = df.pivot_table(index="arg_class", columns="site",
                         values="selective_log_removal", aggfunc="mean").reindex(classes)
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    x = np.arange(len(classes)); w = 0.26
    # fix axis limits first so the shaded "enriched" band fills cleanly to the floor
    vmax = float(np.nanmax(piv.values)); vmin = float(np.nanmin(piv.values))
    ymax = vmax + 0.22
    ymin = min(vmin, 0.0) - 0.12
    ax.set_ylim(ymin, ymax)
    ax.axhspan(ymin, 0, color="#f6e0e0", alpha=0.6, zorder=0)
    for i, s in enumerate(sites):
        bars = ax.bar(x + (i - 1) * w, piv[s].values, w, label=site_lab[s],
                      color=site_col[s], alpha=0.9, edgecolor="white", zorder=2)
        for rect, v in zip(bars, piv[s].values):
            if np.isnan(v):
                continue
            ax.annotate(f"{v:+.2f}", (rect.get_x() + rect.get_width() / 2,
                        v + (0.02 if v >= 0 else -0.02)),
                        ha="center", va="bottom" if v >= 0 else "top", fontsize=6.5, color="#333")
    ax.axhline(0, color="#333", lw=1, zorder=3)
    ax.text(0.015, 0.04, "enriched / no removal", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=8, color="#a33", style="italic")
    ax.set_xticks(x); ax.set_xticklabels(classes)
    ax.set_ylabel(r"Selective log$_{10}$ removal (community-normalized)")
    ax.set_xlim(-0.5, len(classes) - 0.5)
    ax.legend(frameon=False, fontsize=8, title="Site", loc="upper right")
    fig.tight_layout(); fig.savefig(figdir / "fig5_multisite_conventional.png", dpi=200); plt.close(fig)


def main() -> int:
    args = parse_args()
    results = Path(args.results); figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    # remove stale DiD figures
    for old in ["fig1_before_after_box", "fig2_trends_over_time", "fig3_did_estimate", "fig4_change_by_type"]:
        p = figdir / f"{old}.png"
        if p.exists():
            p.unlink()
    clean = pd.read_csv(results / "s10_clean.csv")
    rem = pd.read_csv(results / "removal_by_line_gene.csv")
    mech = pd.read_csv(results / "integron_mechanism.csv")
    unc_path = results / "uncertainty_by_line.csv"
    unc = pd.read_csv(unc_path) if unc_path.exists() else pd.DataFrame(columns=["line", "metric", "t_ci_lo", "t_ci_hi"])
    fig1_cascade(clean, figdir)
    fig2_by_line_gene(rem, figdir)
    fig3_biomass(mech, unc, figdir)
    fig4_integron(mech, unc, figdir)
    fig5_multisite(results, figdir)
    print(f"Wrote figures to {figdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
