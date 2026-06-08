#!/usr/bin/env python3
"""Figures for Paper 4: ARG & class-1 integron removal across quaternary lines.

Reads tables from run_removal_analysis.py / run_uncertainty.py / run_multisite_analysis.py:
  fig1_removal_cascade.png      ARG concentration along IN -> OUT -> CW/GAC/AOP
  fig2_removal_by_line_gene.png absolute vs 16S-normalized log removal, line x gene
  fig3_biomass_vs_selective.png absolute vs normalized ARG removal per line (biomass gap)
  fig4_integron_retention.png   ARG vs intI1 normalized removal per line (the mechanism)
  fig5_multisite_conventional.png cross-site conventional selective removal

Significance stars (Figs 3-4) come from t-tests on n=4 campaigns:
  ns p>=0.05, * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001.
Bars vs zero use a one-sample t-test; brackets use the paired within-campaign test.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- readability: larger default fonts everywhere ----
plt.rcParams.update({
    "font.size": 13,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "legend.title_fontsize": 12,
    "axes.linewidth": 1.0,
})

LINES = ["CW", "GAC", "AOP"]
LINE_LABEL = {"CW": "Constructed\nwetland", "GAC": "GAC", "AOP": "Ozonation\n+ GAC"}
ARGS = ["blaAmpC", "ermB", "sul1", "sul2", "tetW"]
SAMPLE_ORDER = ["IN WWTP", "OUT WWTP", "CW", "GAC", "AOP"]
SAMPLE_LABEL = {"IN WWTP": "Influent", "OUT WWTP": "Conventional\neffluent",
                "CW": "Constructed\nwetland", "GAC": "GAC", "AOP": "Ozonation\n+ GAC"}
COL = {"CW": "#55a868", "GAC": "#dd8452", "AOP": "#4c72b0"}

SIG_KEY = "ns p≥0.05    * p<0.05    ** p<0.01    *** p<0.001    **** p<0.0001"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--results", default="projects/paper4/results")
    p.add_argument("--figdir", default="projects/paper4/analysis/figures")
    return p.parse_args()


def star(p) -> str:
    if pd.isna(p):
        return ""
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _ci_halfwidth(unc: pd.DataFrame, line: str, metric: str) -> float:
    row = unc[(unc.line == line) & (unc.metric == metric)]
    if row.empty or pd.isna(row.iloc[0]["t_ci_lo"]):
        return 0.0
    return float((row.iloc[0]["t_ci_hi"] - row.iloc[0]["t_ci_lo"]) / 2.0)


def _p_vs_zero(unc: pd.DataFrame, line: str, metric: str) -> float:
    row = unc[(unc.line == line) & (unc.metric == metric)]
    return float(row.iloc[0]["p_vs_zero"]) if not row.empty and "p_vs_zero" in unc.columns else np.nan


def _p_paired(pair: pd.DataFrame, contrast: str, metric: str) -> float:
    row = pair[(pair.contrast == contrast) & (pair.metric == metric)]
    return float(row.iloc[0]["p_paired"]) if not row.empty and "p_paired" in pair.columns else np.nan


def _bracket(ax, x1: float, x2: float, y: float, label: str, h: float = 0.07) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.3, color="#222")
    ax.text((x1 + x2) / 2, y + h, label, ha="center", va="bottom", fontsize=14, color="#222")


def fig1_cascade(clean: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    g = (clean[clean.Gene.isin(ARGS)]
         .groupby(["Sample", "Gene"])["conc_gm"]
         .apply(lambda s: np.exp(np.log(s).mean())).reset_index())
    xs = [s for s in SAMPLE_ORDER if s in g.Sample.unique()]
    pos = {s: i for i, s in enumerate(xs)}
    rng = np.random.default_rng(3)
    for s in xs:
        sub = g[g.Sample == s]
        x = pos[s] + rng.uniform(-0.12, 0.12, len(sub))
        ax.scatter(x, np.log10(sub["conc_gm"]), s=42, color="#888", alpha=0.75, zorder=2)
    means = [np.log10(np.exp(np.log(g[g.Sample == s]["conc_gm"]).mean())) for s in xs]
    ax.plot(range(len(xs)), means, "-D", color="#c44e52", lw=2.4, ms=10,
            label="Geometric mean of ARGs", zorder=3)
    ax.set_xticks(range(len(xs))); ax.set_xticklabels([SAMPLE_LABEL[s] for s in xs])
    ax.set_ylabel(r"log$_{10}$ ARG concentration (copies/L)")
    nq = len(xs)
    # shade the two stages distinctly
    ax.axvspan(-0.5, 1.5, color="#e6e6e6", zorder=0)          # conventional plant (grey)
    ax.axvspan(1.5, nq - 0.5, color="#e7f0f7", zorder=0)      # quaternary lines (light blue)
    # clear grouping brackets + bold labels in a header band above the data
    ymin, ymax = ax.get_ylim()
    yr = ymax - ymin
    ax.set_ylim(ymin, ymax + 0.22 * yr)
    band = ymax + 0.05 * yr

    def _group(x1, x2, label):
        ax.plot([x1, x1, x2, x2], [band - 0.02 * yr, band, band, band - 0.02 * yr],
                color="#333", lw=1.6, clip_on=False, zorder=5)
        ax.text((x1 + x2) / 2, band + 0.025 * yr, label, ha="center", va="bottom",
                fontsize=13, fontweight="bold", color="#222")

    _group(-0.35, 1.35, "Conventional plant")
    _group(1.65, nq - 0.65, "Quaternary polishing lines")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout(); fig.savefig(figdir / "fig1_removal_cascade.png", dpi=200); plt.close(fig)


def _grouped(ax, rem, value_col, title):
    genes = ARGS + ["intI1"]
    x = np.arange(len(genes)); w = 0.26
    for i, line in enumerate(LINES):
        sub = rem[rem.line == line].set_index("gene").reindex(genes)
        ax.bar(x + (i - 1) * w, sub[value_col], w, label=LINE_LABEL[line].replace("\n", " "),
               color=COL[line], alpha=0.9, edgecolor="white")
    ax.axhline(0, color="#333", lw=0.9)
    ax.set_xticks(x); ax.set_xticklabels(genes, rotation=25, ha="right", fontsize=12)
    ax.set_title(title, fontsize=14)


def fig2_by_line_gene(rem: pd.DataFrame, figdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.4), sharey=True)
    _grouped(axes[0], rem, "abs_log_removal_mean", "Absolute (total copies/L)")
    _grouped(axes[1], rem, "norm_log_removal_mean", "16S-normalized")
    axes[0].set_ylabel(r"log$_{10}$ removal vs. conventional effluent")
    # anchor the legend outside the axes so it can never overlap the bars
    axes[1].legend(frameon=False, title="Quaternary line", loc="upper left",
                   bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(figdir / "fig2_removal_by_line_gene.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig3_biomass(mech: pd.DataFrame, unc: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    x = np.arange(len(LINES)); w = 0.36
    m = mech.set_index("line").reindex(LINES)
    abs_err = [_ci_halfwidth(unc, l, "ARG_abs") for l in LINES]
    norm_err = [_ci_halfwidth(unc, l, "ARG_norm") for l in LINES]
    ax.bar(x - w / 2, m["arg_abs_removal"], w, yerr=abs_err, capsize=4,
           error_kw={"ecolor": "#444", "lw": 1.2}, label="Absolute (total copies)",
           color="#b0b0b0", edgecolor="white")
    ax.bar(x + w / 2, m["arg_norm_removal"], w, yerr=norm_err, capsize=4,
           error_kw={"ecolor": "#444", "lw": 1.2}, label="Selective (16S-normalized)",
           color="#c44e52", edgecolor="white")
    # significance vs zero on the SELECTIVE bars (the meaningful metric); gap values are
    # given in the caption and Table 1 to keep the plot uncluttered.
    for i, line in enumerate(LINES):
        val = m.loc[line, "arg_norm_removal"]; err = norm_err[i]
        ax.text(i + w / 2, val + err + 0.1, star(_p_vs_zero(unc, line, "ARG_norm")),
                ha="center", va="bottom", fontsize=15, fontweight="bold", color="#222")
    ax.set_xticks(x); ax.set_xticklabels([LINE_LABEL[l] for l in LINES])
    ax.set_ylabel(r"Mean ARG log$_{10}$ removal")
    ax.axhline(0, color="#333", lw=0.9)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(figdir / "fig3_biomass_vs_selective.png", dpi=200); plt.close(fig)


def fig4_integron(mech: pd.DataFrame, unc: pd.DataFrame, pair: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 5.8))
    x = np.arange(len(LINES)); w = 0.36
    m = mech.set_index("line").reindex(LINES)
    arg_err = [_ci_halfwidth(unc, l, "ARG_norm") for l in LINES]
    int_err = [_ci_halfwidth(unc, l, "intI1_norm") for l in LINES]
    ax.bar(x - w / 2, m["arg_norm_removal"], w, yerr=arg_err, capsize=4,
           error_kw={"ecolor": "#444", "lw": 1.2}, label="ARGs (mean)",
           color="#4c72b0", edgecolor="white")
    ax.bar(x + w / 2, m["intI1_norm_removal"], w, yerr=int_err, capsize=4,
           error_kw={"ecolor": "#444", "lw": 1.2}, label="intI1 (class-1 integron)",
           color="#c44e52", edgecolor="white")
    ax.axhline(0, color="#333", lw=0.9)
    # significance vs zero above each bar
    for i, line in enumerate(LINES):
        ax.text(i - w / 2, m.loc[line, "arg_norm_removal"] + arg_err[i] + 0.05,
                star(_p_vs_zero(unc, line, "ARG_norm")), ha="center", va="bottom",
                fontsize=14, fontweight="bold", color="#222")
        ax.text(i + w / 2, m.loc[line, "intI1_norm_removal"] + int_err[i] + 0.05,
                star(_p_vs_zero(unc, line, "intI1_norm")), ha="center", va="bottom",
                fontsize=14, fontweight="bold", color="#222")
    # paired contrast bracket: GAC vs AOP, ARG removal (the resolved between-line difference)
    y_br = max(m.loc["GAC", "arg_norm_removal"] + arg_err[1],
               m.loc["AOP", "arg_norm_removal"] + arg_err[2]) + 0.35
    _bracket(ax, 1 - w / 2, 2 - w / 2, y_br, star(_p_paired(pair, "GAC-AOP", "ARG_norm")))
    ax.set_ylim(top=y_br + 0.55)
    ax.set_xticks(x); ax.set_xticklabels([LINE_LABEL[l] for l in LINES])
    ax.set_ylabel(r"16S-normalized log$_{10}$ removal")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(figdir / "fig4_integron_retention.png", dpi=200); plt.close(fig)


def fig5_multisite(results: Path, figdir: Path) -> None:
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
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    x = np.arange(len(classes)); w = 0.26
    vmax = float(np.nanmax(piv.values)); vmin = float(np.nanmin(piv.values))
    ymax = vmax + 0.24; ymin = min(vmin, 0.0) - 0.14
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
                        ha="center", va="bottom" if v >= 0 else "top", fontsize=8.5, color="#333")
    ax.axhline(0, color="#333", lw=1.1, zorder=3)
    ax.text(0.015, 0.04, "enriched / no removal", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=11, color="#a33", style="italic")
    ax.set_xticks(x); ax.set_xticklabels(classes)
    ax.set_ylabel(r"Selective log$_{10}$ removal (community-normalized)")
    ax.set_xlim(-0.5, len(classes) - 0.5)
    ax.legend(frameon=False, title="Site", loc="upper right")
    fig.tight_layout(); fig.savefig(figdir / "fig5_multisite_conventional.png", dpi=200); plt.close(fig)


def main() -> int:
    args = parse_args()
    results = Path(args.results); figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    for old in ["fig1_before_after_box", "fig2_trends_over_time", "fig3_did_estimate", "fig4_change_by_type"]:
        p = figdir / f"{old}.png"
        if p.exists():
            p.unlink()
    clean = pd.read_csv(results / "s10_clean.csv")
    rem = pd.read_csv(results / "removal_by_line_gene.csv")
    mech = pd.read_csv(results / "integron_mechanism.csv")
    up = results / "uncertainty_by_line.csv"
    pp = results / "pairwise_line_diff.csv"
    unc = pd.read_csv(up) if up.exists() else pd.DataFrame(columns=["line", "metric", "t_ci_lo", "t_ci_hi", "p_vs_zero"])
    pair = pd.read_csv(pp) if pp.exists() else pd.DataFrame(columns=["contrast", "metric", "p_paired"])
    fig1_cascade(clean, figdir)
    fig2_by_line_gene(rem, figdir)
    fig3_biomass(mech, unc, figdir)
    fig4_integron(mech, unc, pair, figdir)
    fig5_multisite(results, figdir)
    print(f"Wrote figures to {figdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
