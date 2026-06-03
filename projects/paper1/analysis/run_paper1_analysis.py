#!/usr/bin/env python3
"""Full statistical battery for Paper 1: Pipeline Uncertainty Decomposition.

Inputs
------
projects/paper1/results/pipeline_long_table.csv
    columns: sample_id, environment, pipeline, ARG_total, ARG_richness

Outputs (projects/paper1/results/)
----------------------------------
variance_decomposition.csv        Type-II ANOVA variance partition of log10(ARG_total).
richness_decomposition.csv        Same partition for ARG_richness.
pipeline_summary.csv              Per-pipeline mean/median ARG_total and richness.
pipeline_env_means.csv            Cell means (pipeline x environment) for the interaction plot.
pipeline_concordance.csv          Pairwise Spearman concordance between pipelines (per sample).
analysis_summary.md               Human-readable summary of all results.

Variance partitioning is performed on log10(ARG_total) because normalized ARG
abundance spans orders of magnitude and the pipeline/environment effects are
multiplicative; on the log scale they are additive, which is the standard scale for
metagenomic abundance variance decomposition.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import spearmanr


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Paper 1 full analysis battery.")
    p.add_argument("--input", default="projects/paper1/results/pipeline_long_table.csv")
    p.add_argument("--outdir", default="projects/paper1/results")
    return p.parse_args()


def variance_partition(df: pd.DataFrame, response: str, outpath: Path) -> pd.DataFrame:
    model = smf.ols(
        f"{response} ~ C(pipeline) + C(environment) + C(pipeline):C(environment)",
        data=df,
    ).fit()
    anova = sm.stats.anova_lm(model, typ=2)
    ss_total = anova["sum_sq"].sum()
    anova = anova.copy()
    anova["variance_percent"] = anova["sum_sq"] / ss_total * 100.0
    rename = {
        "C(pipeline)": "pipeline",
        "C(environment)": "environment",
        "C(pipeline):C(environment)": "pipeline:environment",
        "Residual": "residual",
    }
    out = anova.reset_index().rename(columns={"index": "factor"})
    out["factor"] = out["factor"].map(lambda x: rename.get(x, x))
    out = out[["factor", "sum_sq", "df", "F", "PR(>F)", "variance_percent"]]
    out.to_csv(outpath, index=False)
    return out


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    df["ARG_total"] = pd.to_numeric(df["ARG_total"], errors="coerce")
    df["ARG_richness"] = pd.to_numeric(df["ARG_richness"], errors="coerce")
    df = df.dropna(subset=["ARG_total"]).copy()
    df["log_ARG_total"] = np.log10(df["ARG_total"])

    # 1. Variance partition of log abundance and of richness.
    va = variance_partition(df, "log_ARG_total", outdir / "variance_decomposition.csv")
    vr = variance_partition(df, "ARG_richness", outdir / "richness_decomposition.csv")

    # 2. Per-pipeline summary.
    summ = (
        df.groupby("pipeline")
        .agg(
            n=("ARG_total", "size"),
            mean_ARG_total=("ARG_total", "mean"),
            median_ARG_total=("ARG_total", "median"),
            mean_richness=("ARG_richness", "mean"),
        )
        .reset_index()
    )
    summ.to_csv(outdir / "pipeline_summary.csv", index=False)

    # 3. Pipeline x environment cell means (for interaction figure).
    cell = (
        df.groupby(["pipeline", "environment"])["log_ARG_total"]
        .agg(["mean", "std", "size"])
        .reset_index()
    )
    cell.to_csv(outdir / "pipeline_env_means.csv", index=False)

    # 4. Pairwise pipeline concordance: do pipelines rank samples consistently?
    wide = df.pivot_table(index="sample_id", columns="pipeline", values="ARG_total")
    pipes = ["pipelineA", "pipelineB", "pipelineC"]
    conc_rows = []
    for i in range(len(pipes)):
        for j in range(i + 1, len(pipes)):
            a, b = pipes[i], pipes[j]
            mask = wide[a].notna() & wide[b].notna()
            rho, pval = spearmanr(wide.loc[mask, a], wide.loc[mask, b])
            # fold-change ratio of medians
            ratio = float(wide.loc[mask, b].median() / wide.loc[mask, a].median())
            conc_rows.append({"pipeline_x": a, "pipeline_y": b,
                              "spearman_rho": rho, "p_value": pval,
                              "median_ratio_y_over_x": ratio, "n": int(mask.sum())})
    conc = pd.DataFrame(conc_rows)
    conc.to_csv(outdir / "pipeline_concordance.csv", index=False)

    # 5. Human-readable summary.
    def pct(factor: str, table: pd.DataFrame) -> float:
        return float(table.loc[table["factor"] == factor, "variance_percent"].iloc[0])

    lines = []
    lines.append("# Paper 1 Analysis Summary\n")
    lines.append("## Variance decomposition of log10(ARG_total)\n")
    lines.append("| Factor | SS | df | F | p | Variance % |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for _, r in va.iterrows():
        f = "" if pd.isna(r["F"]) else f"{r['F']:.2f}"
        p = "" if pd.isna(r["PR(>F)"]) else f"{r['PR(>F)']:.2e}"
        lines.append(f"| {r['factor']} | {r['sum_sq']:.4g} | {r['df']:.0f} | {f} | {p} | {r['variance_percent']:.1f} |")
    lines.append("")
    lines.append(f"- Environment explains **{pct('environment', va):.1f}%** of variance in log ARG abundance.")
    lines.append(f"- Pipeline explains **{pct('pipeline', va):.1f}%**.")
    lines.append(f"- Pipeline x environment interaction explains **{pct('pipeline:environment', va):.1f}%**.")
    lines.append(f"- Residual (sample-level) **{pct('residual', va):.1f}%**.\n")

    lines.append("## Richness decomposition (number of ARG classes detected)\n")
    lines.append(f"- Pipeline explains **{pct('pipeline', vr):.1f}%** of richness variance "
                 f"(environment **{pct('environment', vr):.1f}%**).\n")

    lines.append("## Per-pipeline summary\n")
    lines.append("| Pipeline | n | Mean ARG_total | Median ARG_total | Mean richness |")
    lines.append("|---|---:|---:|---:|---:|")
    for _, r in summ.iterrows():
        lines.append(f"| {r['pipeline']} | {r['n']:.0f} | {r['mean_ARG_total']:.3e} | "
                     f"{r['median_ARG_total']:.3e} | {r['mean_richness']:.1f} |")
    lines.append("")

    lines.append("## Pairwise pipeline concordance\n")
    lines.append("| Pair | Spearman rho | p | Median fold-difference |")
    lines.append("|---|---:|---:|---:|")
    for _, r in conc.iterrows():
        lines.append(f"| {r['pipeline_x']} vs {r['pipeline_y']} | {r['spearman_rho']:.3f} | "
                     f"{r['p_value']:.2e} | {r['median_ratio_y_over_x']:.2f}x |")
    lines.append("")

    (outdir / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")

    print("Variance decomposition (log10 ARG_total):")
    print(va.to_string(index=False))
    print("\nRichness decomposition:")
    print(vr.to_string(index=False))
    print("\nPipeline concordance:")
    print(conc.to_string(index=False))
    print(f"\nWrote outputs to {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
