#!/usr/bin/env python3
"""Uncertainty quantification for the German advanced-treatment site (Paper 4).

Per quaternary line, recompute campaign-level removal (n = 4 campaigns) and report
mean, SD, t-based 95% CI, and a campaign bootstrap 95% CI for:
  * ARG_norm   — mean 16S-normalized log10 removal across the 5 ARGs
  * ARG_abs    — mean absolute log10 removal across the 5 ARGs
  * intI1_norm — 16S-normalized log10 removal of the class-1 integron

Also paired (within-campaign) line contrasts (CW vs GAC, CW vs AOP, GAC vs AOP) for
ARG_norm and intI1_norm. With only n = 4 campaigns we report CIs, not p-values.

Outputs: uncertainty_by_line.csv, pairwise_line_diff.csv, uncertainty_summary.md
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RESULTS = Path("projects/paper4/results")
LINES = ["CW", "GAC", "AOP"]
ARGS = ["blaAmpC", "ermB", "sul1", "sul2", "tetW"]
NORM = "16S rRNA"
RNG = np.random.default_rng(42)
NBOOT = 10000


def per_campaign_removals() -> dict:
    clean = pd.read_csv(RESULTS / "s10_clean.csv")
    wide = clean.pivot_table(index=["Sample", "Campaign"], columns="Gene", values="conc_gm")
    norm = wide.div(wide[NORM], axis=0)
    campaigns = sorted({c for (_s, c) in wide.index})

    def removal(table, line, gene, c):
        ko, kx = ("OUT WWTP", c), (line, c)
        if ko in table.index and kx in table.index:
            a, b = table.loc[ko, gene], table.loc[kx, gene]
            if pd.notna(a) and pd.notna(b) and a > 0 and b > 0:
                return np.log10(a / b)
        return np.nan

    out = {}  # (line, metric) -> array over campaigns
    for line in LINES:
        arg_norm, arg_abs, int_norm = [], [], []
        for c in campaigns:
            an = np.nanmean([removal(norm, line, g, c) for g in ARGS])
            aa = np.nanmean([removal(wide, line, g, c) for g in ARGS])
            iN = removal(norm, line, "intI1", c)
            arg_norm.append(an); arg_abs.append(aa); int_norm.append(iN)
        out[(line, "ARG_norm")] = np.array(arg_norm, float)
        out[(line, "ARG_abs")] = np.array(arg_abs, float)
        out[(line, "intI1_norm")] = np.array(int_norm, float)
    return out, campaigns


def ci_row(vals: np.ndarray) -> dict:
    vals = vals[~np.isnan(vals)]
    n = len(vals)
    mean = float(np.mean(vals)); sd = float(np.std(vals, ddof=1)) if n > 1 else np.nan
    if n > 1:
        tcrit = stats.t.ppf(0.975, n - 1)
        half = tcrit * sd / np.sqrt(n)
        t_lo, t_hi = mean - half, mean + half
        boots = [np.mean(RNG.choice(vals, n, replace=True)) for _ in range(NBOOT)]
        b_lo, b_hi = np.percentile(boots, [2.5, 97.5])
        # one-sample t-test vs 0; for a difference array this is the paired test.
        p_vs_zero = float(stats.ttest_1samp(vals, 0.0).pvalue)
    else:
        t_lo = t_hi = b_lo = b_hi = p_vs_zero = np.nan
    return {"mean": mean, "sd": sd, "n": n,
            "t_ci_lo": t_lo, "t_ci_hi": t_hi, "boot_ci_lo": b_lo, "boot_ci_hi": b_hi,
            "p_vs_zero": p_vs_zero}


def main() -> int:
    data, campaigns = per_campaign_removals()

    rows = []
    for line in LINES:
        for metric in ["ARG_abs", "ARG_norm", "intI1_norm"]:
            r = ci_row(data[(line, metric)]); r.update({"line": line, "metric": metric})
            rows.append(r)
    unc = pd.DataFrame(rows)[["line", "metric", "mean", "sd", "n",
                              "t_ci_lo", "t_ci_hi", "boot_ci_lo", "boot_ci_hi", "p_vs_zero"]]
    unc.to_csv(RESULTS / "uncertainty_by_line.csv", index=False)

    # paired within-campaign line contrasts
    pair_rows = []
    for metric in ["ARG_norm", "intI1_norm"]:
        for a, b in [("CW", "GAC"), ("CW", "AOP"), ("GAC", "AOP")]:
            diff = data[(a, metric)] - data[(b, metric)]
            r = ci_row(diff); r.update({"contrast": f"{a}-{b}", "metric": metric})
            pair_rows.append(r)
    pair = pd.DataFrame(pair_rows)[["contrast", "metric", "mean", "sd", "n",
                                    "t_ci_lo", "t_ci_hi", "boot_ci_lo", "boot_ci_hi", "p_vs_zero"]]
    pair = pair.rename(columns={"p_vs_zero": "p_paired"})
    pair.to_csv(RESULTS / "pairwise_line_diff.csv", index=False)

    L = ["# German-site uncertainty (n = 4 campaigns; CIs, not p-values)\n",
         f"Campaigns: {campaigns}\n",
         "## Per-line removal: mean [95% t-CI] (bootstrap CI)\n",
         "| Line | Metric | Mean | SD | 95% t-CI | Bootstrap 95% CI |",
         "|---|---|---:|---:|---|---|"]
    for _, r in unc.iterrows():
        L.append(f"| {r['line']} | {r['metric']} | {r['mean']:.2f} | "
                 f"{r['sd']:.2f} | [{r['t_ci_lo']:.2f}, {r['t_ci_hi']:.2f}] | "
                 f"[{r['boot_ci_lo']:.2f}, {r['boot_ci_hi']:.2f}] |")
    L.append("\n## Paired within-campaign line contrasts (mean difference [95% t-CI])\n")
    L.append("| Contrast | Metric | Mean diff | 95% t-CI |")
    L.append("|---|---|---:|---|")
    for _, r in pair.iterrows():
        L.append(f"| {r['contrast']} | {r['metric']} | {r['mean']:.2f} | "
                 f"[{r['t_ci_lo']:.2f}, {r['t_ci_hi']:.2f}] |")
    L.append("\nNote: with n = 4 campaigns these intervals are wide and indicative; "
             "they bound the line rankings rather than establish significance.")
    (RESULTS / "uncertainty_summary.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nWrote uncertainty_by_line.csv, pairwise_line_diff.csv, uncertainty_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
