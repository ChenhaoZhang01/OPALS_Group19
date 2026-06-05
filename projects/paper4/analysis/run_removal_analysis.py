#!/usr/bin/env python3
"""Paper 4 (corrected): ARG and class-1 integron REMOVAL across quaternary treatment lines.

Replaces the invalid difference-in-differences analysis. The source study (Bliesen, DE;
Supplementary Table S10 of the CEBEDEAU dataset) is a removal-efficiency comparison, not
an upgrade natural experiment: a conventional WWTP effluent (OUT WWTP) is polished by three
parallel quaternary lines — constructed wetland (CW), granular activated carbon (GAC), and
ozonation+GAC (AOP) — each sampled across 4 campaigns with technical replicates, for the
genes 16S rRNA, blaAmpC, ermB, intI1, sul1, sul2, tetW.

We compute, per gene and per quaternary line:
  * absolute log10 removal     = log10( C[OUT WWTP] / C[line] )
  * 16S-normalized log10 removal (gene/16S relative abundance) — separates genuine ARG
    removal from biomass removal.
  * the conventional plant's own removal (IN WWTP -> OUT WWTP) for context.

Mechanism (the "why"):
  * intI1 (class-1 integron-integrase) is the mobile-resistance marker; sul1 sits in the
    conserved 3' end of most class-1 integrons. We test whether ARG removal tracks intI1
    removal (mobile fraction cleared) and whether normalized removal differs from absolute
    (selective vs biomass-driven).

Outputs (projects/paper4/results/):
  s10_clean.csv                 tidy per Sample x Campaign x Gene geometric-mean concentration
  removal_by_line_gene.csv      absolute + normalized log10 removal per line x gene (mean, sd, n)
  removal_summary.csv           per-line mean removal across genes (abs + norm)
  integron_mechanism.csv        intI1 vs ARG removal, abs-vs-norm gaps
  removal_analysis_summary.md   human-readable summary
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

XLSX_DEFAULT = ("projects/paper4/data/external/oa_scan/PMC12735332/PMC12735332/"
                "supp_unzipped/Supplementary-tables-CEBEDEAU.xlsx")
QUAT_LINES = ["CW", "GAC", "AOP"]
LINE_LABEL = {"CW": "Constructed wetland", "GAC": "GAC", "AOP": "Ozonation + GAC"}
ARGS = ["blaAmpC", "ermB", "sul1", "sul2", "tetW"]   # resistance genes
MARKERS = ["intI1"]                                   # mobile-element marker
NORM = "16S rRNA"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--xlsx", default=XLSX_DEFAULT)
    p.add_argument("--sheet", default="Table S10")
    p.add_argument("--outdir", default="projects/paper4/results")
    return p.parse_args()


def load_s10(xlsx: str, sheet: str) -> pd.DataFrame:
    df = pd.read_excel(xlsx, sheet_name=sheet, header=1)
    df.columns = [str(c).strip() for c in df.columns]
    df = df.rename(columns={"Concentration (copies/L)": "conc"})
    df = df.dropna(subset=["Sample", "Gene", "conc"]).copy()
    df["Sample"] = df["Sample"].astype(str).str.strip()
    df["Gene"] = df["Gene"].astype(str).str.strip()
    df["Campaign"] = df["Campaign"].astype(str).str.strip()
    df["conc"] = pd.to_numeric(df["conc"], errors="coerce")
    df = df[df["conc"] > 0].dropna(subset=["conc"])
    return df


def geomean(s: pd.Series) -> float:
    return float(np.exp(np.log(s).mean()))


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_s10(args.xlsx, args.sheet)

    # Geometric-mean concentration per Sample x Campaign x Gene (collapses replicates).
    gm = (df.groupby(["Sample", "Campaign", "Gene"])["conc"]
            .apply(geomean).reset_index(name="conc_gm"))
    gm.to_csv(outdir / "s10_clean.csv", index=False)

    # Wide: index (Sample, Campaign) -> gene columns
    wide = gm.pivot_table(index=["Sample", "Campaign"], columns="Gene", values="conc_gm")

    # 16S-normalized relative abundance per Sample x Campaign.
    norm = wide.div(wide[NORM], axis=0)

    genes = [g for g in ARGS + MARKERS if g in wide.columns]
    campaigns = sorted(gm["Campaign"].unique())

    # ---- Removal of each quaternary line relative to its input (OUT WWTP) ----
    rows = []
    for line in QUAT_LINES:
        for gene in genes:
            abs_vals, norm_vals = [], []
            for c in campaigns:
                key_out, key_line = ("OUT WWTP", c), (line, c)
                if key_out in wide.index and key_line in wide.index:
                    co, cl = wide.loc[key_out, gene], wide.loc[key_line, gene]
                    if pd.notna(co) and pd.notna(cl) and co > 0 and cl > 0:
                        abs_vals.append(np.log10(co / cl))
                    no, nl = norm.loc[key_out, gene], norm.loc[key_line, gene]
                    if pd.notna(no) and pd.notna(nl) and no > 0 and nl > 0:
                        norm_vals.append(np.log10(no / nl))
            rows.append({
                "line": line, "gene": gene,
                "abs_log_removal_mean": np.mean(abs_vals) if abs_vals else np.nan,
                "abs_log_removal_sd": np.std(abs_vals, ddof=1) if len(abs_vals) > 1 else np.nan,
                "norm_log_removal_mean": np.mean(norm_vals) if norm_vals else np.nan,
                "norm_log_removal_sd": np.std(norm_vals, ddof=1) if len(norm_vals) > 1 else np.nan,
                "n_campaigns": len(abs_vals),
            })
    rem = pd.DataFrame(rows)
    rem.to_csv(outdir / "removal_by_line_gene.csv", index=False)

    # Conventional plant removal IN -> OUT (context)
    conv_rows = []
    for gene in genes:
        vals = []
        for c in campaigns:
            ki, ko = ("IN WWTP", c), ("OUT WWTP", c)
            if ki in wide.index and ko in wide.index:
                ci, co = wide.loc[ki, gene], wide.loc[ko, gene]
                if pd.notna(ci) and pd.notna(co) and ci > 0 and co > 0:
                    vals.append(np.log10(ci / co))
        conv_rows.append({"gene": gene, "in_out_log_removal_mean": np.mean(vals) if vals else np.nan,
                          "n": len(vals)})
    conv = pd.DataFrame(conv_rows)

    # ---- Per-line summary across ARGs (markers reported separately) ----
    arg_mask = rem["gene"].isin(ARGS)
    summary = (rem[arg_mask].groupby("line")
               .agg(mean_abs_removal=("abs_log_removal_mean", "mean"),
                    mean_norm_removal=("norm_log_removal_mean", "mean")).reset_index())
    # attach intI1 removal per line
    int_rem = rem[rem.gene == "intI1"].set_index("line")
    summary["intI1_abs_removal"] = summary["line"].map(int_rem["abs_log_removal_mean"])
    summary["intI1_norm_removal"] = summary["line"].map(int_rem["norm_log_removal_mean"])
    summary = summary.sort_values("mean_abs_removal", ascending=False)
    summary.to_csv(outdir / "removal_summary.csv", index=False)

    # ---- Mechanism: does ARG removal track intI1 (mobile) removal? abs vs norm gap ----
    mech_rows = []
    for line in QUAT_LINES:
        sub = rem[(rem.line == line) & (rem.gene.isin(ARGS))]
        arg_abs = sub["abs_log_removal_mean"].mean()
        arg_norm = sub["norm_log_removal_mean"].mean()
        i_abs = float(int_rem.loc[line, "abs_log_removal_mean"]) if line in int_rem.index else np.nan
        i_norm = float(int_rem.loc[line, "norm_log_removal_mean"]) if line in int_rem.index else np.nan
        mech_rows.append({
            "line": line,
            "arg_abs_removal": arg_abs, "intI1_abs_removal": i_abs,
            "arg_norm_removal": arg_norm, "intI1_norm_removal": i_norm,
            "biomass_vs_selective_gap": arg_abs - arg_norm,  # how much removal is just biomass
        })
    mech = pd.DataFrame(mech_rows)
    mech.to_csv(outdir / "integron_mechanism.csv", index=False)

    # ---- Human-readable summary ----
    L = []
    L.append("# Paper 4 (corrected) — ARG & class-1 integron removal across quaternary lines\n")
    L.append(f"Data: {args.xlsx} ({args.sheet}); campaigns {campaigns}; genes {genes}.\n")
    L.append("Positive log10 values = reduction (line has fewer copies than its OUT WWTP input).\n")
    L.append("## Per-line mean removal across ARGs (log10)\n")
    L.append("| Line | Mean ARG removal (abs) | Mean ARG removal (16S-norm) | intI1 removal (abs) | intI1 removal (norm) |")
    L.append("|---|---:|---:|---:|---:|")
    for _, r in summary.iterrows():
        L.append(f"| {LINE_LABEL[r['line']]} | {r['mean_abs_removal']:.2f} | {r['mean_norm_removal']:.2f} | "
                 f"{r['intI1_abs_removal']:.2f} | {r['intI1_norm_removal']:.2f} |")
    L.append("")
    L.append("## Removal by line × gene (absolute log10)\n")
    piv = rem.pivot(index="gene", columns="line", values="abs_log_removal_mean")[QUAT_LINES]
    L.append(piv.round(2).to_markdown())
    L.append("")
    L.append("## 16S-normalized removal by line × gene (log10)\n")
    pivn = rem.pivot(index="gene", columns="line", values="norm_log_removal_mean")[QUAT_LINES]
    L.append(pivn.round(2).to_markdown())
    L.append("")
    L.append("## Conventional plant removal (IN WWTP -> OUT WWTP, absolute log10)\n")
    L.append(conv.round(2).to_markdown(index=False))
    L.append("")
    L.append("## Mechanism (biomass vs selective)\n")
    L.append("biomass_vs_selective_gap = abs removal - 16S-normalized removal; "
             "large positive = removal driven mostly by biomass loss, not selective ARG removal.\n")
    L.append(mech.round(2).to_markdown(index=False))
    (outdir / "removal_analysis_summary.md").write_text("\n".join(L), encoding="utf-8")

    print("\n".join(L))
    print(f"\nWrote outputs to {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
