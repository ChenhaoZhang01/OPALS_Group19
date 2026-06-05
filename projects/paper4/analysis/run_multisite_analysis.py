#!/usr/bin/env python3
"""Multi-site selective (community-normalized) ARG removal through conventional treatment.

Combines:
  * Bliesen, Germany (this study's anchor): qPCR, IN WWTP -> OUT WWTP, 16S-normalized
    removal for sul1, tetW, ermB, blaAmpC  (from s10_clean.csv).
  * Taiwan & Slovakia: community-normalized influent->effluent removal for Sul1, tetM,
    ermB, blaTEM (from external_sites_removal.csv; PMC11471163 Table 4).

Common metric: selective log10 removal = log10( influent_norm / effluent_norm ), where
"norm" = ARG relative to total bacteria (16S for qPCR; relative abundance for metagenomics).
Positive = ARG reduced relative to the community; <=0 = unchanged or ENRICHED by treatment.

ARG classes harmonized across sites: sulfonamide, tetracycline, MLS (macrolide), beta-lactam.

Outputs (projects/paper4/results/):
  multisite_removal.csv          tidy per site x class x method selective log removal
  multisite_summary.md           human-readable cross-site summary
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path("projects/paper4/results")
EXTERNAL = "projects/paper4/data/external_sites_removal.csv"
# Bliesen gene -> harmonized class
BLIESEN_CLASS = {"sul1": "sulfonamide", "tetW": "tetracycline",
                 "ermB": "MLS", "blaAmpC": "beta-lactam"}


def bliesen_conventional_removal() -> list[dict]:
    """IN WWTP -> OUT WWTP, 16S-normalized, per shared class (averaged over campaigns)."""
    clean = pd.read_csv(RESULTS / "s10_clean.csv")
    wide = clean.pivot_table(index=["Sample", "Campaign"], columns="Gene", values="conc_gm")
    norm = wide.div(wide["16S rRNA"], axis=0)
    campaigns = sorted({c for (s, c) in wide.index})
    recs = []
    for gene, cls in BLIESEN_CLASS.items():
        vals = []
        for c in campaigns:
            ki, ko = ("IN WWTP", c), ("OUT WWTP", c)
            if ki in norm.index and ko in norm.index:
                ni, no = norm.loc[ki, gene], norm.loc[ko, gene]
                if pd.notna(ni) and pd.notna(no) and ni > 0 and no > 0:
                    vals.append(np.log10(ni / no))
        if vals:
            recs.append({"country": "Germany", "site": "Bliesen", "method": "qPCR",
                         "arg_class": cls, "gene": gene,
                         "selective_log_removal": float(np.mean(vals)), "n": len(vals)})
    return recs


def external_removal(method: str = "qPCR") -> list[dict]:
    df = pd.read_csv(EXTERNAL)
    df = df[df["method"] == method]
    recs = []
    for _, r in df.iterrows():
        inf, eff = r["influent_relabund"], r["effluent_relabund"]
        if pd.notna(inf) and pd.notna(eff) and inf > 0 and eff > 0:
            recs.append({"country": r["country"], "site": r["site"], "method": method,
                         "arg_class": r["arg_class"], "gene": r["gene"],
                         "selective_log_removal": float(np.log10(inf / eff)), "n": 1})
    return recs


def main() -> int:
    recs = bliesen_conventional_removal() + external_removal("qPCR")
    df = pd.DataFrame(recs)
    df.to_csv(RESULTS / "multisite_removal.csv", index=False)

    # cross-site summary
    classes = ["sulfonamide", "tetracycline", "MLS", "beta-lactam"]
    sites = ["Bliesen", "Taiwan", "Slovakia"]
    piv = (df.pivot_table(index="arg_class", columns="site",
                          values="selective_log_removal", aggfunc="mean")
             .reindex(index=classes, columns=sites))

    L = ["# Multi-site selective (community-normalized) ARG removal — conventional treatment\n",
         "Selective log10 removal = log10(influent_norm / effluent_norm); "
         "positive = reduced relative to community, <=0 = unchanged or ENRICHED.\n",
         "Sites: Bliesen (Germany, qPCR 16S-norm, this study) | Taiwan, Slovakia "
         "(qPCR 16S-norm, PMC11471163).\n",
         "## Selective log10 removal by ARG class x site (qPCR)\n",
         piv.round(2).to_markdown(), "",
         "## Fraction of (class x site) cases with NO selective removal (<=0)\n"]
    nonpos = (df["selective_log_removal"] <= 0).sum()
    L.append(f"- {nonpos}/{len(df)} site×class cases show no selective removal "
             f"(ARG unchanged or enriched relative to the community).")
    L.append(f"- Mean selective removal across all cases: "
             f"{df['selective_log_removal'].mean():.2f} log10.\n")
    L.append("## Interpretation\n")
    L.append("Conventional treatment does not reliably reduce the community-normalized "
             "resistome across sites; several ARG classes are unchanged or enriched in "
             "effluent. This generalizes the single-site finding that bulk ARG removal is "
             "largely biomass-driven, and motivates advanced quaternary treatment — which "
             "(at Bliesen) achieves selective removal but still spares the mobile integron.")
    (RESULTS / "multisite_summary.md").write_text("\n".join(L), encoding="utf-8")

    print("\n".join(L))
    print(f"\nWrote {RESULTS}/multisite_removal.csv and multisite_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
