# Paper 4 Data Provenance and Analysis Correction

## Multi-site data sources (real)

The paper combines two real, published datasets, harmonized only on the unit-cancelling
**selective (community-normalized) log removal** metric:

1. **German site (anchor, detailed):** Brouwir et al. 2025, *Microorganisms* 13(12):2663,
   Supplementary Table S10 — qPCR of 16S, blaAmpC, ermB, intI1, sul1, sul2, tetW across
   IN WWTP / OUT WWTP / CW / GAC / AOP, 4 campaigns. (Details below.)
2. **Taiwan & Slovakia (conventional-treatment breadth):** Chen et al. 2024, *Heliyon* 10:e38723,
   Supplementary Table 4 — community-normalized influent/effluent relative abundance for
   Sul1, tetM, ermB, blaTEM (metagenomic and qPCR), multiple plants per country. Extracted
   programmatically by `../analysis/extract_external_sites.py` →
   `../data/external_sites_removal.csv` (16 records, audited against the source table).

The cross-site comparison (`../analysis/run_multisite_analysis.py` →
`multisite_removal.csv`, `multisite_summary.md`) computes selective log removal
= log10(influent_norm / effluent_norm) for each site, using qPCR community-normalized
values for all three sites so the methods match. Harmonized ARG classes: sulfonamide,
tetracycline, MLS, β-lactam (matched genes are related but not identical across studies —
tetW≈tetM, blaAmpC≈blaTEM — a stated limitation).

**Heterogeneity caveat (important):** these are independent studies with different ARG
panels, sequencing/qPCR protocols, and (for Taiwan/Slovakia) country-pooled reporting.
We therefore compare only the unit-cancelling selective-removal metric and emphasize the
qualitative pattern (weak/inconsistent conventional removal; sulfonamide recalcitrance)
rather than precise magnitudes. No data were simulated or fabricated.

## German-site data source (real)

All results derive from publicly available qPCR data in **Supplementary Table S10** of:

> L. Brouwir, H. KleinJan, C. Balent, G. Quabron, I. Salmerón, S. Venditti, F. Gritten,
> X. Zhao, et al., "Fate and removal of antibiotics and antibiotic resistance genes in a
> rural wastewater treatment plant: a microbial perspective of nature-based versus advanced
> technologies," *Microorganisms* 13(12):2663, 2025. doi:10.3390/microorganisms13122663

Local copy: `../data/external/oa_scan/PMC12735332/PMC12735332/supp_unzipped/Supplementary-tables-CEBEDEAU.xlsx`.

The data are qPCR gene-copy concentrations (copies/L) for 7 targets — 16S rRNA, blaAmpC,
ermB, intI1, sul1, sul2, tetW — at 5 sampling points (IN WWTP, OUT WWTP, and the CW, GAC,
and AOP=ozonation+GAC quaternary line effluents), across 4 campaigns (C2, C5, C10, C12)
with ~13–16 technical replicates each (485 measurements).

## IMPORTANT: analysis correction (difference-in-differences was invalid)

The original draft (`../Group 19 (13).pdf`) applied a **difference-in-differences (DiD)
"treatment-upgrade natural experiment."** That design was a misreading of the source data
and has been **retired**. Three problems:

1. **No upgrade event exists.** C2/C5/C10/C12 are sampling campaigns (dates), not a
   before/after-upgrade timeline. The advanced lines operated throughout; nothing was
   upgraded between campaigns. Coding early vs. late campaigns as "before/after upgrade"
   is meaningless.
2. **"OUT WWTP" is not a control plant.** It is the conventional effluent that *feeds* the
   quaternary lines. The advanced lines treat it, so the correct contrast is removal
   (OUT WWTP → line), not treated-vs-control plants.
3. The resulting "treatment increases ARGs" DiD estimate (β = +1.63×10¹², p = 0.091)
   was an artifact of (1)–(2) plus a 4-of-13-campaign subset and baseline imbalance, and
   it **contradicts the source study's finding** that the treatments remove ARGs.

The earlier draft also misdescribed the cohort as "40 samples, 8 plants, 6 countries,
SRA PRJNA1044402"; the real data are a **single site** (one plant) measured by qPCR.

### Retired DiD files (kept for transparency, not used by the current paper)
- `did_input_table.csv`, `did_coefficients.csv`, `did_model_summary.txt` — outputs of the
  invalid DiD on a 16-row campaign subset. Superseded.
- `../analysis/run_did_model.py`, `build_did_from_cebedeau_s10.py` — DiD build/fit scripts.

## Current (correct) analysis

`../analysis/run_removal_analysis.py` computes, for each quaternary line and gene, the
log₁₀ removal relative to the conventional effluent it polishes, on both an absolute
(copies/L) and a 16S-rRNA-normalized basis (separating selective ARG removal from biomass
loss), plus an intI1/sul1 integron-mechanism analysis. Outputs:
`s10_clean.csv`, `removal_by_line_gene.csv`, `removal_summary.csv`,
`integron_mechanism.csv`, `removal_analysis_summary.md`. Figures:
`../analysis/generate_figures.py`. The write-up is `../paper.md`.

## Reproduce

```bash
pip install pandas numpy scipy openpyxl matplotlib tabulate
python projects/paper4/analysis/extract_external_sites.py     # Taiwan/Slovakia -> external_sites_removal.csv
python projects/paper4/analysis/run_removal_analysis.py       # German site removal + integron mechanism
python projects/paper4/analysis/run_uncertainty.py            # per-line CIs, bootstrap, paired contrasts (n=4)
python projects/paper4/analysis/run_multisite_analysis.py     # cross-site selective removal
python projects/paper4/analysis/generate_figures.py           # 5 figures
```

## Data dictionary (derived tables in this folder)

| File | Contents |
|---|---|
| `s10_clean.csv` | tidy German-site qPCR: Sample, Campaign, Gene, geometric-mean conc (copies/L) |
| `removal_by_line_gene.csv` | per line × gene absolute & 16S-normalized log10 removal (mean, SD, n_campaigns) |
| `removal_summary.csv` | per-line mean ARG removal + intI1 removal (abs/norm) |
| `integron_mechanism.csv` | per line: ARG vs intI1 removal (abs/norm), biomass_vs_selective_gap |
| `uncertainty_by_line.csv` | per line × metric: mean, SD, n, t-CI, bootstrap CI |
| `pairwise_line_diff.csv` | paired within-campaign line contrasts (CW/GAC/AOP), mean diff + CIs |
| `external_sites_removal.csv` (in `../data/`) | Taiwan/Slovakia influent/effluent rel. abundance (sourced) |
| `multisite_removal.csv` | cross-site selective log10 removal by class × site |
| `*_summary.md` | human-readable summaries of the above |

Metric definitions: **absolute removal** = log10(C_OUT/C_line); **selective removal** =
log10(R_OUT/R_line) with R = gene/16S; **biomass gap** = mean absolute − mean selective.
Positive = reduction; ≤ 0 = unchanged or enriched. German-site inference uses n = 4
campaigns (CIs, not p-values); the cross-site table is point estimates only (qualitative).

## Honest limitations

Single site (Bliesen, DE), 6 ARGs + intI1 by qPCR, 4 campaigns (n=4 per removal estimate;
influent only in 2), technical (not biological) replicates, high per-campaign variance,
removal inferred from parallel-line concentration ratios (not paired water-parcel tracking),
intI1 as a mobility proxy (no direct transfer assay). The qualitative pattern (selective vs.
biomass removal; integron retention by ozonation) is robust and internally consistent;
precise magnitudes are indicative. Generalization needs multi-site, multi-campaign,
metagenomic + transfer-level confirmation.
