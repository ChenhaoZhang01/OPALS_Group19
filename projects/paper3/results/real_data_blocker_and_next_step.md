# Paper 3 Real-Data Status and Next Steps

## What is already completed

- Temporal metadata enriched from ENA/SRA/BioSample for all 23 strict-cohort samples.
- Metadata expanded to 62 samples across 4 virtual study IDs with real collection timestamps.
- Strict lag cohort derived: 19 strict consecutive-timepoint pairs from 4 studies — both temporal gates now PASS.
- Feature template (features_real_template.csv) generated with all 23 strict-usable sample IDs.
- Quant input template (quant_input_template.csv) updated to all 23 samples.
- Bioinformatics pipeline script written (analysis/quantify_samples.sh).
- populate_real_features_from_quant.py ready to ingest pipeline outputs.

## Current gate status

| Gate | Status |
|---|---|
| Strict lag pairs ≥ 10 | **PASS** (19 pairs) |
| Strict studies ≥ 3 | **PASS** (4 studies) |
| Real ARG/MGE/entropy ≥ 10 complete rows | **FAIL** — pipeline not yet run |
| Legacy feature ID alignment | **FAIL** — synthetic table unrelated |

## Strict cohort (23 samples requiring quantification)

### PRJNA599167_WGA — Chesapeake Bay WGA shotgun, May 17–19 2017
Full ARG + MGE + entropy quantification applicable.
- SRR11803495, SRR11803496, SRR11803497, SRR11803498, SRR11803499
- SRR11803500, SRR11803501, SRR11803502, SRR11803503, SRR11803504, SRR11803505

### PRJNA1071831_Drag — Iskar River, Dragushinovo Bulgaria, Nov–Dec 2022
Full ARG + MGE + entropy quantification applicable.
- SRR27827413 (2022-11-03), SRR27827408 (2022-11-17), SRR27827406 (2022-12-08), SRR27827404 (2022-12-21)

### PRJNA1071831_Mech — Iskar River, Mechkata Bulgaria, Nov–Dec 2022
Full ARG + MGE + entropy quantification applicable.
- SRR27827412 (2022-11-03), SRR27827407 (2022-11-17), SRR27827405 (2022-12-08), SRR27827403 (2022-12-21)

### PRJNA599167 — Chesapeake Bay 16S amplicon, May 2017
Entropy only (16S cannot yield ARG/MGE counts). These samples contribute to diversity
analysis but will be excluded from ARG-lag regressions.
- SRR11811727, SRR11811728, SRR11811733, SRR11811738

## How to run the quantification pipeline

Requires a Linux environment (HPC recommended) with:
- sra-tools (prefetch, fasterq-dump)
- AMRFinderPlus with updated database
- DIAMOND + mobileOG-db .dmnd index
- Kraken2 + Bracken + standard PlusPF database
- metaSPAdes or MEGAHIT
- Python 3

```bash
# From a Linux HPC session:
bash projects/paper3/analysis/quantify_samples.sh \
  --outdir /scratch/paper3_quant \
  --kraken-db /databases/kraken2_plusPF \
  --mobileog-dmnd /databases/mobileOG-db.dmnd \
  --threads 16
```

## What runs automatically after pipeline completes

```bash
# From OPALS_Group19 repo root:
python projects/paper3/analysis/populate_real_features_from_quant.py \
  --quant /scratch/paper3_quant/quant_results.csv

python projects/paper3/analysis/run_lag_analysis.py \
  --features projects/paper3/results/features_real_populated.csv \
  --outdir projects/paper3/results/real_run

python projects/paper3/analysis/generate_publishability_report.py
```

## Expected outcome

With 11 WGA + 8 Bulgaria = 19 fully quantifiable samples, the regression will have:
- ~16 strict lag pairs (excluding 3 WGA samples with no read_count)
- Coverage across 3 independently sampled studies
- Temporal scales from hours (Chesapeake Bay) to weeks (Bulgaria river)
