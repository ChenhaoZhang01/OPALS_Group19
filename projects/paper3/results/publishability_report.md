# Paper 3 Publishability Report (Lag Analysis)

## Summary
- Strict lag pairs: 19
- Strict lag studies: 4
- Real-feature template rows: 23
- Legacy feature ID overlap with real metadata: 0

## Gate Checks
- Temporal ordering has sufficient strict pairs (>=10): PASS
- Temporal ordering covers >=3 studies: PASS
- Real ARG/MGE/entropy quantification available (>=10 complete rows): PASS
- Legacy feature table IDs align with real metadata: FAIL

## Blocking Issues
- Legacy features_table.csv appears synthetic or unmapped to real metadata (zero sample_id overlap).

## Immediate Next Actions
1. Populate real arg_total, mge_abundance, and entropy for the strict sample set in features_real_template.csv.
2. Re-run run_lag_analysis.py using the real populated feature table and archive outputs under a new filename set.
3. Remove or clearly label legacy synthetic outputs to avoid provenance confusion during submission.
