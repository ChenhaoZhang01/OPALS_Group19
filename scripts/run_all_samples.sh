#!/usr/bin/env bash
set -euo pipefail

# Run ARG pipeline for all samples listed in master_metadata.csv, then build ARG matrix.
# Usage:
# bash scripts/run_all_samples.sh metadata/master_metadata.csv data/raw db/CARD.fasta db/CARD results

METADATA_CSV="${1:-metadata/master_metadata.csv}"
RAW_DIR="${2:-data/raw}"
CARD_FASTA="${3:-db/CARD.fasta}"
CARD_DIAMOND_PREFIX="${4:-db/CARD}"
RESULTS_DIR="${5:-results}"

mkdir -p "$RESULTS_DIR"

tail -n +2 "$METADATA_CSV" | cut -d',' -f1 | while read -r SAMPLE; do
  if [[ -z "$SAMPLE" ]]; then
    continue
  fi

  R1="$RAW_DIR/${SAMPLE}_1.fastq"
  R2="$RAW_DIR/${SAMPLE}_2.fastq"

  if [[ -f "$R1" && -f "$R2" ]]; then
    bash scripts/run_arg_pipeline.sh "$SAMPLE" "$R1" "$R2" "$CARD_FASTA" "$CARD_DIAMOND_PREFIX" "$RESULTS_DIR"
  elif [[ -f "$RAW_DIR/${SAMPLE}.fastq" ]]; then
    bash scripts/run_arg_pipeline.sh "$SAMPLE" "$RAW_DIR/${SAMPLE}.fastq" NA "$CARD_FASTA" "$CARD_DIAMOND_PREFIX" "$RESULTS_DIR"
  else
    echo "WARNING: no FASTQ files found for $SAMPLE, skipping"
  fi

done

python scripts/build_arg_matrix.py --hits-glob "$RESULTS_DIR/*/arg_hits.tsv" --output "$RESULTS_DIR/ARG_matrix.csv"
