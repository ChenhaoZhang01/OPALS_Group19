#!/usr/bin/env bash
set -euo pipefail

# Download SRA runs and convert to FASTQ for downstream assembly.
# Usage: bash scripts/download_reads.sh metadata/master_metadata.csv data/raw

METADATA_CSV="${1:-metadata/master_metadata.csv}"
RAW_DIR="${2:-data/raw}"

mkdir -p "$RAW_DIR"

if ! command -v prefetch >/dev/null 2>&1; then
  echo "ERROR: prefetch not found in PATH (install NCBI SRA Toolkit)."
  exit 1
fi

if ! command -v fasterq-dump >/dev/null 2>&1; then
  echo "ERROR: fasterq-dump not found in PATH (install NCBI SRA Toolkit)."
  exit 1
fi

tail -n +2 "$METADATA_CSV" | cut -d',' -f1 | while read -r SAMPLE; do
  if [[ -z "$SAMPLE" ]]; then
    continue
  fi

  echo "Downloading $SAMPLE"
  prefetch "$SAMPLE"

  echo "Converting $SAMPLE to FASTQ"
  fasterq-dump "$SAMPLE" -O "$RAW_DIR"

done
