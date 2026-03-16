#!/usr/bin/env bash
set -euo pipefail

# Run FastQC on all FASTQ files in a directory.
# Usage: bash scripts/run_fastqc_all.sh data/raw results/fastqc

READS_DIR="${1:-data/raw}"
OUT_DIR="${2:-results/fastqc}"

if ! command -v fastqc >/dev/null 2>&1; then
  echo "ERROR: fastqc not found in PATH." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

shopt -s nullglob
fastq_files=("$READS_DIR"/*.fastq)
shopt -u nullglob

if [[ ${#fastq_files[@]} -eq 0 ]]; then
  echo "ERROR: no .fastq files found in $READS_DIR" >&2
  exit 1
fi

fastqc -o "$OUT_DIR" "${fastq_files[@]}"

echo "Done. FastQC reports are in: $OUT_DIR"
