#!/usr/bin/env bash
set -euo pipefail

# Download SRA reads from accession list and convert to FASTQ.
# Usage: bash scripts/download_from_list.sh download_list.txt data/raw

LIST_FILE="${1:-download_list.txt}"
OUT_DIR="${2:-data/raw}"

if [[ ! -f "$LIST_FILE" ]]; then
  echo "ERROR: list file not found: $LIST_FILE" >&2
  exit 1
fi

if ! command -v prefetch >/dev/null 2>&1; then
  echo "ERROR: prefetch not found in PATH (install SRA Toolkit)." >&2
  exit 1
fi

if ! command -v fastq-dump >/dev/null 2>&1; then
  echo "ERROR: fastq-dump not found in PATH (install SRA Toolkit)." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

while IFS= read -r sample; do
  sample="${sample//$'\r'/}"
  if [[ -z "$sample" ]]; then
    continue
  fi

  echo "Downloading $sample"
  prefetch "$sample"

  echo "Converting $sample to FASTQ"
  fastq-dump --split-files --outdir "$OUT_DIR" "$sample"
done < "$LIST_FILE"

echo "Done. FASTQ files are in: $OUT_DIR"
