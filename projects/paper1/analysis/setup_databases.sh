#!/usr/bin/env bash
# Paper 1 -- download and prepare the CARD database for all three pipelines.
#
#   Pipeline A (DIAMOND blastp)  -> protein homolog model FASTA -> DIAMOND .dmnd
#   Pipeline B (Bowtie2 mapping) -> nucleotide homolog model FASTA -> bowtie2 index
#   Pipeline C (RGI)             -> card.json loaded into RGI's local database
#
# Usage:
#   bash projects/paper1/analysis/setup_databases.sh --dbdir /databases/card
#
# Requires (from environment-linux.yml): diamond, bowtie2, rgi, wget (or curl), tar.
# Run once before run_arg_pipeline_hpc.sh.

set -euo pipefail

DBDIR="/databases/card"
CARD_URL="https://card.mcmaster.ca/latest/data"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dbdir) DBDIR="$2"; shift 2;;
    --url)   CARD_URL="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

mkdir -p "$DBDIR"
cd "$DBDIR"

echo "== Downloading CARD =="
if command -v wget >/dev/null 2>&1; then
  wget -O card-data.tar.bz2 "$CARD_URL"
else
  curl -L -o card-data.tar.bz2 "$CARD_URL"
fi
tar -xjf card-data.tar.bz2

# Expected files after extraction:
#   card.json
#   protein_fasta_protein_homolog_model.fasta
#   nucleotide_fasta_protein_homolog_model.fasta
#   aro_index.tsv
PROT="$DBDIR/protein_fasta_protein_homolog_model.fasta"
NUCL="$DBDIR/nucleotide_fasta_protein_homolog_model.fasta"
JSON="$DBDIR/card.json"
ARO="$DBDIR/aro_index.tsv"

for f in "$PROT" "$NUCL" "$JSON"; do
  [[ -f "$f" ]] || { echo "ERROR: expected CARD file missing: $f" >&2; exit 1; }
done

echo "== Building DIAMOND protein DB (Pipeline A) =="
diamond makedb --in "$PROT" -d "$DBDIR/card_protein"

echo "== Building Bowtie2 nucleotide index (Pipeline B) =="
bowtie2-build "$NUCL" "$DBDIR/card_bt2"

echo "== Loading CARD into RGI (Pipeline C) =="
rgi load --card_json "$JSON" --local

cat > "$DBDIR/DB_PATHS.env" <<EOF
# Source this before running the pipeline.
export CARD_DIAMOND="$DBDIR/card_protein.dmnd"
export CARD_BT2="$DBDIR/card_bt2"
export CARD_JSON="$JSON"
export CARD_ARO_INDEX="$ARO"
EOF

echo "Done. Database paths written to $DBDIR/DB_PATHS.env"
echo "Next: bash projects/paper1/analysis/run_arg_pipeline_hpc.sh --card-dir $DBDIR ..."
