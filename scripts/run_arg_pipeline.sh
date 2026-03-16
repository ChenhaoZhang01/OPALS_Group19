#!/usr/bin/env bash
set -euo pipefail

# End-to-end sample pipeline:
# reads -> assembly -> gene prediction -> ARG detection
#
# Usage:
# bash scripts/run_arg_pipeline.sh SAMPLE_ID READS_1 [READS_2] CARD_DB DIAMOND_DB_PREFIX OUTPUT_DIR
#
# Example (paired):
# bash scripts/run_arg_pipeline.sh SRR12345 data/raw/SRR12345_1.fastq data/raw/SRR12345_2.fastq db/CARD.fasta db/CARD results
#
# Example (single-end):
# bash scripts/run_arg_pipeline.sh SRR12345 data/raw/SRR12345.fastq NA db/CARD.fasta db/CARD results

SAMPLE_ID="$1"
READS_1="$2"
READS_2="$3"
CARD_FASTA="$4"
CARD_DIAMOND_PREFIX="$5"
OUTPUT_DIR="$6"

SAMPLE_DIR="$OUTPUT_DIR/$SAMPLE_ID"
ASSEMBLY_DIR="$SAMPLE_DIR/assembly"

mkdir -p "$SAMPLE_DIR" "$ASSEMBLY_DIR"

if [[ ! -f "$CARD_DIAMOND_PREFIX.dmnd" ]]; then
  echo "Building DIAMOND database from CARD fasta"
  diamond makedb --in "$CARD_FASTA" -d "$CARD_DIAMOND_PREFIX"
fi

echo "[$SAMPLE_ID] Assembling reads with MEGAHIT"
if [[ "$READS_2" == "NA" ]]; then
  megahit -r "$READS_1" -o "$ASSEMBLY_DIR"
else
  megahit -1 "$READS_1" -2 "$READS_2" -o "$ASSEMBLY_DIR"
fi

ASSEMBLY_FA="$ASSEMBLY_DIR/final.contigs.fa"
if [[ ! -f "$ASSEMBLY_FA" ]]; then
  echo "ERROR: expected assembly output not found at $ASSEMBLY_FA"
  exit 1
fi

echo "[$SAMPLE_ID] Predicting genes with Prodigal"
prodigal -i "$ASSEMBLY_FA" -a "$SAMPLE_DIR/proteins.faa" -d "$SAMPLE_DIR/genes.fna" -p meta

echo "[$SAMPLE_ID] Detecting ARGs with DIAMOND against CARD"
diamond blastp \
  -d "$CARD_DIAMOND_PREFIX" \
  -q "$SAMPLE_DIR/proteins.faa" \
  -o "$SAMPLE_DIR/arg_hits.tsv" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --id 80 \
  --query-cover 70 \
  --max-target-seqs 1

echo "[$SAMPLE_ID] Done: $SAMPLE_DIR/arg_hits.tsv"
