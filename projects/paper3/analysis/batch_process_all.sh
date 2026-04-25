#!/usr/bin/env bash
# batch_process_all.sh — process all Paper 3 shotgun samples for entropy + ARG
# Runs Kraken2+Bracken (entropy) and DIAMOND blastx (arg_total) on each sample.
# MEGAHIT assembly for mge_abundance is handled separately (run_single_sample.sh).
#
# Usage (in WSL2 Ubuntu):
#   bash /mnt/c/path/to/batch_process_all.sh
# or copy to WSL and run directly.
#
# Prereqs (already set up in biotools env):
#   /opt/miniforge/envs/biotools/bin/{kraken2,bracken,diamond,prefetch,fasterq-dump}
#   /data/dbs/kraken2_8gb           — Kraken2 standard 8GB DB
#   /data/dbs/amrfinder_diamond/AMRProt.dmnd — DIAMOND ARG database

set -euo pipefail

OUTDIR="/data/paper3_quant"
FASTQ_DIR="$OUTDIR/fastq"
KRAKEN_DIR="$OUTDIR/kraken"
BRACKEN_DIR="$OUTDIR/bracken"
ARG_DIR="$OUTDIR/arg"
QUANT_CSV="$OUTDIR/quant_results.csv"
KRAKEN_DB="/data/dbs/kraken2_8gb"
ARG_DB="/data/dbs/amrfinder_diamond/AMRProt.dmnd"
THREADS=4
ENV_BIN="/opt/miniforge/envs/biotools/bin"

mkdir -p "$FASTQ_DIR" "$KRAKEN_DIR" "$BRACKEN_DIR" "$ARG_DIR"
[[ ! -f "$QUANT_CSV" ]] && echo "sample_id,arg_total,mge_abundance,entropy" > "$QUANT_CSV"

SHOTGUN_SAMPLES=(
  # PRJNA1071831 — Iskar River Bulgaria (biweekly, Nov-Dec 2022)
  SRR27827413 SRR27827408 SRR27827406 SRR27827404   # Dragushinovo
  SRR27827412 SRR27827407 SRR27827405 SRR27827403   # Mechkata
  # PRJNA599167_WGA — Chesapeake Bay WGA shotgun (May 2017, 1.5h intervals)
  SRR11803495 SRR11803496 SRR11803497 SRR11803498
  SRR11803499 SRR11803500 SRR11803501 SRR11803502
  SRR11803503 SRR11803504 SRR11803505
)

compute_entropy() {
  local f="$1"
  "$ENV_BIN/python3" - "$f" <<'PYEOF'
import sys, math
fname = sys.argv[1]
fracs = []
with open(fname) as fh:
    next(fh)
    for line in fh:
        parts = line.rstrip().split('\t')
        if len(parts) >= 7:
            try:
                v = float(parts[6])
                if v > 0:
                    fracs.append(v)
            except ValueError:
                pass
if not fracs:
    print("0.0")
else:
    total = sum(fracs)
    h = -sum((p/total)*math.log2(p/total) for p in fracs if p > 0)
    print(f"{h:.6f}")
PYEOF
}

process_sample() {
  local srr="$1"

  if grep -q "^${srr}," "$QUANT_CSV" 2>/dev/null; then
    echo "=== $srr already in $QUANT_CSV — skipping ==="
    return
  fi

  echo "=== Processing: $srr ==="

  # --- 1. Download ---
  local fq1="$FASTQ_DIR/${srr}_1.fastq"
  local fq2="$FASTQ_DIR/${srr}_2.fastq"
  local fq_se="$FASTQ_DIR/${srr}.fastq"

  if [[ ! -f "$fq1" && ! -f "${fq1}.gz" && ! -f "$fq_se" ]]; then
    echo "  Downloading $srr..."
    "$ENV_BIN/prefetch" --max-size 60G "$srr" -O "$FASTQ_DIR" 2>/dev/null || true
    "$ENV_BIN/fasterq-dump" --split-files --threads "$THREADS" \
      --outdir "$FASTQ_DIR" "$FASTQ_DIR/$srr" 2>/dev/null || \
    "$ENV_BIN/fasterq-dump" --split-files --threads "$THREADS" \
      --outdir "$FASTQ_DIR" "$srr" 2>/dev/null || true
  fi

  local is_paired=false
  [[ -f "$fq1" && -f "$fq2" ]] && is_paired=true

  if [[ "$is_paired" == false && ! -f "$fq_se" ]]; then
    echo "  ERROR: no FASTQ found for $srr"
    echo "${srr},0,0,0" >> "$QUANT_CSV"
    return
  fi

  # --- 2. Kraken2 ---
  local kraken_report="$KRAKEN_DIR/${srr}.report"
  if [[ ! -f "$kraken_report" ]]; then
    echo "  Kraken2..."
    if [[ "$is_paired" == true ]]; then
      "$ENV_BIN/kraken2" --db "$KRAKEN_DB" --paired \
        --threads "$THREADS" --report "$kraken_report" \
        --output /dev/null "$fq1" "$fq2" 2>/dev/null || touch "$kraken_report"
    else
      "$ENV_BIN/kraken2" --db "$KRAKEN_DB" \
        --threads "$THREADS" --report "$kraken_report" \
        --output /dev/null "$fq_se" 2>/dev/null || touch "$kraken_report"
    fi
  fi

  # --- 3. Bracken ---
  local bracken_report="$BRACKEN_DIR/${srr}.bracken"
  if [[ ! -f "$bracken_report" && -s "$kraken_report" ]]; then
    echo "  Bracken..."
    /opt/miniforge/bin/conda run -n biotools bracken \
      -d "$KRAKEN_DB" -i "$kraken_report" -o "$bracken_report" \
      -r 150 -l S -t 10 2>/dev/null || touch "$bracken_report"
  fi

  local entropy=0
  if [[ -f "$bracken_report" && -s "$bracken_report" ]]; then
    entropy=$(compute_entropy "$bracken_report")
  fi
  echo "  Entropy: $entropy"

  # --- 4. DIAMOND ARG (read-based, --block-size 0.5 to stay within RAM) ---
  local arg_out="$ARG_DIR/${srr}.m8"
  if [[ ! -f "$arg_out" || ! -s "$arg_out" ]]; then
    echo "  DIAMOND ARG..."
    local r1_out="$ARG_DIR/${srr}_r1.m8"
    local r2_out="$ARG_DIR/${srr}_r2.m8"
    if [[ "$is_paired" == true ]]; then
      "$ENV_BIN/diamond" blastx \
        --db "$ARG_DB" --query "$fq1" --out "$r1_out" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 1 --evalue 1e-5 --id 70 --query-cover 60 \
        --block-size 0.5 --threads "$THREADS" --quiet 2>/dev/null || touch "$r1_out"
      "$ENV_BIN/diamond" blastx \
        --db "$ARG_DB" --query "$fq2" --out "$r2_out" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 1 --evalue 1e-5 --id 70 --query-cover 60 \
        --block-size 0.5 --threads "$THREADS" --quiet 2>/dev/null || touch "$r2_out"
      cat "$r1_out" "$r2_out" > "$arg_out" 2>/dev/null || touch "$arg_out"
    else
      "$ENV_BIN/diamond" blastx \
        --db "$ARG_DB" --query "$fq_se" --out "$arg_out" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 1 --evalue 1e-5 --id 70 --query-cover 60 \
        --block-size 0.5 --threads "$THREADS" --quiet 2>/dev/null || touch "$arg_out"
    fi
  fi

  local arg_total=0
  if [[ -f "$arg_out" && -s "$arg_out" ]]; then
    arg_total=$(awk '{print $2}' "$arg_out" | sort -u | wc -l)
  fi
  echo "  arg_total: $arg_total"

  # mge_abundance: 0 here (requires MEGAHIT+IntegronFinder; run separately)
  echo "${srr},${arg_total},0,${entropy}" >> "$QUANT_CSV"
  echo "  Done: $srr  ARG=$arg_total  H=$entropy  (MGE=0 pending assembly)"
}

echo "=== Batch pipeline START ==="
for srr in "${SHOTGUN_SAMPLES[@]}"; do
  process_sample "$srr"
done
echo "=== Batch pipeline COMPLETE ==="
echo "Results in: $QUANT_CSV"
