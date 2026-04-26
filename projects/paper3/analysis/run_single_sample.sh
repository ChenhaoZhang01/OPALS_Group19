#!/usr/bin/env bash
# run_single_sample.sh  — process one SRR through the Paper 3 quantification pipeline
# Auto-detects all CPU cores and available RAM for maximum throughput.
#
# Usage: bash run_single_sample.sh <SRR_ID> [--no-assembly]
#
# Steps:
#   1. fasterq-dump  → paired FASTQ
#   2. Kraken2 + Bracken on reads  → Shannon entropy
#   3. DIAMOND blastx reads vs AMRProt  → arg_total
#   4. MEGAHIT assembly (unless --no-assembly) + IntegronFinder  → mge_abundance

set -euo pipefail

SRR="${1:?Usage: $0 <SRR_ID> [--no-assembly]}"
NO_ASSEMBLY=false
[[ "${2:-}" == "--no-assembly" ]] && NO_ASSEMBLY=true

CONDA_BIN="/opt/miniforge/bin/conda run -n biotools"
OUTDIR="/data/paper3_quant"
FASTQ_DIR="$OUTDIR/fastq"
KRAKEN_DIR="$OUTDIR/kraken"
BRACKEN_DIR="$OUTDIR/bracken"
ARG_DIR="$OUTDIR/arg"
MGE_DIR="$OUTDIR/mge"
ASSEMBLE_DIR="$OUTDIR/assembly"
QUANT_CSV="$OUTDIR/quant_results.csv"

KRAKEN_DB="/data/dbs/kraken2_8gb"
ARG_DMND="/data/dbs/amrfinder_diamond/AMRProt.dmnd"

# ---------------------------------------------------------------------------
# Auto-detect hardware
# ---------------------------------------------------------------------------
THREADS=$(nproc)
TOTAL_RAM_MB=$(free -m | awk '/^Mem:/{print $2}')

# DIAMOND block-size: ~6 GB RAM per unit; use 75% of RAM minus 8GB for Kraken2 DB
DIAMOND_BLOCK=$(python3 -c "
ram_gb = $TOTAL_RAM_MB / 1024
block = max(0.5, round((ram_gb * 0.75 - 8) / 6, 1))
print(block)
")

# MEGAHIT memory: 85% of total RAM
MEGAHIT_MEM=$(python3 -c "print(int($TOTAL_RAM_MB * 0.85 * 1024 * 1024))")  # bytes

echo "============================================================"
echo " Sample  : $SRR"
echo " Hardware: ${THREADS} cores | ${TOTAL_RAM_MB} MB RAM"
echo " DIAMOND block-size: ${DIAMOND_BLOCK}"
echo " MEGAHIT memory: $(( MEGAHIT_MEM / 1073741824 )) GB"
echo "============================================================"

mkdir -p "$FASTQ_DIR" "$KRAKEN_DIR" "$BRACKEN_DIR" "$ARG_DIR" "$MGE_DIR" "$ASSEMBLE_DIR"
[[ ! -f "$QUANT_CSV" ]] && echo "sample_id,arg_total,mge_abundance,entropy" > "$QUANT_CSV"

if grep -q "^${SRR}," "$QUANT_CSV" 2>/dev/null; then
  echo "$SRR already in $QUANT_CSV — skipping"
  exit 0
fi

echo "=== $SRR ==="

##############################################################################
# 1. Download FASTQ
##############################################################################
FQ1="$FASTQ_DIR/${SRR}_1.fastq"
FQ2="$FASTQ_DIR/${SRR}_2.fastq"
FQ_SE="$FASTQ_DIR/${SRR}.fastq"

if [[ ! -f "$FQ1" && ! -f "${FQ1}.gz" && ! -f "$FQ_SE" ]]; then
  echo "  Downloading $SRR..."
  $CONDA_BIN prefetch --max-size 60G "$SRR" -O "$FASTQ_DIR" 2>/dev/null || true
  $CONDA_BIN fasterq-dump \
    --split-files --threads "$THREADS" --bufsize 1000MB --curcache 10000MB \
    --outdir "$FASTQ_DIR" "$FASTQ_DIR/$SRR" 2>/dev/null || \
  $CONDA_BIN fasterq-dump \
    --split-files --threads "$THREADS" --bufsize 1000MB --curcache 10000MB \
    --outdir "$FASTQ_DIR" "$SRR" 2>/dev/null || true
fi

IS_PAIRED=false
if [[ -f "$FQ1" && -f "$FQ2" ]]; then
  IS_PAIRED=true
elif [[ -f "${FQ1}.gz" && -f "${FQ2}.gz" ]]; then
  echo "  Decompressing..."
  gunzip -f "${FQ1}.gz" "${FQ2}.gz"
  IS_PAIRED=true
fi

if [[ "$IS_PAIRED" == false && ! -f "$FQ_SE" ]]; then
  echo "  ERROR: no FASTQ found for $SRR"
  echo "${SRR},0,0,0" >> "$QUANT_CSV"
  exit 1
fi

##############################################################################
# 2. Kraken2 + Bracken → Shannon entropy
##############################################################################
KRAKEN_REPORT="$KRAKEN_DIR/${SRR}.report"
if [[ ! -f "$KRAKEN_REPORT" ]]; then
  echo "  Kraken2 (${THREADS} threads)..."
  if [[ "$IS_PAIRED" == true ]]; then
    $CONDA_BIN kraken2 --db "$KRAKEN_DB" --paired \
      --threads "$THREADS" --report "$KRAKEN_REPORT" \
      --output /dev/null "$FQ1" "$FQ2" 2>/dev/null || touch "$KRAKEN_REPORT"
  else
    $CONDA_BIN kraken2 --db "$KRAKEN_DB" \
      --threads "$THREADS" --report "$KRAKEN_REPORT" \
      --output /dev/null "$FQ_SE" 2>/dev/null || touch "$KRAKEN_REPORT"
  fi
fi

BRACKEN_REPORT="$BRACKEN_DIR/${SRR}.bracken"
if [[ ! -f "$BRACKEN_REPORT" && -s "$KRAKEN_REPORT" ]]; then
  echo "  Bracken..."
  $CONDA_BIN bracken -d "$KRAKEN_DB" \
    -i "$KRAKEN_REPORT" -o "$BRACKEN_REPORT" \
    -r 150 -l S -t 10 2>/dev/null || touch "$BRACKEN_REPORT"
fi

ENTROPY=0
if [[ -f "$BRACKEN_REPORT" && -s "$BRACKEN_REPORT" ]]; then
  ENTROPY=$($CONDA_BIN python3 - "$BRACKEN_REPORT" <<'PYEOF'
import sys, math
fname = sys.argv[1]
fracs = []
with open(fname) as fh:
    next(fh)
    for line in fh:
        parts = line.rstrip().split("\t")
        if len(parts) >= 7:
            try:
                fraction = float(parts[6])
                if fraction > 0:
                    fracs.append(fraction)
            except ValueError:
                pass
if not fracs:
    print("0.0")
else:
    total = sum(fracs)
    h = -sum((p/total)*math.log2(p/total) for p in fracs if p > 0)
    print(f"{h:.6f}")
PYEOF
)
fi
echo "  Entropy: $ENTROPY"

##############################################################################
# 3. DIAMOND blastx reads vs AMRProt → arg_total
#    Combine paired reads in one pass (avoids loading DB twice)
##############################################################################
ARG_OUT="$ARG_DIR/${SRR}.m8"
if [[ ! -f "$ARG_OUT" || ! -s "$ARG_OUT" ]]; then
  echo "  DIAMOND (block=${DIAMOND_BLOCK}, ${THREADS} threads)..."
  if [[ "$IS_PAIRED" == true ]]; then
    cat "$FQ1" "$FQ2" | \
      $CONDA_BIN diamond blastx \
        --db "$ARG_DMND" --query - --out "$ARG_OUT" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 1 --evalue 1e-5 --id 70 --query-cover 60 \
        --block-size "$DIAMOND_BLOCK" --threads "$THREADS" --quiet 2>/dev/null || touch "$ARG_OUT"
  else
    $CONDA_BIN diamond blastx \
      --db "$ARG_DMND" --query "$FQ_SE" --out "$ARG_OUT" \
      --outfmt 6 qseqid sseqid pident length evalue bitscore \
      --max-target-seqs 1 --evalue 1e-5 --id 70 --query-cover 60 \
      --block-size "$DIAMOND_BLOCK" --threads "$THREADS" --quiet 2>/dev/null || touch "$ARG_OUT"
  fi
fi

ARG_TOTAL=0
if [[ -f "$ARG_OUT" && -s "$ARG_OUT" ]]; then
  ARG_TOTAL=$(awk '{print $2}' "$ARG_OUT" | sort -u | wc -l)
fi
echo "  arg_total: $ARG_TOTAL"

##############################################################################
# 4. MEGAHIT assembly + IntegronFinder → mge_abundance
##############################################################################
MGE_ABUNDANCE=0
if [[ "$NO_ASSEMBLY" == false ]]; then
  ASSEMBLY="$ASSEMBLE_DIR/${SRR}/final.contigs.fa"
  if [[ ! -f "$ASSEMBLY" ]]; then
    echo "  MEGAHIT (${THREADS} threads, $(( MEGAHIT_MEM / 1073741824 ))GB RAM)..."
    $CONDA_BIN megahit \
      -1 "$FQ1" -2 "$FQ2" \
      -o "$ASSEMBLE_DIR/$SRR" \
      --num-cpu-threads "$THREADS" \
      --memory "$MEGAHIT_MEM" \
      --min-contig-len 500 2>/dev/null || true
  fi

  if [[ -f "$ASSEMBLY" ]]; then
    echo "  Assembly: $(grep -c '>' "$ASSEMBLY") contigs"
  fi

  MGE_DIR_SAMPLE="$MGE_DIR/$SRR"
  if [[ -f "$ASSEMBLY" && ! -d "$MGE_DIR_SAMPLE" ]]; then
    echo "  IntegronFinder (${THREADS} CPUs)..."
    mkdir -p "$MGE_DIR_SAMPLE"
    $CONDA_BIN integron_finder \
      --outdir "$MGE_DIR_SAMPLE" \
      --cpu "$THREADS" \
      --local-max \
      "$ASSEMBLY" 2>/dev/null || true
  fi

  if [[ -d "$MGE_DIR_SAMPLE" ]]; then
    MGE_ABUNDANCE=$(find "$MGE_DIR_SAMPLE" -name "*.integrons" 2>/dev/null | \
      xargs grep -h 'complete\|CALIN' 2>/dev/null | wc -l || echo 0)
    MGE_ABUNDANCE=${MGE_ABUNDANCE:-0}
  fi
fi
echo "  mge_abundance: $MGE_ABUNDANCE"

##############################################################################
# 5. Append result
##############################################################################
echo "${SRR},${ARG_TOTAL},${MGE_ABUNDANCE},${ENTROPY}" >> "$QUANT_CSV"
echo "=== DONE: $SRR → ARG=${ARG_TOTAL} MGE=${MGE_ABUNDANCE} H=${ENTROPY} ==="
