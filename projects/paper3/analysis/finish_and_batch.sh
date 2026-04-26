#!/usr/bin/env bash
# finish_and_batch.sh
# 1. Wait for any running DIAMOND R1 on SRR27827413, then run R2
# 2. Wait for MEGAHIT, then run IntegronFinder
# 3. Write SRR27827413 final row to quant_results.csv
# 4. Run batch_process_all.sh for remaining samples (fully parallelized)
#
# Run in WSL2: bash /mnt/c/.../finish_and_batch.sh

set -uo pipefail

OUTDIR="/data/paper3_quant"
ARG_DIR="$OUTDIR/arg"
KRAKEN_DIR="$OUTDIR/kraken"
BRACKEN_DIR="$OUTDIR/bracken"
MGE_DIR="$OUTDIR/mge"
ASSEMBLE_DIR="$OUTDIR/assembly"
QUANT_CSV="$OUTDIR/quant_results.csv"
ENV_BIN="/opt/miniforge/envs/biotools/bin"
KRAKEN_DB="/data/dbs/kraken2_8gb"
ARG_DB="/data/dbs/amrfinder_diamond/AMRProt.dmnd"
SRR="SRR27827413"

# ---------------------------------------------------------------------------
# Auto-detect hardware
# ---------------------------------------------------------------------------
THREADS=$(nproc)
TOTAL_RAM_MB=$(free -m | awk '/^Mem:/{print $2}')
DIAMOND_BLOCK=$(python3 -c "
ram_gb = $TOTAL_RAM_MB / 1024
block = max(0.5, round((ram_gb * 0.75 - 8) / 6, 1))
print(block)
")

echo "============================================================"
echo " Hardware: ${THREADS} cores | ${TOTAL_RAM_MB} MB RAM"
echo " DIAMOND block-size: ${DIAMOND_BLOCK}"
echo "============================================================"

[[ ! -f "$QUANT_CSV" ]] && echo "sample_id,arg_total,mge_abundance,entropy" > "$QUANT_CSV"

ENTROPY="6.291076"  # already computed from Bracken on previous machine

echo "=== Phase 1: Wait for DIAMOND R1 to finish ==="
while [[ ! -s "$ARG_DIR/${SRR}_r1.m8" ]]; do
  echo "  Waiting for R1... $(date '+%H:%M')"
  sleep 30
done
echo "  R1 done: $(wc -l < "$ARG_DIR/${SRR}_r1.m8") lines"

echo "=== Phase 2: Run DIAMOND R2 ==="
if [[ ! -s "$ARG_DIR/${SRR}_r2.m8" ]]; then
  "$ENV_BIN/diamond" blastx \
    --db "$ARG_DB" \
    --query "$OUTDIR/fastq/${SRR}_2.fastq" \
    --out "$ARG_DIR/${SRR}_r2.m8" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 1 --evalue 1e-5 --id 70 --query-cover 60 \
    --block-size "$DIAMOND_BLOCK" --threads "$THREADS" --quiet 2>/dev/null || true
  echo "  R2 done: $(wc -l < "$ARG_DIR/${SRR}_r2.m8" 2>/dev/null || echo 0) lines"
fi

cat "$ARG_DIR/${SRR}_r1.m8" "$ARG_DIR/${SRR}_r2.m8" > "$ARG_DIR/${SRR}.m8" 2>/dev/null || true
ARG_TOTAL=$(awk '{print $2}' "$ARG_DIR/${SRR}.m8" | sort -u | wc -l)
echo "  Combined ARG_TOTAL=$ARG_TOTAL"

echo "=== Phase 3: Wait for MEGAHIT to finish ==="
while [[ ! -f "$ASSEMBLE_DIR/${SRR}/final.contigs.fa" ]]; do
  echo "  Waiting for MEGAHIT... $(date '+%H:%M')"
  sleep 120
done
echo "  MEGAHIT done: $(grep -c '>' "$ASSEMBLE_DIR/${SRR}/final.contigs.fa") contigs"

echo "=== Phase 4: IntegronFinder (${THREADS} CPUs) ==="
MGE_DIR_SAMPLE="$MGE_DIR/$SRR"
if [[ ! -d "$MGE_DIR_SAMPLE" ]]; then
  mkdir -p "$MGE_DIR_SAMPLE"
  /opt/miniforge/bin/conda run -n biotools integron_finder \
    --outdir "$MGE_DIR_SAMPLE" \
    --cpu "$THREADS" \
    --local-max \
    "$ASSEMBLE_DIR/${SRR}/final.contigs.fa" 2>/dev/null || true
fi

MGE_ABUNDANCE=0
if [[ -d "$MGE_DIR_SAMPLE" ]]; then
  MGE_ABUNDANCE=$(find "$MGE_DIR_SAMPLE" -name "*.integrons" -exec \
    awk '$2 ~ /complete|CALIN/ {count++} END {print count+0}' {} \; 2>/dev/null | \
    awk '{s+=$1} END {print s+0}')
  MGE_ABUNDANCE=${MGE_ABUNDANCE:-0}
fi
echo "  MGE_ABUNDANCE=$MGE_ABUNDANCE"

echo "=== Phase 5: Write SRR27827413 to quant_results.csv ==="
grep -v "^${SRR}," "$QUANT_CSV" > /tmp/quant_tmp.csv 2>/dev/null || cp "$QUANT_CSV" /tmp/quant_tmp.csv
mv /tmp/quant_tmp.csv "$QUANT_CSV"
echo "${SRR},${ARG_TOTAL},${MGE_ABUNDANCE},${ENTROPY}" >> "$QUANT_CSV"
echo "  Written: ${SRR},${ARG_TOTAL},${MGE_ABUNDANCE},${ENTROPY}"

echo ""
echo "=== Phase 6: Run parallel batch pipeline for remaining samples ==="
SCRIPT_DIR="$(dirname "$0")"
bash "$SCRIPT_DIR/batch_process_all.sh"

echo ""
echo "=== All done! ==="
echo "Quant results: $QUANT_CSV"
