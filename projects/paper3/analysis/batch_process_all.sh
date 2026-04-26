#!/usr/bin/env bash
# batch_process_all.sh — Paper 3 shotgun pipeline, FULLY PARALLELIZED
# Auto-detects all CPU cores and RAM; runs multiple samples concurrently.
#
# Usage (in WSL2 Ubuntu):
#   bash /mnt/c/path/to/batch_process_all.sh
#
# Prereqs (biotools conda env):
#   /opt/miniforge/envs/biotools/bin/{kraken2,bracken,diamond,prefetch,fasterq-dump}
#   /data/dbs/kraken2_8gb           — Kraken2 standard 8GB DB
#   /data/dbs/amrfinder_diamond/AMRProt.dmnd — DIAMOND ARG database

set -uo pipefail

OUTDIR="/data/paper3_quant"
FASTQ_DIR="$OUTDIR/fastq"
KRAKEN_DIR="$OUTDIR/kraken"
BRACKEN_DIR="$OUTDIR/bracken"
ARG_DIR="$OUTDIR/arg"
QUANT_CSV="$OUTDIR/quant_results.csv"
KRAKEN_DB="/data/dbs/kraken2_8gb"
ARG_DB="/data/dbs/amrfinder_diamond/AMRProt.dmnd"
ENV_BIN="/opt/miniforge/envs/biotools/bin"
LOCK_FILE="/tmp/paper3_quant_csv.lock"

# ---------------------------------------------------------------------------
# Auto-detect hardware and compute optimal parallelism
# ---------------------------------------------------------------------------
TOTAL_CORES=$(nproc)
TOTAL_RAM_MB=$(free -m | awk '/^Mem:/{print $2}')

# Each sample needs: ~8GB for Kraken2 DB + ~6GB per DIAMOND block unit + 2GB headroom
# We target 80% RAM utilization; split remainder across parallel jobs.
read -r PARALLEL_JOBS THREADS_PER_JOB DIAMOND_BLOCK < <(python3 - "$TOTAL_CORES" "$TOTAL_RAM_MB" <<'PYEOF'
import sys, math
cores = int(sys.argv[1])
ram_gb = int(sys.argv[2]) / 1024

# Kraken2 loads the 8GB DB into shared memory (one copy for all parallel jobs)
usable_ram = ram_gb * 0.80 - 8  # subtract DB overhead, keep 20% for OS
usable_ram = max(usable_ram, 4)

# DIAMOND: ~6 GB RAM per block-size 1.0
# Max parallel jobs limited by RAM and by cores (min 4 threads/job for diamond efficiency)
max_by_ram  = max(1, int(usable_ram / 8))       # ~8 GB per sample slot
max_by_cores= max(1, cores // 4)                 # at least 4 cores per job
jobs = min(max_by_ram, max_by_cores, 8)          # cap at 8 parallel jobs

threads = max(1, cores // jobs)
block   = max(0.5, round((usable_ram / 6) / jobs, 1))

print(jobs, threads, block)
PYEOF
)

echo "============================================================"
echo " Hardware : ${TOTAL_CORES} cores | ${TOTAL_RAM_MB} MB RAM"
echo " Strategy : ${PARALLEL_JOBS} parallel samples"
echo "            ${THREADS_PER_JOB} threads/sample"
echo "            DIAMOND block-size ${DIAMOND_BLOCK} (~$(python3 -c "print(round(float('$DIAMOND_BLOCK')*6,0))" 2>/dev/null || echo '?') GB/sample)"
echo "============================================================"

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

# ---------------------------------------------------------------------------
# Shannon entropy from Bracken report
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Process one sample (runs as a background job)
# ---------------------------------------------------------------------------
process_sample() {
  local srr="$1"
  local threads="$2"
  local block="$3"

  if grep -q "^${srr}," "$QUANT_CSV" 2>/dev/null; then
    echo "[${srr}] already done — skipping"
    return 0
  fi

  echo "[${srr}] START (${threads} threads | block=${block})"

  # ---- 1. Download --------------------------------------------------------
  local fq1="$FASTQ_DIR/${srr}_1.fastq"
  local fq2="$FASTQ_DIR/${srr}_2.fastq"
  local fq_se="$FASTQ_DIR/${srr}.fastq"

  if [[ ! -f "$fq1" && ! -f "${fq1}.gz" && ! -f "$fq_se" ]]; then
    echo "[${srr}] Downloading..."
    "$ENV_BIN/prefetch" --max-size 60G "$srr" -O "$FASTQ_DIR" 2>/dev/null || true
    "$ENV_BIN/fasterq-dump" \
      --split-files --threads "$threads" --bufsize 1000MB --curcache 10000MB \
      --outdir "$FASTQ_DIR" "$FASTQ_DIR/$srr" 2>/dev/null || \
    "$ENV_BIN/fasterq-dump" \
      --split-files --threads "$threads" --bufsize 1000MB --curcache 10000MB \
      --outdir "$FASTQ_DIR" "$srr" 2>/dev/null || true
  fi

  local is_paired=false
  [[ -f "$fq1" && -f "$fq2" ]] && is_paired=true

  if [[ "$is_paired" == false && ! -f "$fq_se" ]]; then
    echo "[${srr}] ERROR: no FASTQ found"
    (flock -x 9; echo "${srr},0,0,0" >> "$QUANT_CSV") 9>"$LOCK_FILE"
    return 1
  fi

  # ---- 2. Kraken2 ---------------------------------------------------------
  local kraken_report="$KRAKEN_DIR/${srr}.report"
  if [[ ! -f "$kraken_report" ]]; then
    echo "[${srr}] Kraken2..."
    if [[ "$is_paired" == true ]]; then
      "$ENV_BIN/kraken2" --db "$KRAKEN_DB" --paired \
        --threads "$threads" --report "$kraken_report" \
        --output /dev/null "$fq1" "$fq2" 2>/dev/null || touch "$kraken_report"
    else
      "$ENV_BIN/kraken2" --db "$KRAKEN_DB" \
        --threads "$threads" --report "$kraken_report" \
        --output /dev/null "$fq_se" 2>/dev/null || touch "$kraken_report"
    fi
  fi

  # ---- 3. Bracken ---------------------------------------------------------
  local bracken_report="$BRACKEN_DIR/${srr}.bracken"
  if [[ ! -f "$bracken_report" && -s "$kraken_report" ]]; then
    echo "[${srr}] Bracken..."
    /opt/miniforge/bin/conda run -n biotools bracken \
      -d "$KRAKEN_DB" -i "$kraken_report" -o "$bracken_report" \
      -r 150 -l S -t 10 2>/dev/null || touch "$bracken_report"
  fi

  local entropy=0
  if [[ -f "$bracken_report" && -s "$bracken_report" ]]; then
    entropy=$(compute_entropy "$bracken_report")
  fi
  echo "[${srr}] Entropy: ${entropy}"

  # ---- 4. DIAMOND ARG (combined paired reads in one pass) -----------------
  local arg_out="$ARG_DIR/${srr}.m8"
  if [[ ! -f "$arg_out" || ! -s "$arg_out" ]]; then
    echo "[${srr}] DIAMOND (block=${block}, ${threads} threads)..."
    if [[ "$is_paired" == true ]]; then
      cat "$fq1" "$fq2" | \
        "$ENV_BIN/diamond" blastx \
          --db "$ARG_DB" --query - --out "$arg_out" \
          --outfmt 6 qseqid sseqid pident length evalue bitscore \
          --max-target-seqs 1 --evalue 1e-5 --id 70 --query-cover 60 \
          --block-size "$block" --threads "$threads" --quiet 2>/dev/null || touch "$arg_out"
    else
      "$ENV_BIN/diamond" blastx \
        --db "$ARG_DB" --query "$fq_se" --out "$arg_out" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 1 --evalue 1e-5 --id 70 --query-cover 60 \
        --block-size "$block" --threads "$threads" --quiet 2>/dev/null || touch "$arg_out"
    fi
  fi

  local arg_total=0
  if [[ -f "$arg_out" && -s "$arg_out" ]]; then
    arg_total=$(awk '{print $2}' "$arg_out" | sort -u | wc -l)
  fi
  echo "[${srr}] arg_total: ${arg_total}"

  # ---- 5. Write result (file-locked to avoid concurrent write corruption) -
  (flock -x 9
    # Remove any partial/stale row then append final result
    tmp="/tmp/quant_tmp_${srr}_$$.csv"
    grep -v "^${srr}," "$QUANT_CSV" > "$tmp" 2>/dev/null || cp "$QUANT_CSV" "$tmp"
    mv "$tmp" "$QUANT_CSV"
    echo "${srr},${arg_total},0,${entropy}" >> "$QUANT_CSV"
  ) 9>"$LOCK_FILE"

  echo "[${srr}] DONE — ARG=${arg_total} H=${entropy}"
}

# ---------------------------------------------------------------------------
# Parallel job pool (bash 4.3+ wait -n)
# ---------------------------------------------------------------------------
echo ""
echo "=== Batch pipeline START ($(date)) ==="

running=0
for srr in "${SHOTGUN_SAMPLES[@]}"; do
  process_sample "$srr" "$THREADS_PER_JOB" "$DIAMOND_BLOCK" &
  running=$((running + 1))
  if [[ $running -ge $PARALLEL_JOBS ]]; then
    wait -n 2>/dev/null || wait
    running=$((running - 1))
  fi
done
wait

echo ""
echo "=== Batch pipeline COMPLETE ($(date)) ==="
echo "Results: $QUANT_CSV"
