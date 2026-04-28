#!/usr/bin/env bash
# ============================================================
# Machine 3 Autonomous Pipeline
# Hardware: Intel i9-14900 | 20GB RAM (WSL: 18GB) | RTX 5070 Ti
# GPU is NOT used — all tools are CPU-only
#
# SETUP (Windows side, do ONCE before running this):
# 1. Enable WSL2, install Ubuntu 22.04 from Microsoft Store
# 2. Create C:\Users\<your-username>\.wslconfig with:
#        [wsl2]
#        memory=18GB
#        processors=20
#        swap=4GB
#    Then restart WSL: wsl --shutdown
# 3. Open Ubuntu terminal and run:
#        sudo bash /path/to/machine3_autonomous.sh
#    OR copy to WSL: cp machine3_autonomous.sh /data/
#    Then: sudo bash /data/machine3_autonomous.sh
#
# The script is fully restartable — skips completed samples.
# Expected runtime: ~24-33 hours (3 parallel jobs × 6 threads)
# Storage needed: ~100GB free (FASTQs auto-deleted after each sample)
# ============================================================
set -uo pipefail

LOGFILE=/data/machine3_pipeline.log
mkdir -p /data
exec > >(tee -a "$LOGFILE") 2>&1
echo "=== Machine 3 Pipeline START $(date) ==="

# ---- Paths ----
OUTDIR="/data/paper3_quant"
FASTQ_DIR="$OUTDIR/fastq"
ASSEMBLE_DIR="$OUTDIR/assembly"
MGE_DIR="$OUTDIR/mge"
KRAKEN2_DB="/data/dbs/kraken2_8gb"
DIAMOND_DB_DIR="/data/dbs/amrprot"
DIAMOND_DB="$DIAMOND_DB_DIR/AMRProt.dmnd"
QUANT_CSV="$OUTDIR/quant_results.csv"
ENV_BIN="/opt/miniforge/envs/biotools/bin"
CSV_LOCK="/tmp/quant_csv.lock"
KRAKEN_LOCK="/tmp/kraken2_run.lock"

# ---- URLs ----
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
KRAKEN2_DB_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz"
AMRPROT_URL="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt.fa"

# ---- Hardware config (fixed for i9-14900, 18GB WSL) ----
# Kraken2 serialized (flock): only 1 runs at a time, so 3rd parallel job stays
# in download/DIAMOND/MEGAHIT while waiting → safe to run 3 concurrent.
# MEGAHIT: at most 2 doing assembly simultaneously (3rd blocked on Kraken lock)
#   2 × 5GB = 10GB → safe within 18GB
# DIAMOND block 1.0 ≈ 6GB/job; 2 concurrent = 12GB (DIAMOND+MEGAHIT don't
#   overlap within a sample, so peak is ~12GB → safe)
PARALLEL_JOBS=3
THREADS_PER_JOB=6          # 3 × 6 = 18 threads (20 WSL cores, 2 spare for OS)
MEM_PER_JOB=5000000000     # 5GB per MEGAHIT job (max 2 concurrent = 10GB)
DIAMOND_BLOCK="1.0"        # ~6GB RAM per DIAMOND job; 2× faster than 0.5

ALL_SAMPLES=(
  SRR27827413 SRR27827408 SRR27827406 SRR27827404
  SRR27827412 SRR27827407 SRR27827405 SRR27827403
  SRR11803495 SRR11803496 SRR11803497 SRR11803498
  SRR11803499 SRR11803500 SRR11803501 SRR11803502
  SRR11803503 SRR11803504 SRR11803505
)

# ============================================================
# PHASE 0: Install Miniforge + biotools conda environment
# ============================================================
echo ""; echo "=== PHASE 0: Environment setup ==="

if [[ ! -f /opt/miniforge/bin/conda ]]; then
  echo "Installing Miniforge to /opt/miniforge..."
  curl -fsSL "$MINIFORGE_URL" -o /tmp/miniforge.sh
  bash /tmp/miniforge.sh -b -p /opt/miniforge
  rm -f /tmp/miniforge.sh
  echo "Miniforge installed."
fi

source /opt/miniforge/etc/profile.d/conda.sh

if ! /opt/miniforge/bin/conda env list | grep -q biotools; then
  echo "Creating biotools environment (installs kraken2, bracken, diamond, megahit, integron_finder, sra-tools)..."
  /opt/miniforge/bin/conda create -n biotools -y \
    -c conda-forge -c bioconda \
    kraken2 bracken diamond megahit integron_finder sra-tools
  echo "biotools environment ready."
fi

mkdir -p "$FASTQ_DIR" "$ASSEMBLE_DIR" "$MGE_DIR" \
  "$KRAKEN2_DB" "$DIAMOND_DB_DIR" \
  "$OUTDIR/kraken" "$OUTDIR/bracken" "$OUTDIR/diamond"
touch "$CSV_LOCK" "$KRAKEN_LOCK"

echo "Cores: $(nproc)  RAM: $(free -h | awk '/^Mem:/{print $2}')"
echo "Config: ${PARALLEL_JOBS} parallel jobs | ${THREADS_PER_JOB} threads/job | $(( MEM_PER_JOB / 1073741824 ))GB MEGAHIT/job | DIAMOND block=${DIAMOND_BLOCK}"

# ============================================================
# PHASE 1: Download databases
# ============================================================
echo ""; echo "=== PHASE 1: Databases ==="

if [[ ! -f "$KRAKEN2_DB/hash.k2d" ]]; then
  echo "Downloading Kraken2 8GB standard DB (~8GB, ~20-40 min)..."
  curl -fSL --retry 3 "$KRAKEN2_DB_URL" | tar -xz -C "$KRAKEN2_DB"
  echo "Kraken2 DB ready."
else
  echo "Kraken2 DB already present."
fi

if [[ ! -f "$DIAMOND_DB" ]]; then
  echo "Downloading AMRProt.fa from NCBI..."
  curl -fSL --retry 3 "$AMRPROT_URL" -o "$DIAMOND_DB_DIR/AMRProt.fa"
  n_proteins=$(grep -c '>' "$DIAMOND_DB_DIR/AMRProt.fa" 2>/dev/null || echo 0)
  echo "AMRProt.fa: ${n_proteins} proteins"
  if [[ "$n_proteins" -lt 1000 ]]; then
    echo "ERROR: AMRProt.fa download failed (only ${n_proteins} sequences). Aborting."
    exit 1
  fi
  echo "Building DIAMOND database..."
  "$ENV_BIN/diamond" makedb \
    --in "$DIAMOND_DB_DIR/AMRProt.fa" \
    --db "$DIAMOND_DB_DIR/AMRProt" \
    --threads "$THREADS_PER_JOB"
  rm -f "$DIAMOND_DB_DIR/AMRProt.fa"
  echo "DIAMOND DB ready: $DIAMOND_DB"
else
  echo "DIAMOND DB already present."
fi

# ============================================================
# PHASE 2: Initialize output CSV
# ============================================================
echo ""; echo "=== PHASE 2: Initialize CSV ==="

if [[ ! -f "$QUANT_CSV" ]]; then
  echo "sample,entropy,arg_total,mge_abundance" > "$QUANT_CSV"
  for srr in "${ALL_SAMPLES[@]}"; do
    echo "$srr,0,0,0" >> "$QUANT_CSV"
  done
  echo "Created $QUANT_CSV with ${#ALL_SAMPLES[@]} samples."
else
  echo "$QUANT_CSV exists — will skip completed samples."
fi

# ============================================================
# PHASE 3: Per-sample processing
# Each sample: download → Kraken2+Bracken → DIAMOND → MEGAHIT → delete FASTQs → IntegronFinder
# ============================================================

process_sample() {
  local srr="$1"
  local fq1="$FASTQ_DIR/${srr}_1.fastq"
  local fq2="$FASTQ_DIR/${srr}_2.fastq"
  local fq_se="$FASTQ_DIR/${srr}.fastq"
  local asm_dir="$ASSEMBLE_DIR/$srr"
  local assembly="$asm_dir/final.contigs.fa"
  local mge_dir="$MGE_DIR/$srr"
  local kraken_report="$OUTDIR/kraken/${srr}.txt"
  local bracken_out="$OUTDIR/bracken/${srr}.bracken"
  local diamond_out="$OUTDIR/diamond/${srr}.m8"

  # Read current CSV values to determine what to skip
  local cur_entropy cur_arg cur_mge
  cur_entropy=$(awk -F',' -v s="$srr" '$1==s{print $2}' "$QUANT_CSV")
  cur_arg=$(awk -F',' -v s="$srr" '$1==s{print $3}' "$QUANT_CSV")
  cur_mge=$(awk -F',' -v s="$srr" '$1==s{print $4}' "$QUANT_CSV")

  local need_entropy=true need_arg=true need_mge=true
  [[ -n "$cur_entropy" && "$cur_entropy" != "0" ]] && need_entropy=false
  [[ -n "$cur_arg" && "$cur_arg" != "0" ]] && need_arg=false
  # mge: skip only if non-zero AND no pending assembly
  [[ -n "$cur_mge" && "$cur_mge" != "0" && ! -d "$asm_dir" ]] && need_mge=false

  if ! $need_entropy && ! $need_arg && ! $need_mge; then
    echo "[$srr] All values computed — skipping."
    return 0
  fi

  echo "[$srr] Starting (need: entropy=$need_entropy arg=$need_arg mge=$need_mge)"

  # If Bracken output already exists, compute entropy from it without re-downloading.
  if $need_entropy && [[ -f "$bracken_out" ]]; then
    local entropy=0
    entropy=$(awk -F'\t' 'NR>1 && $7>0 {p=$7; s+=p*log(p)/log(2)} END{printf "%.6f", -s+0}' \
      "$bracken_out" 2>/dev/null || echo 0)
    [[ -z "$entropy" ]] && entropy=0
    (flock -x 9
      tmp="/tmp/csv_e_${srr}_$$.csv"
      awk -v s="$srr" -v v="$entropy" -F',' \
        'BEGIN{OFS=","} $1==s{$2=v} 1' "$QUANT_CSV" > "$tmp" && mv "$tmp" "$QUANT_CSV"
    ) 9>"$CSV_LOCK"
    echo "[$srr] entropy=${entropy} (reused existing bracken output)"
    need_entropy=false
  fi

  # Download FASTQs only if still needed for Kraken2 or DIAMOND.
  if ($need_entropy || $need_arg) && [[ ! -f "$fq1" && ! -f "$fq_se" ]]; then
    echo "[$srr] Downloading FASTQs..."
    # prefetch downloads .sra to $FASTQ_DIR/$srr/
    "$ENV_BIN/prefetch" --max-size 50G "$srr" -O "$FASTQ_DIR/" 2>/dev/null || true

    # fasterq-dump: try .sra path, fall back to accession
    local sra_file="$FASTQ_DIR/$srr/$srr.sra"
    if [[ -f "$sra_file" ]]; then
      "$ENV_BIN/fasterq-dump" --split-files --outdir "$FASTQ_DIR" \
        --threads "$THREADS_PER_JOB" "$sra_file" 2>/dev/null || \
      "$ENV_BIN/fasterq-dump" --outdir "$FASTQ_DIR" \
        --threads "$THREADS_PER_JOB" "$sra_file" 2>/dev/null || true
    else
      "$ENV_BIN/fasterq-dump" --split-files --outdir "$FASTQ_DIR" \
        --threads "$THREADS_PER_JOB" "$srr" 2>/dev/null || \
      "$ENV_BIN/fasterq-dump" --outdir "$FASTQ_DIR" \
        --threads "$THREADS_PER_JOB" "$srr" 2>/dev/null || true
    fi
    rm -rf "$FASTQ_DIR/$srr/"  # clean up prefetch dir

    echo "[$srr] Download done. Disk free: $(df -h / | awk 'NR==2{print $4}')"
  fi

  local has_pe=false has_se=false
  [[ -f "$fq1" && -f "$fq2" ]] && has_pe=true
  [[ -f "$fq_se" ]] && has_se=true

  if ($need_entropy || $need_arg) && ! $has_pe && ! $has_se; then
    echo "[$srr] WARNING: No FASTQs found after download. Skipping sample."
    return 0
  fi

  # ---- Kraken2 → entropy (serialized with flock: only 1 Kraken2 at a time to prevent OOM) ----
  if $need_entropy; then
    (flock -x 8
      echo "[$srr] Kraken2 (serialized, $(( $(df "$KRAKEN2_DB" | awk 'NR==2{print $4}') / 1024 / 1024 ))GB free)..."
      if $has_pe; then
        "$ENV_BIN/kraken2" --db "$KRAKEN2_DB" \
          --paired "$fq1" "$fq2" \
          --output /dev/null \
          --report "$kraken_report" \
          --threads "$THREADS_PER_JOB" 2>/dev/null || true
      else
        "$ENV_BIN/kraken2" --db "$KRAKEN2_DB" \
          "$fq_se" \
          --output /dev/null \
          --report "$kraken_report" \
          --threads "$THREADS_PER_JOB" 2>/dev/null || true
      fi

      if [[ -f "$kraken_report" ]]; then
        # Use conda run so 'python' resolves correctly inside the bracken wrapper script.
        /opt/miniforge/bin/conda run -n biotools bracken \
          -d "$KRAKEN2_DB" \
          -i "$kraken_report" \
          -o "$bracken_out" \
          -r 150 -l S 2>/dev/null || true
      fi
    ) 8>"$KRAKEN_LOCK"

    local entropy=0
    if [[ -f "$bracken_out" ]]; then
      entropy=$(awk -F'\t' 'NR>1 && $7>0 {p=$7; s+=p*log(p)/log(2)} END{printf "%.6f", -s+0}' \
        "$bracken_out" 2>/dev/null || echo 0)
      [[ -z "$entropy" ]] && entropy=0
    fi

    (flock -x 9
      tmp="/tmp/csv_e_${srr}_$$.csv"
      awk -v s="$srr" -v v="$entropy" -F',' \
        'BEGIN{OFS=","} $1==s{$2=v} 1' "$QUANT_CSV" > "$tmp" && mv "$tmp" "$QUANT_CSV"
    ) 9>"$CSV_LOCK"
    echo "[$srr] entropy=${entropy}"
  fi

  # ---- DIAMOND blastx → arg_total ----
  if $need_arg; then
    echo "[$srr] DIAMOND blastx (block-size=${DIAMOND_BLOCK})..."
    if $has_pe; then
      cat "$fq1" "$fq2" | \
        "$ENV_BIN/diamond" blastx \
          --db "$DIAMOND_DB" \
          --query - \
          --out "$diamond_out" \
          --outfmt 6 \
          --threads "$THREADS_PER_JOB" \
          --block-size "$DIAMOND_BLOCK" \
          --evalue 1e-5 \
          --id 80 \
          --query-cover 80 \
          --max-target-seqs 1 \
          --quiet 2>/dev/null || true
    else
      "$ENV_BIN/diamond" blastx \
        --db "$DIAMOND_DB" \
        --query "$fq_se" \
        --out "$diamond_out" \
        --outfmt 6 \
        --threads "$THREADS_PER_JOB" \
        --block-size "$DIAMOND_BLOCK" \
        --evalue 1e-5 \
        --id 80 \
        --query-cover 80 \
        --max-target-seqs 1 \
        --quiet 2>/dev/null || true
    fi

    local arg_total=0
    [[ -f "$diamond_out" ]] && arg_total=$(wc -l < "$diamond_out" 2>/dev/null || echo 0)
    arg_total=${arg_total:-0}

    (flock -x 9
      tmp="/tmp/csv_a_${srr}_$$.csv"
      awk -v s="$srr" -v v="$arg_total" -F',' \
        'BEGIN{OFS=","} $1==s{$3=v} 1' "$QUANT_CSV" > "$tmp" && mv "$tmp" "$QUANT_CSV"
    ) 9>"$CSV_LOCK"
    echo "[$srr] arg_total=${arg_total}"
  fi

  # ---- MEGAHIT assembly ----
  if $need_mge && [[ ! -f "$assembly" ]]; then
    echo "[$srr] MEGAHIT (${THREADS_PER_JOB} threads, $(( MEM_PER_JOB / 1073741824 ))GB)..."
    rm -rf "$asm_dir"
    if $has_pe; then
      "$ENV_BIN/megahit" -1 "$fq1" -2 "$fq2" \
        -o "$asm_dir" \
        --num-cpu-threads "$THREADS_PER_JOB" \
        --memory "$MEM_PER_JOB" \
        --min-contig-len 500 2>/dev/null || true
    else
      "$ENV_BIN/megahit" -r "$fq_se" \
        -o "$asm_dir" \
        --num-cpu-threads "$THREADS_PER_JOB" \
        --memory "$MEM_PER_JOB" \
        --min-contig-len 500 2>/dev/null || true
    fi
  fi

  # ---- Delete FASTQs now that assembly is done (or failed) ----
  echo "[$srr] Deleting FASTQs to free storage..."
  rm -f "$fq1" "$fq2" "$fq_se"
  echo "[$srr] Disk free after cleanup: $(df -h / | awk 'NR==2{print $4}')"

  # ---- IntegronFinder → mge_abundance ----
  if $need_mge; then
    if [[ -f "$assembly" ]]; then
      local n_contigs
      n_contigs=$(grep -c '>' "$assembly" 2>/dev/null || echo 0)

      # IntegronFinder runs HMMER on every contig regardless of length.
      # Complete integrons require ≥4kb; filtering here cuts runtime ~10-50x.
      local filtered_assembly="$asm_dir/contigs_4kb.fa"
      awk '/^>/{h=$0; next} length($0)>=4000{print h; print}' \
        "$assembly" > "$filtered_assembly"
      local n_filtered
      n_filtered=$(grep -c '>' "$filtered_assembly" 2>/dev/null || echo 0)
      echo "[$srr] Assembly: ${n_contigs} total contigs, ${n_filtered} ≥4kb — running IntegronFinder..."

      if [[ "$n_filtered" -eq 0 ]]; then
        echo "[$srr] No contigs ≥4kb — mge_abundance=0"
      elif ! find "$mge_dir" -name "*.integrons" 2>/dev/null | grep -q .; then
        rm -rf "$mge_dir"
        mkdir -p "$mge_dir"
        /opt/miniforge/bin/conda run -n biotools integron_finder \
          --outdir "$mge_dir" \
          --cpu "$THREADS_PER_JOB" \
          --local-max \
          "$filtered_assembly" 2>/dev/null || true
      fi

      rm -rf "$asm_dir/intermediate_contigs" "$asm_dir/tmp"
    else
      echo "[$srr] MEGAHIT produced no assembly — mge_abundance=0"
    fi
  fi

  local mge_abundance=0
  mge_abundance=$(find "$mge_dir" -name "*.integrons" 2>/dev/null | \
    xargs grep -h 'complete\|CALIN' 2>/dev/null | wc -l || echo 0)
  mge_abundance=${mge_abundance:-0}

  (flock -x 9
    tmp="/tmp/csv_m_${srr}_$$.csv"
    awk -v s="$srr" -v v="$mge_abundance" -F',' \
      'BEGIN{OFS=","} $1==s{$4=v} 1' "$QUANT_CSV" > "$tmp" && mv "$tmp" "$QUANT_CSV"
  ) 9>"$CSV_LOCK"

  local final_entropy final_arg
  final_entropy=$(awk -F',' -v s="$srr" '$1==s{print $2}' "$QUANT_CSV")
  final_arg=$(awk -F',' -v s="$srr" '$1==s{print $3}' "$QUANT_CSV")
  echo "[$srr] DONE — entropy=${final_entropy} | arg_total=${final_arg} | mge=${mge_abundance}"
  echo "[$srr] Disk free: $(df -h / | awk 'NR==2{print $4}')"
}

# ============================================================
# Launch parallel job pool
# ============================================================
echo ""; echo "=== PHASE 3: Processing ${#ALL_SAMPLES[@]} samples (${PARALLEL_JOBS} parallel, ${THREADS_PER_JOB} threads each) ==="

running=0
for srr in "${ALL_SAMPLES[@]}"; do
  process_sample "$srr" &
  running=$((running + 1))
  if [[ $running -ge $PARALLEL_JOBS ]]; then
    wait -n 2>/dev/null || true
    running=$((running - 1))
  fi
done
wait || true

# ============================================================
# Done
# ============================================================
echo ""
echo "=== Pipeline COMPLETE $(date) ==="
echo ""
echo "Final quant_results.csv:"
cat "$QUANT_CSV"
echo ""
echo "Full log: $LOGFILE"
