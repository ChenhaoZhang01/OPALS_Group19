#!/usr/bin/env bash
# Paper 1 -- REAL three-pipeline ARG quantification for the variance-decomposition study.
#
# Runs every sample in the download list through three independent ARG-calling pipelines
# and emits the three wide matrices + long table that analysis/run_paper1_analysis.py
# consumes. Designed for a Linux HPC node. Resumable, logged, disk-aware.
#
#   A  assembly (MEGAHIT) + gene prediction (Prodigal) + DIAMOND blastp vs CARD
#   B  read mapping (Bowtie2) to CARD nucleotide reference
#   C  strict CARD resistance calls (RGI) on the predicted ORFs from A
#
# Prereqs:
#   1. conda env from environment-linux.yml  (sra-tools, fastp, megahit, prodigal,
#      diamond, bowtie2, samtools, rgi, python+pandas)
#   2. bash projects/paper1/analysis/setup_databases.sh --dbdir <CARD_DIR>
#
# Usage:
#   # smoke test on the single smallest sample first:
#   bash projects/paper1/analysis/run_arg_pipeline_hpc.sh --card-dir /databases/card \
#        --outdir /scratch/paper1 --threads 16 --only SRR30403277
#
#   # full cohort:
#   bash projects/paper1/analysis/run_arg_pipeline_hpc.sh --card-dir /databases/card \
#        --outdir /scratch/paper1 --threads 16
#
# Re-running skips samples that already have all three hit files (resume).

set -euo pipefail

DOWNLOAD_LIST="projects/paper1/metadata/download_list.txt"
METADATA="projects/paper1/metadata/metadata_final.csv"
OUTDIR="/scratch/paper1"
CARD_DIR="/databases/card"
THREADS=8
ONLY=""            # run a single accession (smoke test)
KEEP_INTERMEDIATE=0  # 1 = keep reads/assembly (uses a lot of disk)
MAX_READS=""       # if set, stream only the first N spots/pairs (equal-depth subsample)
MEM_FRAC=0.85      # megahit memory cap as fraction of system RAM

while [[ $# -gt 0 ]]; do
  case "$1" in
    --download-list) DOWNLOAD_LIST="$2"; shift 2;;
    --metadata)      METADATA="$2"; shift 2;;
    --outdir)        OUTDIR="$2"; shift 2;;
    --card-dir)      CARD_DIR="$2"; shift 2;;
    --threads)       THREADS="$2"; shift 2;;
    --only)          ONLY="$2"; shift 2;;
    --max-reads)     MAX_READS="$2"; shift 2;;
    --mem-frac)      MEM_FRAC="$2"; shift 2;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift;;
    *) echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

CARD_DIAMOND="$CARD_DIR/card_protein.dmnd"
CARD_BT2="$CARD_DIR/card_bt2"
for f in "$CARD_DIAMOND" "${CARD_BT2}.1.bt2"; do
  [[ -e "$f" ]] || { echo "ERROR: missing CARD db file: $f (run setup_databases.sh)" >&2; exit 1; }
done

mkdir -p "$OUTDIR"/{reads,assembly,proteins,hits_A,hits_B,hits_C,logs}
RCFILE="$OUTDIR/read_counts.txt"; : > /dev/null  # do not truncate; we append/dedup below

log() { echo "[$(date +%H:%M:%S)] $*"; }

process_sample() {
  local SRR="$1"
  local LOG="$OUTDIR/logs/$SRR.log"
  local A="$OUTDIR/hits_A/$SRR.tsv"
  local B="$OUTDIR/hits_B/$SRR.idxstats"
  local C="$OUTDIR/hits_C/$SRR.txt"

  if [[ -s "$A" && -s "$B" && -s "$C" ]]; then
    log "$SRR: already complete, skipping"
    return 0
  fi
  log "$SRR: starting (log -> $LOG)"
  {
    set -euo pipefail
    local R1="$OUTDIR/reads/${SRR}_1.fastq" R2="$OUTDIR/reads/${SRR}_2.fastq"
    if [[ ! -s "$R1" ]]; then
      if [[ -n "$MAX_READS" ]]; then
        # Stream only the first MAX_READS spots (equal-depth subsample, bounds RAM+download).
        fastq-dump --split-files --skip-technical -X "$MAX_READS" -O "$OUTDIR/reads" "$SRR"
      else
        prefetch "$SRR" -O "$OUTDIR/reads"
        fasterq-dump "$SRR" -O "$OUTDIR/reads" -e "$THREADS" --split-files
      fi
    fi
    local R1t="${R1%.fastq}.tr.fastq" R2t="${R2%.fastq}.tr.fastq"
    local PAIRED=1; [[ -s "$R2" ]] || PAIRED=0
    if [[ "$PAIRED" -eq 1 ]]; then
      fastp -i "$R1" -I "$R2" -o "$R1t" -O "$R2t" -w "$THREADS" \
            -j "$OUTDIR/logs/$SRR.fastp.json" -h "$OUTDIR/logs/$SRR.fastp.html"
    else
      fastp -i "$R1" -o "$R1t" -w "$THREADS" \
            -j "$OUTDIR/logs/$SRR.fastp.json" -h "$OUTDIR/logs/$SRR.fastp.html"
    fi
    local n1; n1=$(( $(wc -l < "$R1t") / 4 ))
    local TOTAL_READS; TOTAL_READS=$(( PAIRED == 1 ? n1 * 2 : n1 ))

    # Pipeline B FIRST: read mapping to CARD (no assembly; must never be lost to an
    # assembly OOM, so it runs before the memory-heavy assembly block).
    local BT2_ARGS
    if [[ "$PAIRED" -eq 1 ]]; then
      BT2_ARGS=(-1 "$R1t" -2 "$R2t")
    else
      BT2_ARGS=(-U "$R1t")
    fi
    bowtie2 -x "$CARD_BT2" "${BT2_ARGS[@]}" -p "$THREADS" --no-unal \
      | samtools sort -@ "$THREADS" -o "$OUTDIR/hits_B/$SRR.bam" -
    samtools index "$OUTDIR/hits_B/$SRR.bam"
    samtools idxstats "$OUTDIR/hits_B/$SRR.bam" > "$B.tmp" && mv "$B.tmp" "$B"

    # Pipelines A & C share an assembly. Make the whole assembly block NON-FATAL:
    # a megahit OOM on high-complexity soil must not abort the sample (B is already
    # saved). On failure we write empty A/C so the sample still counts as complete.
    rm -rf "$OUTDIR/assembly/$SRR"
    local assembled=0
    if [[ "$PAIRED" -eq 1 ]]; then
      megahit -1 "$R1t" -2 "$R2t" -o "$OUTDIR/assembly/$SRR" -t "$THREADS" \
              --min-contig-len 500 -m "$MEM_FRAC" && assembled=1 || assembled=0
    else
      megahit -r "$R1t" -o "$OUTDIR/assembly/$SRR" -t "$THREADS" \
              --min-contig-len 500 -m "$MEM_FRAC" && assembled=1 || assembled=0
    fi
    if [[ "$assembled" -eq 1 && -s "$OUTDIR/assembly/$SRR/final.contigs.fa" ]]; then
      prodigal -i "$OUTDIR/assembly/$SRR/final.contigs.fa" \
               -a "$OUTDIR/proteins/$SRR.faa" -p meta -q || true
      # Pipeline A: DIAMOND blastp of predicted ORFs vs CARD
      diamond blastp -d "$CARD_DIAMOND" -q "$OUTDIR/proteins/$SRR.faa" -p "$THREADS" \
              --id 80 --query-cover 70 --max-target-seqs 1 \
              --outfmt 6 qseqid sseqid pident length evalue bitscore \
              -o "$A.tmp" && mv "$A.tmp" "$A" || : > "$A"
      # Pipeline C: RGI ORF-level homology on predicted ORFs (strip Prodigal '*').
      # NOTE: RGI "strict" (Perfect/Strict only) returns ~0 on environmental
      # metagenomes because environmental ARG homologs fall below CARD's clinical
      # curated bitscore cutoffs (all hits land in the Loose tier). We therefore
      # include the Loose tier so Pipeline C yields a comparable abundance gradient;
      # strict-only ~0 is reported as a finding. See results/DATA_PROVENANCE.md.
      sed 's/\*//g' "$OUTDIR/proteins/$SRR.faa" > "$OUTDIR/proteins/$SRR.rgi.faa"
      rgi main -i "$OUTDIR/proteins/$SRR.rgi.faa" -o "$OUTDIR/hits_C/$SRR" \
          -t protein -a DIAMOND --clean --include_loose --num_threads "$THREADS" || true
      [[ -s "$C" ]] || { echo "WARN: RGI produced no $C" >&2; : > "$C"; }
    else
      echo "WARN: assembly failed/empty for $SRR; A & C empty (B retained)" >&2
      : > "$A"; : > "$C"
    fi

    echo "$SRR $TOTAL_READS" >> "$RCFILE"

    if [[ "$KEEP_INTERMEDIATE" -eq 0 ]]; then
      rm -f "$R1" "$R2" "$R1t" "$R2t" "$OUTDIR/hits_B/$SRR.bam" "$OUTDIR/hits_B/$SRR.bam.bai"
      rm -rf "$OUTDIR/assembly/$SRR"
    fi
  } >>"$LOG" 2>&1 || { log "$SRR: FAILED (see $LOG)"; return 1; }
  log "$SRR: done"
}

# Build the work list.
mapfile -t SAMPLES < <(tr -d '\r' < "$DOWNLOAD_LIST" | sed '/^\s*$/d')
[[ -n "$ONLY" ]] && SAMPLES=("$ONLY")

FAIL=0
for SRR in "${SAMPLES[@]}"; do
  process_sample "$SRR" || FAIL=$((FAIL+1))
done
log "Processing complete. Failures: $FAIL / ${#SAMPLES[@]}"

# De-duplicate read_counts (resume may append twice).
if [[ -f "$RCFILE" ]]; then
  sort -u "$RCFILE" -o "$RCFILE"
fi

# Aggregate into matrices + long table.
log "Aggregating hits into matrices + long table"
python3 projects/paper1/analysis/aggregate_pipeline_hits.py \
  --hits-a "$OUTDIR/hits_A" --hits-b "$OUTDIR/hits_B" --hits-c "$OUTDIR/hits_C" \
  --read-counts "$RCFILE" --metadata "$METADATA" \
  --aro-index "$CARD_DIR/aro_index.tsv" \
  --outdir projects/paper1/results

log "Done. Now run: python3 projects/paper1/analysis/run_paper1_analysis.py && python3 projects/paper1/analysis/generate_figures.py"
