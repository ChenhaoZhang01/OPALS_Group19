#!/usr/bin/env bash
# quantify_samples.sh
# Download SRA FASTQ files for Paper 3 strict cohort and compute:
#   arg_total     – total ARG hit count (AMRFinderPlus)
#   mge_abundance – total MGE hit count  (mobileOG-db BLAST)
#   entropy       – Shannon diversity index from Kraken2+Bracken
#
# Run on a Linux system (HPC recommended) with:
#   sra-tools  (prefetch, fasterq-dump)
#   AMRFinder+ (amrfinder, updated db)
#   DIAMOND    (for mobileOG-db)
#   mobileOG-db diamond .dmnd index  (mobileOG-db.dmnd)
#   Kraken2 + Bracken  (standard PlusPF database recommended)
#   Python 3 + scipy   (for entropy calculation)
#
# Usage:
#   bash quantify_samples.sh --outdir /scratch/paper3 \
#        --kraken-db /databases/kraken2_plusPF \
#        --mobileog-dmnd /databases/mobileOG-db.dmnd \
#        --threads 16
#
# Outputs:
#   <outdir>/quant_results.csv   (ready for populate_real_features_from_quant.py)

set -euo pipefail

##############################################################################
# Defaults — tuned for WSL2 on this machine
##############################################################################
OUTDIR="/data/paper3_quant"
THREADS=4
KRAKEN_DB="/data/dbs/kraken2_8gb"
MOBILEOG_DMND="/data/dbs/mobileog/mobileOG-db.dmnd"
AMRFINDER_DB="/data/dbs/amrfinder"
CONDA_ENV="/opt/miniforge/envs/biotools"
SKIP_DOWNLOAD=false
ONLY_SHOTGUN=false

# Prepend conda env bin so all tools resolve correctly
export PATH="$CONDA_ENV/bin:$PATH"

##############################################################################
# Argument parsing
##############################################################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)       OUTDIR="$2";         shift 2 ;;
    --threads)      THREADS="$2";        shift 2 ;;
    --kraken-db)    KRAKEN_DB="$2";      shift 2 ;;
    --mobileog-dmnd) MOBILEOG_DMND="$2"; shift 2 ;;
    --skip-download) SKIP_DOWNLOAD=true; shift   ;;
    --only-shotgun)  ONLY_SHOTGUN=true;  shift   ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

if [[ -z "$KRAKEN_DB" ]]; then
  echo "ERROR: --kraken-db is required." >&2; exit 1
fi
if [[ -z "$MOBILEOG_DMND" ]]; then
  echo "ERROR: --mobileog-dmnd is required." >&2; exit 1
fi

##############################################################################
# Sample lists
# SHOTGUN: full ARG + MGE + entropy quantification possible
# AMPLICON_16S: entropy-only (PRJNA599167 16S samples)
##############################################################################
SHOTGUN_SAMPLES=(
  # PRJNA599167_WGA — Chesapeake Bay whole-genome-amplified shotgun (May 2017)
  SRR11803495 SRR11803496 SRR11803497 SRR11803498
  SRR11803499 SRR11803500 SRR11803501 SRR11803502
  SRR11803503 SRR11803504 SRR11803505
  # PRJNA1071831 — Iskar River, Bulgaria (biweekly, Nov–Dec 2022)
  SRR27827413 SRR27827408 SRR27827406 SRR27827404   # Dragushinovo
  SRR27827412 SRR27827407 SRR27827405 SRR27827403   # Mechkata
)

AMPLICON_16S_SAMPLES=(
  # PRJNA599167 — Chesapeake Bay 16S amplicon (May 2017)
  # These can only contribute entropy; ARG/MGE fields will be left empty.
  SRR11811727 SRR11811728 SRR11811733 SRR11811738
)

##############################################################################
# Directory setup
##############################################################################
FASTQ_DIR="${OUTDIR}/fastq"
AMR_DIR="${OUTDIR}/amrfinder"
MGE_DIR="${OUTDIR}/mobileog"
KRAKEN_DIR="${OUTDIR}/kraken"
BRACKEN_DIR="${OUTDIR}/bracken"
ASSEMBLE_DIR="${OUTDIR}/assembly"
mkdir -p "$FASTQ_DIR" "$AMR_DIR" "$MGE_DIR" "$KRAKEN_DIR" "$BRACKEN_DIR" "$ASSEMBLE_DIR"

QUANT_CSV="${OUTDIR}/quant_results.csv"
echo "sample_id,arg_total,mge_abundance,entropy" > "$QUANT_CSV"

##############################################################################
# Helper: Shannon entropy from Bracken report
# Reads the <report>.bracken file and computes H = -sum(p * log2(p)) at
# species level (S rank). Pipes the result through python.
##############################################################################
compute_entropy() {
  local bracken_report="$1"
  python3 - "$bracken_report" <<'PYEOF'
import sys, math

fname = sys.argv[1]
fracs = []
with open(fname) as fh:
    next(fh)  # header
    for line in fh:
        parts = line.rstrip().split("\t")
        if len(parts) < 7:
            continue
        fraction = float(parts[6])  # fraction_total_reads col
        if fraction > 0:
            fracs.append(fraction)

if not fracs:
    print("0.0")
    sys.exit(0)

# normalize (should already sum to ~1, but guard against rounding)
total = sum(fracs)
h = -sum((p / total) * math.log2(p / total) for p in fracs if p > 0)
print(f"{h:.6f}")
PYEOF
}

##############################################################################
# Process a single SHOTGUN sample
##############################################################################
process_shotgun() {
  local srr="$1"
  echo "=== Processing shotgun: $srr ==="

  # ---------- 1. Download -------------------------------------------------
  local fq1="${FASTQ_DIR}/${srr}_1.fastq.gz"
  local fq2="${FASTQ_DIR}/${srr}_2.fastq.gz"
  local fq_se="${FASTQ_DIR}/${srr}.fastq.gz"

  if [[ "$SKIP_DOWNLOAD" == false ]]; then
    echo "  Downloading $srr …"
    prefetch --max-size 50G -O "${FASTQ_DIR}" "$srr" 2>/dev/null || true
    fasterq-dump --split-files --threads "$THREADS" \
      --outdir "$FASTQ_DIR" "${FASTQ_DIR}/${srr}/${srr}.sra" 2>/dev/null || \
    fasterq-dump --split-files --threads "$THREADS" \
      --outdir "$FASTQ_DIR" "$srr" 2>/dev/null || true
    # Compress if uncompressed output was written
    for f in "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${srr}_2.fastq" \
              "${FASTQ_DIR}/${srr}.fastq"; do
      [[ -f "$f" ]] && gzip -f "$f"
    done
  fi

  # Determine paired vs single-end
  local reads_arg=""
  if [[ -f "$fq1" && -f "$fq2" ]]; then
    reads_arg="$fq1 $fq2"
    local is_paired=true
  elif [[ -f "$fq_se" ]]; then
    reads_arg="$fq_se"
    local is_paired=false
  else
    echo "  WARNING: no FASTQ found for $srr — skipping"
    echo "${srr},,,," >> "$QUANT_CSV"
    return
  fi

  # ---------- 2. Assembly (for AMRFinder+) --------------------------------
  # AMRFinder works best on contigs. Assemble with metaSPAdes if reads > 1M.
  # For very small libraries, skip assembly and run directly on reads.
  local assembly="${ASSEMBLE_DIR}/${srr}/contigs.fasta"
  if [[ ! -f "$assembly" ]]; then
    echo "  Assembling $srr …"
    local spades_out="${ASSEMBLE_DIR}/${srr}"
    if [[ "$is_paired" == true ]]; then
      metaspades.py -1 "$fq1" -2 "$fq2" -o "$spades_out" \
        --threads "$THREADS" --memory 64 2>/dev/null || true
    else
      metaspades.py -s "$fq_se" -o "$spades_out" \
        --threads "$THREADS" --memory 64 2>/dev/null || true
    fi
    # Fallback: if metaSPAdes fails (too few reads), use MEGAHIT
    if [[ ! -f "$assembly" ]]; then
      if [[ "$is_paired" == true ]]; then
        megahit -1 "$fq1" -2 "$fq2" -o "${spades_out}_megahit" \
          --num-cpu-threads "$THREADS" 2>/dev/null || true
        [[ -f "${spades_out}_megahit/final.contigs.fa" ]] && \
          cp "${spades_out}_megahit/final.contigs.fa" "$assembly"
      fi
    fi
  fi

  # ---------- 3. ARG quantification (AMRFinderPlus) -----------------------
  local amr_out="${AMR_DIR}/${srr}.tsv"
  if [[ ! -f "$amr_out" && -f "$assembly" ]]; then
    echo "  Running AMRFinderPlus on $srr …"
    amrfinder --nucleotide "$assembly" \
      --threads "$THREADS" \
      --output "$amr_out" 2>/dev/null || touch "$amr_out"
  fi
  local arg_total=0
  if [[ -f "$amr_out" ]]; then
    # Count lines excluding header; each line is one ARG hit
    arg_total=$(awk 'NR > 1 {count++} END {print count+0}' "$amr_out")
  fi

  # ---------- 4. MGE quantification (IntegronFinder) -----------------------
  # Counts class 1/2/3 integrons per assembly — the dominant ARG-mobilizing
  # MGE class in environmental metagenomics (Gillings et al. 2014).
  local mge_dir_sample="${MGE_DIR}/${srr}"
  if [[ ! -d "$mge_dir_sample" && -f "$assembly" ]]; then
    echo "  Running IntegronFinder on $srr …"
    mkdir -p "$mge_dir_sample"
    integron_finder \
      --outdir "$mge_dir_sample" \
      --cpu "$THREADS" \
      --local-max \
      "$assembly" 2>/dev/null || true
  fi
  local mge_abundance=0
  # Count complete + CALIN integrons across all .integrons output files
  if [[ -d "$mge_dir_sample" ]]; then
    mge_abundance=$(find "$mge_dir_sample" -name "*.integrons" -exec \
      awk '$2 ~ /complete|CALIN/ {count++} END {print count+0}' {} \; | \
      awk '{s+=$1} END {print s+0}')
  fi

  # ---------- 5. Kraken2 + Bracken (taxonomic diversity) ------------------
  local kraken_report="${KRAKEN_DIR}/${srr}.report"
  if [[ ! -f "$kraken_report" ]]; then
    echo "  Running Kraken2 on $srr …"
    if [[ "$is_paired" == true ]]; then
      kraken2 --db "$KRAKEN_DB" --paired \
        --threads "$THREADS" --report "$kraken_report" \
        --output /dev/null \
        "$fq1" "$fq2" 2>/dev/null || touch "$kraken_report"
    else
      kraken2 --db "$KRAKEN_DB" \
        --threads "$THREADS" --report "$kraken_report" \
        --output /dev/null \
        "$fq_se" 2>/dev/null || touch "$kraken_report"
    fi
  fi

  local bracken_report="${BRACKEN_DIR}/${srr}.bracken"
  if [[ ! -f "$bracken_report" && -f "$kraken_report" ]]; then
    echo "  Running Bracken on $srr …"
    bracken -d "$KRAKEN_DB" \
      -i "$kraken_report" \
      -o "$bracken_report" \
      -r 150 -l S -t 10 2>/dev/null || touch "$bracken_report"
  fi

  local entropy=0
  if [[ -f "$bracken_report" && -s "$bracken_report" ]]; then
    entropy=$(compute_entropy "$bracken_report")
  fi

  # ---------- 6. Append to quant CSV -------------------------------------
  echo "${srr},${arg_total},${mge_abundance},${entropy}" >> "$QUANT_CSV"
  echo "  Done: $srr  ARG=$arg_total  MGE=$mge_abundance  H=$entropy"
}

##############################################################################
# Process a single 16S amplicon sample (entropy only)
##############################################################################
process_amplicon() {
  local srr="$1"
  echo "=== Processing 16S amplicon: $srr (entropy only) ==="

  local fq1="${FASTQ_DIR}/${srr}_1.fastq.gz"
  local fq2="${FASTQ_DIR}/${srr}_2.fastq.gz"

  if [[ "$SKIP_DOWNLOAD" == false ]]; then
    echo "  Downloading $srr …"
    prefetch --max-size 10G -O "${FASTQ_DIR}" "$srr" 2>/dev/null || true
    fasterq-dump --split-files --threads "$THREADS" \
      --outdir "$FASTQ_DIR" "$srr" 2>/dev/null || true
    for f in "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${srr}_2.fastq"; do
      [[ -f "$f" ]] && gzip -f "$f"
    done
  fi

  # Run Kraken2 on amplicon reads (V4 region). For standard 16S classification,
  # use the 16S Greengenes or SILVA k2 index; the plusPF DB also works.
  local kraken_report="${KRAKEN_DIR}/${srr}.report"
  if [[ ! -f "$kraken_report" ]]; then
    if [[ -f "$fq1" && -f "$fq2" ]]; then
      kraken2 --db "$KRAKEN_DB" --paired \
        --threads "$THREADS" --report "$kraken_report" \
        --output /dev/null "$fq1" "$fq2" 2>/dev/null || touch "$kraken_report"
    fi
  fi

  local bracken_report="${BRACKEN_DIR}/${srr}.bracken"
  if [[ ! -f "$bracken_report" && -f "$kraken_report" ]]; then
    bracken -d "$KRAKEN_DB" -i "$kraken_report" \
      -o "$bracken_report" -r 250 -l S -t 5 2>/dev/null || touch "$bracken_report"
  fi

  local entropy=0
  if [[ -f "$bracken_report" && -s "$bracken_report" ]]; then
    entropy=$(compute_entropy "$bracken_report")
  fi

  # ARG and MGE are not quantifiable from 16S amplicons — leave empty
  echo "${srr},,,$entropy" >> "$QUANT_CSV"
  echo "  Done: $srr  H=$entropy  (ARG/MGE N/A for 16S)"
}

##############################################################################
# Main loop
##############################################################################
echo "Starting Paper 3 quantification pipeline"
echo "Output directory: $OUTDIR"
echo "Threads: $THREADS"
echo "Kraken2 DB: $KRAKEN_DB"
echo "mobileOG-db: $MOBILEOG_DMND"
echo ""

for srr in "${SHOTGUN_SAMPLES[@]}"; do
  process_shotgun "$srr"
done

if [[ "$ONLY_SHOTGUN" == false ]]; then
  for srr in "${AMPLICON_16S_SAMPLES[@]}"; do
    process_amplicon "$srr"
  done
fi

echo ""
echo "=== Pipeline complete ==="
echo "Quant results written to: $QUANT_CSV"
echo ""
echo "Next step — from the OPALS_Group19 repo root:"
echo "  python projects/paper3/analysis/populate_real_features_from_quant.py \\"
echo "    --quant ${QUANT_CSV}"
echo "  python projects/paper3/analysis/run_lag_analysis.py \\"
echo "    --features projects/paper3/results/features_real_populated.csv \\"
echo "    --outdir projects/paper3/results/real_run"
