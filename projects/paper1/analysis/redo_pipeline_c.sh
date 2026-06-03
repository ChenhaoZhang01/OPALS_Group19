#!/usr/bin/env bash
# Regenerate Pipeline C (RGI) for all samples that already have predicted proteins,
# using the Loose-inclusive tier. Needed because RGI "strict" returns ~0 on
# environmental metagenomes (all ARG homologs fall below CARD clinical cutoffs).
# Safe to run after run_arg_pipeline_hpc.sh; only rewrites hits_C/.
#
# Usage: bash redo_pipeline_c.sh --outdir /home/chenh/paper1_scratch --threads 16

set -uo pipefail
OUTDIR="/home/chenh/paper1_scratch"
THREADS=8
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

shopt -s nullglob
for faa in "$OUTDIR"/proteins/*.faa; do
  base="$(basename "$faa")"
  [[ "$base" == *.rgi.faa ]] && continue
  SRR="${base%.faa}"
  rgifaa="$OUTDIR/proteins/$SRR.rgi.faa"
  sed 's/\*//g' "$faa" > "$rgifaa"
  echo "[C] RGI (loose) for $SRR ($(grep -c '>' "$rgifaa") ORFs)"
  rgi main -i "$rgifaa" -o "$OUTDIR/hits_C/$SRR" \
      -t protein -a DIAMOND --clean --include_loose --num_threads "$THREADS" \
      > "$OUTDIR/logs/$SRR.rgi_loose.log" 2>&1 || echo "  RGI failed for $SRR (see log)"
  rows=$(($(wc -l < "$OUTDIR/hits_C/$SRR.txt" 2>/dev/null) - 1))
  echo "  -> $rows hits"
done
echo "Pipeline C regeneration complete."
