#!/usr/bin/env bash
set -euo pipefail

# Fast phase-1 bootstrap: create venv if missing and generate metadata only.
OUTPUT="${1:-metadata/master_metadata.csv}"
PER_QUERY_LIMIT="${2:-5000}"
SLEEP_SECONDS="${3:-0.1}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

mkdir -p metadata

if [[ ! -d .venv ]]; then
  if command -v python3 >/dev/null 2>&1; then
    python3 -m venv .venv
  elif command -v python >/dev/null 2>&1; then
    python -m venv .venv
  else
    echo "ERROR: Python 3.10+ was not found." >&2
    exit 1
  fi
fi

if [[ -x .venv/bin/python ]]; then
  PYTHON_BIN=".venv/bin/python"
elif [[ -x .venv/Scripts/python.exe ]]; then
  PYTHON_BIN=".venv/Scripts/python.exe"
else
  echo "ERROR: Could not find virtual environment python executable." >&2
  exit 1
fi

"$PYTHON_BIN" scripts/build_master_metadata.py \
  --output "$OUTPUT" \
  --per-query-limit "$PER_QUERY_LIMIT" \
  --sleep "$SLEEP_SECONDS"

echo "Phase 1 complete: $OUTPUT"
