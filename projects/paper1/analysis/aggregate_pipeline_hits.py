#!/usr/bin/env python3
"""Aggregate raw three-pipeline ARG hits into the matrices + long table for Paper 1.

Consumes the per-sample outputs written by run_arg_pipeline_hpc.sh:
  hits_A/<SRR>.tsv        DIAMOND blastp tabular: col 1 qseqid, col 2 sseqid (CARD ref)
  hits_B/<SRR>.idxstats   samtools idxstats: ref, length, mapped, unmapped
  hits_C/<SRR>.txt        RGI main tabular (has Best_Hit_ARO / ARO / Drug Class columns)

Counts are normalized by library size (read_counts.txt) and written as:
  ARG_matrix_pipeline{A,B,C}.csv   wide, one column per ARG feature (drug class if an
                                   --aro-index is given and maps, else raw CARD ref)
  pipeline_long_table.csv          sample_id, environment, pipeline, ARG_total, ARG_richness

This is the EMPIRICAL replacement for the calibrated build_pipeline_matrices.py; the
downstream run_paper1_analysis.py / generate_figures.py are unchanged. See
results/DATA_PROVENANCE.md.

ARG_total  = sum of normalized hits (pipeline-comparable abundance).
ARG_richness = number of distinct ARG features detected (drug classes if mapped, else refs).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

ARO_RE = re.compile(r"ARO:\d+")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--hits-a", required=True)
    p.add_argument("--hits-b", required=True)
    p.add_argument("--hits-c", required=True)
    p.add_argument("--read-counts", required=True, help="whitespace: SRR total_reads")
    p.add_argument("--metadata", required=True)
    p.add_argument("--aro-index", default="", help="CARD aro_index.tsv for drug-class collapse (optional)")
    p.add_argument("--outdir", default="projects/paper1/results")
    return p.parse_args()


def load_read_counts(path: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[1].isdigit():
            counts[parts[0]] = int(parts[1])
    return counts


def load_env(path: str) -> dict[str, str]:
    with open(path, encoding="utf-8", newline="") as fh:
        return {r["sample_id"]: r["environment"] for r in csv.DictReader(fh)}


def load_aro_to_class(path: str) -> dict[str, str]:
    """Map any CARD identifier token (ARO accession, protein/DNA accession) -> Drug Class."""
    mapping: dict[str, str] = {}
    if not path:
        return mapping
    p = Path(path)
    if not p.exists():
        print(f"WARN: aro-index not found ({path}); using raw CARD references.", file=sys.stderr)
        return mapping
    with p.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            dclass = (row.get("Drug Class") or "").strip() or "unclassified"
            for col in ("ARO Accession", "Protein Accession", "DNA Accession", "ARO Name"):
                tok = (row.get(col) or "").strip()
                if tok:
                    mapping[tok] = dclass
    return mapping


def feature_of(subject: str, aro_map: dict[str, str]) -> tuple[str, bool]:
    """Return (feature_name, mapped?). feature = drug class if resolvable, else raw subject."""
    if not aro_map:
        return subject, False
    # try whole token, then any ARO:xxxx inside, then pipe-split tokens
    if subject in aro_map:
        return aro_map[subject], True
    m = ARO_RE.search(subject)
    if m and m.group(0) in aro_map:
        return aro_map[m.group(0)], True
    for tok in re.split(r"[|\s]", subject):
        if tok and tok in aro_map:
            return aro_map[tok], True
    return subject, False


def hits_from_diamond(tsv: Path) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    if not tsv.exists():
        return counts
    with tsv.open(encoding="utf-8") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                counts[cols[1]] += 1
    return counts


def hits_from_idxstats(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    if not path.exists():
        return counts
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            parts = (line.rstrip("\n").split("\t") + ["0", "0", "0", "0"])[:4]
            ref, _length, mapped, _unmapped = parts
            if ref != "*" and mapped.isdigit() and int(mapped) > 0:
                counts[ref] = int(mapped)
    return counts


def hits_from_rgi(path: Path) -> tuple[dict[str, int], bool]:
    """RGI tabular -> counts keyed by Drug Class if present (mapped=True), else Best_Hit_ARO."""
    counts: dict[str, int] = defaultdict(int)
    if not path.exists() or path.stat().st_size == 0:
        return counts, False
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        if "Drug Class" in fields:
            key, mapped = "Drug Class", True
        elif "Best_Hit_ARO" in fields:
            key, mapped = "Best_Hit_ARO", False
        elif "ARO" in fields:
            key, mapped = "ARO", False
        else:
            return counts, False
        for row in reader:
            val = (row.get(key) or "").strip()
            if not val:
                continue
            # RGI "Drug Class" can be "; "-separated; count each class once per hit.
            for piece in (val.split(";") if mapped else [val]):
                piece = piece.strip()
                if piece:
                    counts[piece] += 1
    return counts, mapped


def write_matrix(outdir: Path, pipe: str, reads: dict[str, int],
                 per_sample: dict[str, dict[str, float]], features: list[str]) -> None:
    with (outdir / f"ARG_matrix_{pipe}.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", *features])
        for srr in reads:
            row = per_sample.get(srr, {})
            w.writerow([srr, *[f"{row.get(g, 0.0):.8e}" for g in features]])
    print(f"Wrote ARG_matrix_{pipe}.csv ({len(reads)} samples x {len(features)} features)")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    reads = load_read_counts(args.read_counts)
    env = load_env(args.metadata)
    aro_map = load_aro_to_class(args.aro_index)
    if not reads:
        print("ERROR: no read counts loaded.", file=sys.stderr)
        return 1

    long_rows = []

    for pipe, hdir, ext in [
        ("pipelineA", Path(args.hits_a), ".tsv"),
        ("pipelineB", Path(args.hits_b), ".idxstats"),
        ("pipelineC", Path(args.hits_c), ".txt"),
    ]:
        per_sample: dict[str, dict[str, float]] = {}
        all_feats: set[str] = set()
        matched = total = 0
        for srr, lib in reads.items():
            if pipe == "pipelineA":
                raw = hits_from_diamond(hdir / f"{srr}{ext}")
            elif pipe == "pipelineB":
                raw = hits_from_idxstats(hdir / f"{srr}{ext}")
            else:
                raw, _ = hits_from_rgi(hdir / f"{srr}{ext}")

            # Collapse to features (drug class where possible).
            feat_counts: dict[str, float] = defaultdict(float)
            for subj, c in raw.items():
                if pipe == "pipelineC":
                    feat = subj  # already a class or ARO from RGI
                else:
                    feat, ok = feature_of(subj, aro_map)
                    total += 1
                    matched += int(ok)
                feat_counts[feat] += c / lib if lib else 0.0
            per_sample[srr] = dict(feat_counts)
            all_feats.update(feat_counts)
            long_rows.append({
                "sample_id": srr,
                "environment": env.get(srr, "NA"),
                "pipeline": pipe,
                "ARG_total": sum(feat_counts.values()),
                "ARG_richness": len(feat_counts),
            })
        if pipe in ("pipelineA", "pipelineB") and aro_map and total:
            rate = 100.0 * matched / total
            if rate < 80:
                print(f"WARN: {pipe} drug-class mapping only {rate:.0f}% "
                      f"(check CARD header format vs aro_index).", file=sys.stderr)
        write_matrix(outdir, pipe, reads, per_sample, sorted(all_feats))

    with (outdir / "pipeline_long_table.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample_id", "environment", "pipeline", "ARG_total", "ARG_richness"])
        w.writeheader()
        for lr in long_rows:
            w.writerow({**lr, "ARG_total": f"{lr['ARG_total']:.8e}"})
    print(f"Wrote pipeline_long_table.csv ({len(long_rows)} rows)")

    # This is REAL pipeline output -> clear the simulation marker if present.
    marker = outdir / "SIMULATED_PLACEHOLDER.txt"
    if marker.exists():
        marker.unlink()
        print("Removed SIMULATED_PLACEHOLDER.txt (results are now real pipeline output).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
