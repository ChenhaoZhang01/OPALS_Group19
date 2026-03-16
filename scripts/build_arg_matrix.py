#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import glob
import os
from collections import Counter, defaultdict


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build ARG abundance matrix from arg_hits.tsv files")
    parser.add_argument(
        "--hits-glob",
        default="results/*/arg_hits.tsv",
        help="Glob pattern to find arg_hits.tsv files (default: results/*/arg_hits.tsv)",
    )
    parser.add_argument(
        "--output",
        default="results/ARG_matrix.csv",
        help="Output matrix CSV path (default: results/ARG_matrix.csv)",
    )
    return parser.parse_args()


def sample_id_from_path(path: str) -> str:
    return os.path.basename(os.path.dirname(path))


def parse_arg_hits(path: str) -> Counter:

    counts: Counter = Counter()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            arg_id = cols[1]
            counts[arg_id] += 1
    return counts


def write_matrix(sample_to_counts: dict[str, Counter], output_path: str) -> None:
    all_args = sorted({arg for counts in sample_to_counts.values() for arg in counts.keys()})
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Sample", *all_args])

        for sample in sorted(sample_to_counts.keys()):
            row = [sample]
            counts = sample_to_counts[sample]
            row.extend(str(counts.get(arg, 0)) for arg in all_args)
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    hit_files = sorted(glob.glob(args.hits_glob))

    if not hit_files:
        print(f"No hit files found with pattern: {args.hits_glob}")
        return 1

    sample_to_counts: dict[str, Counter] = defaultdict(Counter)
    for hit_file in hit_files:
        sample = sample_id_from_path(hit_file)
        sample_to_counts[sample] = parse_arg_hits(hit_file)

    write_matrix(sample_to_counts=sample_to_counts, output_path=args.output)
    print(f"Wrote ARG matrix to: {args.output}")
    print(f"Samples: {len(sample_to_counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
