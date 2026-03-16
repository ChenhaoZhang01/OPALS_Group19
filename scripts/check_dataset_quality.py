#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate ARG dataset quality thresholds")
    parser.add_argument(
        "--dataset",
        default="results/ARG_dataset.csv",
        help="Path to merged ARG dataset (default: results/ARG_dataset.csv)",
    )
    parser.add_argument(
        "--arg-matrix",
        default="results/ARG_matrix.csv",
        help="Path to ARG matrix (default: results/ARG_matrix.csv)",
    )
    parser.add_argument(
        "--report",
        default="results/dataset_quality_report.txt",
        help="Path to quality report output (default: results/dataset_quality_report.txt)",
    )
    parser.add_argument("--min-samples", type=int, default=60)
    parser.add_argument("--max-samples", type=int, default=90)
    parser.add_argument("--min-detections", type=int, default=50)
    parser.add_argument("--min-environments", type=int, default=2)
    return parser.parse_args()


def parse_float(value: str) -> float:
    if value is None:
        return 0.0
    text = value.strip()
    if not text:
        return 0.0
    try:
        return float(text)
    except ValueError:
        return 0.0


def evaluate(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


def main() -> int:
    args = parse_args()

    with open(args.dataset, "r", encoding="utf-8", newline="") as handle:
        dreader = csv.DictReader(handle)
        dataset_rows = [row for row in dreader]

    sample_count = len(dataset_rows)
    envs = {(row.get("environment") or "NA").strip() or "NA" for row in dataset_rows}

    with open(args.arg_matrix, "r", encoding="utf-8", newline="") as handle:
        mreader = csv.DictReader(handle)
        matrix_rows = [row for row in mreader]
        fieldnames = mreader.fieldnames or []

    sample_col = "Sample" if "Sample" in fieldnames else ("sample_id" if "sample_id" in fieldnames else (fieldnames[0] if fieldnames else "sample_id"))
    arg_cols = [c for c in fieldnames if c != sample_col]

    detections = 0
    for row in matrix_rows:
        for col in arg_cols:
            if parse_float(row.get(col, "0")) > 0:
                detections += 1

    checks = [
        (
            f"Samples between {args.min_samples} and {args.max_samples}",
            evaluate(args.min_samples <= sample_count <= args.max_samples),
            f"observed={sample_count}",
        ),
        (
            f"ARG detections >= {args.min_detections}",
            evaluate(detections >= args.min_detections),
            f"observed={detections}",
        ),
        (
            f"Environments >= {args.min_environments}",
            evaluate(len(envs) >= args.min_environments),
            f"observed={len(envs)}",
        ),
    ]

    os.makedirs(os.path.dirname(args.report), exist_ok=True)
    with open(args.report, "w", encoding="utf-8") as handle:
        handle.write("Dataset Quality Report\n")
        handle.write("======================\n")
        for label, status, observed in checks:
            handle.write(f"{status}: {label} ({observed})\n")

    for label, status, observed in checks:
        print(f"{status}: {label} ({observed})")
    print(f"Wrote report: {args.report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
