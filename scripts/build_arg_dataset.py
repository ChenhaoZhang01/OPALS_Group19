#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from typing import Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build normalized ARG matrix and merged ARG dataset"
    )
    parser.add_argument(
        "--metadata",
        default="metadata_final.csv",
        help="Path to cleaned metadata CSV (default: metadata_final.csv)",
    )
    parser.add_argument(
        "--arg-matrix",
        default="results/ARG_matrix.csv",
        help="Path to raw ARG matrix CSV (default: results/ARG_matrix.csv)",
    )
    parser.add_argument(
        "--normalized-out",
        default="results/ARG_matrix_normalized.csv",
        help="Path to write normalized ARG matrix (default: results/ARG_matrix_normalized.csv)",
    )
    parser.add_argument(
        "--dataset-out",
        default="results/ARG_dataset.csv",
        help="Path to write merged dataset (default: results/ARG_dataset.csv)",
    )
    return parser.parse_args()


def parse_read_count(value: str) -> int:
    if value is None:
        return 0
    raw = value.strip().upper()
    if not raw or raw == "NA":
        return 0
    if re.fullmatch(r"\d+", raw):
        return int(raw)
    match = re.fullmatch(r"(\d+(?:\.\d+)?)([KMB])", raw)
    if not match:
        return 0

    number = float(match.group(1))
    suffix = match.group(2)
    multiplier = {"K": 1_000, "M": 1_000_000, "B": 1_000_000_000}[suffix]
    return int(round(number * multiplier))


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


def shannon_diversity(values: Iterable[float]) -> float:
    positive = [v for v in values if v > 0]
    total = sum(positive)
    if total <= 0:
        return 0.0

    entropy = 0.0
    for value in positive:
        p = value / total
        entropy -= p * math.log(p)
    return entropy


def detect_sample_column(fieldnames: list[str]) -> str:
    candidates = ["sample_id", "Sample", "sample", "SampleID"]
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    return fieldnames[0]


def load_metadata(path: str) -> dict[str, dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = {}
        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                continue
            row["read_count"] = str(parse_read_count(row.get("read_count", "0")))
            rows[sample_id] = row
        return rows


def load_arg_matrix(path: str) -> tuple[str, list[str], list[dict[str, str]]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("ARG matrix has no header")

        sample_col = detect_sample_column(reader.fieldnames)
        arg_cols = [c for c in reader.fieldnames if c != sample_col]
        rows = [row for row in reader]
        return sample_col, arg_cols, rows


def write_normalized_matrix(
    rows: list[dict[str, str]],
    sample_col: str,
    arg_cols: list[str],
    metadata: dict[str, dict[str, str]],
    output_path: str,
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", *arg_cols])

        for row in rows:
            sample = (row.get(sample_col) or "").strip()
            if not sample:
                continue
            read_count = parse_read_count(metadata.get(sample, {}).get("read_count", "0"))
            normalized = []
            for arg in arg_cols:
                raw_count = parse_float(row.get(arg, "0"))
                value = (raw_count / read_count) if read_count > 0 else 0.0
                normalized.append(f"{value:.10f}")
            writer.writerow([sample, *normalized])


def write_dataset(
    rows: list[dict[str, str]],
    sample_col: str,
    arg_cols: list[str],
    metadata: dict[str, dict[str, str]],
    output_path: str,
) -> tuple[int, int]:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fieldnames = [
        "sample_id",
        "environment",
        "country",
        "year",
        "study",
        "timepoint",
        "treatment",
        "read_count",
        "ARG_total",
        "ARG_richness",
        "ARG_diversity",
        "ARG_total_normalized",
    ]

    written = 0
    missing_metadata = 0

    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            sample = (row.get(sample_col) or "").strip()
            if not sample:
                continue

            meta = metadata.get(sample)
            if meta is None:
                missing_metadata += 1
                continue

            values = [parse_float(row.get(arg, "0")) for arg in arg_cols]
            arg_total = float(sum(values))
            arg_richness = sum(1 for v in values if v > 0)
            arg_diversity = shannon_diversity(values)

            read_count = parse_read_count(meta.get("read_count", "0"))
            arg_total_normalized = (arg_total / read_count) if read_count > 0 else 0.0

            writer.writerow(
                {
                    "sample_id": sample,
                    "environment": meta.get("environment", "NA"),
                    "country": meta.get("country", "NA"),
                    "year": meta.get("year", "NA"),
                    "study": meta.get("study", "NA"),
                    "timepoint": meta.get("timepoint", "NA"),
                    "treatment": meta.get("treatment", "NA"),
                    "read_count": str(read_count),
                    "ARG_total": f"{arg_total:.6f}",
                    "ARG_richness": str(arg_richness),
                    "ARG_diversity": f"{arg_diversity:.6f}",
                    "ARG_total_normalized": f"{arg_total_normalized:.10f}",
                }
            )
            written += 1

    return written, missing_metadata


def main() -> int:
    args = parse_args()

    metadata = load_metadata(args.metadata)
    sample_col, arg_cols, rows = load_arg_matrix(args.arg_matrix)

    write_normalized_matrix(
        rows=rows,
        sample_col=sample_col,
        arg_cols=arg_cols,
        metadata=metadata,
        output_path=args.normalized_out,
    )

    written, missing_metadata = write_dataset(
        rows=rows,
        sample_col=sample_col,
        arg_cols=arg_cols,
        metadata=metadata,
        output_path=args.dataset_out,
    )

    print(f"Input metadata samples: {len(metadata)}")
    print(f"ARG matrix samples: {len(rows)}")
    print(f"ARG features: {len(arg_cols)}")
    print(f"Rows written to dataset: {written}")
    if missing_metadata > 0:
        print(f"WARNING: samples skipped due to missing metadata: {missing_metadata}")
    print(f"Wrote normalized matrix: {args.normalized_out}")
    print(f"Wrote merged dataset: {args.dataset_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
