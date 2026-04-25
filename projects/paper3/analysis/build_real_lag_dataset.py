#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build real lag sample/pair tables from enriched metadata")
    parser.add_argument(
        "--metadata",
        default="projects/paper3/metadata/metadata_enriched.csv",
        help="Enriched metadata CSV",
    )
    parser.add_argument(
        "--samples-out",
        default="projects/paper3/results/lag_samples_actual.csv",
        help="Output CSV with per-sample lag usability",
    )
    parser.add_argument(
        "--pairs-relaxed-out",
        default="projects/paper3/results/lag_pairs_relaxed_actual.csv",
        help="Output CSV with relaxed lag pairs",
    )
    parser.add_argument(
        "--pairs-strict-out",
        default="projects/paper3/results/lag_pairs_strict_actual.csv",
        help="Output CSV with strict lag pairs",
    )
    parser.add_argument(
        "--features-template-out",
        default="projects/paper3/results/features_real_template.csv",
        help="Output CSV template for lag regression with real sample IDs",
    )
    return parser.parse_args()


def _clean(value: str | None) -> str:
    return (value or "").strip()


def _parse_dt(value: str) -> datetime | None:
    text = _clean(value)
    if not text or text.upper() == "NA":
        return None
    formats = [
        "%Y-%m-%d",
        "%Y-%m",
        "%Y",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%S",
    ]
    for fmt in formats:
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00").replace("+00:00", ""))
    except ValueError:
        return None


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()

    metadata_path = Path(args.metadata)
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    required = {"sample_id", "study", "read_count", "collection_date", "temporal_source"}
    missing = [col for col in required if col not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"Missing columns in metadata: {', '.join(missing)}")

    by_study: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        study = _clean(row.get("study")) or "unknown_study"
        by_study[study].append(row)

    sample_rows: list[dict[str, str]] = []
    relaxed_pairs: list[dict[str, str]] = []
    strict_pairs: list[dict[str, str]] = []
    features_template_rows: list[dict[str, str]] = []

    total_relaxed_samples = 0
    total_strict_samples = 0

    for study, study_rows in sorted(by_study.items(), key=lambda kv: kv[0]):
        parsed_items: list[tuple[dict[str, str], datetime]] = []
        for row in study_rows:
            dt = _parse_dt(row.get("collection_date", ""))
            if dt is not None:
                parsed_items.append((row, dt))

        parsed_items.sort(key=lambda x: (x[1], _clean(x[0].get("sample_id"))))
        n_parsed = len(parsed_items)

        ts_counter = Counter(dt.isoformat() for _, dt in parsed_items)

        relaxed_order: dict[str, int] = {}
        strict_order: dict[str, int] = {}

        for i, (row, dt) in enumerate(parsed_items, start=1):
            sample_id = _clean(row.get("sample_id"))
            relaxed_order[sample_id] = i

        strict_candidates: list[tuple[dict[str, str], datetime]] = []
        for row, dt in parsed_items:
            if ts_counter[dt.isoformat()] == 1:
                strict_candidates.append((row, dt))

        for i, (row, dt) in enumerate(strict_candidates, start=1):
            sample_id = _clean(row.get("sample_id"))
            strict_order[sample_id] = i

        total_relaxed_samples += len(parsed_items)
        total_strict_samples += len(strict_candidates)

        # Build sample rows.
        for row in study_rows:
            sample_id = _clean(row.get("sample_id"))
            read_count = _clean(row.get("read_count"))
            temporal_source = _clean(row.get("temporal_source"))
            collection_date = _clean(row.get("collection_date"))
            dt = _parse_dt(collection_date)

            relaxed_ok = "yes" if sample_id in relaxed_order and n_parsed >= 2 else "no"
            strict_ok = "yes" if sample_id in strict_order and len(strict_candidates) >= 3 else "no"
            ambiguous = "yes" if dt is not None and ts_counter[dt.isoformat()] > 1 else "no"

            note = (
                f"parsed_in_study={n_parsed}; strict_unique_dates={len(strict_candidates)}; "
                f"temporal_source={temporal_source}"
            )

            sample_rows.append(
                {
                    "sample_id": sample_id,
                    "study": study,
                    "collection_date": collection_date,
                    "read_count": read_count,
                    "temporal_source": temporal_source,
                    "order_relaxed": str(relaxed_order.get(sample_id, "")),
                    "order_strict": str(strict_order.get(sample_id, "")),
                    "ambiguous_date_tie": ambiguous,
                    "usable_relaxed": relaxed_ok,
                    "usable_strict": strict_ok,
                    "notes": note,
                }
            )

            if strict_ok == "yes":
                features_template_rows.append(
                    {
                        "study": study,
                        "sample_id": sample_id,
                        "order": str(strict_order[sample_id]),
                        "mge_abundance": "",
                        "entropy": "",
                        "arg_total": "",
                        "sequencing_depth": read_count,
                        "metric_status": "missing_quantification",
                    }
                )

        # Build relaxed pairs.
        if n_parsed >= 2:
            for i in range(n_parsed - 1):
                left_row, left_dt = parsed_items[i]
                right_row, right_dt = parsed_items[i + 1]
                relaxed_pairs.append(
                    {
                        "study": study,
                        "t_sample_id": _clean(left_row.get("sample_id")),
                        "t1_sample_id": _clean(right_row.get("sample_id")),
                        "t_collection_date": left_dt.isoformat(),
                        "t1_collection_date": right_dt.isoformat(),
                        "delta_days": f"{(right_dt - left_dt).total_seconds() / 86400:.3f}",
                        "pair_quality": "relaxed",
                    }
                )

        # Build strict pairs.
        if len(strict_candidates) >= 2:
            for i in range(len(strict_candidates) - 1):
                left_row, left_dt = strict_candidates[i]
                right_row, right_dt = strict_candidates[i + 1]
                strict_pairs.append(
                    {
                        "study": study,
                        "t_sample_id": _clean(left_row.get("sample_id")),
                        "t1_sample_id": _clean(right_row.get("sample_id")),
                        "t_collection_date": left_dt.isoformat(),
                        "t1_collection_date": right_dt.isoformat(),
                        "delta_days": f"{(right_dt - left_dt).total_seconds() / 86400:.3f}",
                        "pair_quality": "strict",
                    }
                )

    _write_csv(
        Path(args.samples_out),
        [
            "sample_id",
            "study",
            "collection_date",
            "read_count",
            "temporal_source",
            "order_relaxed",
            "order_strict",
            "ambiguous_date_tie",
            "usable_relaxed",
            "usable_strict",
            "notes",
        ],
        sample_rows,
    )

    _write_csv(
        Path(args.pairs_relaxed_out),
        [
            "study",
            "t_sample_id",
            "t1_sample_id",
            "t_collection_date",
            "t1_collection_date",
            "delta_days",
            "pair_quality",
        ],
        relaxed_pairs,
    )

    _write_csv(
        Path(args.pairs_strict_out),
        [
            "study",
            "t_sample_id",
            "t1_sample_id",
            "t_collection_date",
            "t1_collection_date",
            "delta_days",
            "pair_quality",
        ],
        strict_pairs,
    )

    _write_csv(
        Path(args.features_template_out),
        [
            "study",
            "sample_id",
            "order",
            "mge_abundance",
            "entropy",
            "arg_total",
            "sequencing_depth",
            "metric_status",
        ],
        features_template_rows,
    )

    print(f"Metadata rows processed: {len(rows)}")
    print(f"Relaxed lag-usable samples: {total_relaxed_samples}")
    print(f"Strict lag-usable samples: {total_strict_samples}")
    print(f"Relaxed t->t+1 pairs: {len(relaxed_pairs)}")
    print(f"Strict t->t+1 pairs: {len(strict_pairs)}")
    print(f"Wrote: {args.samples_out}")
    print(f"Wrote: {args.pairs_relaxed_out}")
    print(f"Wrote: {args.pairs_strict_out}")
    print(f"Wrote: {args.features_template_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
