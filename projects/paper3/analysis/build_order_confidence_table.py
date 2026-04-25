#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable


TIMEPOINT_PATTERNS: list[tuple[re.Pattern[str], Callable[[re.Match[str]], float]]] = [
    (re.compile(r"\bt\s*(\d+)\b", flags=re.IGNORECASE), lambda m: float(m.group(1))),
    (re.compile(r"\btime\s*point\s*(\d+)\b", flags=re.IGNORECASE), lambda m: float(m.group(1))),
    (re.compile(r"\bday\s*(\d+)\b", flags=re.IGNORECASE), lambda m: float(m.group(1))),
    (re.compile(r"\bweek\s*(\d+)\b", flags=re.IGNORECASE), lambda m: float(m.group(1) * 7)),
    (re.compile(r"\bmonth\s*(\d+)\b", flags=re.IGNORECASE), lambda m: float(m.group(1) * 30)),
]


@dataclass
class RankedSample:
    index: int
    sample_id: str
    study: str
    parsed_value: float | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build order-confidence table for lag analysis from metadata",
    )
    parser.add_argument(
        "--metadata",
        default="projects/paper3/metadata/metadata_final.csv",
        help="Input metadata CSV",
    )
    parser.add_argument(
        "--out",
        default="projects/paper3/results/order_confidence.csv",
        help="Output order-confidence CSV",
    )
    return parser.parse_args()


def _clean(value: str | None) -> str:
    return (value or "").strip()


def _parse_timepoint(value: str) -> float | None:
    text = _clean(value)
    if not text or text.upper() == "NA":
        return None
    for pattern, transform in TIMEPOINT_PATTERNS:
        match = pattern.search(text)
        if match:
            try:
                return float(transform(match))
            except (TypeError, ValueError):
                return None
    return None


def _parse_date_like(value: str) -> float | None:
    text = _clean(value)
    if not text or text.upper() == "NA":
        return None

    # Support full date, date-time, and year-month strings.
    candidates = [
        "%Y-%m-%d",
        "%Y/%m/%d",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m",
        "%Y",
    ]
    for fmt in candidates:
        try:
            dt = datetime.strptime(text, fmt)
            return dt.timestamp()
        except ValueError:
            continue

    # Handle ISO timestamps that include timezone suffixes.
    iso_text = text.replace("Z", "+00:00")
    try:
        return datetime.fromisoformat(iso_text).timestamp()
    except ValueError:
        return None


def _confidence_from_source(source: str, coverage: float, unique_ratio: float) -> str:
    if source == "timepoint":
        if coverage >= 0.8 and unique_ratio >= 0.8:
            return "high"
        if coverage >= 0.5:
            return "medium"
        return "low"

    if source == "collection_date":
        if coverage >= 0.8 and unique_ratio >= 0.8:
            return "high"
        if coverage >= 0.5:
            return "medium"
        return "low"

    if source in {"first_public", "last_updated"}:
        if coverage >= 0.8:
            return "medium"
        if coverage >= 0.5:
            return "low"
        return "none"

    return "none"


def _choose_source(rows: list[dict[str, str]]) -> tuple[str, list[RankedSample], float, float]:
    source_specs = [
        ("timepoint", _parse_timepoint),
        ("collection_date", _parse_date_like),
        ("first_public", _parse_date_like),
        ("last_updated", _parse_date_like),
    ]

    best_source = "none"
    best_ranked: list[RankedSample] = []
    best_coverage = 0.0
    best_unique_ratio = 0.0

    for source, parser in source_specs:
        ranked: list[RankedSample] = []
        parsed_count = 0
        for idx, row in enumerate(rows):
            parsed = parser(row.get(source, ""))
            if parsed is not None:
                parsed_count += 1
            ranked.append(
                RankedSample(
                    index=idx,
                    sample_id=_clean(row.get("sample_id", "")),
                    study=_clean(row.get("study", "")),
                    parsed_value=parsed,
                )
            )

        if not rows:
            continue

        coverage = parsed_count / len(rows)
        if parsed_count > 0:
            unique_values = {r.parsed_value for r in ranked if r.parsed_value is not None}
            unique_ratio = len(unique_values) / parsed_count
        else:
            unique_ratio = 0.0

        if parsed_count >= 2:
            best_source = source
            best_ranked = ranked
            best_coverage = coverage
            best_unique_ratio = unique_ratio
            break

        if coverage > best_coverage:
            best_source = source
            best_ranked = ranked
            best_coverage = coverage
            best_unique_ratio = unique_ratio

    return best_source, best_ranked, best_coverage, best_unique_ratio


def _assign_order(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_study: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        study = _clean(row.get("study", "")) or "unknown_study"
        by_study.setdefault(study, []).append(row)

    out_rows: list[dict[str, str]] = []

    for study, study_rows in sorted(by_study.items(), key=lambda item: item[0]):
        source, ranked, coverage, unique_ratio = _choose_source(study_rows)
        confidence = _confidence_from_source(source, coverage, unique_ratio)

        # Rank all rows. Missing parsed values are placed after known ones.
        sortable: list[tuple[int, float, str, int]] = []
        for i, item in enumerate(ranked):
            missing_flag = 1 if item.parsed_value is None else 0
            numeric_value = float("inf") if item.parsed_value is None else item.parsed_value
            sortable.append((missing_flag, numeric_value, item.sample_id, i))

        sortable.sort()
        assigned_order: dict[int, int] = {}
        for rank_pos, (_, _, _, original_index) in enumerate(sortable, start=1):
            assigned_order[original_index] = rank_pos

        for idx, row in enumerate(study_rows):
            sample_id = _clean(row.get("sample_id", ""))
            parsed_value = ranked[idx].parsed_value if ranked else None
            has_parsed = parsed_value is not None

            if source == "none" or (source != "none" and not has_parsed and coverage < 0.5):
                proposed_order = ""
                usable = "no"
            else:
                proposed_order = str(assigned_order.get(idx, ""))
                usable = "yes" if confidence in {"high", "medium"} and has_parsed else "no"

            notes = (
                f"source={source}; coverage={coverage:.2f}; unique_ratio={unique_ratio:.2f}"
                if source != "none"
                else "no usable ordering field found"
            )

            out_rows.append(
                {
                    "sample_id": sample_id,
                    "study": study,
                    "proposed_order": proposed_order,
                    "order_source": source,
                    "confidence": confidence,
                    "is_usable_for_lag": usable,
                    "notes": notes,
                }
            )

    return out_rows


def main() -> int:
    args = parse_args()
    metadata_path = Path(args.metadata)
    out_path = Path(args.out)

    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    required = {"sample_id", "study"}
    missing = [col for col in required if col not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"Missing required columns in metadata: {', '.join(missing)}")

    out_rows = _assign_order(rows)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "sample_id",
            "study",
            "proposed_order",
            "order_source",
            "confidence",
            "is_usable_for_lag",
            "notes",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    usable_count = sum(1 for row in out_rows if row["is_usable_for_lag"] == "yes")
    print(f"Wrote {len(out_rows)} rows to {out_path}")
    print(f"Usable lag-order rows: {usable_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
