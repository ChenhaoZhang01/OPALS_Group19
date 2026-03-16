#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
import urllib.error
from dataclasses import dataclass

ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FIELDS = [
    "run_accession",
    "study_accession",
    "country",
    "collection_date",
    "description",
    "sample_description",
    "sample_alias",
    "study_title",
    "read_count",
    "base_count",
]

TARGETS = {
    "wastewater": 40,
    "soil": 25,
    "river": 20,
    "irrigation": 15,
}

KEYWORDS = {
    "wastewater": ["wastewater", "wwtp", "sewage", "effluent", "sludge"],
    "soil": ["soil", "rhizosphere", "topsoil", "farmland soil"],
    "river": ["river", "stream", "freshwater", "surface water"],
    "irrigation": ["irrigation", "irrigated", "canal water", "agricultural water"],
}

ENVIRONMENT_TERMS = {
    "wastewater": ["wastewater", "sewage", "effluent", "sludge", "wwtp"],
    "soil": ["soil", "topsoil", "rhizosphere", "farmland soil"],
    "river": ["river", "stream", "freshwater", "surface water"],
    "irrigation": ["irrigation", "irrigated", "canal water", "agricultural water"],
}

TREATMENT_HINTS = [
    "upgraded",
    "treated",
    "untreated",
    "control",
    "amended",
    "fertilized",
    "manure",
]


@dataclass
class Record:
    sample_id: str
    environment: str
    country: str
    year: str
    study: str
    timepoint: str
    treatment: str
    read_count: str


def classify_environment(text: str) -> str | None:
    text_l = text.lower()
    for env, words in KEYWORDS.items():
        if any(word in text_l for word in words):
            return env
    return None


def extract_year(collection_date: str) -> str:
    if not collection_date:
        return "NA"
    match = re.search(r"(19|20)\\d{2}", collection_date)
    return match.group(0) if match else "NA"


def extract_timepoint(text: str) -> str:
    patterns = [
        r"\\bt\\s*\\d+\\b",
        r"\\btime\\s*point\\s*\\d+\\b",
        r"\\bday\\s*\\d+\\b",
        r"\\bweek\\s*\\d+\\b",
        r"\\bmonth\\s*\\d+\\b",
    ]
    text_l = text.lower()
    for pattern in patterns:
        match = re.search(pattern, text_l)
        if match:
            return match.group(0).replace(" ", "")
    return "NA"


def extract_treatment(text: str) -> str:
    text_l = text.lower()
    for hint in TREATMENT_HINTS:
        if hint in text_l:
            return hint
    return "NA"


def human_read_count(raw_read_count: str) -> str:
    if not raw_read_count or not raw_read_count.isdigit():
        return "NA"
    reads = int(raw_read_count)
    if reads >= 1_000_000:
        return f"{round(reads / 1_000_000)}M"
    if reads >= 1_000:
        return f"{round(reads / 1_000)}K"
    return str(reads)


def fetch_rows(query: str, limit: int) -> list[dict[str, str]]:
    params = {
        "result": "read_run",
        "query": query,
        "fields": ",".join(ENA_FIELDS),
        "format": "json",
        "limit": str(limit),
    }
    url = f"{ENA_API_URL}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(url, headers={"User-Agent": "ARG-metadata-builder/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read().decode("utf-8")
        return json.loads(payload)
    except urllib.error.HTTPError as exc:
        detail = ""
        try:
            detail = exc.read().decode("utf-8")
        except Exception:
            detail = ""
        print(f"WARNING: ENA request failed ({exc.code})")
        if detail:
            print(detail)
        return []


def make_query(term: str) -> str:
    token = term.replace('"', "")
    return (
        'library_source="METAGENOMIC" AND '
        f'(description="*{token}*" OR '
        f'sample_description="*{token}*" OR '
        f'study_title="*{token}*" OR '
        f'sample_alias="*{token}*")'
    )


def build_records(per_query_limit: int, sleep_seconds: float) -> list[Record]:
    seen_runs: set[str] = set()
    records: list[Record] = []

    for environment, target_count in TARGETS.items():
        needed = target_count
        terms = ENVIRONMENT_TERMS.get(environment, [environment])

        for term in terms:
            if needed <= 0:
                break

            query = make_query(term)
            rows = fetch_rows(query=query, limit=per_query_limit)
            if not rows:
                if sleep_seconds > 0:
                    time.sleep(sleep_seconds)
                continue

            for row in rows:
                run_id = (row.get("run_accession") or "").strip()
                if not run_id or run_id in seen_runs:
                    continue

                text_fields = [
                    row.get("description") or "",
                    row.get("sample_description") or "",
                    row.get("sample_alias") or "",
                    row.get("study_title") or "",
                ]
                merged_text = " ".join(text_fields)

                inferred_environment = classify_environment(merged_text)
                if inferred_environment and inferred_environment != environment:
                    continue

                record = Record(
                    sample_id=run_id,
                    environment=environment,
                    country=(row.get("country") or "NA").strip() or "NA",
                    year=extract_year(row.get("collection_date") or ""),
                    study=(row.get("study_accession") or "NA").strip() or "NA",
                    timepoint=extract_timepoint(merged_text),
                    treatment=extract_treatment(merged_text),
                    read_count=human_read_count((row.get("read_count") or "").strip()),
                )

                records.append(record)
                seen_runs.add(run_id)
                needed -= 1

                if needed <= 0:
                    break

            if sleep_seconds > 0:
                time.sleep(sleep_seconds)

        if needed > 0:
            print(f"WARNING: {environment} shortfall: missing {needed} samples")

    return records


def write_csv(records: list[Record], output_csv: str) -> None:
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "sample_id",
                "environment",
                "country",
                "year",
                "study",
                "timepoint",
                "treatment",
                "read_count",
            ]
        )
        for record in records:
            writer.writerow(
                [
                    record.sample_id,
                    record.environment,
                    record.country,
                    record.year,
                    record.study,
                    record.timepoint,
                    record.treatment,
                    record.read_count,
                ]
            )


def summarize(records: list[Record]) -> None:
    counts = {env: 0 for env in TARGETS}
    for record in records:
        counts[record.environment] = counts.get(record.environment, 0) + 1

    print("Collected records:")
    for env in ("wastewater", "soil", "river", "irrigation"):
        print(f"  - {env}: {counts.get(env, 0)} / {TARGETS[env]}")
    print(f"Total: {len(records)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build master_metadata.csv from ENA metadata")
    parser.add_argument(
        "--output",
        default="metadata/master_metadata.csv",
        help="Output CSV path (default: metadata/master_metadata.csv)",
    )
    parser.add_argument(
        "--per-query-limit",
        type=int,
        default=5000,
        help="Maximum ENA rows returned per term query (default: 5000)",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.2,
        help="Seconds to sleep between ENA requests (default: 0.2)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    records = build_records(
        per_query_limit=args.per_query_limit,
        sleep_seconds=args.sleep,
    )
    write_csv(records=records, output_csv=args.output)
    summarize(records)

    if len(records) < sum(TARGETS.values()):
        print("WARNING: Did not reach full target counts. Consider increasing --per-query-limit.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
