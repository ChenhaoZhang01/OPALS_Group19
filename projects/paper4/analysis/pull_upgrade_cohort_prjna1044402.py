#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import io
import re
import urllib.parse
import urllib.request
from pathlib import Path


ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/filereport"
BIOPROJECT = "PRJNA1044402"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build Paper 4 replacement cohort from PRJNA1044402"
    )
    parser.add_argument(
        "--metadata-out",
        default="projects/paper4/results/prjna1044402_metadata_pulled.csv",
        help="Output CSV with pulled metadata and suggested coding",
    )
    parser.add_argument(
        "--did-out",
        default="projects/paper4/results/did_input_table.csv",
        help="Output DiD table path",
    )
    parser.add_argument(
        "--notes-out",
        default="projects/paper4/results/prjna1044402_assumptions.txt",
        help="Output text note describing coding assumptions",
    )
    return parser.parse_args()


def pull_rows() -> list[dict[str, str]]:
    fields = [
        "run_accession",
        "study_accession",
        "sample_accession",
        "sample_alias",
        "sample_title",
        "collection_date",
        "location",
        "country",
        "first_public",
        "last_updated",
    ]
    params = {
        "accession": BIOPROJECT,
        "result": "read_run",
        "fields": ",".join(fields),
        "format": "tsv",
    }
    url = f"{ENA_BASE}?{urllib.parse.urlencode(params)}"
    text = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "ignore")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def alias_prefix(alias: str) -> str:
    match = re.match(r"([A-Za-z]+)", alias or "")
    return match.group(1) if match else ""


def suggest_time(alias: str) -> str:
    prefix = alias_prefix(alias)
    if prefix == "WA":
        return "0"
    if prefix == "WB":
        return "1"
    return ""


def suggest_treatment(country: str) -> str:
    c = (country or "").lower()
    if "wascana creek" in c:
        return "1"
    if "qu'appelle river" in c:
        return "0"
    return ""


def derive_plant_id(country: str, alias: str) -> str:
    location = (country or "unknown_site").strip()
    prefix = alias_prefix(alias)
    core = re.sub(r"[^A-Za-z0-9]+", "_", location).strip("_")
    core = re.sub(r"_+", "_", core)
    if prefix:
        return f"{core}_{prefix}"
    return core or "UNKNOWN_PLANT"


def write_metadata(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "study_accession",
        "sample_accession",
        "sample_alias",
        "sample_title",
        "collection_date",
        "location",
        "country",
        "first_public",
        "last_updated",
        "plant_id",
        "time_suggested",
        "treatment_suggested",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            alias = row.get("sample_alias", "")
            country = row.get("country", "")
            writer.writerow(
                {
                    "sample_id": row.get("run_accession", ""),
                    "study_accession": row.get("study_accession", ""),
                    "sample_accession": row.get("sample_accession", ""),
                    "sample_alias": alias,
                    "sample_title": row.get("sample_title", ""),
                    "collection_date": row.get("collection_date", ""),
                    "location": row.get("location", ""),
                    "country": country,
                    "first_public": row.get("first_public", ""),
                    "last_updated": row.get("last_updated", ""),
                    "plant_id": derive_plant_id(country, alias),
                    "time_suggested": suggest_time(alias),
                    "treatment_suggested": suggest_treatment(country),
                }
            )


def write_did(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "ARG_total", "treatment", "time", "plant_id"])
        for row in rows:
            alias = row.get("sample_alias", "")
            country = row.get("country", "")
            writer.writerow(
                [
                    row.get("run_accession", ""),
                    "",
                    suggest_treatment(country),
                    suggest_time(alias),
                    derive_plant_id(country, alias),
                ]
            )


def write_notes(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "Replacement cohort source: PRJNA1044402",
        "Project title contains wastewater treatment upgrade context.",
        "",
        "Suggested coding used to prefill did_input_table.csv:",
        "1) time: WA sample_alias prefix -> 0 (earlier period, 2012); WB prefix -> 1 (later period, 2018).",
        "2) treatment: country containing 'Wascana Creek' -> 1; 'Qu'Appelle River' -> 0.",
        "",
        "IMPORTANT: Verify this treatment mapping against the study manuscript/supplement before final analysis.",
        "ARG_total remains blank and must be merged from ARG quantification outputs.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()

    rows = pull_rows()
    rows = [r for r in rows if (r.get("run_accession") or "").strip()]
    rows.sort(key=lambda r: (r.get("sample_alias", ""), r.get("run_accession", "")))

    metadata_out = Path(args.metadata_out)
    did_out = Path(args.did_out)
    notes_out = Path(args.notes_out)

    write_metadata(metadata_out, rows)
    write_did(did_out, rows)
    write_notes(notes_out)

    prefix_counts: dict[str, int] = {}
    for row in rows:
        prefix = alias_prefix(row.get("sample_alias", ""))
        prefix_counts[prefix] = prefix_counts.get(prefix, 0) + 1

    print(f"BioProject: {BIOPROJECT}")
    print(f"Rows pulled: {len(rows)}")
    print(f"Alias prefix counts: {prefix_counts}")
    print(f"Wrote metadata: {metadata_out}")
    print(f"Wrote DiD table: {did_out}")
    print(f"Wrote assumptions note: {notes_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
