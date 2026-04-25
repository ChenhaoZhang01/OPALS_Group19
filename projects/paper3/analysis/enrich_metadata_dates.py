#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import time
import xml.etree.ElementTree as ET
import urllib.parse
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError


ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/filereport"
SRA_RUNINFO_BASE = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"
BIOSAMPLE_EFETCH_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Enrich Paper 3 metadata with temporal fields from ENA/SRA")
    parser.add_argument(
        "--metadata",
        default="projects/paper3/metadata/metadata_final.csv",
        help="Input metadata CSV with sample_id and study",
    )
    parser.add_argument(
        "--out",
        default="projects/paper3/metadata/metadata_enriched.csv",
        help="Output enriched metadata CSV",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=40,
        help="Run accessions per API request",
    )
    return parser.parse_args()


def _clean(value: str | None) -> str:
    return (value or "").strip()


def _chunks(items: list[str], size: int) -> list[list[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def fetch_ena_rows(run_ids: list[str]) -> dict[str, dict[str, str]]:
    if not run_ids:
        return {}

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
        "accession": ",".join(run_ids),
        "result": "read_run",
        "fields": ",".join(fields),
        "format": "json",
    }
    url = f"{ENA_BASE}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=60) as response:
        parsed = json.loads(response.read().decode("utf-8", "ignore"))

    by_run: dict[str, dict[str, str]] = {}
    for row in parsed:
        run = _clean(row.get("run_accession"))
        if run:
            by_run[run] = {k: _clean(v) for k, v in row.items()}
    return by_run


def fetch_sra_rows(run_ids: list[str]) -> dict[str, dict[str, str]]:
    if not run_ids:
        return {}

    params = {"acc": ",".join(run_ids)}
    url = f"{SRA_RUNINFO_BASE}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8", "ignore")

    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return {}

    reader = csv.DictReader(lines)
    by_run: dict[str, dict[str, str]] = {}
    for row in reader:
        run = _clean(row.get("Run"))
        if run:
            by_run[run] = {k: _clean(v) for k, v in row.items()}
    return by_run


def fetch_biosample_attributes(biosample_acc: str) -> dict[str, str]:
    biosample_acc = _clean(biosample_acc)
    if not biosample_acc:
        return {}

    params = {
        "db": "biosample",
        "id": biosample_acc,
        "retmode": "xml",
    }
    url = f"{BIOSAMPLE_EFETCH_BASE}?{urllib.parse.urlencode(params)}"
    xml_text = ""
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                xml_text = response.read().decode("utf-8", "ignore")
            break
        except HTTPError as exc:
            if exc.code == 429 and attempt < 3:
                time.sleep(1.0 * (attempt + 1))
                continue
            return {}
        except URLError:
            return {}

    attrs: dict[str, str] = {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return attrs

    for node in root.findall(".//Attribute"):
        key = _clean(node.attrib.get("attribute_name", "")).lower()
        val = _clean(node.text)
        if key and val and key not in attrs:
            attrs[key] = val
    return attrs


def choose_temporal_source(
    ena_row: dict[str, str],
    sra_row: dict[str, str],
    biosample_attrs: dict[str, str],
) -> tuple[str, str, str]:
    biosample_collection = _clean(
        biosample_attrs.get("collection date") or biosample_attrs.get("collection_date")
    )
    ena_collection = _clean(ena_row.get("collection_date"))
    first_public = _clean(ena_row.get("first_public"))
    sra_release = _clean(sra_row.get("ReleaseDate"))
    last_updated = _clean(ena_row.get("last_updated"))

    candidates = [
        ("biosample_collection_date", biosample_collection),
        ("ena_collection_date", ena_collection),
        ("ena_first_public", first_public),
        ("sra_release_date", sra_release),
        ("ena_last_updated", last_updated),
    ]
    for source, value in candidates:
        if value:
            # Keep collection_date strict: only true collection date fields.
            collection_date = biosample_collection or ena_collection
            return collection_date, source, value
    return "", "none", ""


def main() -> int:
    args = parse_args()
    in_path = Path(args.metadata)
    out_path = Path(args.out)

    with in_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        original_fields = list(reader.fieldnames or [])

    if "sample_id" not in original_fields:
        raise ValueError("Input metadata must contain sample_id")

    run_ids = sorted({_clean(r.get("sample_id")) for r in rows if _clean(r.get("sample_id"))})

    ena_by_run: dict[str, dict[str, str]] = {}
    sra_by_run: dict[str, dict[str, str]] = {}
    biosample_cache: dict[str, dict[str, str]] = {}
    for batch in _chunks(run_ids, max(1, args.batch_size)):
        ena_by_run.update(fetch_ena_rows(batch))
        sra_by_run.update(fetch_sra_rows(batch))

    enriched_fields = [
        "collection_date",
        "first_public",
        "last_updated",
        "sample_alias",
        "sample_title",
        "location",
        "country_ena",
        "study_accession",
        "sample_accession",
        "biosample_accession",
        "biosample_collection_date",
        "sra_release_date",
        "temporal_source",
    ]

    fieldnames = original_fields + [f for f in enriched_fields if f not in original_fields]

    matched_ena = 0
    matched_sra = 0
    filled_collection = 0
    biosample_hits = 0

    for row in rows:
        run = _clean(row.get("sample_id"))
        ena = ena_by_run.get(run, {})
        sra = sra_by_run.get(run, {})
        biosample_acc = _clean(sra.get("BioSample"))
        biosample_attrs = {}
        if biosample_acc:
            if biosample_acc not in biosample_cache:
                biosample_cache[biosample_acc] = fetch_biosample_attributes(biosample_acc)
            biosample_attrs = biosample_cache[biosample_acc]
            if biosample_attrs:
                biosample_hits += 1

        if ena:
            matched_ena += 1
        if sra:
            matched_sra += 1

        collection_date, temporal_source, _temporal_value = choose_temporal_source(ena, sra, biosample_attrs)
        first_public = _clean(ena.get("first_public"))
        last_updated = _clean(ena.get("last_updated"))
        sra_release = _clean(sra.get("ReleaseDate"))
        biosample_collection = _clean(biosample_attrs.get("collection date") or biosample_attrs.get("collection_date"))
        if collection_date:
            filled_collection += 1

        row["collection_date"] = collection_date or _clean(row.get("collection_date"))
        row["first_public"] = first_public
        row["last_updated"] = last_updated
        row["sample_alias"] = _clean(ena.get("sample_alias"))
        row["sample_title"] = _clean(ena.get("sample_title"))
        row["location"] = _clean(ena.get("location"))
        row["country_ena"] = _clean(ena.get("country"))
        row["study_accession"] = _clean(ena.get("study_accession"))
        row["sample_accession"] = _clean(ena.get("sample_accession"))
        row["biosample_accession"] = biosample_acc
        row["biosample_collection_date"] = biosample_collection
        row["sra_release_date"] = sra_release
        row["temporal_source"] = temporal_source

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Input rows: {len(rows)}")
    print(f"ENA matched runs: {matched_ena}")
    print(f"SRA matched runs: {matched_sra}")
    print(f"BioSample rows with attributes: {biosample_hits}")
    print(f"Rows with collection date (ENA or BioSample): {filled_collection}")
    print(f"Wrote enriched metadata: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
