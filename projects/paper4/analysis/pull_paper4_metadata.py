#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import time
import xml.etree.ElementTree as ET
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError
from pathlib import Path


ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/filereport"
SRA_RUNINFO_BASE = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"
BIOSAMPLE_EFETCH_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pull ENA metadata for Paper 4 wastewater samples and prefill DiD table"
    )
    parser.add_argument(
        "--metadata",
        default="projects/paper4/metadata/metadata_final.csv",
        help="Input metadata file with sample_id/environment",
    )
    parser.add_argument(
        "--pulled-out",
        default="projects/paper4/results/paper4_sample_metadata_pulled.csv",
        help="Output CSV with enriched sample metadata",
    )
    parser.add_argument(
        "--did-out",
        default="projects/paper4/results/did_input_table.csv",
        help="Output CSV for DiD input table",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=40,
        help="Number of run accessions per ENA request",
    )
    return parser.parse_args()


def read_wastewater_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            if (row.get("environment") or "").strip().lower() != "wastewater":
                continue
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                continue
            rows.append(row)
    return rows


def chunked(items: list[str], size: int) -> list[list[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def fetch_ena_rows(run_ids: list[str]) -> dict[str, dict[str, str]]:
    if not run_ids:
        return {}

    fields = [
        "run_accession",
        "study_accession",
        "sample_accession",
        "experiment_accession",
        "sample_alias",
        "sample_title",
        "collection_date",
        "location",
        "country",
        "submitted_ftp",
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
        payload = response.read().decode("utf-8")

    # Avoid importing json at module load in environments where this endpoint is unavailable.
    import json

    parsed = json.loads(payload)

    result: dict[str, dict[str, str]] = {}
    for row in parsed:
        run = (row.get("run_accession") or "").strip()
        if run:
            result[run] = {k: (v or "") for k, v in row.items()}
    return result


def fetch_sra_runinfo_rows(run_ids: list[str]) -> dict[str, dict[str, str]]:
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
    result: dict[str, dict[str, str]] = {}
    for row in reader:
        run = (row.get("Run") or "").strip()
        if run:
            result[run] = {k: (v or "") for k, v in row.items()}
    return result


def fetch_biosample_attributes(biosample_acc: str) -> dict[str, str]:
    biosample_acc = (biosample_acc or "").strip()
    if not biosample_acc:
        return {}

    params = {
        "db": "biosample",
        "id": biosample_acc,
        "retmode": "xml",
    }
    url = f"{BIOSAMPLE_EFETCH_BASE}?{urllib.parse.urlencode(params)}"
    xml_text = ""
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                xml_text = response.read().decode("utf-8", "ignore")
            break
        except HTTPError as exc:
            if exc.code == 429 and attempt < 4:
                time.sleep(1.5 * (attempt + 1))
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
        key = (node.attrib.get("attribute_name") or "").strip().lower()
        val = (node.text or "").strip()
        if key and val and key not in attrs:
            attrs[key] = val
    return attrs


def derive_plant_id(country_fallback: str, ena_row: dict[str, str]) -> str:
    # Prefer site-like free-text fields that distinguish plants, then clean to a stable ID.
    candidate = (
        ena_row.get("location")
        or ena_row.get("sample_alias")
        or ena_row.get("sample_title")
        or ena_row.get("country")
        or country_fallback
    )
    candidate = candidate.strip()
    if not candidate:
        return "UNKNOWN_PLANT"

    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", candidate).strip("_")
    cleaned = re.sub(r"_+", "_", cleaned)
    if not cleaned:
        return "UNKNOWN_PLANT"
    return cleaned[:80]


def choose_collection_date(sra_row: dict[str, str], biosample_attrs: dict[str, str], ena_row: dict[str, str]) -> str:
    candidates = [
        biosample_attrs.get("collection date", ""),
        biosample_attrs.get("collection_date", ""),
        ena_row.get("collection_date", ""),
        sra_row.get("ReleaseDate", ""),
    ]
    for value in candidates:
        value = value.strip()
        if value:
            return value
    return ""


def choose_location(country_fallback: str, biosample_attrs: dict[str, str], ena_row: dict[str, str]) -> str:
    candidates = [
        biosample_attrs.get("geo_loc_name", ""),
        biosample_attrs.get("geographic location (country and/or sea)", ""),
        biosample_attrs.get("geographic location", ""),
        ena_row.get("location", ""),
        ena_row.get("country", ""),
        country_fallback,
    ]
    for value in candidates:
        value = value.strip()
        if value:
            return value
    return ""


def write_pulled_output(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "study",
        "country_repo",
        "country_ena",
        "collection_date",
        "location",
        "sample_alias",
        "sample_title",
        "first_public",
        "last_updated",
        "plant_id_derived",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_did_prefill(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "ARG_total", "treatment", "time", "plant_id"])
        for row in rows:
            writer.writerow([row["sample_id"], "", "", "", row["plant_id_derived"]])


def main() -> int:
    args = parse_args()
    metadata_path = Path(args.metadata)
    pulled_out = Path(args.pulled_out)
    did_out = Path(args.did_out)

    wastewater_rows = read_wastewater_rows(metadata_path)
    run_ids = sorted({(r.get("sample_id") or "").strip() for r in wastewater_rows if r.get("sample_id")})

    pulled_ena_by_run: dict[str, dict[str, str]] = {}
    pulled_sra_by_run: dict[str, dict[str, str]] = {}
    for batch in chunked(run_ids, max(1, args.batch_size)):
        pulled_ena_by_run.update(fetch_ena_rows(batch))
        pulled_sra_by_run.update(fetch_sra_runinfo_rows(batch))

    biosample_by_run: dict[str, dict[str, str]] = {}
    for run in run_ids:
        sra_row = pulled_sra_by_run.get(run, {})
        biosample_acc = sra_row.get("BioSample", "")
        biosample_by_run[run] = fetch_biosample_attributes(biosample_acc)

    out_rows: list[dict[str, str]] = []
    for row in wastewater_rows:
        run = (row.get("sample_id") or "").strip()
        ena_row = pulled_ena_by_run.get(run, {})
        sra_row = pulled_sra_by_run.get(run, {})
        biosample_attrs = biosample_by_run.get(run, {})
        picked_location = choose_location(row.get("country", ""), biosample_attrs, ena_row)
        out_rows.append(
            {
                "sample_id": run,
                "study": row.get("study", ""),
                "country_repo": row.get("country", ""),
                "country_ena": ena_row.get("country", ""),
                "collection_date": choose_collection_date(sra_row, biosample_attrs, ena_row),
                "location": picked_location,
                "sample_alias": ena_row.get("sample_alias", ""),
                "sample_title": ena_row.get("sample_title", ""),
                "first_public": ena_row.get("first_public", ""),
                "last_updated": sra_row.get("LoadDate", "") or ena_row.get("last_updated", ""),
                "plant_id_derived": derive_plant_id(
                    row.get("country", ""),
                    {"location": picked_location, "sample_alias": "", "sample_title": "", "country": row.get("country", "")},
                ),
            }
        )

    write_pulled_output(pulled_out, out_rows)
    write_did_prefill(did_out, out_rows)

    missing_ena = sum(1 for r in out_rows if not r["collection_date"] and not r["location"] and not r["sample_alias"])
    print(f"Wastewater sample rows: {len(out_rows)}")
    print(f"Rows with sparse ENA metadata: {missing_ena}")
    print(f"Wrote: {pulled_out}")
    print(f"Wrote: {did_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
