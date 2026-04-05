#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict

NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build Paper 4 DiD table from CEBEDEAU Supplementary Table S10 qPCR ARG data"
    )
    parser.add_argument(
        "--xlsx",
        required=True,
        help="Path to Supplementary-tables-CEBEDEAU.xlsx",
    )
    parser.add_argument(
        "--did-out",
        required=True,
        help="Output CSV for DiD input table",
    )
    parser.add_argument(
        "--long-out",
        required=True,
        help="Output CSV for long-form extracted S10 rows",
    )
    parser.add_argument(
        "--notes-out",
        required=True,
        help="Output TXT describing coding assumptions",
    )
    return parser.parse_args()


def _shared_strings(zf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    out: list[str] = []
    for si in root.findall(".//a:si", NS):
        txt = "".join(t.text or "" for t in si.findall(".//a:t", NS))
        out.append(txt)
    return out


def _sheet_path_for_name(zf: zipfile.ZipFile, sheet_name: str) -> str:
    wb = ET.fromstring(zf.read("xl/workbook.xml"))
    rel = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rid_map = {r.attrib["Id"]: r.attrib["Target"] for r in rel}

    for sh in wb.findall(".//a:sheets/a:sheet", NS):
        if sh.attrib.get("name") != sheet_name:
            continue
        rid = sh.attrib.get(REL_NS)
        if not rid:
            break
        path = "xl/" + rid_map[rid]
        path = path.replace("xl/xl/", "xl/")
        return path
    raise ValueError(f"Sheet not found: {sheet_name}")


def _cell_value(cell: ET.Element, shared: list[str]) -> str:
    t = cell.attrib.get("t")
    v = cell.find("a:v", NS)
    if v is None:
        return ""
    txt = v.text or ""
    if t == "s" and txt:
        try:
            return shared[int(txt)]
        except (ValueError, IndexError):
            return txt
    return txt


def load_s10_rows(xlsx_path: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with zipfile.ZipFile(xlsx_path, "r") as zf:
        shared = _shared_strings(zf)
        sheet_path = _sheet_path_for_name(zf, "Table S10")
        sh = ET.fromstring(zf.read(sheet_path))

        header: list[str] | None = None
        for row in sh.findall(".//a:sheetData/a:row", NS):
            vals = [_cell_value(c, shared).strip() for c in row.findall("a:c", NS)]
            if not vals:
                continue
            # Header row has canonical columns used below.
            if vals[:4] == ["Sample", "Campaign", "Gene", "Concentration (copies/L)"]:
                header = vals[:4]
                continue
            if header is None:
                continue
            if len(vals) < 4:
                continue
            sample, campaign, gene, concentration = vals[:4]
            if not sample or not campaign or not gene or not concentration:
                continue
            rows.append(
                {
                    "Sample": sample,
                    "Campaign": campaign,
                    "Gene": gene,
                    "Concentration (copies/L)": concentration,
                }
            )
    return rows


def build_did_rows(long_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    # Keep treatment-line comparison with OUT WWTP as control.
    treatment_map = {
        "OUT WWTP": 0,
        "AOP": 1,
        "CW": 1,
        "GAC": 1,
    }
    time_map = {
        "C2-I": 0,
        "C5-I": 0,
        "C10-I": 1,
        "C12-I": 1,
    }

    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in long_rows:
        sample = row["Sample"].strip()
        campaign = row["Campaign"].strip()
        gene = row["Gene"].strip()
        if sample not in treatment_map or campaign not in time_map:
            continue
        try:
            value = float(row["Concentration (copies/L)"])
        except ValueError:
            continue
        grouped[(sample, campaign, gene)].append(value)

    # Aggregate replicates per gene first, then sum gene means into ARG_total.
    by_sample_campaign: dict[tuple[str, str], list[float]] = defaultdict(list)
    for (sample, campaign, _gene), values in grouped.items():
        by_sample_campaign[(sample, campaign)].append(sum(values) / len(values))

    out: list[dict[str, str]] = []
    for (sample, campaign), gene_means in sorted(by_sample_campaign.items()):
        arg_total = sum(gene_means)
        out.append(
            {
                "sample_id": f"{sample}_{campaign}".replace(" ", "_"),
                "ARG_total": f"{arg_total:.6f}",
                "treatment": str(treatment_map[sample]),
                "time": str(time_map[campaign]),
                "plant_id": sample.replace(" ", "_"),
                "campaign": campaign,
            }
        )
    return out


def write_csv(path: str, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()

    long_rows = load_s10_rows(args.xlsx)
    write_csv(
        args.long_out,
        long_rows,
        ["Sample", "Campaign", "Gene", "Concentration (copies/L)"],
    )

    did_rows = build_did_rows(long_rows)
    write_csv(
        args.did_out,
        did_rows,
        ["sample_id", "ARG_total", "treatment", "time", "plant_id", "campaign"],
    )

    os.makedirs(os.path.dirname(args.notes_out), exist_ok=True)
    with open(args.notes_out, "w", encoding="utf-8") as f:
        f.write("Source: Supplementary Table S10 from Microorganisms 2025 (DOI: 10.3390/microorganisms13122663).\n")
        f.write("Outcome: qPCR ARG concentration (copies/L), aggregated as sum of per-gene mean concentrations by sample line and campaign.\n")
        f.write("Control group: OUT WWTP.\n")
        f.write("Treatment group: AOP, CW, GAC quaternary lines.\n")
        f.write("Pre period: C2-I, C5-I. Post period: C10-I, C12-I.\n")
        f.write("Excluded from DiD: IN WWTP rows (not a post-WWTP treatment line).\n")

    print(f"Wrote long table: {args.long_out}")
    print(f"Wrote DiD input: {args.did_out}")
    print(f"Wrote assumptions: {args.notes_out}")
    print(f"Rows (long): {len(long_rows)}")
    print(f"Rows (DiD): {len(did_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
