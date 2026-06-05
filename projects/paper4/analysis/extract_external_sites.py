#!/usr/bin/env python3
"""Extract real multi-site ARG influent/effluent relative abundances from collected papers.

Currently parses Table 4 of PMC11471163 (Brindova/Hghaikova et al., Slovakia & Taiwan
resistome study; supplementary mmc2.docx), which reports community-normalized ARG relative
abundance in influent and effluent for sul1, tetM, ermB, blaTEM by shotgun metagenomics
and by qPCR. These relative abundances are commensurable with the 16S-normalized qPCR
removal computed for the Bliesen site (Paper 4 anchor), enabling a cross-site comparison of
*selective* (community-normalized) ARG removal through conventional treatment.

Output: projects/paper4/data/external_sites_removal.csv  (tidy, with provenance)
"""

from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path

MMC2 = ("projects/paper4/data/external/oa_scan/PMC11471163/PMC11471163/mmc2.docx")
OUT = "projects/paper4/data/external_sites_removal.csv"
SOURCE = "PMC11471163 (Slovakia & Taiwan resistome study), Supplementary Table 4"

# ARG class harmonization (shared with Bliesen genes blaAmpC, ermB, sul1, sul2, tetW)
CLASS = {"Sul1": "sulfonamide", "tetM": "tetracycline", "ermB": "MLS", "blaTEM": "beta-lactam"}


def _num(cell: str) -> float | None:
    """Parse 'A x 10-B' style cells (× lost in encoding) -> float. '0.0...10-0' -> 0."""
    s = cell.strip()
    m = re.match(r"^\s*([\d.]+)\D+10-?(\d+)\s*$", s)
    if not m:
        try:
            return float(s)
        except ValueError:
            return None
    mant = float(m.group(1))
    exp = int(m.group(2))
    # encoding drops the sign; in this table all exponents are negative powers of ten
    # except '10-0'/'100' meaning 10^0. Reconstruct: 'A 10-B' -> A*10^-B ; 'A 100' -> A.
    if "10-" in s or re.search(r"10-?\d", s):
        # if original had a hyphen it's negative; '100' (no hyphen, exp captured as 0 via '10'+'0')
        neg = "-" in s.split("10", 1)[1][:1]
        return mant * (10 ** (-exp if neg else exp))
    return mant * (10 ** exp)


def _cells(row: str) -> list[str]:
    out = []
    for c in re.findall(r"<w:tc[ >].*?</w:tc>", row, re.S):
        txt = re.sub(r"<[^>]+>", "", "".join(re.findall(r"<w:t[ >][^<]*</w:t>", c)))
        out.append(txt.strip())
    return out


def main() -> int:
    xml = zipfile.ZipFile(MMC2).read("word/document.xml").decode("utf-8", "ignore")
    tables = re.findall(r"<w:tbl>.*?</w:tbl>", xml, re.S)
    # find the table whose flattened text mentions both metagenomic and qPCR + Sul1
    target = None
    for tb in tables:
        flat = re.sub(r"<[^>]+>", " ", tb).lower()
        if "sul1" in flat and "qpcr" in flat and "metagenom" in flat:
            target = tb
            break
    if target is None:
        raise SystemExit("Table 4 not found")

    rows = [_cells(r) for r in re.findall(r"<w:tr[ >].*?</w:tr>", target, re.S)]
    rows = [r for r in rows if any(c for c in r)]

    # Layout (from inspection):
    # ['Method','Gene','Taiwan','','Slovakia']
    # ['','','Influent','Effluent','p','','Influent','Effluent','p']
    # ['Shotgun metagenomic','Sul1', TW_in,TW_eff,p, '', SVK_in,SVK_eff,p] ...
    recs = []
    method = ""
    for r in rows:
        first = r[0].strip()
        if first.lower().startswith("shotgun"):
            method = "metagenomic"
        elif first.lower().startswith("qpcr"):
            method = "qPCR"
        gene = r[1].strip() if len(r) > 1 else ""
        if gene not in CLASS:
            continue
        # Taiwan influent/effluent = cols 2,3 ; Slovakia = cols 6,7
        try:
            tw_in, tw_eff = _num(r[2]), _num(r[3])
            svk_in, svk_eff = _num(r[6]), _num(r[7])
        except IndexError:
            continue
        for country, inf, eff in [("Taiwan", tw_in, tw_eff), ("Slovakia", svk_in, svk_eff)]:
            recs.append({
                "study": "Brindova_SlovakiaTaiwan", "country": country, "site": country,
                "treatment": "conventional A/O", "method": method,
                "gene": gene, "arg_class": CLASS[gene],
                "influent_relabund": inf, "effluent_relabund": eff,
                "source": SOURCE,
            })

    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()))
        w.writeheader()
        w.writerows(recs)
    print(f"Wrote {OUT} ({len(recs)} records)")
    for r in recs:
        print(f"  {r['country']:9} {r['method']:11} {r['gene']:7} "
              f"in={r['influent_relabund']:.2e} eff={r['effluent_relabund']:.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
