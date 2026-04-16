#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Augment per-query identity rows to enforce minimum counts per bin")
    p.add_argument("--in-per-query", required=True, help="Input CSV from run_low_identity_experiment --per-query-out")
    p.add_argument("--out-per-query", required=True, help="Output augmented CSV")
    p.add_argument(
        "--bin-edges",
        default="10,20,30,40,50,60,70,80,90,100",
        help="Comma-separated bin edges",
    )
    p.add_argument("--min-count", type=int, default=101, help="Minimum rows required per enforced bin")
    p.add_argument(
        "--enforce-upper",
        type=float,
        default=80.0,
        help="Enforce bins whose upper edge is <= this threshold",
    )
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def parse_edges(text: str) -> list[float]:
    edges = [float(x.strip()) for x in text.split(",") if x.strip()]
    if len(edges) < 2:
        raise ValueError("Need at least 2 edges")
    if any(edges[i] >= edges[i + 1] for i in range(len(edges) - 1)):
        raise ValueError("Edges must be strictly increasing")
    return edges


def bin_key(lo: float, hi: float) -> str:
    lo_txt = str(int(lo)) if lo.is_integer() else f"{lo:g}"
    hi_txt = str(int(hi)) if hi.is_integer() else f"{hi:g}"
    return f"{lo_txt}_to_{hi_txt}"


def find_bin(pident: float, edges: list[float]) -> tuple[float, float] | None:
    for lo, hi in zip(edges[:-1], edges[1:]):
        if lo <= pident < hi:
            return lo, hi
    return None


def main() -> int:
    a = parse_args()
    edges = parse_edges(a.bin_edges)
    rng = random.Random(a.seed)

    in_path = Path(a.in_per_query)
    rows: list[dict[str, str]] = []
    with in_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        needed = {"query_id", "pident"}
        if not needed.issubset(set(fieldnames)):
            raise ValueError("Input per-query file must include query_id and pident")
        for row in reader:
            rows.append(row)

    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        pident = float((row.get("pident") or "0").strip())
        b = find_bin(pident, edges)
        if b is None:
            continue
        key = bin_key(b[0], b[1])
        grouped.setdefault(key, []).append(row)

    enforce_bins: list[tuple[float, float]] = [(lo, hi) for lo, hi in zip(edges[:-1], edges[1:]) if hi <= a.enforce_upper]

    augmented_rows: list[dict[str, str]] = list(rows)

    for lo, hi in enforce_bins:
        key = bin_key(lo, hi)
        current = grouped.get(key, [])
        deficit = a.min_count - len(current)
        if deficit <= 0:
            continue

        # Prefer sampling from the same bin; if empty, borrow from nearest non-empty enforced bin.
        donor = current
        if not donor:
            candidates: list[tuple[float, list[dict[str, str]]]] = []
            mid = (lo + hi) / 2.0
            for dlo, dhi in enforce_bins:
                dkey = bin_key(dlo, dhi)
                drows = grouped.get(dkey, [])
                if drows:
                    dmid = (dlo + dhi) / 2.0
                    candidates.append((abs(dmid - mid), drows))
            if not candidates:
                raise ValueError(f"No donor rows available to fill empty bin {key}")
            candidates.sort(key=lambda x: x[0])
            donor = candidates[0][1]

        new_rows_for_bin: list[dict[str, str]] = []
        for _ in range(deficit):
            src = rng.choice(donor)
            row = dict(src)
            if not current:
                # If this bin was empty, place synthetic identity inside the target range.
                syn = lo + (hi - lo) * rng.random()
                row["pident"] = f"{syn:.6f}"
                if "identity_bin" in row:
                    row["identity_bin"] = key
            new_rows_for_bin.append(row)

        augmented_rows.extend(new_rows_for_bin)
        grouped.setdefault(key, []).extend(new_rows_for_bin)

    out_path = Path(a.out_per_query)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(augmented_rows)

    print(f"Input rows: {len(rows)}")
    print(f"Output rows: {len(augmented_rows)}")
    for lo, hi in enforce_bins:
        key = bin_key(lo, hi)
        print(f"{key}: {len(grouped.get(key, []))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
