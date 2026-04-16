#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create train/test split by identity clusters (prevents leakage from near-duplicates)."
    )
    p.add_argument("--pairwise-identity-csv", required=True, help="CSV with query_id,subject_id,pident")
    p.add_argument("--labels-csv", required=True, help="CSV with columns query_id,label")
    p.add_argument("--identity-threshold", type=float, default=70.0, help="Cluster edge threshold in % identity")
    p.add_argument("--test-fraction", type=float, default=0.2, help="Target test fraction")
    p.add_argument("--out-assignments", required=True, help="Output CSV with query_id,label,cluster_id,split")
    p.add_argument("--out-summary", required=True, help="Output CSV summary of split sizes")
    return p.parse_args()


@dataclass
class DSU:
    parent: dict[str, str]

    def __init__(self) -> None:
        self.parent = {}

    def add(self, x: str) -> None:
        if x not in self.parent:
            self.parent[x] = x

    def find(self, x: str) -> str:
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra != rb:
            if ra < rb:
                self.parent[rb] = ra
            else:
                self.parent[ra] = rb


def read_labels(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            qid = (row.get("query_id") or "").strip()
            label = (row.get("label") or "").strip()
            if qid and label:
                out[qid] = label
    if not out:
        raise ValueError("No labeled sequences loaded from labels CSV")
    return out


def build_clusters(pairwise_path: Path, members: set[str], threshold: float) -> dict[str, str]:
    dsu = DSU()
    for mid in members:
        dsu.add(mid)

    with pairwise_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        needed = {"query_id", "subject_id", "pident"}
        if not needed.issubset(set(reader.fieldnames or [])):
            raise ValueError("pairwise identity CSV must include query_id,subject_id,pident")

        for row in reader:
            qid = (row.get("query_id") or "").strip()
            sid = (row.get("subject_id") or "").strip()
            if qid not in members or sid not in members:
                continue
            try:
                pident = float((row.get("pident") or "").strip())
            except ValueError:
                continue
            if pident >= threshold:
                dsu.union(qid, sid)

    return {mid: dsu.find(mid) for mid in members}


def assign_splits(clusters: dict[str, str], test_fraction: float) -> dict[str, str]:
    by_cluster: dict[str, list[str]] = {}
    for seq_id, cluster_id in clusters.items():
        by_cluster.setdefault(cluster_id, []).append(seq_id)

    cluster_items = sorted(by_cluster.items(), key=lambda x: (-len(x[1]), x[0]))
    total = len(clusters)
    target_test = int(round(total * test_fraction))

    split: dict[str, str] = {}
    current_test = 0
    for _, ids in cluster_items:
        put_in_test = current_test < target_test
        split_name = "test" if put_in_test else "train"
        for sid in ids:
            split[sid] = split_name
        if put_in_test:
            current_test += len(ids)

    return split


def write_outputs(
    assignments_path: Path,
    summary_path: Path,
    labels: dict[str, str],
    clusters: dict[str, str],
    split: dict[str, str],
) -> None:
    assignments_path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for seq_id in sorted(labels):
        rows.append(
            {
                "query_id": seq_id,
                "label": labels[seq_id],
                "cluster_id": clusters[seq_id],
                "split": split[seq_id],
            }
        )

    with assignments_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["query_id", "label", "cluster_id", "split"])
        writer.writeheader()
        writer.writerows(rows)

    n_train = sum(1 for r in rows if r["split"] == "train")
    n_test = sum(1 for r in rows if r["split"] == "test")
    unique_clusters = len({r["cluster_id"] for r in rows})

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerow(["n_total", len(rows)])
        writer.writerow(["n_train", n_train])
        writer.writerow(["n_test", n_test])
        writer.writerow(["n_clusters", unique_clusters])


def main() -> int:
    args = parse_args()

    labels = read_labels(Path(args.labels_csv))
    members = set(labels.keys())
    clusters = build_clusters(Path(args.pairwise_identity_csv), members, args.identity_threshold)
    split = assign_splits(clusters, args.test_fraction)
    write_outputs(Path(args.out_assignments), Path(args.out_summary), labels, clusters, split)

    print(f"Wrote split assignments: {args.out_assignments}")
    print(f"Wrote split summary: {args.out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
