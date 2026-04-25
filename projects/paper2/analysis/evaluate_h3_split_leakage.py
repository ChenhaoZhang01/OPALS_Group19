#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Evaluate H3 by comparing random split vs identity-clustered split."
    )
    p.add_argument("--embeddings", required=True, help="Path to .npy embeddings for card_query_homolog.faa")
    p.add_argument("--fasta", required=True, help="FASTA used to generate embeddings")
    p.add_argument("--labels-csv", required=True, help="CSV with query_id,label")
    p.add_argument("--pairwise-identity-csv", required=True, help="CSV with query_id,subject_id,pident")
    p.add_argument("--identity-threshold", type=float, default=70.0, help="Leakage threshold in percent identity")
    p.add_argument("--test-fraction", type=float, default=0.2, help="Target test fraction for clustered split")
    p.add_argument("--random-seed", type=int, default=42, help="Seed for one-shot random split")
    p.add_argument("--repeat-random", type=int, default=20, help="Number of repeated random splits")
    p.add_argument("--out-metrics-csv", required=True, help="Output CSV with split comparison metrics")
    p.add_argument("--out-summary-json", required=True, help="Output JSON summary for manuscript reporting")
    return p.parse_args()


def read_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def read_labels(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            qid = (row.get("query_id") or "").strip()
            label = (row.get("label") or "").strip()
            if qid and label:
                out[qid] = label
    return out


def read_cluster_assignments(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            qid = (row.get("query_id") or "").strip()
            split = (row.get("split") or "").strip()
            if qid and split:
                out[qid] = split
    return out


class DSU:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

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


def load_adjacency(path: Path, threshold: float, members: set[str]) -> dict[str, set[str]]:
    adj: dict[str, set[str]] = {m: set() for m in members}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            qid = (row.get("query_id") or "").strip()
            sid = (row.get("subject_id") or "").strip()
            if qid not in members or sid not in members or qid == sid:
                continue
            try:
                pident = float((row.get("pident") or "").strip())
            except ValueError:
                continue
            if pident >= threshold:
                adj[qid].add(sid)
                adj[sid].add(qid)
    return adj


def build_clusters(path: Path, threshold: float, members: set[str]) -> dict[str, str]:
    dsu = DSU()
    for m in members:
        dsu.add(m)

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
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
    return {m: dsu.find(m) for m in members}


def assign_cluster_split(clusters: dict[str, str], test_fraction: float) -> dict[str, str]:
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


def leakage_rate(train_ids: set[str], test_ids: set[str], adj: dict[str, set[str]]) -> float:
    if not test_ids:
        return 0.0
    leak = 0
    for tid in test_ids:
        if any(n in train_ids for n in adj.get(tid, set())):
            leak += 1
    return leak / len(test_ids)


def class_counts(labels: list[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for y in labels:
        out[y] = out.get(y, 0) + 1
    return out


def evaluate_split(X_train, y_train, X_test, y_test):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score

    clf = RandomForestClassifier(n_estimators=300, random_state=42)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    return {
        "precision": float(precision_score(y_test, y_pred, average="weighted", zero_division=0)),
        "recall": float(recall_score(y_test, y_pred, average="weighted", zero_division=0)),
        "f1": float(f1_score(y_test, y_pred, average="weighted", zero_division=0)),
        "accuracy": float(accuracy_score(y_test, y_pred)),
    }


def main() -> int:
    args = parse_args()

    import numpy as np
    from sklearn.model_selection import train_test_split

    embeddings = np.load(args.embeddings)
    fasta_ids = read_fasta_ids(Path(args.fasta))
    labels_by_id = read_labels(Path(args.labels_csv))

    if len(fasta_ids) != len(embeddings):
        raise ValueError(
            f"FASTA IDs ({len(fasta_ids)}) and embedding rows ({len(embeddings)}) do not match"
        )

    usable_ids = [qid for qid in fasta_ids if qid in labels_by_id]
    if len(usable_ids) < 20:
        raise ValueError("Too few usable labeled IDs")

    id_to_index = {qid: i for i, qid in enumerate(fasta_ids)}
    X = np.stack([embeddings[id_to_index[qid]] for qid in usable_ids], axis=0)
    y = [labels_by_id[qid] for qid in usable_ids]

    members = set(usable_ids)
    adj = load_adjacency(Path(args.pairwise_identity_csv), args.identity_threshold, members)

    # Random split (seeded, stratified)
    X_train_r, X_test_r, y_train_r, y_test_r, ids_train_r, ids_test_r = train_test_split(
        X,
        y,
        usable_ids,
        test_size=0.2,
        random_state=args.random_seed,
        stratify=y,
    )
    random_metrics = evaluate_split(X_train_r, y_train_r, X_test_r, y_test_r)
    random_leak = leakage_rate(set(ids_train_r), set(ids_test_r), adj)

    # Repeated random splits for stability
    repeated: list[float] = []
    for seed in range(args.repeat_random):
        X_tr, X_te, y_tr, y_te = train_test_split(
            X,
            y,
            test_size=0.2,
            random_state=seed,
            stratify=y,
        )
        m = evaluate_split(X_tr, y_tr, X_te, y_te)
        repeated.append(m["recall"])
    repeated_mean = float(np.mean(repeated)) if repeated else 0.0
    repeated_std = float(np.std(repeated)) if repeated else 0.0

    # Clustered split (stratified by label with cluster groups)
    clusters = build_clusters(Path(args.pairwise_identity_csv), args.identity_threshold, members)
    groups = [clusters[qid] for qid in usable_ids]

    n_splits = max(2, int(round(1.0 / args.test_fraction)))
    train_idx_c = None
    test_idx_c = None
    try:
        from sklearn.model_selection import StratifiedGroupKFold

        sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
        for tr_idx, te_idx in sgkf.split(X, y, groups):
            y_tr_try = [y[i] for i in tr_idx]
            y_te_try = [y[i] for i in te_idx]
            if len(set(y_tr_try)) >= 2 and len(set(y_te_try)) >= 2:
                train_idx_c = tr_idx
                test_idx_c = te_idx
                break
    except Exception:
        train_idx_c = None
        test_idx_c = None

    if train_idx_c is None or test_idx_c is None:
        split_by_id = assign_cluster_split(clusters, args.test_fraction)
        train_ids_c = [qid for qid in usable_ids if split_by_id[qid] == "train"]
        test_ids_c = [qid for qid in usable_ids if split_by_id[qid] == "test"]
    else:
        train_ids_c = [usable_ids[i] for i in train_idx_c]
        test_ids_c = [usable_ids[i] for i in test_idx_c]
    if not train_ids_c or not test_ids_c:
        raise ValueError("Clustered split produced empty train or test set")

    X_train_c = np.stack([embeddings[id_to_index[qid]] for qid in train_ids_c], axis=0)
    y_train_c = [labels_by_id[qid] for qid in train_ids_c]
    X_test_c = np.stack([embeddings[id_to_index[qid]] for qid in test_ids_c], axis=0)
    y_test_c = [labels_by_id[qid] for qid in test_ids_c]

    if len(set(y_train_c)) < 2 or len(set(y_test_c)) < 2:
        raise ValueError("Clustered split does not contain both classes in train/test")

    clustered_metrics = evaluate_split(X_train_c, y_train_c, X_test_c, y_test_c)
    clustered_leak = leakage_rate(set(train_ids_c), set(test_ids_c), adj)

    out_rows = [
        {
            "split": "random_stratified_seed42",
            "n_train": len(ids_train_r),
            "n_test": len(ids_test_r),
            "precision": f"{random_metrics['precision']:.6f}",
            "recall": f"{random_metrics['recall']:.6f}",
            "f1": f"{random_metrics['f1']:.6f}",
            "accuracy": f"{random_metrics['accuracy']:.6f}",
            "leakage_rate": f"{random_leak:.6f}",
            "threshold": f"{args.identity_threshold:.1f}",
        },
        {
            "split": "identity_clustered",
            "n_train": len(train_ids_c),
            "n_test": len(test_ids_c),
            "precision": f"{clustered_metrics['precision']:.6f}",
            "recall": f"{clustered_metrics['recall']:.6f}",
            "f1": f"{clustered_metrics['f1']:.6f}",
            "accuracy": f"{clustered_metrics['accuracy']:.6f}",
            "leakage_rate": f"{clustered_leak:.6f}",
            "threshold": f"{args.identity_threshold:.1f}",
        },
    ]

    out_metrics = Path(args.out_metrics_csv)
    out_metrics.parent.mkdir(parents=True, exist_ok=True)
    with out_metrics.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "split",
                "n_train",
                "n_test",
                "precision",
                "recall",
                "f1",
                "accuracy",
                "leakage_rate",
                "threshold",
            ],
        )
        writer.writeheader()
        writer.writerows(out_rows)

    summary = {
        "identity_threshold": args.identity_threshold,
        "n_total": len(usable_ids),
        "class_counts_total": class_counts(y),
        "random_seed42": {
            **random_metrics,
            "n_train": len(ids_train_r),
            "n_test": len(ids_test_r),
            "class_counts_train": class_counts(y_train_r),
            "class_counts_test": class_counts(y_test_r),
            "leakage_rate": random_leak,
        },
        "random_repeats": {
            "n_splits": args.repeat_random,
            "recall_mean": repeated_mean,
            "recall_std": repeated_std,
        },
        "identity_clustered": {
            **clustered_metrics,
            "n_train": len(train_ids_c),
            "n_test": len(test_ids_c),
            "n_clusters": len(set(clusters.values())),
            "class_counts_train": class_counts(y_train_c),
            "class_counts_test": class_counts(y_test_c),
            "leakage_rate": clustered_leak,
        },
        "delta_random_minus_clustered": {
            "precision": random_metrics["precision"] - clustered_metrics["precision"],
            "recall": random_metrics["recall"] - clustered_metrics["recall"],
            "f1": random_metrics["f1"] - clustered_metrics["f1"],
            "accuracy": random_metrics["accuracy"] - clustered_metrics["accuracy"],
        },
    }

    out_summary = Path(args.out_summary_json)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"Wrote metrics: {out_metrics}")
    print(f"Wrote summary: {out_summary}")
    print(
        "H3 snapshot: "
        f"random recall={random_metrics['recall']:.3f} (leak={random_leak:.3f}) vs "
        f"clustered recall={clustered_metrics['recall']:.3f} (leak={clustered_leak:.3f})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
