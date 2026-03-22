#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run BLAST baseline and write metrics")
    parser.add_argument("--query-fasta", required=True, help="Query proteins FASTA")
    parser.add_argument("--db-fasta", required=True, help="Reference proteins FASTA")
    parser.add_argument("--query-labels", required=True, help="CSV with columns query_id,label")
    parser.add_argument("--db-labels", required=True, help="CSV with columns db_id,label")
    parser.add_argument("--metrics-out", required=True, help="Shared metrics CSV")
    parser.add_argument("--blast-bin", default="blastp", help="blastp executable path")
    parser.add_argument("--makeblastdb-bin", default="makeblastdb", help="makeblastdb executable path")
    parser.add_argument("--evalue", type=float, default=1e-5, help="BLAST e-value threshold")
    return parser.parse_args()


def read_label_csv(path: Path, id_col: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            seq_id = (row.get(id_col) or "").strip()
            label = (row.get("label") or "").strip()
            if seq_id and label:
                out[seq_id] = label
    return out


def compute_weighted_metrics(y_true: list[str], y_pred: list[str]) -> tuple[float, float, float]:
    classes = sorted(set(y_true) | set(y_pred))
    support = Counter(y_true)
    n = len(y_true)
    if n == 0:
        return 0.0, 0.0, 0.0

    weighted_precision = 0.0
    weighted_recall = 0.0
    weighted_f1 = 0.0

    for cls in classes:
        tp = sum(1 for t, p in zip(y_true, y_pred) if t == cls and p == cls)
        fp = sum(1 for t, p in zip(y_true, y_pred) if t != cls and p == cls)
        fn = sum(1 for t, p in zip(y_true, y_pred) if t == cls and p != cls)

        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0

        w = support.get(cls, 0)
        weighted_precision += precision * w
        weighted_recall += recall * w
        weighted_f1 += f1 * w

    return weighted_precision / n, weighted_recall / n, weighted_f1 / n


def upsert_metrics_row(path: Path, method: str, precision: str, recall: str, f1: str, roc_auc: str) -> None:
    fieldnames = ["method", "precision", "recall", "f1", "roc_auc"]
    rows: list[dict[str, str]] = []

    if path.exists():
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                name = (row.get("method") or "").strip()
                if not name:
                    continue
                rows.append(
                    {
                        "method": name,
                        "precision": (row.get("precision") or "").strip(),
                        "recall": (row.get("recall") or "").strip(),
                        "f1": (row.get("f1") or "").strip(),
                        "roc_auc": (row.get("roc_auc") or "").strip(),
                    }
                )

    target = {
        "method": method,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "roc_auc": roc_auc,
    }

    replaced = False
    for i, row in enumerate(rows):
        if row["method"] == method:
            rows[i] = target
            replaced = True
            break
    if not replaced:
        rows.append(target)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()

    blastp = shutil.which(args.blast_bin) if Path(args.blast_bin).name == args.blast_bin else args.blast_bin
    makeblastdb = (
        shutil.which(args.makeblastdb_bin)
        if Path(args.makeblastdb_bin).name == args.makeblastdb_bin
        else args.makeblastdb_bin
    )

    if not blastp or not Path(blastp).exists():
        raise SystemExit(f"blastp not found: {args.blast_bin}")
    if not makeblastdb or not Path(makeblastdb).exists():
        raise SystemExit(f"makeblastdb not found: {args.makeblastdb_bin}")

    query_fasta = Path(args.query_fasta)
    db_fasta = Path(args.db_fasta)
    query_labels = read_label_csv(Path(args.query_labels), "query_id")
    db_labels = read_label_csv(Path(args.db_labels), "db_id")

    if not query_labels:
        raise SystemExit("No query labels found")
    if not db_labels:
        raise SystemExit("No db labels found")

    with tempfile.TemporaryDirectory(prefix="paper2_blast_") as tmpdir:
        tmp = Path(tmpdir)
        staged_query = tmp / "query.faa"
        staged_db = tmp / "db.faa"
        db_prefix = tmp / "arg_db"
        out_file = tmp / "blast_hits.tsv"

        shutil.copyfile(query_fasta, staged_query)
        shutil.copyfile(db_fasta, staged_db)

        subprocess.run(
            [
                makeblastdb,
                "-in",
                str(staged_db),
                "-dbtype",
                "prot",
                "-out",
                str(db_prefix),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        subprocess.run(
            [
                blastp,
                "-query",
                str(staged_query),
                "-db",
                str(db_prefix),
                "-outfmt",
                "6 qseqid sseqid bitscore evalue",
                "-max_target_seqs",
                "1",
                "-max_hsps",
                "1",
                "-evalue",
                str(args.evalue),
                "-out",
                str(out_file),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        pred_by_query: dict[str, str] = {}
        if out_file.exists():
            with out_file.open("r", encoding="utf-8") as handle:
                for line in handle:
                    parts = line.strip().split("\t")
                    if len(parts) < 2:
                        continue
                    qid, sid = parts[0], parts[1]
                    if qid not in query_labels:
                        continue
                    if sid in db_labels and qid not in pred_by_query:
                        pred_by_query[qid] = db_labels[sid]

    y_true: list[str] = []
    y_pred: list[str] = []
    for qid, true_label in query_labels.items():
        if qid in pred_by_query:
            y_true.append(true_label)
            y_pred.append(pred_by_query[qid])

    if not y_true:
        raise SystemExit("No labeled BLAST predictions matched query labels")

    precision, recall, f1 = compute_weighted_metrics(y_true, y_pred)

    upsert_metrics_row(
        path=Path(args.metrics_out),
        method="BLAST_alignment",
        precision=f"{precision:.6f}",
        recall=f"{recall:.6f}",
        f1=f"{f1:.6f}",
        roc_auc="NA",
    )

    print(f"BLAST rows scored: {len(y_true)}")
    print(f"BLAST weighted precision: {precision:.6f}")
    print(f"BLAST weighted recall: {recall:.6f}")
    print(f"BLAST weighted f1: {f1:.6f}")
    print(f"Updated metrics: {args.metrics_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
