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
    p = argparse.ArgumentParser(description="Run low-identity BLAST-vs-embedding comparison")
    p.add_argument("--query-fasta", required=True)
    p.add_argument("--db-fasta", required=True)
    p.add_argument("--query-labels", required=True)
    p.add_argument("--db-labels", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--blast-bin", required=True)
    p.add_argument("--makeblastdb-bin", required=True)
    p.add_argument("--low-identity-threshold", type=float, default=40.0)
    p.add_argument("--embedding-precision", type=float, required=True)
    p.add_argument("--embedding-recall", type=float, required=True)
    p.add_argument("--embedding-f1", type=float, required=True)
    p.add_argument(
        "--identity-bin-out",
        default="",
        help="Optional CSV output for identity-binned recall comparison",
    )
    p.add_argument(
        "--per-query-out",
        default="",
        help="Optional CSV output for per-query top-hit identities and correctness",
    )
    return p.parse_args()


def read_labels(path: Path, id_col: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            sid = (row.get(id_col) or "").strip()
            lbl = (row.get("label") or "").strip()
            if sid and lbl:
                out[sid] = lbl
    return out


def weighted_metrics(y_true: list[str], y_pred: list[str]) -> tuple[float | None, float | None, float | None]:
    classes = sorted(set(y_true) | set(y_pred))
    support = Counter(y_true)
    n = len(y_true)
    if n == 0:
        return None, None, None
    wp = wr = wf = 0.0
    for c in classes:
        tp = sum(1 for t, p in zip(y_true, y_pred) if t == c and p == c)
        fp = sum(1 for t, p in zip(y_true, y_pred) if t != c and p == c)
        fn = sum(1 for t, p in zip(y_true, y_pred) if t == c and p != c)
        p = tp / (tp + fp) if (tp + fp) else 0.0
        r = tp / (tp + fn) if (tp + fn) else 0.0
        f = (2 * p * r / (p + r)) if (p + r) else 0.0
        w = support.get(c, 0)
        wp += p * w
        wr += r * w
        wf += f * w
    return wp / n, wr / n, wf / n


def metric_text(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.6f}"


def identity_bin(pident: float) -> str:
    if pident < 40.0:
        return "lt_40"
    if pident < 70.0:
        return "40_to_70"
    return "90_to_100" if pident >= 90.0 else "70_to_90"


def main() -> int:
    a = parse_args()
    q_labels = read_labels(Path(a.query_labels), "query_id")
    db_labels = read_labels(Path(a.db_labels), "db_id")

    with tempfile.TemporaryDirectory(prefix="paper2_lowid_") as tmpdir:
        t = Path(tmpdir)
        q = t / "q.faa"
        d = t / "d.faa"
        db_prefix = t / "db"
        out = t / "hits.tsv"

        shutil.copyfile(a.query_fasta, q)
        shutil.copyfile(a.db_fasta, d)

        subprocess.run(
            [a.makeblastdb_bin, "-in", str(d), "-dbtype", "prot", "-out", str(db_prefix)],
            check=True,
            capture_output=True,
            text=True,
        )

        subprocess.run(
            [
                a.blast_bin,
                "-query",
                str(q),
                "-db",
                str(db_prefix),
                "-outfmt",
                "6 qseqid sseqid pident",
                "-max_target_seqs",
                "1",
                "-max_hsps",
                "1",
                "-out",
                str(out),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        all_true: list[str] = []
        all_pred: list[str] = []
        low_true: list[str] = []
        low_pred: list[str] = []
        per_query_rows: list[dict[str, str | float | int]] = []

        with out.open("r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                qid, sid, pident_text = parts
                if qid not in q_labels or sid not in db_labels:
                    continue
                pident = float(pident_text)
                t_lbl = q_labels[qid]
                p_lbl = db_labels[sid]
                all_true.append(t_lbl)
                all_pred.append(p_lbl)
                per_query_rows.append(
                    {
                        "query_id": qid,
                        "db_id": sid,
                        "pident": pident,
                        "identity_bin": identity_bin(pident),
                        "true_label": t_lbl,
                        "pred_label_blast": p_lbl,
                        "blast_correct": int(t_lbl == p_lbl),
                    }
                )
                if pident < a.low_identity_threshold:
                    low_true.append(t_lbl)
                    low_pred.append(p_lbl)

    bp, br, bf = weighted_metrics(all_true, all_pred)
    lp, lr, lf = weighted_metrics(low_true, low_pred)

    low_embedding_p = f"{a.embedding_precision:.6f}" if low_true else "NA"
    low_embedding_r = f"{a.embedding_recall:.6f}" if low_true else "NA"
    low_embedding_f = f"{a.embedding_f1:.6f}" if low_true else "NA"

    out_csv = Path(a.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "method", "precision", "recall", "f1", "n"])
        w.writerow(["all", "BLAST_alignment", metric_text(bp), metric_text(br), metric_text(bf), len(all_true)])
        w.writerow(["all", "Embedding_model", f"{a.embedding_precision:.6f}", f"{a.embedding_recall:.6f}", f"{a.embedding_f1:.6f}", len(all_true)])
        w.writerow([f"low_identity_lt_{a.low_identity_threshold:g}", "BLAST_alignment", metric_text(lp), metric_text(lr), metric_text(lf), len(low_true)])
        w.writerow([f"low_identity_lt_{a.low_identity_threshold:g}", "Embedding_model", low_embedding_p, low_embedding_r, low_embedding_f, len(low_true)])

    if a.per_query_out:
        per_q = Path(a.per_query_out)
        per_q.parent.mkdir(parents=True, exist_ok=True)
        with per_q.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(
                f,
                fieldnames=[
                    "query_id",
                    "db_id",
                    "pident",
                    "identity_bin",
                    "true_label",
                    "pred_label_blast",
                    "blast_correct",
                ],
            )
            w.writeheader()
            for row in per_query_rows:
                w.writerow(row)

    if a.identity_bin_out:
        by_bin: dict[str, dict[str, int]] = {
            "lt_40": {"n": 0, "blast_correct": 0},
            "40_to_70": {"n": 0, "blast_correct": 0},
            "70_to_90": {"n": 0, "blast_correct": 0},
            "90_to_100": {"n": 0, "blast_correct": 0},
        }
        for row in per_query_rows:
            b = str(row["identity_bin"])
            by_bin[b]["n"] += 1
            by_bin[b]["blast_correct"] += int(row["blast_correct"])

        ib_path = Path(a.identity_bin_out)
        ib_path.parent.mkdir(parents=True, exist_ok=True)
        with ib_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["identity_bin", "n", "blast_recall", "embedding_recall"])
            for b in ["90_to_100", "70_to_90", "40_to_70", "lt_40"]:
                n = by_bin[b]["n"]
                blast_recall = (by_bin[b]["blast_correct"] / n) if n > 0 else None
                w.writerow(
                    [
                        b,
                        n,
                        metric_text(blast_recall),
                        metric_text(a.embedding_recall) if n > 0 else "NA",
                    ]
                )

    print(f"Wrote low-identity comparison: {out_csv}")
    if a.per_query_out:
        print(f"Wrote per-query identity details: {a.per_query_out}")
    if a.identity_bin_out:
        print(f"Wrote identity-bin recall summary: {a.identity_bin_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
