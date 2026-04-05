#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build Paper 3 divergence benchmark outputs")
    p.add_argument("--blast-per-query", required=True, help="CSV with query_id,pident,blast_correct")
    p.add_argument("--query-labels", required=True, help="CSV with query_id,label")
    p.add_argument("--db-labels", required=True, help="CSV with db_id,label")
    p.add_argument("--query-fasta", required=True, help="Query FASTA")
    p.add_argument("--db-fasta", required=True, help="DB FASTA")
    p.add_argument(
        "--embedding-model",
        choices=["hashed", "protbert", "esm2"],
        default="hashed",
        help="Embedding backend. Use protbert/esm2 only if transformers+torch are available.",
    )
    p.add_argument(
        "--hf-model-id",
        default="",
        help="Optional Hugging Face model ID override for protbert/esm2 backends.",
    )
    p.add_argument("--batch-size", type=int, default=16, help="Batch size for transformer embedding")
    p.add_argument(
        "--max-length",
        type=int,
        default=256,
        help="Maximum token length for transformer models (CPU-friendly speed control).",
    )
    p.add_argument("--identity-bin-out", required=True, help="Output CSV for identity-binned recall")
    p.add_argument("--benchmark-table-out", required=True, help="Output CSV for benchmark table")
    p.add_argument("--figure-out", required=True, help="Output PNG for recall-vs-identity figure")
    p.add_argument("--embedding-per-query-out", required=True, help="Output CSV for per-query embedding prediction")
    return p.parse_args()


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    cur_id = ""
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id:
                    out[cur_id] = "".join(seq_parts)
                cur_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if cur_id:
            out[cur_id] = "".join(seq_parts)
    return out


def read_label_csv(path: Path, id_col: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            sid = (row.get(id_col) or "").strip()
            label = (row.get("label") or "").strip()
            if sid and label:
                out[sid] = label
    return out


def cosine_similarity(a: list[float], b: list[float]) -> float:
    dot = 0.0
    na = 0.0
    nb = 0.0
    for x, y in zip(a, b):
        dot += x * y
        na += x * x
        nb += y * y
    if na == 0.0 or nb == 0.0:
        return -1.0
    return dot / math.sqrt(na * nb)


def seq_to_hashed_embedding(seq: str, dim: int = 256, k: int = 3) -> list[float]:
    vec = [0.0] * dim
    clean = "".join(ch for ch in seq.upper() if "A" <= ch <= "Z")
    if len(clean) < k:
        return vec
    n = len(clean) - k + 1
    for i in range(n):
        token = clean[i : i + k]
        idx = (hash(token) & 0xFFFFFFFF) % dim
        vec[idx] += 1.0
    inv = 1.0 / n
    return [v * inv for v in vec]


def transformer_embeddings(
    sequences: list[str],
    model_name: str,
    batch_size: int,
    max_length: int,
) -> list[list[float]]:
    try:
        import importlib

        torch = importlib.import_module("torch")
        transformers = importlib.import_module("transformers")
        AutoModel = transformers.AutoModel
        AutoTokenizer = transformers.AutoTokenizer
    except ImportError as exc:
        raise RuntimeError("transformers and torch are required for protbert/esm2 mode") from exc

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()

    all_vecs: list[list[float]] = []
    with torch.no_grad():
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i : i + batch_size]
            # Protein models generally expect space-delimited amino acids.
            prepared = [" ".join(list(seq)) for seq in batch]
            toks = tokenizer(
                prepared,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=max_length,
            )
            out = model(**toks)
            hidden = out.last_hidden_state
            mask = toks["attention_mask"].unsqueeze(-1)
            summed = (hidden * mask).sum(dim=1)
            counts = mask.sum(dim=1).clamp(min=1)
            pooled = summed / counts
            for row in pooled:
                all_vecs.append(row.detach().cpu().tolist())
    return all_vecs


def embed_sequences(seqs: dict[str, str], model: str, batch_size: int) -> dict[str, list[float]]:
    ids = sorted(seqs)
    sequences = [seqs[i] for i in ids]

    if model == "hashed":
        vecs = [seq_to_hashed_embedding(s) for s in sequences]
    else:
        raise ValueError("embed_sequences requires a model_name for transformer backends")

    return {sid: vec for sid, vec in zip(ids, vecs)}


def embed_sequences_transformer(
    seqs: dict[str, str],
    model_name: str,
    batch_size: int,
    max_length: int,
) -> dict[str, list[float]]:
    ids = sorted(seqs)
    sequences = [seqs[i] for i in ids]
    vecs = transformer_embeddings(sequences, model_name, batch_size, max_length)
    return {sid: vec for sid, vec in zip(ids, vecs)}


def predict_by_nearest_db(
    query_emb: dict[str, list[float]],
    db_emb: dict[str, list[float]],
    db_labels: dict[str, str],
) -> dict[str, str]:
    db_items = [(db_id, db_labels[db_id], db_emb[db_id]) for db_id in sorted(db_emb) if db_id in db_labels]
    if not db_items:
        raise ValueError("No db embeddings with labels available")

    out: dict[str, str] = {}
    for qid, qvec in query_emb.items():
        best_label = ""
        best_score = -2.0
        for _dbid, lbl, dvec in db_items:
            score = cosine_similarity(qvec, dvec)
            if score > best_score:
                best_score = score
                best_label = lbl
        out[qid] = best_label
    return out


def id_bin_3way(pident: float) -> str | None:
    if 30.0 <= pident < 50.0:
        return "30_to_50"
    if 50.0 <= pident < 70.0:
        return "50_to_70"
    if 70.0 <= pident < 90.0:
        return "70_to_90"
    return None


def fmt(x: float | None) -> str:
    if x is None:
        return "NA"
    return f"{x:.6f}"


def safe_ratio(num: int, den: int) -> float | None:
    if den == 0:
        return None
    return num / den


def main() -> int:
    args = parse_args()

    blast_rows: list[dict[str, str]] = []
    with Path(args.blast_per_query).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        needed = {"query_id", "pident", "blast_correct"}
        if not needed.issubset(set(reader.fieldnames or [])):
            raise ValueError("blast-per-query CSV must include query_id,pident,blast_correct")
        blast_rows = [row for row in reader]

    query_labels = read_label_csv(Path(args.query_labels), "query_id")
    db_labels = read_label_csv(Path(args.db_labels), "db_id")
    query_fasta = read_fasta(Path(args.query_fasta))
    db_fasta = read_fasta(Path(args.db_fasta))

    usable_query_ids = sorted({(r.get("query_id") or "").strip() for r in blast_rows if (r.get("query_id") or "").strip() in query_labels and (r.get("query_id") or "").strip() in query_fasta})
    query_seqs = {qid: query_fasta[qid] for qid in usable_query_ids}
    db_seqs = {dbid: seq for dbid, seq in db_fasta.items() if dbid in db_labels}

    if args.embedding_model == "hashed":
        query_emb = embed_sequences(query_seqs, args.embedding_model, args.batch_size)
        db_emb = embed_sequences(db_seqs, args.embedding_model, args.batch_size)
    else:
        if args.hf_model_id:
            model_name = args.hf_model_id
        elif args.embedding_model == "protbert":
            model_name = "Rostlab/prot_bert"
        else:
            model_name = "facebook/esm2_t6_8M_UR50D"

        query_emb = embed_sequences_transformer(query_seqs, model_name, args.batch_size, args.max_length)
        db_emb = embed_sequences_transformer(db_seqs, model_name, args.batch_size, args.max_length)
    emb_pred = predict_by_nearest_db(query_emb, db_emb, db_labels)

    embedding_per_query_rows: list[dict[str, str]] = []
    blast_correct_total = 0
    embedding_correct_total = 0

    by_bin_counts: dict[str, dict[str, int]] = {
        "30_to_50": {"n": 0, "blast_correct": 0, "embedding_correct": 0},
        "50_to_70": {"n": 0, "blast_correct": 0, "embedding_correct": 0},
        "70_to_90": {"n": 0, "blast_correct": 0, "embedding_correct": 0},
    }

    for row in blast_rows:
        qid = (row.get("query_id") or "").strip()
        if qid not in query_labels or qid not in emb_pred:
            continue

        pident = float((row.get("pident") or "0").strip())
        true_label = query_labels[qid]
        pred_label = emb_pred[qid]

        emb_correct = int(pred_label == true_label)
        blast_correct = int((row.get("blast_correct") or "0").strip())

        embedding_per_query_rows.append(
            {
                "query_id": qid,
                "pident": f"{pident:.6f}",
                "true_label": true_label,
                "pred_label_embedding": pred_label,
                "embedding_correct": str(emb_correct),
                "blast_correct": str(blast_correct),
            }
        )

        blast_correct_total += blast_correct
        embedding_correct_total += emb_correct

        b = id_bin_3way(pident)
        if b is not None:
            by_bin_counts[b]["n"] += 1
            by_bin_counts[b]["blast_correct"] += blast_correct
            by_bin_counts[b]["embedding_correct"] += emb_correct

    n_total = len(embedding_per_query_rows)
    blast_overall = safe_ratio(blast_correct_total, n_total)
    emb_overall = safe_ratio(embedding_correct_total, n_total)

    low_bins = ["30_to_50", "50_to_70"]
    high_bins = ["70_to_90"]

    low_n = sum(by_bin_counts[b]["n"] for b in low_bins)
    high_n = sum(by_bin_counts[b]["n"] for b in high_bins)

    blast_low = safe_ratio(sum(by_bin_counts[b]["blast_correct"] for b in low_bins), low_n)
    emb_low = safe_ratio(sum(by_bin_counts[b]["embedding_correct"] for b in low_bins), low_n)
    blast_high = safe_ratio(sum(by_bin_counts[b]["blast_correct"] for b in high_bins), high_n)
    emb_high = safe_ratio(sum(by_bin_counts[b]["embedding_correct"] for b in high_bins), high_n)

    per_query_path = Path(args.embedding_per_query_out)
    per_query_path.parent.mkdir(parents=True, exist_ok=True)
    with per_query_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "query_id",
                "pident",
                "true_label",
                "pred_label_embedding",
                "embedding_correct",
                "blast_correct",
            ],
        )
        writer.writeheader()
        writer.writerows(embedding_per_query_rows)

    identity_path = Path(args.identity_bin_out)
    identity_path.parent.mkdir(parents=True, exist_ok=True)
    with identity_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["identity_bin", "n", "blast_recall", "embedding_recall"])
        for b in ["30_to_50", "50_to_70", "70_to_90"]:
            counts = by_bin_counts[b]
            writer.writerow(
                [
                    b,
                    counts["n"],
                    fmt(safe_ratio(counts["blast_correct"], counts["n"])),
                    fmt(safe_ratio(counts["embedding_correct"], counts["n"])),
                ]
            )

    benchmark_path = Path(args.benchmark_table_out)
    benchmark_path.parent.mkdir(parents=True, exist_ok=True)
    with benchmark_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "overall_recall", "low_id_recall", "high_id_recall"])
        writer.writerow(["BLAST_alignment", fmt(blast_overall), fmt(blast_low), fmt(blast_high)])
        writer.writerow([f"Embedding_{args.embedding_model}", fmt(emb_overall), fmt(emb_low), fmt(emb_high)])

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to generate the figure") from exc

    x_labels = ["30-50", "50-70", "70-90"]
    x = [40, 60, 80]

    blast_y = []
    emb_y = []
    for b in ["30_to_50", "50_to_70", "70_to_90"]:
        c = by_bin_counts[b]
        blast_y.append(safe_ratio(c["blast_correct"], c["n"]) or 0.0)
        emb_y.append(safe_ratio(c["embedding_correct"], c["n"]) or 0.0)

    fig_path = Path(args.figure_out)
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7, 4.5))
    plt.plot(x, blast_y, marker="o", linewidth=2, label="BLAST")
    plt.plot(x, emb_y, marker="s", linewidth=2, label=f"Embedding ({args.embedding_model})")
    plt.xticks(x, x_labels)
    plt.ylim(0.0, 1.05)
    plt.xlabel("Sequence identity bin (%)")
    plt.ylabel("Recall")
    plt.title("Recall versus sequence identity")
    plt.grid(alpha=0.25)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=220)
    plt.close()

    print(f"Wrote: {args.identity_bin_out}")
    print(f"Wrote: {args.benchmark_table_out}")
    print(f"Wrote: {args.figure_out}")
    print(f"Wrote: {args.embedding_per_query_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
