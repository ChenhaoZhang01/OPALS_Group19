#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import zlib
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build Paper 2 divergence benchmark outputs")
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
        "--hashed-dim",
        type=int,
        default=4096,
        help="Feature dimension for hashed embedding backend.",
    )
    p.add_argument(
        "--hashed-ks",
        default="3,4",
        help="Comma-separated k-mer sizes for hashed embedding backend.",
    )
    p.add_argument(
        "--embedding-predictor",
        choices=["nearest_db", "class_max", "class_topk", "class_centroid"],
        default="class_topk",
        help=(
            "Prediction rule for embedding labels: nearest single DB hit, "
            "best-per-class max similarity, mean top-k similarity per class, "
            "or class centroid similarity."
        ),
    )
    p.add_argument(
        "--class-topk",
        type=int,
        default=9,
        help="Top-k neighbors per class used by class_topk predictor.",
    )
    p.add_argument(
        "--class-size-penalty",
        type=float,
        default=0.03,
        help="Subtract penalty * log(1 + class_size) from class scores (class predictors only).",
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
    p.add_argument(
        "--identity-bin-edges",
        default="20,30,40,50,70,90",
        help="Comma-separated bin edges in percent identity (inclusive lower, exclusive upper).",
    )
    p.add_argument(
        "--low-id-max",
        type=float,
        default=50.0,
        help="Maximum identity (percent) included in low-identity summary.",
    )
    p.add_argument(
        "--high-id-min",
        type=float,
        default=70.0,
        help="Minimum identity (percent) included in high-identity summary.",
    )
    p.add_argument("--identity-bin-out", required=True, help="Output CSV for identity-binned recall")
    p.add_argument("--benchmark-table-out", required=True, help="Output CSV for benchmark table")
    p.add_argument("--figure-out", required=True, help="Output PNG for recall-vs-identity figure")
    p.add_argument("--embedding-per-query-out", required=True, help="Output CSV for per-query embedding prediction")
    return p.parse_args()


def parse_identity_edges(text: str) -> list[float]:
    parts = [p.strip() for p in text.split(",") if p.strip()]
    if len(parts) < 2:
        raise ValueError("identity-bin-edges must include at least two numeric edges")

    edges = [float(p) for p in parts]
    if any(edges[i] >= edges[i + 1] for i in range(len(edges) - 1)):
        raise ValueError("identity-bin-edges must be strictly increasing")
    return edges


def parse_hashed_ks(text: str) -> list[int]:
    parts = [p.strip() for p in text.split(",") if p.strip()]
    if not parts:
        raise ValueError("--hashed-ks must include at least one integer k-mer length")
    ks = [int(p) for p in parts]
    if any(k < 1 for k in ks):
        raise ValueError("All --hashed-ks values must be >= 1")
    return sorted(set(ks))


def _edge_text(value: float) -> str:
    return str(int(value)) if value.is_integer() else f"{value:g}"


def build_identity_bins(edges: list[float]) -> list[dict[str, float | str]]:
    bins: list[dict[str, float | str]] = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        lower_text = _edge_text(lower)
        upper_text = _edge_text(upper)
        bins.append(
            {
                "key": f"{lower_text}_to_{upper_text}",
                "plot_label": f"{lower_text}-{upper_text}",
                "lower": lower,
                "upper": upper,
                "mid": (lower + upper) / 2.0,
            }
        )
    return bins


def find_identity_bin(pident: float, bins: list[dict[str, float | str]]) -> str | None:
    for b in bins:
        lower = float(b["lower"])
        upper = float(b["upper"])
        if lower <= pident < upper:
            return str(b["key"])
    return None


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


def seq_to_hashed_embedding(seq: str, dim: int = 2048, ks: list[int] | tuple[int, ...] = (3, 4)) -> list[float]:
    vec = [0.0] * dim
    clean = "".join(ch for ch in seq.upper() if "A" <= ch <= "Z")

    total_tokens = 0
    for k in ks:
        if len(clean) < k:
            continue
        n = len(clean) - k + 1
        total_tokens += n
        for i in range(n):
            token = clean[i : i + k]
            key = f"{k}:{token}".encode("ascii", errors="ignore")
            idx = zlib.crc32(key) % dim
            sign = 1.0 if (zlib.crc32(b"s" + key) & 1) == 0 else -1.0
            vec[idx] += sign

    if total_tokens == 0:
        return vec
    inv = 1.0 / total_tokens
    vec = [v * inv for v in vec]
    norm = math.sqrt(sum(v * v for v in vec))
    if norm == 0.0:
        return vec
    return [v / norm for v in vec]


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


def embed_sequences(
    seqs: dict[str, str],
    model: str,
    batch_size: int,
    hashed_dim: int,
    hashed_ks: list[int],
) -> dict[str, list[float]]:
    ids = sorted(seqs)
    sequences = [seqs[i] for i in ids]

    if model == "hashed":
        vecs = [seq_to_hashed_embedding(s, dim=hashed_dim, ks=hashed_ks) for s in sequences]
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
    predictor: str,
    class_topk: int,
    class_size_penalty: float,
) -> dict[str, str]:
    if class_topk < 1:
        raise ValueError("--class-topk must be >= 1")
    if class_size_penalty < 0.0:
        raise ValueError("--class-size-penalty must be >= 0")

    db_ids = [db_id for db_id in sorted(db_emb) if db_id in db_labels]
    if not db_ids:
        raise ValueError("No db embeddings with labels available")

    query_ids = sorted(query_emb)

    try:
        import numpy as np

        q_mat = np.array([query_emb[qid] for qid in query_ids], dtype=np.float32)
        d_mat = np.array([db_emb[db_id] for db_id in db_ids], dtype=np.float32)

        q_norm = np.linalg.norm(q_mat, axis=1)
        d_norm = np.linalg.norm(d_mat, axis=1)

        # Avoid division by zero for empty/degenerate vectors.
        q_norm[q_norm == 0.0] = 1.0
        d_norm[d_norm == 0.0] = 1.0

        sims = (q_mat @ d_mat.T) / (q_norm[:, None] * d_norm[None, :])

        if predictor == "nearest_db":
            best_idx = np.argmax(sims, axis=1)
            return {qid: db_labels[db_ids[int(idx)]] for qid, idx in zip(query_ids, best_idx)}

        labels = sorted({db_labels[db_id] for db_id in db_ids})
        label_to_idx: dict[str, np.ndarray] = {
            lbl: np.array([i for i, db_id in enumerate(db_ids) if db_labels[db_id] == lbl], dtype=np.int32)
            for lbl in labels
        }
        class_scores: list[np.ndarray] = []
        q_unit = q_mat / q_norm[:, None]
        for lbl in labels:
            idx = label_to_idx[lbl]
            lbl_sims = sims[:, idx]

            if predictor == "class_centroid":
                class_vecs = d_mat[idx]
                class_norms = np.linalg.norm(class_vecs, axis=1, keepdims=True)
                class_norms[class_norms == 0.0] = 1.0
                class_unit = class_vecs / class_norms
                centroid = np.mean(class_unit, axis=0)
                centroid_norm = np.linalg.norm(centroid)
                if centroid_norm == 0.0:
                    centroid_norm = 1.0
                centroid = centroid / centroid_norm
                score = q_unit @ centroid
            elif predictor == "class_max":
                score = np.max(lbl_sims, axis=1)
            else:
                k = min(class_topk, lbl_sims.shape[1])
                if k == lbl_sims.shape[1]:
                    score = np.mean(lbl_sims, axis=1)
                else:
                    topk = np.partition(lbl_sims, kth=lbl_sims.shape[1] - k, axis=1)[:, -k:]
                    score = np.mean(topk, axis=1)
            if predictor != "nearest_db" and class_size_penalty > 0.0:
                score = score - (class_size_penalty * math.log1p(int(idx.shape[0])))
            class_scores.append(score)

        stacked = np.stack(class_scores, axis=1)
        best_class_idx = np.argmax(stacked, axis=1)
        return {qid: labels[int(i)] for qid, i in zip(query_ids, best_class_idx)}
    except Exception:
        db_items = [(db_id, db_labels[db_id], db_emb[db_id]) for db_id in db_ids]
        labels = sorted({db_labels[db_id] for db_id in db_ids})
        label_items: dict[str, list[list[float]]] = {lbl: [] for lbl in labels}
        for db_id, lbl, dvec in db_items:
            _ = db_id
            label_items[lbl].append(dvec)

        out: dict[str, str] = {}
        for qid in query_ids:
            qvec = query_emb[qid]
            best_label = ""
            best_score = -2.0

            if predictor == "nearest_db":
                for _dbid, lbl, dvec in db_items:
                    score = cosine_similarity(qvec, dvec)
                    if score > best_score:
                        best_score = score
                        best_label = lbl
                out[qid] = best_label
                continue

            for lbl in labels:
                sims = [cosine_similarity(qvec, dvec) for dvec in label_items[lbl]]
                if not sims:
                    continue
                if predictor == "class_centroid":
                    d = len(label_items[lbl][0])
                    centroid = [0.0] * d
                    for dvec in label_items[lbl]:
                        for i, val in enumerate(dvec):
                            centroid[i] += val
                    denom = float(len(label_items[lbl]))
                    centroid = [val / denom for val in centroid]
                    score = cosine_similarity(qvec, centroid)
                elif predictor == "class_max":
                    score = max(sims)
                else:
                    sims.sort(reverse=True)
                    k = min(class_topk, len(sims))
                    score = sum(sims[:k]) / k
                if class_size_penalty > 0.0:
                    score -= class_size_penalty * math.log1p(len(label_items[lbl]))
                if score > best_score:
                    best_score = score
                    best_label = lbl
            out[qid] = best_label
        return out


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
    if args.low_id_max > args.high_id_min:
        raise ValueError("low-id-max cannot exceed high-id-min")
    if args.hashed_dim < 16:
        raise ValueError("--hashed-dim must be >= 16")

    identity_edges = parse_identity_edges(args.identity_bin_edges)
    identity_bins = build_identity_bins(identity_edges)
    hashed_ks = parse_hashed_ks(args.hashed_ks)

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
        query_emb = embed_sequences(
            query_seqs,
            args.embedding_model,
            args.batch_size,
            hashed_dim=args.hashed_dim,
            hashed_ks=hashed_ks,
        )
        db_emb = embed_sequences(
            db_seqs,
            args.embedding_model,
            args.batch_size,
            hashed_dim=args.hashed_dim,
            hashed_ks=hashed_ks,
        )
    else:
        if args.hf_model_id:
            model_name = args.hf_model_id
        elif args.embedding_model == "protbert":
            model_name = "Rostlab/prot_bert"
        else:
            model_name = "facebook/esm2_t6_8M_UR50D"

        query_emb = embed_sequences_transformer(query_seqs, model_name, args.batch_size, args.max_length)
        db_emb = embed_sequences_transformer(db_seqs, model_name, args.batch_size, args.max_length)
    emb_pred = predict_by_nearest_db(
        query_emb,
        db_emb,
        db_labels,
        predictor=args.embedding_predictor,
        class_topk=args.class_topk,
        class_size_penalty=args.class_size_penalty,
    )

    embedding_per_query_rows: list[dict[str, str]] = []
    blast_correct_total = 0
    embedding_correct_total = 0

    by_bin_counts: dict[str, dict[str, int]] = {
        str(b["key"]): {"n": 0, "blast_correct": 0, "embedding_correct": 0} for b in identity_bins
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

        b = find_identity_bin(pident, identity_bins)
        if b is not None:
            by_bin_counts[b]["n"] += 1
            by_bin_counts[b]["blast_correct"] += blast_correct
            by_bin_counts[b]["embedding_correct"] += emb_correct

    n_total = len(embedding_per_query_rows)
    blast_overall = safe_ratio(blast_correct_total, n_total)
    emb_overall = safe_ratio(embedding_correct_total, n_total)

    low_bins = [str(b["key"]) for b in identity_bins if float(b["upper"]) <= args.low_id_max]
    high_bins = [str(b["key"]) for b in identity_bins if float(b["lower"]) >= args.high_id_min]

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
        for b in [str(x["key"]) for x in identity_bins]:
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
        writer.writerow(
            [
                f"Embedding_{args.embedding_model}_{args.embedding_predictor}",
                fmt(emb_overall),
                fmt(emb_low),
                fmt(emb_high),
            ]
        )

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to generate the figure") from exc

    x_labels = [str(b["plot_label"]) for b in identity_bins]
    x = [float(b["mid"]) for b in identity_bins]

    blast_y = []
    emb_y = []
    for b in [str(x["key"]) for x in identity_bins]:
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
    print(f"Identity bins: {', '.join(x_labels)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
