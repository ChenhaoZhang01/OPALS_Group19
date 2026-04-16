#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import random
import tarfile
from collections import Counter, defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare Paper 2 benchmark inputs from CARD data tar.bz2")
    p.add_argument("--card-archive", required=True, help="Path to CARD archive (latest/data payload)")
    p.add_argument(
        "--query-member",
        default="./protein_fasta_protein_variant_model.fasta",
        help="TAR member path for query FASTA",
    )
    p.add_argument(
        "--query-member-2",
        default="",
        help="Optional second TAR member path for query FASTA (merged with query-member)",
    )
    p.add_argument(
        "--db-member",
        default="./protein_fasta_protein_homolog_model.fasta",
        help="TAR member path for DB FASTA",
    )
    p.add_argument(
        "--label-category",
        default="Resistance Mechanism",
        choices=["Resistance Mechanism", "Drug Class", "AMR Gene Family"],
        help="ARO category to use as class label",
    )
    p.add_argument("--max-query", type=int, default=800, help="Maximum query sequences to keep")
    p.add_argument("--max-db", type=int, default=3000, help="Maximum DB sequences to keep")
    p.add_argument("--min-class-count", type=int, default=8, help="Minimum query samples per class")
    p.add_argument(
        "--min-db-class-count",
        type=int,
        default=1,
        help="Minimum DB samples per class",
    )
    p.add_argument("--embedding-dim", type=int, default=128, help="Hashed k-mer embedding dimension")
    p.add_argument("--kmer", type=int, default=3, help="K-mer size for hashed embeddings")
    p.add_argument(
        "--augment-mutation-rates",
        default="",
        help="Optional comma-separated mutation rates for synthetic query augmentation, e.g. 0.2,0.4",
    )
    p.add_argument(
        "--augment-copies-per-rate",
        type=int,
        default=0,
        help="Synthetic copies per sequence for each mutation rate (0 disables augmentation)",
    )
    p.add_argument("--augmentation-seed", type=int, default=42, help="Seed for synthetic augmentation")
    p.add_argument("--out-query-fasta", required=True)
    p.add_argument("--out-db-fasta", required=True)
    p.add_argument("--out-query-labels", required=True)
    p.add_argument("--out-db-labels", required=True)
    p.add_argument("--out-embeddings", required=True)
    p.add_argument("--out-training-labels", required=True)
    return p.parse_args()


def parse_fasta(text: str) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    cur_id = ""
    cur_header = ""
    seq_parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_id:
                records.append((cur_id, cur_header, "".join(seq_parts)))
            cur_header = line[1:]
            cur_id = cur_header.split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)
    if cur_id:
        records.append((cur_id, cur_header, "".join(seq_parts)))
    return records


def dedupe_records(records: list[tuple[str, str, str]]) -> list[tuple[str, str, str]]:
    # Keep the longest sequence per ID when the same ID appears across merged query sources.
    by_id: dict[str, tuple[str, str, str]] = {}
    for seq_id, header, seq in records:
        prev = by_id.get(seq_id)
        if prev is None or len(seq) > len(prev[2]):
            by_id[seq_id] = (seq_id, header, seq)
    return list(by_id.values())


def extract_aro_id(header: str) -> str | None:
    marker = "ARO:"
    idx = header.find(marker)
    if idx < 0:
        return None
    end = idx + len(marker)
    while end < len(header) and header[end].isdigit():
        end += 1
    aro = header[idx:end]
    if aro == marker:
        return None
    return aro


def load_label_map(tf: tarfile.TarFile, category: str) -> dict[str, str]:
    member = tf.extractfile("./aro_index.tsv")
    if member is None:
        raise ValueError("Missing ./aro_index.tsv in CARD archive")
    text = member.read().decode("utf-8", errors="ignore")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    out: dict[str, str] = {}

    category_column = {
        "Resistance Mechanism": "Resistance Mechanism",
        "Drug Class": "Drug Class",
        "AMR Gene Family": "AMR Gene Family",
    }[category]

    for row in reader:
        aro = (row.get("ARO Accession") or "").strip()
        name = (row.get(category_column) or "").strip()
        if category in {"Drug Class", "Resistance Mechanism"} and ";" in name:
            name = name.split(";", 1)[0].strip()
        if aro and name:
            out[aro] = name
    if not out:
        raise ValueError(f"No labels loaded for category: {category}")
    return out


def stable_hash_int(text: str) -> int:
    h = hashlib.md5(text.encode("utf-8")).digest()
    return int.from_bytes(h[:8], "big", signed=False)


def parse_rates(text: str) -> list[float]:
    if not text.strip():
        return []
    rates: list[float] = []
    for part in text.split(","):
        piece = part.strip()
        if not piece:
            continue
        val = float(piece)
        if not (0.0 < val < 1.0):
            raise ValueError(f"Mutation rates must be in (0,1): {val}")
        rates.append(val)
    return rates


def mutate_sequence(seq: str, rate: float, rng: random.Random) -> str:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    chars = list(seq.upper())
    for i, ch in enumerate(chars):
        if ch not in alphabet:
            continue
        if rng.random() < rate:
            choices = [aa for aa in alphabet if aa != ch]
            chars[i] = choices[rng.randrange(len(choices))]
    return "".join(chars)


def seq_to_embedding(seq: str, dim: int, k: int) -> list[float]:
    vec = [0.0] * dim
    clean = "".join(ch for ch in seq.upper() if "A" <= ch <= "Z")
    if len(clean) < k:
        return vec
    n = len(clean) - k + 1
    for i in range(n):
        token = clean[i : i + k]
        j = stable_hash_int(token) % dim
        vec[j] += 1.0
    inv = 1.0 / n
    return [v * inv for v in vec]


def write_fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as h:
        for seq_id, seq in rows:
            h.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                h.write(seq[i : i + 80] + "\n")


def write_label_csv(path: Path, id_col: str, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h)
        w.writerow([id_col, "label"])
        for seq_id, label in rows:
            w.writerow([seq_id, label])


def main() -> int:
    args = parse_args()

    import numpy as np

    with tarfile.open(args.card_archive, "r:bz2") as tf:
        label_map = load_label_map(tf, args.label_category)

        q_member = tf.extractfile(args.query_member)
        if q_member is None:
            raise ValueError(f"Missing query member: {args.query_member}")
        d_member = tf.extractfile(args.db_member)
        if d_member is None:
            raise ValueError(f"Missing db member: {args.db_member}")

        query_records = parse_fasta(q_member.read().decode("utf-8", errors="ignore"))
        if args.query_member_2:
            q_member2 = tf.extractfile(args.query_member_2)
            if q_member2 is None:
                raise ValueError(f"Missing query member 2: {args.query_member_2}")
            query_records.extend(parse_fasta(q_member2.read().decode("utf-8", errors="ignore")))
        query_records = dedupe_records(query_records)
        db_records = parse_fasta(d_member.read().decode("utf-8", errors="ignore"))

    def filter_and_label(records: list[tuple[str, str, str]]) -> list[tuple[str, str, str]]:
        out: list[tuple[str, str, str]] = []
        for seq_id, header, seq in records:
            aro = extract_aro_id(header)
            if not aro:
                continue
            label = label_map.get(aro)
            if not label:
                continue
            if len(seq) < args.kmer:
                continue
            out.append((seq_id, label, seq))
        return out

    query_labeled = filter_and_label(query_records)
    db_labeled = filter_and_label(db_records)

    db_count_by_label = Counter(label for _, label, _ in db_labeled)
    query_count_by_label = Counter(label for _, label, _ in query_labeled)
    shared_labels = {
        lab
        for lab, qn in query_count_by_label.items()
        if qn >= args.min_class_count and db_count_by_label.get(lab, 0) >= args.min_db_class_count
    }
    query_labeled = [r for r in query_labeled if r[1] in shared_labels]
    db_labeled = [r for r in db_labeled if r[1] in shared_labels]

    by_label_query: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for rec in query_labeled:
        by_label_query[rec[1]].append(rec)

    by_label_db: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for rec in db_labeled:
        by_label_db[rec[1]].append(rec)

    labels_sorted = sorted(shared_labels)

    selected_query: list[tuple[str, str, str]] = []
    selected_db: list[tuple[str, str, str]] = []

    if not labels_sorted:
        raise ValueError("No shared labels left after filtering; try lower --min-class-count")

    per_label_query_cap = max(1, args.max_query // len(labels_sorted))
    per_label_db_cap = max(1, args.max_db // len(labels_sorted))

    same_source = args.query_member == args.db_member or args.query_member_2 == args.db_member

    for lab in labels_sorted:
        q_rows = by_label_query[lab]
        d_rows = by_label_db[lab]

        if same_source:
            # Prefer query IDs not present in DB to avoid eroding DB class support.
            db_ids = {seq_id for seq_id, _, _ in d_rows}
            q_non_overlap = [r for r in q_rows if r[0] not in db_ids]
            q_overlap = [r for r in q_rows if r[0] in db_ids]

            take_non_overlap = min(per_label_query_cap, len(q_non_overlap))
            remaining = per_label_query_cap - take_non_overlap

            # Keep at least min_db_class_count examples per class in DB after overlap removal.
            max_overlap_take = max(0, len(d_rows) - args.min_db_class_count)
            take_overlap = min(remaining, max_overlap_take, len(q_overlap))

            q_pick = q_non_overlap[:take_non_overlap] + q_overlap[:take_overlap]
        else:
            q_pick = q_rows[:per_label_query_cap]
        selected_query.extend(q_pick)

        if same_source:
            selected_ids = {seq_id for seq_id, _, _ in q_pick}
            d_rows = [r for r in d_rows if r[0] not in selected_ids]

        selected_db.extend(d_rows[:per_label_db_cap])

    selected_query = selected_query[: args.max_query]
    selected_db = selected_db[: args.max_db]

    rates = parse_rates(args.augment_mutation_rates)
    if rates and args.augment_copies_per_rate > 0:
        augmented_rows: list[tuple[str, str, str]] = []
        for seq_id, label, seq in selected_query:
            for rate in rates:
                for copy_idx in range(args.augment_copies_per_rate):
                    seed_key = f"{args.augmentation_seed}|{seq_id}|{rate}|{copy_idx}"
                    rng = random.Random(stable_hash_int(seed_key))
                    mut_seq = mutate_sequence(seq, rate, rng)
                    mut_id = f"{seq_id}_mut{int(round(rate * 100)):02d}_{copy_idx + 1}"
                    augmented_rows.append((mut_id, label, mut_seq))
        selected_query.extend(augmented_rows)

    if len(selected_query) < 50:
        raise ValueError(f"Too few query samples after filtering: {len(selected_query)}")
    if len(selected_db) < 50:
        raise ValueError(f"Too few db samples after filtering: {len(selected_db)}")

    query_fasta_rows = [(seq_id, seq) for seq_id, _, seq in selected_query]
    db_fasta_rows = [(seq_id, seq) for seq_id, _, seq in selected_db]
    query_label_rows = [(seq_id, label) for seq_id, label, _ in selected_query]
    db_label_rows = [(seq_id, label) for seq_id, label, _ in selected_db]

    write_fasta(Path(args.out_query_fasta), query_fasta_rows)
    write_fasta(Path(args.out_db_fasta), db_fasta_rows)
    write_label_csv(Path(args.out_query_labels), "query_id", query_label_rows)
    write_label_csv(Path(args.out_db_labels), "db_id", db_label_rows)

    emb = np.array([seq_to_embedding(seq, dim=args.embedding_dim, k=args.kmer) for _, _, seq in selected_query], dtype=np.float32)
    out_emb = Path(args.out_embeddings)
    out_emb.parent.mkdir(parents=True, exist_ok=True)
    np.save(out_emb, emb)

    training_rows = [(idx, label) for idx, (_, label, _) in enumerate(selected_query)]
    out_train = Path(args.out_training_labels)
    out_train.parent.mkdir(parents=True, exist_ok=True)
    with out_train.open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h)
        w.writerow(["row_index", "label"])
        for idx, label in training_rows:
            w.writerow([idx, label])

    q_counts = Counter(label for _, label, _ in selected_query)
    d_counts = Counter(label for _, label, _ in selected_db)

    print(f"Prepared query proteins: {len(selected_query)}")
    print(f"Prepared db proteins: {len(selected_db)}")
    print(f"Label category: {args.label_category}")
    print(f"Classes used: {len(q_counts)}")
    for lab in sorted(q_counts):
        print(f"  {lab}: query={q_counts[lab]} db={d_counts.get(lab, 0)}")
    print(f"Embeddings shape: {emb.shape[0]} x {emb.shape[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
