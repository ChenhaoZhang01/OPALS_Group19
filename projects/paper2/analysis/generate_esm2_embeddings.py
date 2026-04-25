#!/usr/bin/env python3
"""Generate ESM2 embeddings for paper2 protein sequences.

Uses esm2_t6_8M_UR50D (320-dim, fastest CPU model).
Sequences longer than 1022 residues are truncated to 1022 (ESM2 position limit).
Output: results/embeddings/protein_embeddings.npy  (shape: n_seqs x 320, float32)
Row order matches the order of sequences in the input FASTA.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import torch


MAX_SEQ_LEN = 1022  # ESM2 hard limit (excludes BOS/EOS tokens)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    seqs: list[tuple[str, str]] = []
    header: str | None = None
    buf: list[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(buf)))
                header = line[1:]
                buf = []
            elif line:
                buf.append(line)
    if header is not None:
        seqs.append((header, "".join(buf)))
    return seqs


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate ESM2 embeddings")
    p.add_argument(
        "--fasta",
        default=str(Path(__file__).resolve().parents[1] / "proteins" / "card_query_homolog.faa"),
    )
    p.add_argument(
        "--out",
        default=str(Path(__file__).resolve().parents[1] / "results" / "embeddings" / "protein_embeddings.npy"),
    )
    p.add_argument(
        "--model", default="esm2_t6_8M_UR50D",
        help="ESM2 model name (esm2_t6_8M_UR50D=320d, esm2_t12_35M_UR50D=480d)",
    )
    p.add_argument("--batch-size", type=int, default=8)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    import esm  # noqa: PLC0415

    print(f"Loading model: {args.model}")
    loader = getattr(esm.pretrained, args.model)
    model, alphabet = loader()
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    # Identify the final representation layer
    num_layers = model.num_layers

    fasta_path = Path(args.fasta)
    sequences = read_fasta(fasta_path)
    print(f"Sequences loaded: {len(sequences)}")

    truncated = sum(1 for _, s in sequences if len(s) > MAX_SEQ_LEN)
    if truncated:
        print(f"Sequences truncated to {MAX_SEQ_LEN} aa: {truncated}")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_embeddings: list[np.ndarray] = []

    batch_size = args.batch_size
    for start in range(0, len(sequences), batch_size):
        batch = sequences[start : start + batch_size]
        batch_data = [
            (hdr, seq[:MAX_SEQ_LEN]) for hdr, seq in batch
        ]

        _, _, tokens = batch_converter(batch_data)

        with torch.no_grad():
            results = model(tokens, repr_layers=[num_layers], return_contacts=False)

        token_reps = results["representations"][num_layers]  # (B, L+2, D)

        for i, (_, seq) in enumerate(batch_data):
            seq_len = min(len(seq), tokens.shape[1] - 2)
            # mean-pool over residue positions, exclude BOS (0) and EOS (-1)
            emb = token_reps[i, 1 : seq_len + 1].mean(dim=0).float().numpy()
            all_embeddings.append(emb)

        done = min(start + batch_size, len(sequences))
        print(f"  {done}/{len(sequences)} sequences embedded", end="\r")

    print()
    embeddings = np.stack(all_embeddings, axis=0).astype(np.float32)
    print(f"Embedding matrix shape: {embeddings.shape}")
    np.save(out_path, embeddings)
    print(f"Saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
