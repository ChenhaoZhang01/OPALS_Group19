#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Paper 3 lag publishability report")
    parser.add_argument(
        "--lag-samples",
        default="projects/paper3/results/lag_samples_actual.csv",
        help="Lag sample usability table",
    )
    parser.add_argument(
        "--lag-pairs-strict",
        default="projects/paper3/results/lag_pairs_strict_actual.csv",
        help="Strict lag pair table",
    )
    parser.add_argument(
        "--features-real-template",
        default="projects/paper3/results/features_real_template.csv",
        help="Real-ID feature template",
    )
    parser.add_argument(
        "--legacy-features",
        default="projects/paper3/results/features_table.csv",
        help="Legacy feature table to check for synthetic IDs",
    )
    parser.add_argument(
        "--metadata-enriched",
        default="projects/paper3/metadata/metadata_enriched.csv",
        help="Enriched metadata",
    )
    parser.add_argument(
        "--out",
        default="projects/paper3/results/publishability_report.md",
        help="Markdown report output",
    )
    return parser.parse_args()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _status(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def main() -> int:
    args = parse_args()

    lag_samples = _read_csv(Path(args.lag_samples))
    lag_pairs_strict = _read_csv(Path(args.lag_pairs_strict))
    features_real = _read_csv(Path(args.features_real_template))
    legacy_features = _read_csv(Path(args.legacy_features))
    metadata = _read_csv(Path(args.metadata_enriched))

    strict_samples = [r for r in lag_samples if r.get("usable_strict") == "yes"]
    strict_studies = Counter(r.get("study", "") for r in strict_samples)

    strict_pairs_ok = len(lag_pairs_strict) >= 10
    strict_studies_ok = len(strict_studies) >= 3

    # Real quantification check.
    arg_present = sum(1 for r in features_real if (r.get("arg_total") or "").strip())
    mge_present = sum(1 for r in features_real if (r.get("mge_abundance") or "").strip())
    entropy_present = sum(1 for r in features_real if (r.get("entropy") or "").strip())
    quant_ok = min(arg_present, mge_present, entropy_present) >= 10

    metadata_ids = {(r.get("sample_id") or "").strip() for r in metadata if (r.get("sample_id") or "").strip()}
    legacy_ids = {(r.get("sample_id") or "").strip() for r in legacy_features if (r.get("sample_id") or "").strip()}
    legacy_overlap = len(metadata_ids & legacy_ids)
    legacy_ok = legacy_overlap > 0

    lines: list[str] = []
    lines.append("# Paper 3 Publishability Report (Lag Analysis)")
    lines.append("")
    lines.append("## Summary")
    lines.append(f"- Strict lag pairs: {len(lag_pairs_strict)}")
    lines.append(f"- Strict lag studies: {len(strict_studies)}")
    lines.append(f"- Real-feature template rows: {len(features_real)}")
    lines.append(f"- Legacy feature ID overlap with real metadata: {legacy_overlap}")
    lines.append("")

    lines.append("## Gate Checks")
    lines.append(f"- Temporal ordering has sufficient strict pairs (>=10): {_status(strict_pairs_ok)}")
    lines.append(f"- Temporal ordering covers >=3 studies: {_status(strict_studies_ok)}")
    lines.append(f"- Real ARG/MGE/entropy quantification available (>=10 complete rows): {_status(quant_ok)}")
    lines.append(f"- Legacy feature table IDs align with real metadata: {_status(legacy_ok)}")
    lines.append("")

    lines.append("## Blocking Issues")
    if not strict_pairs_ok:
        lines.append(f"- Not enough strict t->t+1 lag pairs ({len(lag_pairs_strict)} found).")
    if not strict_studies_ok:
        lines.append(f"- Strict lag coverage too narrow ({len(strict_studies)} studies found).")
    if not quant_ok:
        lines.append(
            "- Real lag feature quantification missing: at least one of arg_total, mge_abundance, entropy is largely empty in features_real_template.csv."
        )
    if not legacy_ok:
        lines.append(
            "- Legacy features_table.csv appears synthetic or unmapped to real metadata (zero sample_id overlap)."
        )
    if strict_pairs_ok and strict_studies_ok and quant_ok and legacy_ok:
        lines.append("- None. All gate checks passed.")
    lines.append("")

    lines.append("## Immediate Next Actions")
    lines.append("1. Populate real arg_total, mge_abundance, and entropy for the strict sample set in features_real_template.csv.")
    lines.append("2. Re-run run_lag_analysis.py using the real populated feature table and archive outputs under a new filename set.")
    lines.append("3. Remove or clearly label legacy synthetic outputs to avoid provenance confusion during submission.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote publishability report: {out_path}")
    print(f"Strict pairs: {len(lag_pairs_strict)}")
    print(f"Strict studies: {len(strict_studies)}")
    print(f"Quantified rows (arg/mge/entropy): {arg_present}/{mge_present}/{entropy_present}")
    print(f"Legacy overlap count: {legacy_overlap}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
