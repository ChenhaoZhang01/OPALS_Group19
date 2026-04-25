#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Populate Paper 3 real feature template from a quantification table",
    )
    parser.add_argument(
        "--template",
        default="projects/paper3/results/features_real_template.csv",
        help="Input real-ID feature template CSV",
    )
    parser.add_argument(
        "--quant",
        required=True,
        help="Quantification CSV containing sample IDs and metric columns",
    )
    parser.add_argument(
        "--sample-col",
        default="sample_id",
        help="Sample ID column name in quant CSV",
    )
    parser.add_argument(
        "--arg-col",
        default="arg_total",
        help="ARG total abundance column in quant CSV",
    )
    parser.add_argument(
        "--mge-col",
        default="mge_abundance",
        help="MGE abundance column in quant CSV",
    )
    parser.add_argument(
        "--entropy-col",
        default="entropy",
        help="Entropy column in quant CSV",
    )
    parser.add_argument(
        "--out",
        default="projects/paper3/results/features_real_populated.csv",
        help="Output populated feature CSV",
    )
    parser.add_argument(
        "--missing-out",
        default="projects/paper3/results/features_real_missing_report.csv",
        help="Output missing-data report CSV",
    )
    return parser.parse_args()


def _clean(value: str | None) -> str:
    return (value or "").strip()


def _safe_float_text(value: str | None) -> str:
    text = _clean(value)
    if not text:
        return ""
    try:
        return str(float(text))
    except ValueError:
        return ""


def _read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return rows, list(reader.fieldnames or [])


def main() -> int:
    args = parse_args()

    template_path = Path(args.template)
    quant_path = Path(args.quant)
    out_path = Path(args.out)
    missing_path = Path(args.missing_out)

    template_rows, template_fields = _read_csv(template_path)
    quant_rows, quant_fields = _read_csv(quant_path)

    required_quant_cols = [args.sample_col, args.arg_col, args.mge_col, args.entropy_col]
    missing_quant_cols = [c for c in required_quant_cols if c not in quant_fields]
    if missing_quant_cols:
        raise ValueError(
            "Quant CSV missing required columns: " + ", ".join(missing_quant_cols)
        )

    by_sample: dict[str, dict[str, str]] = {}
    for row in quant_rows:
        sample_id = _clean(row.get(args.sample_col))
        if sample_id:
            by_sample[sample_id] = row

    populated = 0
    missing_rows: list[dict[str, str]] = []

    for row in template_rows:
        sample_id = _clean(row.get("sample_id"))
        q = by_sample.get(sample_id)
        if not q:
            row["metric_status"] = "missing_quant_row"
            missing_rows.append(
                {
                    "sample_id": sample_id,
                    "study": _clean(row.get("study")),
                    "reason": "sample_id_not_found_in_quant",
                }
            )
            continue

        arg_value = _safe_float_text(q.get(args.arg_col))
        mge_value = _safe_float_text(q.get(args.mge_col))
        entropy_value = _safe_float_text(q.get(args.entropy_col))

        row["arg_total"] = arg_value
        row["mge_abundance"] = mge_value
        row["entropy"] = entropy_value

        if arg_value and mge_value and entropy_value:
            row["metric_status"] = "ok"
            populated += 1
        else:
            row["metric_status"] = "partial_missing_metric"
            missing = []
            if not arg_value:
                missing.append("arg_total")
            if not mge_value:
                missing.append("mge_abundance")
            if not entropy_value:
                missing.append("entropy")
            missing_rows.append(
                {
                    "sample_id": sample_id,
                    "study": _clean(row.get("study")),
                    "reason": "missing_metric_fields:" + "|".join(missing),
                }
            )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=template_fields)
        writer.writeheader()
        writer.writerows(template_rows)

    missing_path.parent.mkdir(parents=True, exist_ok=True)
    with missing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "study", "reason"])
        writer.writeheader()
        writer.writerows(missing_rows)

    print(f"Template rows: {len(template_rows)}")
    print(f"Quant rows: {len(quant_rows)}")
    print(f"Rows fully populated: {populated}")
    print(f"Rows with missing data: {len(missing_rows)}")
    print(f"Wrote populated features: {out_path}")
    print(f"Wrote missing report: {missing_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
