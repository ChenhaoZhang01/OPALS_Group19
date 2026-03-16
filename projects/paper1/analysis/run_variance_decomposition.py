#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run variance decomposition for Paper 1")
    parser.add_argument("--input", required=True, help="Long table CSV with columns sample_id, environment, pipeline, ARG_total")
    parser.add_argument("--output", required=True, help="Output CSV for variance decomposition")
    return parser.parse_args()


def ensure_columns(fieldnames: list[str], required: list[str]) -> None:
    missing = [c for c in required if c not in fieldnames]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")


def main() -> int:
    args = parse_args()

    with open(args.input, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        if not reader.fieldnames:
            raise ValueError("Input file has no header")
        ensure_columns(reader.fieldnames, ["sample_id", "environment", "pipeline", "ARG_total"])

    if not rows:
        raise ValueError("Input file has no data rows")

    try:
        import pandas as pd
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
    except ImportError as exc:
        raise SystemExit(
            "This script needs pandas and statsmodels. Install with: pip install pandas statsmodels"
        ) from exc

    df = pd.DataFrame(rows)
    df["ARG_total"] = pd.to_numeric(df["ARG_total"], errors="coerce")
    df = df.dropna(subset=["ARG_total"]).copy()

    if df.empty:
        raise ValueError("No valid numeric ARG_total values found")

    model = smf.ols(
        "ARG_total ~ C(pipeline) + C(environment) + C(pipeline):C(environment)",
        data=df,
    ).fit()

    anova = sm.stats.anova_lm(model, typ=2)
    ss_total = anova["sum_sq"].sum()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["factor", "sum_sq", "variance_percent"])
        for factor, row in anova.iterrows():
            sum_sq = float(row["sum_sq"])
            pct = (sum_sq / ss_total * 100.0) if ss_total > 0 else 0.0
            writer.writerow([factor, f"{sum_sq:.6f}", f"{pct:.4f}"])

    print(f"Rows analyzed: {len(df)}")
    print(f"Wrote variance decomposition: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
