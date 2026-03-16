#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run difference-in-differences model for Paper 4")
    parser.add_argument("--input", required=True, help="CSV with sample_id,ARG_total,treatment,time,plant_id")
    parser.add_argument("--coef-out", required=True, help="CSV output for model coefficients")
    parser.add_argument("--summary-out", required=True, help="Text output for model summary")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        import pandas as pd
        import statsmodels.formula.api as smf
    except ImportError as exc:
        raise SystemExit(
            "This script needs pandas and statsmodels. Install with: pip install pandas statsmodels"
        ) from exc

    df = pd.read_csv(args.input)
    required = ["ARG_total", "treatment", "time"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")

    df = df.copy()
    df["ARG_total"] = pd.to_numeric(df["ARG_total"], errors="coerce")
    df["treatment"] = pd.to_numeric(df["treatment"], errors="coerce")
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df = df.dropna(subset=["ARG_total", "treatment", "time"])

    if df.empty:
        raise ValueError("No valid rows available for model fitting")

    model = smf.ols("ARG_total ~ treatment + time + treatment:time", data=df).fit()

    os.makedirs(os.path.dirname(args.coef_out), exist_ok=True)
    with open(args.coef_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["term", "coefficient", "p_value"])
        for term in model.params.index:
            coef = float(model.params[term])
            pval = float(model.pvalues[term])
            writer.writerow([term, f"{coef:.6f}", f"{pval:.6g}"])

    os.makedirs(os.path.dirname(args.summary_out), exist_ok=True)
    with open(args.summary_out, "w", encoding="utf-8") as handle:
        handle.write(model.summary().as_text())

    print(f"Rows analyzed: {len(df)}")
    print(f"Wrote coefficients: {args.coef_out}")
    print(f"Wrote summary: {args.summary_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
