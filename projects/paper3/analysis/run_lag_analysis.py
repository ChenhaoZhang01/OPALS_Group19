#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run lag analysis for Paper 3")
    parser.add_argument("--input", required=True, help="CSV with study,sample_id,order,mge_abundance,entropy,arg_total")
    parser.add_argument("--model-out", required=True, help="Output CSV for lag regression")
    parser.add_argument("--corr-out", required=True, help="Output CSV for lag correlation")
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
    required = ["study", "sample_id", "order", "mge_abundance", "arg_total"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")

    df = df.copy()
    df["order"] = pd.to_numeric(df["order"], errors="coerce")
    df["mge_abundance"] = pd.to_numeric(df["mge_abundance"], errors="coerce")
    df["arg_total"] = pd.to_numeric(df["arg_total"], errors="coerce")
    df = df.dropna(subset=["order", "mge_abundance", "arg_total"])

    df = df.sort_values(["study", "order"])
    df["arg_t1"] = df.groupby("study")["arg_total"].shift(-1)
    lag_df = df.dropna(subset=["arg_t1"]).copy()

    if lag_df.empty:
        raise ValueError("No lag pairs available after shift. Check study/order values.")

    model = smf.ols("arg_t1 ~ mge_abundance", data=lag_df).fit()
    slope = float(model.params.get("mge_abundance", 0.0))
    p_value = float(model.pvalues.get("mge_abundance", 1.0))
    r_squared = float(model.rsquared)

    corr = float(lag_df["mge_abundance"].corr(lag_df["arg_t1"]))

    os.makedirs(os.path.dirname(args.model_out), exist_ok=True)
    with open(args.model_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerow(["slope_mge_to_arg_t1", f"{slope:.6f}"])
        writer.writerow(["p_value", f"{p_value:.6g}"])
        writer.writerow(["r_squared", f"{r_squared:.6f}"])
        writer.writerow(["rows_used", str(len(lag_df))])

    os.makedirs(os.path.dirname(args.corr_out), exist_ok=True)
    with open(args.corr_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["lag", "correlation"])
        writer.writerow(["t_to_t_plus_1", f"{corr:.6f}"])

    print(f"Rows used: {len(lag_df)}")
    print(f"Wrote model summary: {args.model_out}")
    print(f"Wrote correlation: {args.corr_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
