#!/usr/bin/env python3
"""Run a Difference-in-Differences OLS model for Paper 4.

Expected input columns:
- ARG_total (numeric outcome)
- treatment (0/1)
- time (0/1)

Outputs:
- Coefficient table CSV
- Human-readable model summary TXT
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


REQUIRED_COLUMNS = ["ARG_total", "treatment", "time"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run DiD model: ARG_total ~ treatment + time + treatment:time"
    )
    parser.add_argument("--input", required=True, help="Path to did_input_table.csv")
    parser.add_argument(
        "--coef-out",
        required=True,
        help="Path to write coefficient table CSV",
    )
    parser.add_argument(
        "--summary-out",
        required=True,
        help="Path to write model summary TXT",
    )
    parser.add_argument(
        "--robust",
        action="store_true",
        help="Use HC3 robust standard errors",
    )
    parser.add_argument(
        "--log1p",
        action="store_true",
        help="Fit model on log1p(ARG_total) instead of ARG_total",
    )
    return parser.parse_args()


def validate_and_prepare(df: pd.DataFrame, use_log1p: bool) -> pd.DataFrame:
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    work = df.copy()

    # Enforce numeric coding and fail early if values are invalid.
    for col in REQUIRED_COLUMNS:
        work[col] = pd.to_numeric(work[col], errors="coerce")

    null_mask = work[REQUIRED_COLUMNS].isnull().any(axis=1)
    if null_mask.any():
        bad_rows = work.loc[null_mask, ["sample_id"] + REQUIRED_COLUMNS if "sample_id" in work.columns else REQUIRED_COLUMNS]
        raise ValueError(
            "Input contains non-numeric or missing values in required columns. "
            f"Problem rows:\n{bad_rows.to_string(index=False)}"
        )

    if use_log1p:
        if (work["ARG_total"] < 0).any():
            raise ValueError("Cannot use --log1p with negative ARG_total values.")
        work["ARG_model"] = np.log1p(work["ARG_total"])
        formula = "ARG_model ~ treatment + time + treatment:time"
    else:
        formula = "ARG_total ~ treatment + time + treatment:time"

    work["_formula"] = formula
    return work


def fit_model(df: pd.DataFrame, robust: bool):
    formula = df["_formula"].iloc[0]
    model = smf.ols(formula, data=df).fit()
    if robust:
        model = model.get_robustcov_results(cov_type="HC3")
    return model


def coefficient_table(model) -> pd.DataFrame:
    terms = list(model.model.exog_names)
    params = np.asarray(model.params)
    bse = np.asarray(model.bse)
    tvals = np.asarray(model.tvalues)
    pvals = np.asarray(model.pvalues)
    conf = np.asarray(model.conf_int())

    table = pd.DataFrame(
        {
            "term": terms,
            "coef": params,
            "std_err": bse,
            "t": tvals,
            "p_value": pvals,
            "ci_low": conf[:, 0],
            "ci_high": conf[:, 1],
        }
    )

    table["nobs"] = float(model.nobs)
    table["r2"] = float(model.rsquared)
    table["adj_r2"] = float(model.rsquared_adj)
    return table


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    coef_out = Path(args.coef_out)
    summary_out = Path(args.summary_out)

    df = pd.read_csv(input_path)
    prepared = validate_and_prepare(df, use_log1p=args.log1p)
    model = fit_model(prepared, robust=args.robust)

    coef_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.parent.mkdir(parents=True, exist_ok=True)

    coef_df = coefficient_table(model)
    coef_df.to_csv(coef_out, index=False)

    model_header = [
        "Paper 4 Difference-in-Differences Model",
        f"Input: {input_path}",
        f"Formula: {prepared['_formula'].iloc[0]}",
        f"Robust HC3: {args.robust}",
        f"log1p outcome: {args.log1p}",
        "",
    ]
    summary_text = "\n".join(model_header) + model.summary().as_text() + "\n"
    summary_out.write_text(summary_text, encoding="utf-8")

    print(f"Wrote coefficients: {coef_out}")
    print(f"Wrote summary: {summary_out}")


if __name__ == "__main__":
    main()
