#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run lag analysis for Paper 3")
    parser.add_argument("--input", required=True, help="CSV with study,sample_id,order,mge_abundance,entropy,arg_total")
    parser.add_argument("--model-out", required=True, help="Output CSV for lag regression")
    parser.add_argument("--corr-out", required=True, help="Output CSV for lag correlation")
    parser.add_argument(
        "--comparison-out",
        default="",
        help="Output CSV for model comparison (default: alongside model-out)",
    )
    parser.add_argument(
        "--direction-out",
        default="",
        help="Output CSV for forward vs reverse direction test (default: alongside model-out)",
    )
    parser.add_argument(
        "--difference-out",
        default="",
        help="Output CSV for first-difference model (default: alongside model-out)",
    )
    parser.add_argument(
        "--granger-out",
        default="",
        help="Output CSV for Granger-style added-value test (default: alongside model-out)",
    )
    parser.add_argument(
        "--metadata",
        default="",
        help="Optional metadata CSV with sample_id and read_count for sequencing_depth",
    )
    parser.add_argument(
        "--figures-dir",
        default="",
        help="Output directory for figures (default: ../analysis/figures)",
    )
    return parser.parse_args()


def _format_num(value: float | None, digits: int = 6) -> str:
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"


def _safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _lag_correlation(df, x_col: str, y_col: str) -> tuple[float | None, int]:
    sub = df[[x_col, y_col]].dropna()
    n = len(sub)
    if n < 2:
        return None, n
    return _safe_float(sub[x_col].corr(sub[y_col])), n


def _build_formula(response: str, predictor: str, include_entropy: bool, include_depth: bool) -> str:
    terms: list[str] = [predictor]
    if include_entropy:
        terms.append("entropy")
    terms.append("C(study)")
    if include_depth:
        terms.append("sequencing_depth")
    return f"{response} ~ " + " + ".join(terms)


def _extract_model_stats(model, predictor: str) -> dict[str, float | int | str | None]:
    ci = model.conf_int()
    ci_low = _safe_float(ci.loc[predictor, 0]) if predictor in ci.index else None
    ci_high = _safe_float(ci.loc[predictor, 1]) if predictor in ci.index else None
    return {
        "coefficient": _safe_float(model.params.get(predictor)),
        "std_error": _safe_float(model.bse.get(predictor)),
        "p_value": _safe_float(model.pvalues.get(predictor)),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "r_squared": _safe_float(model.rsquared),
        "adj_r_squared": _safe_float(model.rsquared_adj),
        "rows_used": int(model.nobs),
    }


def _empty_model_stats(rows_used: int = 0) -> dict[str, float | int | str | None]:
    return {
        "coefficient": None,
        "std_error": None,
        "p_value": None,
        "ci_low": None,
        "ci_high": None,
        "r_squared": None,
        "adj_r_squared": None,
        "rows_used": rows_used,
    }


def _write_key_value_csv(path: Path, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def _make_figures(figures_dir: Path, lag1_df, diff_df, corr_rows: list[tuple[str, float | None, int]]) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib is not installed; skipping figure generation.")
        return

    figures_dir.mkdir(parents=True, exist_ok=True)

    # 1) Scatter plot: MGE(t) vs ARG(t+1)
    scatter_df = lag1_df[["mge_abundance", "arg_t1"]].dropna()
    if len(scatter_df) >= 2:
        x = scatter_df["mge_abundance"].to_numpy()
        y = scatter_df["arg_t1"].to_numpy()
        m, b = np.polyfit(x, y, 1)

        plt.figure(figsize=(7, 5))
        plt.scatter(x, y, alpha=0.85)
        x_line = np.linspace(float(x.min()), float(x.max()), 100)
        plt.plot(x_line, m * x_line + b, linewidth=2)
        plt.xlabel("MGE(t)")
        plt.ylabel("ARG(t+1)")
        plt.title("Lagged Association: MGE(t) vs ARG(t+1)")
        plt.tight_layout()
        plt.savefig(figures_dir / "scatter_mge_t_vs_arg_t1.png", dpi=180)
        plt.close()

    # 2) Time-series plot per study
    studies = sorted(lag1_df["study"].dropna().unique())
    if studies:
        n = len(studies)
        fig, axes = plt.subplots(n, 1, figsize=(9, max(3 * n, 4)), squeeze=False)
        for i, study in enumerate(studies):
            ax = axes[i][0]
            sub = lag1_df[lag1_df["study"] == study].sort_values("order")
            ax.plot(sub["order"], sub["mge_abundance"], marker="o", label="MGE")
            ax.plot(sub["order"], sub["arg_total"], marker="s", label="ARG")
            ax.set_title(f"Study {study}")
            ax.set_xlabel("Time order")
            ax.set_ylabel("Abundance")
            ax.legend(loc="best")
        fig.suptitle("Time-Series Trajectories by Study")
        fig.tight_layout()
        fig.savefig(figures_dir / "timeseries_mge_arg_by_study.png", dpi=180)
        plt.close(fig)

    # 3) Lag comparison plot
    labels: list[str] = []
    values: list[float] = []
    for lag, corr, _n in corr_rows:
        if corr is not None:
            labels.append(lag)
            values.append(corr)
    if values:
        plt.figure(figsize=(7, 4))
        plt.bar(labels, values)
        plt.axhline(0.0, color="black", linewidth=1)
        plt.ylabel("Correlation")
        plt.title("Lag Comparison of MGE-ARG Correlation")
        plt.figtext(
            0.5,
            0.01,
            "High correlation across all lags indicates temporal coupling, not directional prediction",
            ha="center",
            fontsize=9,
        )
        plt.tight_layout()
        plt.savefig(figures_dir / "lag_comparison_correlation.png", dpi=180)
        plt.close()

    # 4) Differenced scatter: dMGE(t) vs dARG(t+1)
    if diff_df is not None and len(diff_df) >= 2:
        scatter_diff = diff_df[["d_mge_t", "d_arg_t1"]].dropna()
        if len(scatter_diff) >= 2:
            x = scatter_diff["d_mge_t"].to_numpy()
            y = scatter_diff["d_arg_t1"].to_numpy()
            m, b = np.polyfit(x, y, 1)

            plt.figure(figsize=(7, 5))
            plt.scatter(x, y, alpha=0.85)
            x_line = np.linspace(float(x.min()), float(x.max()), 100)
            plt.plot(x_line, m * x_line + b, linewidth=2)
            plt.xlabel("dMGE(t) = MGE(t) - MGE(t-1)")
            plt.ylabel("dARG(t+1) = ARG(t+1) - ARG(t)")
            plt.title("Differenced Association: dMGE(t) vs dARG(t+1)")
            plt.tight_layout()
            plt.savefig(figures_dir / "differenced_scatter_dmge_t_vs_darg_t1.png", dpi=180)
            plt.close()


def main() -> int:
    args = parse_args()

    try:
        import pandas as pd
        import statsmodels.formula.api as smf
        from statsmodels.stats.anova import anova_lm
    except ImportError as exc:
        raise SystemExit(
            "This script needs pandas and statsmodels. Install with: pip install pandas statsmodels"
        ) from exc

    df = pd.read_csv(args.input)
    required = ["study", "sample_id", "order", "mge_abundance", "entropy", "arg_total"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")

    df = df.copy()
    df["order"] = pd.to_numeric(df["order"], errors="coerce")
    df["mge_abundance"] = pd.to_numeric(df["mge_abundance"], errors="coerce")
    df["entropy"] = pd.to_numeric(df["entropy"], errors="coerce")
    df["arg_total"] = pd.to_numeric(df["arg_total"], errors="coerce")

    if "sequencing_depth" not in df.columns:
        df["sequencing_depth"] = None

    metadata_path = Path(args.metadata) if args.metadata else Path(__file__).resolve().parents[1] / "metadata" / "metadata_final.csv"
    if metadata_path.exists():
        meta = pd.read_csv(metadata_path)
        if "sample_id" in meta.columns and "read_count" in meta.columns:
            merged = df.merge(meta[["sample_id", "read_count"]], on="sample_id", how="left")
            if "sequencing_depth" in merged.columns:
                merged["sequencing_depth"] = merged["sequencing_depth"].fillna(merged["read_count"])
            else:
                merged["sequencing_depth"] = merged["read_count"]
            df = merged.drop(columns=["read_count"], errors="ignore")

    # Sequencing depth is used as a raw covariate (read count scale), not log-transformed.
    df["sequencing_depth"] = pd.to_numeric(df["sequencing_depth"], errors="coerce")
    sequencing_depth_transform = "raw_read_count"

    df = df.dropna(subset=["order", "mge_abundance", "entropy", "arg_total"])

    df = df.sort_values(["study", "order"])
    df["arg_t1"] = df.groupby("study")["arg_total"].shift(-1)
    df["arg_t2"] = df.groupby("study")["arg_total"].shift(-2)
    df["mge_t1"] = df.groupby("study")["mge_abundance"].shift(-1)
    df["arg_tminus1"] = df.groupby("study")["arg_total"].shift(1)
    df["mge_tminus1"] = df.groupby("study")["mge_abundance"].shift(1)
    df["d_mge_t"] = df["mge_abundance"] - df["mge_tminus1"]
    df["d_arg_t1"] = df["arg_t1"] - df["arg_total"]
    lag_df = df.dropna(subset=["arg_t1"]).copy()
    lag2_df = df.dropna(subset=["arg_t2"]).copy()
    diff_df = df.dropna(subset=["d_mge_t", "d_arg_t1"]).copy()

    if lag_df.empty:
        raise ValueError("No lag pairs available after shift. Check study/order values.")

    include_depth = lag_df["sequencing_depth"].notna().sum() >= 3
    forward_formula = _build_formula(
        response="arg_t1",
        predictor="mge_abundance",
        include_entropy=False,
        include_depth=include_depth,
    )
    forward_model = smf.ols(forward_formula, data=lag_df).fit()
    forward_stats = _extract_model_stats(forward_model, "mge_abundance")

    enriched_formula = _build_formula(
        response="arg_t1",
        predictor="mge_abundance",
        include_entropy=True,
        include_depth=include_depth,
    )
    enriched_model = smf.ols(enriched_formula, data=lag_df).fit()
    enriched_stats = _extract_model_stats(enriched_model, "mge_abundance")

    reverse_formula = _build_formula(
        response="mge_t1",
        predictor="arg_total",
        include_entropy=False,
        include_depth=include_depth,
    )
    reverse_df = df.dropna(subset=["mge_t1", "arg_total"]).copy()
    reverse_model = smf.ols(reverse_formula, data=reverse_df).fit()
    reverse_stats = _extract_model_stats(reverse_model, "arg_total")

    # First-difference model to reduce trend/autocorrelation confounding.
    diff_formula = "d_arg_t1 ~ d_mge_t"
    if len(diff_df) >= 3:
        diff_model = smf.ols(diff_formula, data=diff_df).fit()
        diff_stats = _extract_model_stats(diff_model, "d_mge_t")
    else:
        diff_stats = _empty_model_stats(rows_used=len(diff_df))

    # Granger-style added-value test:
    # Base: ARG(t+1) ~ ARG(t) (+ controls)
    # Full: ARG(t+1) ~ ARG(t) + MGE(t) (+ controls)
    granger_df = df.dropna(subset=["arg_t1", "arg_total", "mge_abundance"]).copy()
    if include_depth:
        granger_base_formula = "arg_t1 ~ arg_total + C(study) + sequencing_depth"
        granger_full_formula = "arg_t1 ~ arg_total + mge_abundance + C(study) + sequencing_depth"
    else:
        granger_base_formula = "arg_t1 ~ arg_total + C(study)"
        granger_full_formula = "arg_t1 ~ arg_total + mge_abundance + C(study)"

    if len(granger_df) >= 5:
        granger_base_model = smf.ols(granger_base_formula, data=granger_df).fit()
        granger_full_model = smf.ols(granger_full_formula, data=granger_df).fit()
        nested = anova_lm(granger_base_model, granger_full_model)
        granger_f = _safe_float(nested.loc[1, "F"]) if 1 in nested.index else None
        granger_p = _safe_float(nested.loc[1, "Pr(>F)"]) if 1 in nested.index else None
        granger_delta_r2 = None
        base_r2 = _safe_float(granger_base_model.rsquared)
        full_r2 = _safe_float(granger_full_model.rsquared)
        if base_r2 is not None and full_r2 is not None:
            granger_delta_r2 = full_r2 - base_r2
        granger_mge_stats = _extract_model_stats(granger_full_model, "mge_abundance")
        granger_n = int(granger_full_model.nobs)
    else:
        granger_f = None
        granger_p = None
        granger_delta_r2 = None
        granger_mge_stats = _empty_model_stats(rows_used=len(granger_df))
        granger_n = len(granger_df)

    corr_t_t, n_t_t = _lag_correlation(df, "mge_abundance", "arg_total")
    corr_t_t1, n_t_t1 = _lag_correlation(df, "mge_abundance", "arg_t1")
    corr_t_t2, n_t_t2 = _lag_correlation(df, "mge_abundance", "arg_t2")
    corr_rows = [
        ("t_to_t", corr_t_t, n_t_t),
        ("t_to_t_plus_1", corr_t_t1, n_t_t1),
        ("t_to_t_plus_2", corr_t_t2, n_t_t2),
    ]

    model_out = Path(args.model_out)
    corr_out = Path(args.corr_out)
    comparison_out = Path(args.comparison_out) if args.comparison_out else model_out.with_name("model_comparison.csv")
    direction_out = Path(args.direction_out) if args.direction_out else model_out.with_name("directionality_test.csv")
    difference_out = Path(args.difference_out) if args.difference_out else model_out.with_name("difference_model.csv")
    granger_out = Path(args.granger_out) if args.granger_out else model_out.with_name("granger_test.csv")
    figures_dir = (
        Path(args.figures_dir)
        if args.figures_dir
        else Path(__file__).resolve().parent / "figures"
    )

    rows_used = int(forward_stats["rows_used"] or 0)
    proof_of_concept_mode = rows_used < 30

    _write_key_value_csv(
        model_out,
        [
            ("formula_forward", forward_formula),
            ("slope_mge_to_arg_t1", _format_num(forward_stats["coefficient"])),
            ("std_error_mge_to_arg_t1", _format_num(forward_stats["std_error"])),
            ("p_value_mge_to_arg_t1", _format_num(forward_stats["p_value"], digits=8)),
            ("ci_low_mge_to_arg_t1", _format_num(forward_stats["ci_low"])),
            ("ci_high_mge_to_arg_t1", _format_num(forward_stats["ci_high"])),
            ("r_squared", _format_num(forward_stats["r_squared"])),
            ("adjusted_r_squared", _format_num(forward_stats["adj_r_squared"])),
            ("rows_used", str(rows_used)),
            ("sequencing_depth_transform", sequencing_depth_transform),
            ("recommended_interpretation", "proof_of_concept" if proof_of_concept_mode else "predictive_signal_test"),
            ("sample_size_flag", "small_n_lt_30" if proof_of_concept_mode else "adequate_n"),
        ],
    )

    _write_csv(
        comparison_out,
        ["model", "formula", "r_squared", "adjusted_r_squared", "n"],
        [
            [
                "MGE_only_with_controls",
                forward_formula,
                _format_num(forward_stats["r_squared"]),
                _format_num(forward_stats["adj_r_squared"]),
                str(rows_used),
            ],
            [
                "MGE_plus_entropy_with_controls",
                enriched_formula,
                _format_num(enriched_stats["r_squared"]),
                _format_num(enriched_stats["adj_r_squared"]),
                str(int(enriched_stats["rows_used"] or 0)),
            ],
            [
                "difference_model",
                diff_formula,
                _format_num(diff_stats["r_squared"]),
                _format_num(diff_stats["adj_r_squared"]),
                str(int(diff_stats["rows_used"] or 0)),
            ],
            [
                "granger_base_arg_history",
                granger_base_formula,
                _format_num(_safe_float(granger_base_model.rsquared) if len(granger_df) >= 5 else None),
                _format_num(_safe_float(granger_base_model.rsquared_adj) if len(granger_df) >= 5 else None),
                str(granger_n),
            ],
            [
                "granger_full_arg_history_plus_mge",
                granger_full_formula,
                _format_num(_safe_float(granger_full_model.rsquared) if len(granger_df) >= 5 else None),
                _format_num(_safe_float(granger_full_model.rsquared_adj) if len(granger_df) >= 5 else None),
                str(granger_n),
            ],
        ],
    )

    _write_csv(
        difference_out,
        ["formula", "coefficient", "std_error", "p_value", "ci_low", "ci_high", "r_squared", "adjusted_r_squared", "n"],
        [
            [
                diff_formula,
                _format_num(diff_stats["coefficient"]),
                _format_num(diff_stats["std_error"]),
                _format_num(diff_stats["p_value"], digits=8),
                _format_num(diff_stats["ci_low"]),
                _format_num(diff_stats["ci_high"]),
                _format_num(diff_stats["r_squared"]),
                _format_num(diff_stats["adj_r_squared"]),
                str(int(diff_stats["rows_used"] or 0)),
            ]
        ],
    )

    _write_csv(
        granger_out,
        [
            "base_formula",
            "full_formula",
            "base_r_squared",
            "base_adjusted_r_squared",
            "full_r_squared",
            "full_adjusted_r_squared",
            "delta_r_squared",
            "added_value_f_stat",
            "added_value_p_value",
            "full_model_mge_coefficient",
            "full_model_mge_p_value",
            "n",
        ],
        [
            [
                granger_base_formula,
                granger_full_formula,
                _format_num(_safe_float(granger_base_model.rsquared) if len(granger_df) >= 5 else None),
                _format_num(_safe_float(granger_base_model.rsquared_adj) if len(granger_df) >= 5 else None),
                _format_num(_safe_float(granger_full_model.rsquared) if len(granger_df) >= 5 else None),
                _format_num(_safe_float(granger_full_model.rsquared_adj) if len(granger_df) >= 5 else None),
                _format_num(granger_delta_r2),
                _format_num(granger_f),
                _format_num(granger_p, digits=8),
                _format_num(granger_mge_stats["coefficient"]),
                _format_num(granger_mge_stats["p_value"], digits=8),
                str(granger_n),
            ]
        ],
    )

    _write_csv(
        direction_out,
        ["direction", "formula", "coefficient", "std_error", "p_value", "ci_low", "ci_high", "r_squared", "adjusted_r_squared", "n"],
        [
            [
                "forward_mge_t_to_arg_t1",
                forward_formula,
                _format_num(forward_stats["coefficient"]),
                _format_num(forward_stats["std_error"]),
                _format_num(forward_stats["p_value"], digits=8),
                _format_num(forward_stats["ci_low"]),
                _format_num(forward_stats["ci_high"]),
                _format_num(forward_stats["r_squared"]),
                _format_num(forward_stats["adj_r_squared"]),
                str(rows_used),
            ],
            [
                "reverse_arg_t_to_mge_t1",
                reverse_formula,
                _format_num(reverse_stats["coefficient"]),
                _format_num(reverse_stats["std_error"]),
                _format_num(reverse_stats["p_value"], digits=8),
                _format_num(reverse_stats["ci_low"]),
                _format_num(reverse_stats["ci_high"]),
                _format_num(reverse_stats["r_squared"]),
                _format_num(reverse_stats["adj_r_squared"]),
                str(int(reverse_stats["rows_used"] or 0)),
            ],
        ],
    )

    _write_csv(
        corr_out,
        ["lag", "correlation", "n"],
        [[lag, _format_num(corr), str(n)] for lag, corr, n in corr_rows],
    )

    _make_figures(figures_dir, lag_df, diff_df, corr_rows)

    print(f"Rows used: {rows_used}")
    if proof_of_concept_mode:
        print("Warning: rows_used < 30, interpret as proof-of-concept only.")
    print(f"Wrote model summary: {model_out}")
    print(f"Wrote model comparison: {comparison_out}")
    print(f"Wrote directionality test: {direction_out}")
    print(f"Wrote difference model: {difference_out}")
    print(f"Wrote Granger-style test: {granger_out}")
    print(f"Wrote correlation: {corr_out}")
    print(f"Figures dir: {figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
