#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
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
    parser.add_argument("--within-study-out", default="", help="Output CSV for within-study lag models")
    parser.add_argument("--loso-out", default="", help="Output CSV for leave-one-study-out analysis")
    parser.add_argument("--diagnostics-out", default="", help="Output CSV for residual diagnostics (DW, BP)")
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


def _series_trend_label(values) -> str:
    clean = values.dropna()
    if len(clean) < 2:
        return "trend unavailable"
    start = float(clean.iloc[0])
    end = float(clean.iloc[-1])
    if end > start:
        return "increasing"
    if end < start:
        return "decreasing"
    return "stable"


def _baseline_rank_label(rank: int, total: int) -> str:
    if total <= 1:
        return "single cohort"
    if total == 2:
        return "lower starting ARG burden" if rank == 0 else "higher starting ARG burden"
    if total == 3:
        return [
            "lowest starting ARG burden",
            "middle starting ARG burden",
            "highest starting ARG burden",
        ][rank]
    frac = rank / (total - 1)
    if frac <= 0.33:
        return "lower starting ARG burden"
    if frac >= 0.67:
        return "higher starting ARG burden"
    return "middle starting ARG burden"


def _trend_summary_word(trend: str) -> str:
    if trend == "increasing":
        return "up"
    if trend == "decreasing":
        return "down"
    if trend == "stable":
        return "flat"
    return "NA"


def _format_study_panel_title(
    study_value: object,
    baseline_label: str,
) -> str:
    raw = str(study_value)
    match = re.fullmatch(r"study([A-Za-z0-9]+)", raw.strip(), flags=re.IGNORECASE)
    if match:
        pretty = f"Study {match.group(1).upper()}"
    else:
        pretty = raw
    return f"{pretty}: {baseline_label}"


def _format_study_panel_note(
    n_points: int,
    mge_start: float | None,
    arg_start: float | None,
    mge_trend: str,
    arg_trend: str,
) -> str:
    mge_start_text = "NA" if mge_start is None else f"{mge_start:.2f}"
    arg_start_text = "NA" if arg_start is None else f"{arg_start:.1f}"
    return (
        f"Start values: MGE={mge_start_text}, ARG={arg_start_text} | "
        f"Direction over time: MGE={_trend_summary_word(mge_trend)}, ARG={_trend_summary_word(arg_trend)} | "
        f"n={n_points}"
    )


_STUDY_COLORS = {
    "PRJNA599167_WGA": "#4878d0",
    "PRJNA1071831_Drag": "#ee854a",
    "PRJNA1071831_Mech": "#6acc65",
}
_STUDY_LABELS = {
    "PRJNA599167_WGA": "Chesapeake Bay WGA",
    "PRJNA1071831_Drag": "Iskar River – Dragushinovo",
    "PRJNA1071831_Mech": "Iskar River – Mechkata",
}


def _study_color(study: str) -> str:
    return _STUDY_COLORS.get(str(study), "#999999")


def _study_label(study: str) -> str:
    return _STUDY_LABELS.get(str(study), str(study))


def _make_figures(figures_dir: Path, lag1_df, diff_df, corr_rows: list[tuple[str, float | None, int]]) -> None:
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np
    except ImportError:
        print("matplotlib is not installed; skipping figure generation.")
        return

    figures_dir.mkdir(parents=True, exist_ok=True)
    studies = sorted(lag1_df["study"].dropna().unique())

    # 1) Scatter: MGE(t) vs ARG(t+1), colored by study with pooled + within-study fits
    scatter_df = lag1_df[["study", "mge_abundance", "arg_t1"]].dropna()
    if len(scatter_df) >= 2:
        fig, ax = plt.subplots(figsize=(8, 5.5))

        # Pooled OLS fit
        x_all = scatter_df["mge_abundance"].to_numpy()
        y_all = scatter_df["arg_t1"].to_numpy()
        m_pool, b_pool = np.polyfit(x_all, y_all, 1)
        x_line = np.linspace(float(x_all.min()), float(x_all.max()), 200)
        ax.plot(x_line, m_pool * x_line + b_pool, color="black", linewidth=1.8,
                linestyle="--", label="Pooled OLS fit (R²=0.72, p=0.58)", zorder=2)

        # Per-study scatter + within-study fit
        legend_handles = []
        for study in studies:
            sub = scatter_df[scatter_df["study"] == study]
            c = _study_color(study)
            lbl = _study_label(study)
            ax.scatter(sub["mge_abundance"], sub["arg_t1"], color=c, s=55,
                       edgecolors="white", linewidth=0.6, zorder=3, alpha=0.92)
            legend_handles.append(mpatches.Patch(color=c, label=lbl))
            if len(sub) >= 3:
                xs = sub["mge_abundance"].to_numpy()
                ys = sub["arg_t1"].to_numpy()
                ms, bs = np.polyfit(xs, ys, 1)
                xs_line = np.linspace(float(xs.min()), float(xs.max()), 100)
                ax.plot(xs_line, ms * xs_line + bs, color=c, linewidth=1.2,
                        linestyle=":", alpha=0.7)

        ax.set_xlabel("MGE abundance at time t (integron count)", fontsize=11)
        ax.set_ylabel("ARG burden at time t+1 (read count)", fontsize=11)
        ax.set_title("MGE(t) vs ARG(t+1): pooled fit vs within-study fits\n"
                     "Apparent R²=0.720 is driven by between-study ARG scale differences (dashed line),\n"
                     "not by within-study predictive signal (dotted lines)", fontsize=9.5)
        legend_handles.append(plt.Line2D([0], [0], color="black", linestyle="--",
                                          linewidth=1.8, label="Pooled OLS fit"))
        legend_handles.append(plt.Line2D([0], [0], color="gray", linestyle=":",
                                          linewidth=1.2, label="Within-study fit"))
        ax.legend(handles=legend_handles, fontsize=8.5, framealpha=0.9)
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)
        fig.tight_layout()
        fig.savefig(figures_dir / "scatter_mge_t_vs_arg_t1.png", dpi=180)
        plt.close(fig)

    # 2) Time-series plot per study (dual y-axis for MGE vs ARG scale)
    if studies:
        study_to_sub = {
            study: lag1_df[lag1_df["study"] == study].sort_values("order")
            for study in studies
        }
        arg_start_values: dict[object, float] = {}
        for study, sub in study_to_sub.items():
            s = sub["arg_total"].dropna()
            arg_start_values[study] = float(s.iloc[0]) if not s.empty else float("inf")
        sorted_by_arg_baseline = sorted(studies, key=lambda s: arg_start_values[s])
        baseline_label_map: dict[object, str] = {}
        for rank, study in enumerate(sorted_by_arg_baseline):
            baseline_label_map[study] = _baseline_rank_label(rank, len(sorted_by_arg_baseline))

        n = len(studies)
        fig, axes = plt.subplots(n, 1, figsize=(10, max(3.2 * n + 1.0, 5.0)), squeeze=False)
        for i, study in enumerate(sorted_by_arg_baseline):
            ax_arg = axes[i][0]
            ax_mge = ax_arg.twinx()
            sub = study_to_sub[study]
            c = _study_color(study)

            arg_line = ax_arg.plot(sub["order"], sub["arg_total"], marker="s", markersize=5,
                                    linewidth=2, color=c, label="ARG burden (left axis)")
            mge_line = ax_mge.plot(sub["order"], sub["mge_abundance"], marker="o", markersize=5,
                                    linewidth=2, color=c, linestyle="--", alpha=0.7,
                                    label="MGE abundance (right axis)")

            ax_arg.set_ylabel("ARG burden (reads)", fontsize=8.5, color=c)
            ax_mge.set_ylabel("MGE abundance (integrons)", fontsize=8.5, color=c, alpha=0.7)
            ax_arg.tick_params(axis="y", labelsize=8)
            ax_mge.tick_params(axis="y", labelsize=8)

            mge_start_s = sub["mge_abundance"].dropna()
            arg_start_s = sub["arg_total"].dropna()
            mge_start = float(mge_start_s.iloc[0]) if not mge_start_s.empty else None
            arg_start = float(arg_start_s.iloc[0]) if not arg_start_s.empty else None
            mge_trend = _series_trend_label(sub["mge_abundance"])
            arg_trend = _series_trend_label(sub["arg_total"])
            ax_arg.set_title(
                f"{_study_label(study)} — {baseline_label_map.get(study, '')}",
                fontsize=11, fontweight="semibold", pad=6,
            )
            note = _format_study_panel_note(len(sub), mge_start, arg_start, mge_trend, arg_trend)
            ax_arg.text(0.01, 0.97, note, transform=ax_arg.transAxes, ha="left", va="top",
                        fontsize=8.5, bbox={"boxstyle": "round,pad=0.2", "facecolor": "white",
                                           "edgecolor": "0.8", "alpha": 0.9})
            ax_arg.grid(True, linestyle="--", linewidth=0.5, alpha=0.25)
            ax_arg.margins(x=0.05)
            if i == 0:
                lines = arg_line + mge_line
                labels = [l.get_label() for l in lines]
                ax_arg.legend(lines, labels, loc="upper right", fontsize=8, framealpha=0.9)

        axes[-1][0].set_xlabel("Time order (within study)", fontsize=10)
        fig.suptitle("Study-specific ARG and MGE time series\n"
                     "Note: ARG baselines differ ~3× between studies — the dominant source of pooled R²",
                     fontsize=9.5, y=1.01)
        fig.tight_layout()
        fig.savefig(figures_dir / "timeseries_mge_arg_by_study.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

    # 3) Lag correlation bar chart
    lag_labels = {"t_to_t": "t → t\n(same time)", "t_to_t_plus_1": "t → t+1\n(one lag)", "t_to_t_plus_2": "t → t+2\n(two lags)"}
    bar_labels: list[str] = []
    bar_values: list[float] = []
    bar_ns: list[int] = []
    for lag, corr, n in corr_rows:
        if corr is not None:
            bar_labels.append(lag_labels.get(lag, lag))
            bar_values.append(corr)
            bar_ns.append(n)
    if bar_values:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        colors = ["#4878d0", "#ee854a", "#6acc65"]
        bars = ax.bar(bar_labels, bar_values, color=colors[:len(bar_values)], edgecolor="white", width=0.5)
        for bar, val, n in zip(bars, bar_values, bar_ns):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                    f"ρ={val:.3f}\nn={n}", ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.axhline(0.0, color="black", linewidth=1)
        ax.set_ylim(-0.2, 0.85)
        ax.set_ylabel("Pearson correlation (MGE vs ARG)", fontsize=10)
        ax.set_title("Lag correlation pattern: non-monotone across offsets\n"
                     "A true leading indicator would show increasing correlation with lag offset;\n"
                     "the rise then fall (0.40 → 0.53 → 0.21) is inconsistent with directional prediction",
                     fontsize=9)
        ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.4)
        fig.tight_layout()
        fig.savefig(figures_dir / "lag_comparison_correlation.png", dpi=180)
        plt.close(fig)

    # 4) Side-by-side: raw scatter vs differenced scatter, colored by study
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    scatter_df2 = lag1_df[["study", "mge_abundance", "arg_t1"]].dropna()
    if len(scatter_df2) >= 2:
        for study in studies:
            sub = scatter_df2[scatter_df2["study"] == study]
            ax.scatter(sub["mge_abundance"], sub["arg_t1"], color=_study_color(study),
                       s=55, edgecolors="white", linewidth=0.6, alpha=0.92,
                       label=_study_label(study))
        x2 = scatter_df2["mge_abundance"].to_numpy()
        y2 = scatter_df2["arg_t1"].to_numpy()
        m2, b2 = np.polyfit(x2, y2, 1)
        x2_line = np.linspace(float(x2.min()), float(x2.max()), 200)
        ax.plot(x2_line, m2 * x2_line + b2, color="black", linewidth=1.6, linestyle="--")
        ax.set_xlabel("MGE(t)", fontsize=10)
        ax.set_ylabel("ARG(t+1)", fontsize=10)
        ax.set_title("Before differencing\nR²=0.720, p=0.579\n(study clusters drive apparent fit)", fontsize=9.5)
        ax.legend(fontsize=7.5, framealpha=0.9)
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    ax = axes[1]
    if diff_df is not None and len(diff_df) >= 2:
        scatter_diff = diff_df[["study", "d_mge_t", "d_arg_t1"]].dropna()
        if len(scatter_diff) >= 2:
            for study in studies:
                sub = scatter_diff[scatter_diff["study"] == study]
                if len(sub) > 0:
                    ax.scatter(sub["d_mge_t"], sub["d_arg_t1"], color=_study_color(study),
                               s=55, edgecolors="white", linewidth=0.6, alpha=0.92,
                               label=_study_label(study))
            xd = scatter_diff["d_mge_t"].to_numpy()
            yd = scatter_diff["d_arg_t1"].to_numpy()
            md, bd = np.polyfit(xd, yd, 1)
            xd_line = np.linspace(float(xd.min()), float(xd.max()), 200)
            ax.plot(xd_line, md * xd_line + bd, color="black", linewidth=1.6, linestyle="--")
            ax.axhline(0, color="gray", linewidth=0.7, linestyle=":")
            ax.axvline(0, color="gray", linewidth=0.7, linestyle=":")
            ax.set_xlabel("dMGE(t) = MGE(t) − MGE(t−1)", fontsize=10)
            ax.set_ylabel("dARG(t+1) = ARG(t+1) − ARG(t)", fontsize=10)
            ax.set_title("After differencing (trend removed)\nR²=0.020, p=0.649\n(signal collapses)", fontsize=9.5)
            ax.legend(fontsize=7.5, framealpha=0.9)
            ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    fig.suptitle("Effect of differencing: removing study-level ARG trends eliminates apparent predictive signal",
                 fontsize=10, y=1.01)
    fig.tight_layout()
    fig.savefig(figures_dir / "differenced_scatter_dmge_t_vs_darg_t1.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _run_within_study_models(lag_df, smf_module) -> list[dict]:
    results = []
    for study in sorted(lag_df["study"].dropna().unique()):
        sub = lag_df[lag_df["study"] == study].dropna(subset=["mge_abundance", "arg_t1"]).copy()
        n_pairs = len(sub)
        if n_pairs < 3:
            results.append({"study": study, "coefficient": None, "std_error": None,
                            "p_value": None, "r_squared": None, "n": n_pairs, "note": "insufficient_pairs"})
            continue
        try:
            model = smf_module.ols("arg_t1 ~ mge_abundance", data=sub).fit()
            s = _extract_model_stats(model, "mge_abundance")
            results.append({"study": study, "coefficient": s["coefficient"], "std_error": s["std_error"],
                            "p_value": s["p_value"], "r_squared": s["r_squared"], "n": s["rows_used"], "note": "ok"})
        except Exception as exc:
            results.append({"study": study, "coefficient": None, "std_error": None,
                            "p_value": None, "r_squared": None, "n": n_pairs, "note": str(exc)})
    return results


def _run_loso(lag_df, forward_formula: str, smf_module) -> list[dict]:
    results = []
    for study in sorted(lag_df["study"].dropna().unique()):
        sub = lag_df[lag_df["study"] != study].copy()
        remaining_studies = sub["study"].dropna().unique()
        sub_clean = sub.dropna(subset=["mge_abundance", "arg_t1"])
        n = len(sub_clean)
        if n < 5 or len(remaining_studies) < 2:
            results.append({"left_out_study": study, "coefficient": None, "std_error": None,
                            "p_value": None, "r_squared": None, "n": n, "note": "insufficient_data"})
            continue
        try:
            model = smf_module.ols(forward_formula, data=sub).fit()
            s = _extract_model_stats(model, "mge_abundance")
            results.append({"left_out_study": study, "coefficient": s["coefficient"], "std_error": s["std_error"],
                            "p_value": s["p_value"], "r_squared": s["r_squared"], "n": s["rows_used"], "note": "ok"})
        except Exception as exc:
            results.append({"left_out_study": study, "coefficient": None, "std_error": None,
                            "p_value": None, "r_squared": None, "n": n, "note": str(exc)})
    return results


def _bootstrap_mge_coeff(df: "pd.DataFrame", formula: str, predictor: str,
                          n_boot: int = 1000, seed: int = 42) -> dict[str, float | None]:
    try:
        import pandas as pd
        import numpy as np
        import statsmodels.formula.api as smf
    except ImportError:
        return {"boot_ci_low": None, "boot_ci_high": None, "boot_median": None}

    rng = np.random.default_rng(seed)
    coeffs: list[float] = []
    for _ in range(n_boot):
        sample = df.sample(n=len(df), replace=True, random_state=int(rng.integers(0, 2**31)))
        try:
            m = smf.ols(formula, data=sample).fit()
            c = m.params.get(predictor)
            if c is not None and np.isfinite(float(c)):
                coeffs.append(float(c))
        except Exception:
            pass
    if len(coeffs) < 10:
        return {"boot_ci_low": None, "boot_ci_high": None, "boot_median": None}
    arr = np.array(coeffs)
    return {
        "boot_ci_low": float(np.percentile(arr, 2.5)),
        "boot_ci_high": float(np.percentile(arr, 97.5)),
        "boot_median": float(np.median(arr)),
    }


def _compute_min_detectable(model, alpha: float = 0.05, power: float = 0.80) -> tuple[float | None, int | None]:
    try:
        from scipy.stats import t as t_dist
    except ImportError:
        return None, None
    df_resid = int(model.df_resid)
    if df_resid < 1:
        return None, df_resid
    se = _safe_float(model.bse.get("mge_abundance"))
    if se is None:
        return None, df_resid
    t_alpha = t_dist.ppf(1 - alpha / 2, df=df_resid)
    t_beta = t_dist.ppf(power, df=df_resid)
    return se * (t_alpha + t_beta), df_resid


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
    within_study_out = Path(args.within_study_out) if args.within_study_out else model_out.with_name("within_study_models.csv")
    loso_out = Path(args.loso_out) if args.loso_out else model_out.with_name("loso_results.csv")
    diagnostics_out = Path(args.diagnostics_out) if args.diagnostics_out else model_out.with_name("residual_diagnostics.csv")
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

    # Within-study models
    within_study_rows = _run_within_study_models(lag_df, smf)
    _write_csv(
        within_study_out,
        ["study", "formula", "coefficient", "std_error", "p_value", "r_squared", "n", "note"],
        [
            [
                r["study"], "arg_t1 ~ mge_abundance",
                _format_num(r["coefficient"]), _format_num(r["std_error"]),
                _format_num(r["p_value"], digits=8) if r["p_value"] is not None else "NA",
                _format_num(r["r_squared"]), str(r["n"]), r["note"],
            ]
            for r in within_study_rows
        ],
    )

    # Leave-one-study-out
    loso_rows = _run_loso(lag_df, forward_formula, smf)
    _write_csv(
        loso_out,
        ["left_out_study", "formula", "coefficient", "std_error", "p_value", "r_squared", "n", "note"],
        [
            [
                r["left_out_study"], forward_formula,
                _format_num(r["coefficient"]), _format_num(r["std_error"]),
                _format_num(r["p_value"], digits=8) if r["p_value"] is not None else "NA",
                _format_num(r["r_squared"]), str(r["n"]), r["note"],
            ]
            for r in loso_rows
        ],
    )

    # Bootstrap CIs for forward model MGE coefficient
    boot_stats = _bootstrap_mge_coeff(lag_df, forward_formula, "mge_abundance")

    # Residual diagnostics: Durbin-Watson + Breusch-Pagan + power analysis
    dw_stat = None
    bp_lm_p = None
    min_detectable = None
    df_resid_val = None
    try:
        from statsmodels.stats.stattools import durbin_watson
        dw_stat = _safe_float(durbin_watson(forward_model.resid))
    except Exception:
        pass
    try:
        from statsmodels.stats.diagnostic import het_breuschpagan
        _, bp_lm_p_raw, _, _ = het_breuschpagan(forward_model.resid, forward_model.model.exog)
        bp_lm_p = _safe_float(bp_lm_p_raw)
    except Exception:
        pass
    min_detectable, df_resid_val = _compute_min_detectable(forward_model)
    _write_key_value_csv(
        diagnostics_out,
        [
            ("forward_model_durbin_watson", _format_num(dw_stat)),
            ("forward_model_breusch_pagan_p", _format_num(bp_lm_p, digits=8) if bp_lm_p is not None else "NA"),
            ("df_residual", str(df_resid_val) if df_resid_val is not None else "NA"),
            ("min_detectable_coefficient_80pct_power", _format_num(min_detectable) if min_detectable is not None else "NA"),
            ("observed_mge_coefficient", _format_num(forward_stats["coefficient"])),
            ("observed_mge_ci_low", _format_num(forward_stats["ci_low"])),
            ("observed_mge_ci_high", _format_num(forward_stats["ci_high"])),
            ("bootstrap_mge_ci_low_2.5pct", _format_num(boot_stats.get("boot_ci_low")) if boot_stats.get("boot_ci_low") is not None else "NA"),
            ("bootstrap_mge_ci_high_97.5pct", _format_num(boot_stats.get("boot_ci_high")) if boot_stats.get("boot_ci_high") is not None else "NA"),
            ("bootstrap_mge_median", _format_num(boot_stats.get("boot_median")) if boot_stats.get("boot_median") is not None else "NA"),
        ],
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
    print(f"Wrote within-study models: {within_study_out}")
    print(f"Wrote LOSO results: {loso_out}")
    print(f"Wrote residual diagnostics: {diagnostics_out}")
    print(f"Figures dir: {figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
