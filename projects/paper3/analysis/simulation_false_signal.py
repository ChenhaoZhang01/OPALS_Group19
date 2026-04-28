#!/usr/bin/env python3
"""
Simulation demonstrating how between-study ARG scale differences inflate
lag model R² even when MGE has zero true predictive signal.

Calibrated to match the real data structure (3 shotgun studies, same n).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


# Study parameters calibrated to real data distributions
# (n_pairs, arg_mean, arg_sd, mge_scale)
STUDY_PARAMS = {
    "WGA":  (10, 15_000,  6_000,  40),
    "Drag": ( 3, 47_000, 10_000, 200),
    "Mech": ( 3, 47_000, 13_000,  70),
}

N_SIMS = 2_000
SEED = 42


def _simulate_once(rng: np.random.Generator) -> pd.DataFrame:
    rows: list[dict] = []
    for study, (n_pairs, arg_mean, arg_sd, mge_scale) in STUDY_PARAMS.items():
        n_samples = n_pairs + 1
        # ARG: random fluctuations around study baseline — no MGE dependence
        arg = np.clip(rng.normal(arg_mean, arg_sd, n_samples), 0, None)
        mge = np.clip(rng.exponential(mge_scale, n_samples), 0, None)

        for i in range(n_pairs):
            rows.append({
                "study": study,
                "mge_t": mge[i],
                "mge_tminus1": mge[i - 1] if i > 0 else np.nan,
                "arg_t": arg[i],
                "arg_t1": arg[i + 1],
            })

    df = pd.DataFrame(rows)
    df["d_mge_t"] = df["mge_t"] - df["mge_tminus1"]
    df["d_arg_t1"] = df["arg_t1"] - df["arg_t"]
    return df


def main() -> int:
    out_dir = Path(__file__).resolve().parents[1] / "results"
    fig_dir = Path(__file__).resolve().parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)

    naive_r2_list: list[float] = []
    diff_r2_list: list[float] = []
    naive_mge_p_list: list[float] = []

    for _ in range(N_SIMS):
        df = _simulate_once(rng)

        try:
            m = smf.ols("arg_t1 ~ mge_t + C(study)", data=df).fit()
            naive_r2_list.append(float(m.rsquared))
            naive_mge_p_list.append(float(m.pvalues.get("mge_t", np.nan)))
        except Exception:
            naive_r2_list.append(np.nan)
            naive_mge_p_list.append(np.nan)

        try:
            diff_df = df.dropna(subset=["d_mge_t", "d_arg_t1"])
            if len(diff_df) >= 3:
                m_d = smf.ols("d_arg_t1 ~ d_mge_t", data=diff_df).fit()
                diff_r2_list.append(float(m_d.rsquared))
            else:
                diff_r2_list.append(np.nan)
        except Exception:
            diff_r2_list.append(np.nan)

    naive_r2 = np.array(naive_r2_list)
    diff_r2 = np.array(diff_r2_list)
    naive_p = np.array(naive_mge_p_list)

    # Summary statistics
    def _pct_above(arr: np.ndarray, threshold: float) -> float:
        return float(np.mean(arr[~np.isnan(arr)] > threshold) * 100)

    summary = [
        ("n_simulations", str(N_SIMS)),
        ("naive_r2_median", f"{np.nanmedian(naive_r2):.4f}"),
        ("naive_r2_mean", f"{np.nanmean(naive_r2):.4f}"),
        ("naive_r2_p25", f"{np.nanpercentile(naive_r2, 25):.4f}"),
        ("naive_r2_p75", f"{np.nanpercentile(naive_r2, 75):.4f}"),
        ("naive_r2_pct_above_0.5", f"{_pct_above(naive_r2, 0.5):.1f}"),
        ("naive_r2_pct_above_0.7", f"{_pct_above(naive_r2, 0.7):.1f}"),
        ("naive_mge_p_pct_sig_0.05", f"{float(np.mean(naive_p[~np.isnan(naive_p)] < 0.05) * 100):.1f}"),
        ("diff_r2_median", f"{np.nanmedian(diff_r2):.4f}"),
        ("diff_r2_mean", f"{np.nanmean(diff_r2):.4f}"),
        ("diff_r2_p75", f"{np.nanpercentile(diff_r2, 75):.4f}"),
        ("observed_naive_r2", "0.720"),
        ("observed_diff_r2", "0.020"),
    ]

    sim_out = out_dir / "simulation_results.csv"
    with sim_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        w.writerows(summary)

    print("=== Simulation results ===")
    for k, v in summary:
        print(f"  {k}: {v}")

    # Figure
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

        ax = axes[0]
        ax.hist(naive_r2[~np.isnan(naive_r2)], bins=50, color="#4878d0", edgecolor="white", linewidth=0.4)
        ax.axvline(np.nanmedian(naive_r2), color="#d62728", linestyle="--", linewidth=1.8,
                   label=f"Null median = {np.nanmedian(naive_r2):.2f}")
        ax.axvline(0.720, color="#ff7f0e", linestyle=":", linewidth=2.0,
                   label="Observed R² = 0.720")
        ax.set_xlabel("R² — naive model: ARG(t+1) ~ MGE(t) + C(study)", fontsize=10)
        ax.set_ylabel(f"Count (out of {N_SIMS} simulations)", fontsize=10)
        ax.set_title("Naive lag model R² distribution under null\n(zero true MGE→ARG signal)", fontsize=10)
        ax.legend(fontsize=9)

        ax = axes[1]
        ax.hist(diff_r2[~np.isnan(diff_r2)], bins=50, color="#6acc65", edgecolor="white", linewidth=0.4)
        ax.axvline(np.nanmedian(diff_r2), color="#d62728", linestyle="--", linewidth=1.8,
                   label=f"Null median = {np.nanmedian(diff_r2):.2f}")
        ax.axvline(0.020, color="#ff7f0e", linestyle=":", linewidth=2.0,
                   label="Observed R² = 0.020")
        ax.set_xlabel("R² — differenced model: dARG(t+1) ~ dMGE(t)", fontsize=10)
        ax.set_ylabel(f"Count (out of {N_SIMS} simulations)", fontsize=10)
        ax.set_title("After differencing: R² distribution under null\n(shared trends removed)", fontsize=10)
        ax.legend(fontsize=9)

        fig.suptitle(
            "Simulation: study-level ARG scale differences inflate naive lag model R²\n"
            "without any true predictive signal (N=2,000 synthetic datasets)",
            fontsize=10,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        fig_path = fig_dir / "simulation_false_signal.png"
        fig.savefig(fig_path, dpi=180)
        plt.close(fig)
        print(f"Saved figure: {fig_path}")
    except ImportError:
        print("matplotlib not available; skipping figure.")

    print(f"Wrote: {sim_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
