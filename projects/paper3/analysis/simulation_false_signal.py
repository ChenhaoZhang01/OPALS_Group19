#!/usr/bin/env python3
"""
Two-scenario simulation for diagnostic framework calibration:

  Scenario A (null):   no true MGE→ARG predictive signal
  Scenario B (signal): MGE truly causes ARG(t+1) with beta=150

Both scenarios use the same study structure as the real data
(3 studies, same n per study, ARG baselines calibrated to real data).
The comparison shows the framework correctly flags Scenario A while
detecting Scenario B — i.e., it has power when signal exists.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


STUDY_PARAMS = {
    "WGA":  (10, 15_000,  6_000,  40),   # (n_pairs, arg_mean, arg_sd, mge_scale)
    "Drag": ( 3, 47_000, 10_000, 200),
    "Mech": ( 3, 47_000, 13_000,  70),
}

N_SIMS = 2_000
SEED   = 42
TRUE_BETA = 150   # MGE→ARG coefficient used in Scenario B


def _simulate_once(rng: np.random.Generator, beta: float) -> pd.DataFrame:
    rows: list[dict] = []
    for study, (n_pairs, arg_mean, arg_sd, mge_scale) in STUDY_PARAMS.items():
        n_samples = n_pairs + 1
        mge = np.clip(rng.exponential(mge_scale, n_samples), 0, None)
        # ARG(t+1) = baseline + beta*MGE(t) + noise
        noise = rng.normal(0, arg_sd, n_samples)
        arg = np.zeros(n_samples)
        arg[0] = arg_mean + noise[0]
        for i in range(1, n_samples):
            arg[i] = arg_mean + beta * mge[i - 1] + noise[i]
        arg = np.clip(arg, 0, None)

        for i in range(n_pairs):
            rows.append({
                "study":       study,
                "mge_t":       mge[i],
                "mge_tminus1": mge[i - 1] if i > 0 else np.nan,
                "arg_t":       arg[i],
                "arg_t1":      arg[i + 1],
            })

    df = pd.DataFrame(rows)
    df["d_mge_t"] = df["mge_t"] - df["mge_tminus1"]
    df["d_arg_t1"] = df["arg_t1"] - df["arg_t"]
    return df


def _run_scenario(rng: np.random.Generator, beta: float, n_sims: int) -> dict[str, np.ndarray]:
    naive_r2, diff_r2, naive_p, within_wga_coeff = [], [], [], []

    for _ in range(n_sims):
        df = _simulate_once(rng, beta)

        try:
            m = smf.ols("arg_t1 ~ mge_t + C(study)", data=df).fit()
            naive_r2.append(float(m.rsquared))
            naive_p.append(float(m.pvalues.get("mge_t", np.nan)))
        except Exception:
            naive_r2.append(np.nan); naive_p.append(np.nan)

        try:
            ddf = df.dropna(subset=["d_mge_t", "d_arg_t1"])
            m_d = smf.ols("d_arg_t1 ~ d_mge_t", data=ddf).fit() if len(ddf) >= 3 else None
            diff_r2.append(float(m_d.rsquared) if m_d else np.nan)
        except Exception:
            diff_r2.append(np.nan)

        try:
            wga = df[df["study"] == "WGA"].dropna(subset=["mge_t", "arg_t1"])
            m_w = smf.ols("arg_t1 ~ mge_t", data=wga).fit() if len(wga) >= 3 else None
            within_wga_coeff.append(float(m_w.params.get("mge_t", np.nan)) if m_w else np.nan)
        except Exception:
            within_wga_coeff.append(np.nan)

    return {
        "naive_r2":         np.array(naive_r2),
        "diff_r2":          np.array(diff_r2),
        "naive_p":          np.array(naive_p),
        "within_wga_coeff": np.array(within_wga_coeff),
    }


def main() -> int:
    out_dir = Path(__file__).resolve().parents[1] / "results"
    fig_dir = Path(__file__).resolve().parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    rng_a = np.random.default_rng(SEED)
    rng_b = np.random.default_rng(SEED + 1)

    print(f"Running Scenario A (null, beta=0) — {N_SIMS} simulations …")
    res_a = _run_scenario(rng_a, beta=0.0, n_sims=N_SIMS)

    print(f"Running Scenario B (true signal, beta={TRUE_BETA}) — {N_SIMS} simulations …")
    res_b = _run_scenario(rng_b, beta=float(TRUE_BETA), n_sims=N_SIMS)

    def _pct_above(arr, thr):
        a = arr[~np.isnan(arr)]
        return float(np.mean(a > thr) * 100)

    def _pct_sig(p_arr):
        a = p_arr[~np.isnan(p_arr)]
        return float(np.mean(a < 0.05) * 100)

    summary = [
        # Scenario A
        ("A_null_naive_r2_median",          f"{np.nanmedian(res_a['naive_r2']):.4f}"),
        ("A_null_naive_r2_pct_above_0.70",  f"{_pct_above(res_a['naive_r2'], 0.70):.1f}"),
        ("A_null_naive_mge_pct_sig",        f"{_pct_sig(res_a['naive_p']):.1f}"),
        ("A_null_diff_r2_median",           f"{np.nanmedian(res_a['diff_r2']):.4f}"),
        ("A_null_within_wga_coeff_median",  f"{np.nanmedian(res_a['within_wga_coeff']):.1f}"),
        ("A_null_within_wga_pct_positive",  f"{_pct_above(res_a['within_wga_coeff'], 0):.1f}"),
        # Scenario B
        ("B_signal_naive_r2_median",             f"{np.nanmedian(res_b['naive_r2']):.4f}"),
        ("B_signal_naive_r2_pct_above_0.70",     f"{_pct_above(res_b['naive_r2'], 0.70):.1f}"),
        ("B_signal_naive_mge_pct_sig",           f"{_pct_sig(res_b['naive_p']):.1f}"),
        ("B_signal_diff_r2_median",              f"{np.nanmedian(res_b['diff_r2']):.4f}"),
        ("B_signal_within_wga_coeff_median",     f"{np.nanmedian(res_b['within_wga_coeff']):.1f}"),
        ("B_signal_within_wga_pct_positive",     f"{_pct_above(res_b['within_wga_coeff'], 0):.1f}"),
        # Observed real data
        ("observed_naive_r2",    "0.720"),
        ("observed_diff_r2",     "0.020"),
        ("observed_wga_coeff",   "-33.0"),
        ("true_beta_scenario_b", str(TRUE_BETA)),
    ]

    sim_out = out_dir / "simulation_results.csv"
    with sim_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        w.writerows(summary)

    print("\n=== Summary ===")
    for k, v in summary:
        print(f"  {k}: {v}")

    # Figure: 2x2 — rows = naive / differenced, columns = null / true signal
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey="row")

        BLUE   = "#4878d0"
        GREEN  = "#6acc65"
        RED    = "#d62728"
        ORANGE = "#ff7f0e"

        def _annotate_observed(ax, val, label, color=ORANGE):
            ax.axvline(val, color=color, linestyle=":", linewidth=2.0, label=label)

        # Row 0: naive model R²
        for col, (res, title, obs_label) in enumerate([
            (res_a, "Scenario A — null (no true MGE signal)", None),
            (res_b, f"Scenario B — true signal (β={TRUE_BETA})", None),
        ]):
            ax = axes[0][col]
            ax.hist(res["naive_r2"][~np.isnan(res["naive_r2"])], bins=50,
                    color=BLUE, edgecolor="white", linewidth=0.3)
            med = np.nanmedian(res["naive_r2"])
            ax.axvline(med, color=RED, linestyle="--", linewidth=1.8,
                       label=f"Median = {med:.2f}")
            if col == 0:
                _annotate_observed(ax, 0.720, "Observed R² = 0.720")
            pct = _pct_above(res["naive_r2"], 0.70)
            ax.set_title(f"{title}", fontsize=9.5)
            ax.set_xlabel("R² — naive model: ARG(t+1) ~ MGE(t) + C(study)", fontsize=9)
            ax.set_ylabel("Count" if col == 0 else "", fontsize=9)
            ax.legend(fontsize=8.5)
            ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.4)

        # Row 1: differenced model R²
        for col, (res, title) in enumerate([
            (res_a, "Scenario A — null"),
            (res_b, f"Scenario B — true signal (β={TRUE_BETA})"),
        ]):
            ax = axes[1][col]
            ax.hist(res["diff_r2"][~np.isnan(res["diff_r2"])], bins=50,
                    color=GREEN, edgecolor="white", linewidth=0.3)
            med = np.nanmedian(res["diff_r2"])
            ax.axvline(med, color=RED, linestyle="--", linewidth=1.8,
                       label=f"Median = {med:.2f}")
            if col == 0:
                _annotate_observed(ax, 0.020, "Observed R² = 0.020")
            ax.set_title(f"{title}", fontsize=9.5)
            ax.set_xlabel("R² — differenced: dARG(t+1) ~ dMGE(t)", fontsize=9)
            ax.set_ylabel("Count" if col == 0 else "", fontsize=9)
            ax.legend(fontsize=8.5)
            ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.4)

        fig.tight_layout()
        fig_path = fig_dir / "simulation_false_signal.png"
        fig.savefig(fig_path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"\nSaved figure: {fig_path}")
    except ImportError:
        print("matplotlib not available; skipping figure.")

    print(f"Wrote: {sim_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
