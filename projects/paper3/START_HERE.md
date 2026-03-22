# Paper 3 Intern Start Guide (Detailed)

Read this file top to bottom. Copy and paste commands exactly.

## Project goal

Test whether MGE abundance contains predictive information for future ARG change after accounting for temporal coupling.

Core warning to test: near-perfect lag-model fit can arise from temporal coupling/autocorrelation alone.

## Where to run commands

Always run commands from the repository root folder (`OPALS_Group19`) on your own computer.

Activate your local virtual environment first.

Windows PowerShell:

```powershell
./.venv/Scripts/Activate.ps1
```

Linux/WSL/macOS:

```bash
source .venv/bin/activate
```

## Before you start

1. Use this starter file:
	- `projects/paper3/results/features_template.csv`
2. Make a copy named:
	- `projects/paper3/results/features_table.csv`

```powershell
Copy-Item projects/paper3/results/features_template.csv projects/paper3/results/features_table.csv -Force
```

## Step 1: Fill the feature table

Open `projects/paper3/results/features_table.csv` and fill all rows.

Required columns:

- `study`
- `sample_id`
- `order` (time order, 1 then 2 then 3)
- `mge_abundance`
- `entropy`
- `arg_total`

Important:

- `order`, `mge_abundance`, `entropy`, and `arg_total` must be numeric.
- Keep each study sorted by sampling order.

## Step 2: Run lag model

```powershell
python projects/paper3/analysis/run_lag_analysis.py --input projects/paper3/results/features_table.csv --model-out projects/paper3/results/lag_regression_results.csv --corr-out projects/paper3/results/cross_correlation.csv
```

Success check:

- `projects/paper3/results/lag_regression_results.csv` exists.
- `projects/paper3/results/cross_correlation.csv` exists.
- `projects/paper3/results/difference_model.csv` exists.
- `projects/paper3/results/granger_test.csv` exists.

Check that the model outputs include `adjusted_r_squared` and `sequencing_depth_transform`.

## Step 3: Make figures

Save under `projects/paper3/analysis/figures/`:

1. ARG time series.
2. MGE time series.
3. Lag correlation plot.
4. Differenced scatter: `dMGE(t)` vs `dARG(t+1)`.

## Files to submit

- `projects/paper3/results/features_table.csv`
- `projects/paper3/results/lag_regression_results.csv`
- `projects/paper3/results/cross_correlation.csv`
- `projects/paper3/results/difference_model.csv`
- `projects/paper3/results/granger_test.csv`

## Common mistakes

- Missing or non-numeric values in `order`, `mge_abundance`, `arg_total`.
- Mixing sample order across studies.
- Running model before filling enough rows (needs at least two time points per study).
- Interpreting very high R-squared as directional prediction without checking differenced and Granger-style outputs.
- Ignoring uniformly high cross-lag correlation as evidence of temporal coupling rather than directional prediction.
