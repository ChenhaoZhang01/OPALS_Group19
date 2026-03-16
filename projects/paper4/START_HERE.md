# Paper 4 Intern Start Guide (Detailed)

Read this file top to bottom. Copy and paste commands exactly.

## Project goal

Estimate whether wastewater treatment upgrades reduce ARG levels using a difference-in-differences (DiD) model.

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

Use this file:

- `projects/paper4/results/did_input_table.csv`

## Step 1: Fill DiD input table

Required columns:

- `sample_id`
- `ARG_total`
- `treatment`
- `time`
- `plant_id`

Coding rules:

- `treatment = 1` for plants that received the upgrade
- `treatment = 0` for control plants
- `time = 0` for before-upgrade period
- `time = 1` for after-upgrade period
- `ARG_total` must be numeric

Add multiple rows for each plant across time.

## Step 2: Run DiD model

```powershell
python projects/paper4/analysis/run_did_model.py --input projects/paper4/results/did_input_table.csv --coef-out projects/paper4/results/did_coefficients.csv --summary-out projects/paper4/results/did_model_summary.txt
```

Success check:

- `projects/paper4/results/did_coefficients.csv` exists.
- `projects/paper4/results/did_model_summary.txt` exists.

Main coefficient to report:

- `treatment:time`

## Step 3: Make figures

Save under `projects/paper4/analysis/figures/`:

1. ARG before vs after by group.
2. Treated vs control trend lines.
3. DiD coefficient plot with confidence interval.
4. Country-level comparison.

## Files to submit

- `projects/paper4/results/did_input_table.csv`
- `projects/paper4/results/did_coefficients.csv`
- `projects/paper4/results/did_model_summary.txt`

## Common mistakes

- Swapping meaning of `treatment` and `time`.
- Using text values instead of 0/1 for `treatment` and `time`.
- Including missing `ARG_total` values.
