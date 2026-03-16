# Paper 4: Natural Experiment on Wastewater Treatment

For interns with limited coding experience: start with `START_HERE.md` in this folder.

## Objective

Estimate whether wastewater treatment upgrades reduce ARG abundance using causal inference.

Primary design:

Difference-in-differences (DiD)

## Inputs

- metadata/metadata_final.csv (wastewater subset)
- results/ARG_dataset.csv
- plant-level timeline notes (upgrade date by plant)

## Outputs

- results/causal/wastewater_subset.csv
- results/causal/did_input_table.csv
- results/causal/did_model_summary.txt
- results/causal/robustness_checks.csv
- analysis/figures/

## Workflow

1. Build wastewater analysis table.

Add columns:

- plant_id
- time
- treatment_before_after (0 before, 1 after)
- upgraded_plant (1 treated plant, 0 control plant)

2. Select controls.

Choose plants with no treatment upgrade in same time window.

3. Fit DiD model.

```python
import statsmodels.formula.api as smf

model = smf.ols(
	"ARG_total ~ treatment_before_after + time + treatment_before_after:time",
	data=df,
).fit()
```

Key effect: `treatment_before_after:time`

4. Run robustness checks.

Control for:

- read depth
- seasonality
- pipeline differences

## Figures

1. ARG levels before versus after upgrades.
2. Treated versus control trend plot.
3. DiD coefficient with confidence interval.
4. Country-level comparison.

## Expected signal

If upgrades work, treated plants should show larger ARG reductions than controls after intervention.

