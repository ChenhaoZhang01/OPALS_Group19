# Paper 1 — Running the REAL Three-Pipeline Quantification

This is the procedure that turns Paper 1 from a simulation template into a real,
publishable study. It must run on **Linux** (an HPC allocation or a cloud VM); it
cannot run on the Windows laptop, and several samples are too large for a laptop
regardless (up to ~59M read pairs).

When you finish this, the fabricated `results/*.csv` are overwritten with real
measurements, the `SIMULATED_PLACEHOLDER.txt` marker is auto-removed, and you update
the numbers in `paper.md` from the regenerated tables. **Whatever the real result is —
even if it contradicts the current draft — that is the finding you report.**

---

## 0. What you need

| Resource | Estimate |
|---|---|
| OS | Linux (HPC node or cloud VM) |
| CPU | 16+ cores recommended |
| RAM | ≥ 64 GB (MEGAHIT assembly of the 40–59M-read soil/irrigation samples is the bottleneck) |
| Transient disk | ~0.5–1 TB free (reads + assembly are deleted per-sample by default) |
| Wall time | ~100–300 CPU-hours total; a few hours per large sample, minutes for small ones |
| Network | Downloads ~100s of GB of FASTQ from SRA over the run |

The 41 accessions are in `projects/paper1/metadata/download_list.txt`; their depths are
in `metadata/metadata_final.csv` (read_count column) so you can prioritize/triage.

## 1. Create the environment

```bash
conda env create -f environment-linux.yml
conda activate resistome
# sanity check all tools resolve:
which prefetch fasterq-dump fastp megahit prodigal diamond bowtie2 samtools rgi
```

## 2. Build the CARD databases (once)

```bash
bash projects/paper1/analysis/setup_databases.sh --dbdir /databases/card
```

This downloads CARD, builds the DIAMOND protein DB (Pipeline A), the Bowtie2 nucleotide
index (Pipeline B), and loads `card.json` into RGI (Pipeline C). It writes
`/databases/card/DB_PATHS.env`.

## 3. Smoke-test on ONE small sample first

Always validate the whole chain on the smallest library before launching the cohort.
`SRR30403277` (~211k reads, France wastewater) is a good choice:

```bash
bash projects/paper1/analysis/run_arg_pipeline_hpc.sh \
  --card-dir /databases/card --outdir /scratch/paper1 \
  --threads 16 --only SRR30403277 --keep-intermediate
```

Confirm these exist and are non-empty:
- `/scratch/paper1/hits_A/SRR30403277.tsv`   (DIAMOND hits)
- `/scratch/paper1/hits_B/SRR30403277.idxstats`
- `/scratch/paper1/hits_C/SRR30403277.txt`   (RGI table)

Inspect `/scratch/paper1/logs/SRR30403277.log` if anything is empty. (RGI sometimes needs
`rgi load ... --local` to have been run in step 2 — confirm it was.)

## 4. Run the full cohort

```bash
bash projects/paper1/analysis/run_arg_pipeline_hpc.sh \
  --card-dir /databases/card --outdir /scratch/paper1 --threads 16
```

- **Resumable:** re-running skips any sample that already has all three hit files, so it
  is safe to relaunch after a timeout or node failure (e.g. inside a SLURM array job).
- Per-sample reads/assembly are deleted after aggregation (`--keep-intermediate` to keep).
- The script auto-runs the aggregator at the end, writing the real
  `ARG_matrix_pipeline{A,B,C}.csv` and `pipeline_long_table.csv` into
  `projects/paper1/results/` and removing `SIMULATED_PLACEHOLDER.txt`.

### Optional: SLURM array (one sample per task)

```bash
# submit one task per accession; each task runs --only "$ACC", then aggregate once at the end
sbatch --array=1-41 --cpus-per-task=16 --mem=64G --wrap '
  ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" projects/paper1/metadata/download_list.txt)
  bash projects/paper1/analysis/run_arg_pipeline_hpc.sh --card-dir /databases/card \
       --outdir /scratch/paper1 --threads 16 --only "$ACC"'
```

## 5. Re-run the analysis and figures (unchanged)

```bash
python3 projects/paper1/analysis/run_paper1_analysis.py
python3 projects/paper1/analysis/generate_figures.py
```

These read the real `pipeline_long_table.csv` and regenerate
`variance_decomposition.csv`, `richness_decomposition.csv`, `pipeline_summary.csv`,
`pipeline_concordance.csv`, and the four figures — no code changes needed.

## 6. Update the paper with the REAL numbers

1. Remove the `⚠ DRAFT — NOT FOR SUBMISSION` banner at the top of `paper.md`.
2. Replace every numeric value (abstract, Tables 1–4, Figure captions, Discussion,
   Conclusion) with the values from the regenerated CSVs. The narrative structure and
   figure layout do not change.
3. Update `results/DATA_PROVENANCE.md`: flip the "Raw three-pipeline quantification"
   gate to PASS and delete the "what is modeled" section.

## 7. Quality checks before claiming results

- Did all 41 samples produce hits? Check `wc -l projects/paper1/results/pipeline_long_table.csv`
  is `124` (41×3 + header). Investigate any sample missing from a pipeline (look in `logs/`).
- Are amplicon/16S or otherwise non-shotgun libraries in the cohort? ARG read-mapping on
  amplicon data is meaningless — exclude such samples and note it (cross-check the
  BioProject library strategy on SRA).
- Sanity-check magnitudes against literature (wastewater should be the most ARG-rich).
- Report the CARD version and all tool versions in Methods.

---

**Bottom line:** the code, cohort, and paper scaffold are ready. The only thing standing
between this draft and a real submission is executing steps 2–6 on real compute. Until
then the results are a simulation and must not be submitted.
