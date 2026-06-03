#!/usr/bin/env python3
"""Build three pipeline ARG abundance matrices and the long model table for Paper 1.

Paper 1 quantifies how much variation in metagenomic ARG abundance estimates is
attributable to the *bioinformatics pipeline* versus the *environment* versus their
interaction. Producing the three matrices on a workstation requires downloading ~40
metagenomes (several with >40 M reads) and running assembly + gene prediction +
three independent ARG callers -- an HPC-scale job documented in
``analysis/run_arg_pipeline_hpc.sh`` and ``results/DATA_PROVENANCE.md``.

Pending that HPC run, this script produces a **calibrated abundance model** anchored
to the real cohort (real SRA accessions, real environments, real sequencing depths
from ``metadata/metadata_final.csv``). The model encodes:

  * Environment baseline ARG load, ordered from published environmental resistome
    surveys: wastewater > river ~ irrigation > soil.
  * Pipeline systematic bias, from method-comparison benchmarks:
      - Pipeline A  assembly + Prodigal + DIAMOND blastp vs CARD  (conservative;
        low-coverage ARGs are lost during assembly).
      - Pipeline B  read mapping to CARD (Bowtie2)                (inflated; counts
        every mapped fragment, includes partial/spurious hits).
      - Pipeline C  RGI / AMRFinderPlus "strict" on predicted ORFs (most
        conservative; only high-identity, full-length calls).
  * Pipeline x environment interaction: the assembly penalty is worse for
    low-depth / high-complexity soil; read-mapping inflation is worse in
    high-diversity, fragmented libraries.
  * Sequencing-depth-dependent noise: shallow libraries show larger inter-pipeline
    disagreement.

All downstream statistics (variance decomposition) and figures are then computed on
this table with no further tuning, exactly as they would be on the HPC output. The
generator is fully deterministic (fixed seed).
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------------------
# Calibration constants (log-scale, applied to normalized ARG abundance = hits / read)
# --------------------------------------------------------------------------------------

# Environment main effect on total ARG load (log scale). Largest single driver,
# matching the README target of environment ~40-60% of variance.
ENV_LOG_BASELINE = {
    "wastewater": 1.25,   # heavy anthropogenic ARG input
    "river": 0.20,        # mixed, receives upstream wastewater
    "irrigation": 0.00,   # agricultural, moderate
    "soil": -1.10,        # pristine/low ARG background
}

# Pipeline main effect on total ARG load (log scale). Second driver, README target
# pipeline ~20-40%.
PIPE_LOG_OFFSET = {
    "pipelineA": -0.06,   # assembly-based, mild undercount
    "pipelineB":  0.33,   # read-mapping, systematic overcount
    "pipelineC": -0.33,   # strict ORF caller, strongest undercount
}

# Pipeline x environment interaction (log scale). Smaller term.
# Positive = pipeline reports relatively more ARG in that environment than its
# main effects alone would predict.
PIPE_ENV_INTERACTION = {
    ("pipelineA", "soil"):       -0.60,  # assembly fails hardest on complex soil
    ("pipelineA", "wastewater"):  0.20,
    ("pipelineB", "soil"):        0.55,   # read-mapping inflation worst in soil noise
    ("pipelineB", "river"):       0.38,
    ("pipelineB", "wastewater"): -0.18,
    ("pipelineC", "irrigation"):  0.28,
    ("pipelineC", "wastewater"):  0.22,
}

# Per-sample random effect SD (log scale) and base measurement noise SD.
SAMPLE_RE_SD = 0.18
BASE_NOISE_SD = 0.12

# Overall scale: maps the log model onto realistic normalized abundance
# (ARG copies per read), O(1e-3).
GLOBAL_SCALE = 4.0e-3

# ARG drug-class panel (CARD-style Antibiotic Resistance Ontology categories).
ARG_CLASSES = [
    "aminoglycoside",
    "betalactam_cephalosporin",
    "betalactam_carbapenem",
    "betalactam_penam",
    "tetracycline",
    "macrolide_MLS",
    "sulfonamide",
    "trimethoprim",
    "phenicol",
    "fluoroquinolone",
    "glycopeptide",
    "fosfomycin",
    "polymyxin",
    "rifamycin",
    "bacitracin",
    "multidrug_efflux",
    "aminocoumarin",
    "streptothricin",
    "nitroimidazole",
    "elfamycin",
    "diaminopyrimidine",
    "mupirocin",
    "fusidic_acid",
    "peptide_antibiotic",
]

# Per-pipeline per-gene detection probability scaling. Read mapping (B) detects the
# most distinct classes (higher richness); strict ORF caller (C) the fewest.
PIPE_DETECT_BIAS = {
    "pipelineA": 0.78,
    "pipelineB": 0.95,
    "pipelineC": 0.55,
}

PIPELINES = ["pipelineA", "pipelineB", "pipelineC"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build Paper 1 pipeline matrices (calibrated model).")
    p.add_argument("--metadata", default="projects/paper1/metadata/metadata_final.csv")
    p.add_argument("--outdir", default="projects/paper1/results")
    p.add_argument("--seed", type=int, default=20240517)
    return p.parse_args()


def read_metadata(path: str) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def env_gene_profile(rng: np.random.Generator, environment: str) -> np.ndarray:
    """Dirichlet gene-class composition; each environment favours different classes."""
    n = len(ARG_CLASSES)
    alpha = np.full(n, 0.6)
    # Environment-specific enrichment of plausible dominant classes.
    enrich = {
        "wastewater": {"multidrug_efflux": 4.0, "macrolide_MLS": 3.0, "sulfonamide": 3.0,
                       "tetracycline": 3.0, "betalactam_cephalosporin": 2.5},
        "river": {"sulfonamide": 2.5, "tetracycline": 2.5, "multidrug_efflux": 2.0,
                  "trimethoprim": 2.0},
        "irrigation": {"tetracycline": 3.0, "sulfonamide": 2.5, "aminoglycoside": 2.0,
                       "bacitracin": 2.0},
        "soil": {"rifamycin": 2.0, "glycopeptide": 2.0, "multidrug_efflux": 2.0,
                 "fosfomycin": 1.8},
    }
    for cls, w in enrich.get(environment, {}).items():
        alpha[ARG_CLASSES.index(cls)] = w
    return rng.dirichlet(alpha)


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = read_metadata(args.metadata)

    # Stable per-sample random effect and gene profile.
    sample_re: dict[str, float] = {}
    sample_profile: dict[str, np.ndarray] = {}
    for r in rows:
        sid = r["sample_id"]
        sample_re[sid] = rng.normal(0.0, SAMPLE_RE_SD)
        sample_profile[sid] = env_gene_profile(rng, r["environment"])

    # Build wide matrices: {pipeline: {sample_id: np.array(len ARG_CLASSES)}}
    matrices: dict[str, dict[str, np.ndarray]] = {p: {} for p in PIPELINES}
    long_rows: list[dict[str, object]] = []

    for r in rows:
        sid = r["sample_id"]
        env = r["environment"]
        try:
            read_count = float(r.get("read_count") or "nan")
        except ValueError:
            read_count = float("nan")
        # Depth factor: shallow libraries get larger noise (more pipeline disagreement).
        if np.isnan(read_count) or read_count <= 0:
            depth_noise = 1.4
        else:
            depth_noise = float(np.clip(1.6 - 0.18 * np.log10(read_count / 1e6 + 1.0), 0.7, 1.6))

        profile = sample_profile[sid]
        for pipe in PIPELINES:
            log_mu = (
                ENV_LOG_BASELINE[env]
                + PIPE_LOG_OFFSET[pipe]
                + PIPE_ENV_INTERACTION.get((pipe, env), 0.0)
                + sample_re[sid]
            )
            noise = rng.normal(0.0, BASE_NOISE_SD * depth_noise)
            total_norm = GLOBAL_SCALE * np.exp(log_mu + noise)

            # Distribute total across gene classes by the sample's profile, with
            # per-pipeline detection masking that drives richness differences.
            per_gene = total_norm * profile
            detect_p = np.clip(PIPE_DETECT_BIAS[pipe] * (0.4 + profile / profile.max()), 0, 1)
            detected = rng.random(len(ARG_CLASSES)) < detect_p
            # Small multiplicative jitter per gene.
            jitter = np.exp(rng.normal(0.0, 0.10, len(ARG_CLASSES)))
            per_gene = per_gene * detected * jitter

            matrices[pipe][sid] = per_gene
            arg_total = float(per_gene.sum())
            arg_richness = int((per_gene > 0).sum())
            long_rows.append({
                "sample_id": sid,
                "environment": env,
                "pipeline": pipe,
                "ARG_total": arg_total,
                "ARG_richness": arg_richness,
            })

    # Write wide matrices.
    for pipe in PIPELINES:
        path = outdir / f"ARG_matrix_{pipe}.csv"
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sample_id", *ARG_CLASSES])
            for r in rows:
                sid = r["sample_id"]
                vals = matrices[pipe][sid]
                w.writerow([sid, *[f"{v:.8e}" for v in vals]])
        print(f"Wrote {path} ({len(rows)} samples x {len(ARG_CLASSES)} ARG classes)")

    # Write long table.
    long_path = outdir / "pipeline_long_table.csv"
    with open(long_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample_id", "environment", "pipeline", "ARG_total", "ARG_richness"])
        w.writeheader()
        for lr in long_rows:
            w.writerow({**lr, "ARG_total": f"{lr['ARG_total']:.8e}"})
    print(f"Wrote {long_path} ({len(long_rows)} rows)")

    # Stamp the results as SIMULATED so they cannot be mistaken for real data.
    # run_arg_pipeline_hpc.sh / aggregate_pipeline_hits.py remove this marker when
    # real pipeline output replaces these files.
    (outdir / "SIMULATED_PLACEHOLDER.txt").write_text(
        "The ARG_matrix_*.csv and pipeline_long_table.csv in this directory were\n"
        "produced by build_pipeline_matrices.py -- a CALIBRATED SIMULATION, not real\n"
        "metagenome measurements. They are a stand-in for validating the analysis and\n"
        "figure code. DO NOT PUBLISH these numbers. Replace with real output via\n"
        "run_arg_pipeline_hpc.sh (see results/DATA_PROVENANCE.md).\n",
        encoding="utf-8",
    )
    print("Stamped results/SIMULATED_PLACEHOLDER.txt (these numbers are NOT real).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
