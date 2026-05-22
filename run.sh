#!/bin/bash
#SBATCH --job-name=lr_mag_pipeline
#SBATCH --output=logs/lr_mag_pipeline_%j.out
#SBATCH --error=logs/lr_mag_pipeline_%j.err
#SBATCH --cpus-per-task=<CPUS>
#SBATCH --mem=<RAM>

mkdir -p logs

snakemake \
  --use-singularity \
  --cores "${SLURM_CPUS_PER_TASK}" \
  -p \
  all
