#!/bin/bash
#SBATCH --job-name=snakemake_single
#SBATCH --partition=general
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=168:00:00
#SBATCH --output=snakemake_%j.out

# Activate env
source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake-env

# Run locally using the allocated resources
snakemake --configfile config/config.yaml --snakefile workflow/Snakefile --use-conda --cores 16 

