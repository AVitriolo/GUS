#!/usr/bin/env bash
#SBATCH --job-name=meteor_submit
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/submit_all_%j.out
# Selects the cohorts and submits the whole pipeline from a compute node. From the pipeline folder:
#   mkdir -p logs && sbatch submit_all.sbatch
#   tail -f logs/submit_all_<jobid>.out        (progress, one line per step / cohort)
cd "$SLURM_SUBMIT_DIR"
bash ./submit_all.sh --yes
