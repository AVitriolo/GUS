#!/usr/bin/env bash
# Submit the whole pipeline for one or more TCGA cohorts, e.g.:  ./submit_cohort.sh BRCA LUAD COAD
# Each cohort: prepare job (download + preprocessing) -> it submits the gene job array itself.
# Rerunning is safe: finished downloads, preprocessing and genes are skipped.
set -euo pipefail
# folder of this script; when submitted with sbatch, SLURM runs a copy from its spool folder, so use the submit folder
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -f "$HERE/config/pipeline.conf" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then HERE=$SLURM_SUBMIT_DIR; fi
[ -f "$HERE/config/pipeline.conf" ] || { echo "config/pipeline.conf not found next to this script (run it from the pipeline folder)"; exit 1; }
export PIPE_DIR=$HERE
source "$PIPE_DIR/config/pipeline.conf"
EXTRA=""; [ -n "${SLURM_PARTITION:-}" ] && EXTRA="$EXTRA --partition=$SLURM_PARTITION"
[ -n "${SLURM_ACCOUNT:-}" ] && EXTRA="$EXTRA --account=$SLURM_ACCOUNT"; [ -n "${SLURM_QOS:-}" ] && EXTRA="$EXTRA --qos=$SLURM_QOS"
for COHORT in "$@"; do
  mkdir -p "$WORK_DIR/$COHORT/logs"
  sbatch --job-name=prep_$COHORT $EXTRA --cpus-per-task=$PREP_CPUS --mem=$PREP_MEM --time=$PREP_TIME \
    --output="$WORK_DIR/$COHORT/logs/prepare_%j.out" --export=ALL,PIPE_DIR=$PIPE_DIR,WORK_DIR=$WORK_DIR "$PIPE_DIR/slurm/prepare.sbatch" $COHORT
done
