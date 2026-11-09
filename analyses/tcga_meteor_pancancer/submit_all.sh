#!/usr/bin/env bash
# Select the TCGA cohorts with >= MIN_SAMPLES samples having both 450K methylation and HiSeqV2 expression (Xena),
# then submit the pipeline for each of them. Run on a login node with internet access:  ./submit_all.sh  (add --yes to skip the question)
set -euo pipefail
# folder of this script; when submitted with sbatch, SLURM runs a copy from its spool folder, so use the submit folder
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -f "$HERE/config/pipeline.conf" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then HERE=$SLURM_SUBMIT_DIR; fi
[ -f "$HERE/config/pipeline.conf" ] || { echo "config/pipeline.conf not found next to this script (run it from the pipeline folder)"; exit 1; }
export PIPE_DIR=$HERE
source "$PIPE_DIR/config/pipeline.conf"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
mkdir -p "$WORK_DIR"
echo "[$(date +%T)] pipeline: $PIPE_DIR | work dir: $WORK_DIR"
if ! curl -s -m 30 -o /dev/null https://tcga.xenahubs.net; then
  echo "[$(date +%T)] ERROR: this node cannot reach the Xena server (https://tcga.xenahubs.net). The downloads need internet."; exit 1
fi
echo "[$(date +%T)] Xena server reachable; selecting cohorts (one line per cohort follows)"
Rscript "$PIPE_DIR/scripts/00_select_cohorts.R" --min_samples=$MIN_SAMPLES --out="$WORK_DIR/cohorts_selected.txt"
echo "[$(date +%T)] Selected: $(tr '\n' ' ' < "$WORK_DIR/cohorts_selected.txt")"
if [ "${1:-}" = "--yes" ]; then ok=y; elif [ -t 0 ]; then read -r -p "Submit these cohorts? [y/N] " ok; else echo "not interactive: rerun with --yes to submit"; ok=n; fi
[ "$ok" = y ] && "$PIPE_DIR/submit_cohort.sh" $(cat "$WORK_DIR/cohorts_selected.txt") || true
