#!/usr/bin/env bash
#SBATCH --job-name=analyses
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# All the WGBS vs TCGA analyses on a set of cohorts, written to new folders ending in SUFFIX (nothing is overwritten):
#   saturation (target-only and footprint WGBS files), probe overlap (both), never-selected probes, confirm_analyses.
# Usage (from the pipeline folder): sbatch --output=analyses_%j.out slurm/rerun_analyses.sbatch <SUFFIX> COHORT [COHORT ...]
set -euo pipefail
SUFFIX=$1; shift
PIPE_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}; source "$PIPE_DIR/config/pipeline.conf"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
cd "$PIPE_DIR"
G=$(dirname "$WORK_DIR")    # where the WGBS files were copied (GUS_PANCANCER/)
for f in wgbs_prostate_functional_allcpgs wgbs_prostate_functional_footprint; do
  echo "[$(date)] saturation + overlap: $f"
  python scripts/saturation_plot.py "$WORK_DIR" --out_suffix=$SUFFIX "--first=Prostate WGBS=$G/$f.tsv.gz:342" "$@"
  python scripts/probe_overlap.py "$WORK_DIR" --out_suffix=$SUFFIX "--first=Prostate WGBS=$G/$f.tsv.gz" "$@"
done
echo "[$(date)] never selected"; python scripts/never_selected.py "$WORK_DIR" --out_suffix=$SUFFIX "$@"
echo "[$(date)] confirm";        python scripts/confirm_analyses.py "$WORK_DIR" --wgbs=$G/wgbs_prostate_functional_footprint.tsv.gz --out_suffix=$SUFFIX "$@"
echo "[$(date)] done. Folders: $WORK_DIR/{saturation_*,probe_overlap_*,never_selected,confirm}$SUFFIX"
